"""One continuous weighted MIP search with observed flow-proof timing.

The first-phase incumbent and bound are the last observations at or before its
runtime cap. A later callback never replaces that reported incumbent. The
search tree, objective, and focus remain unchanged throughout the extension.
"""

import math

from gurobipy import GRB

from static_lexicographic import CERTIFICATE_GAP_MARGIN, _integer_value
from static_weighted_certification import (
    RELIABLE_FINISHED_STATUSES, WEIGHTED_BOUND_MARGIN, _safe_weight_instance,
    gap_flow_certificate, safe_weight_parameters,
)


class SafeWeightedContinuation:
    def __init__(self, flow_expr, movement_expr, flow_weight, context,
                 extract_callback_metrics=None):
        self.context = dict(context)
        self.flow_weight = flow_weight
        self.phase_limit = float(context["weighted_time_limit"])
        self.extension_limit = float(context["extension_time_limit"])
        self.stop_on_flow_proof = context.get("stop_on_flow_proof", True)
        if not isinstance(self.stop_on_flow_proof, bool):
            raise ValueError("stop_on_flow_proof must be boolean")
        self.stop_at_flow_proof = context.get("stop_at_flow_proof", False)
        if not isinstance(self.stop_at_flow_proof, bool):
            raise ValueError("stop_at_flow_proof must be boolean")
        self.load_count, self.distance_sum, self.distance_max = _safe_weight_instance(
            self.context["targets"], self.context["outputs"],
            self.context["cell_count"], self.context["escort_count"])
        self.flow_criteria = {}
        self.flow_vars = [flow_expr.getVar(i) for i in range(flow_expr.size())]
        self.flow_coefficients = [flow_expr.getCoeff(i) for i in range(flow_expr.size())]
        self.flow_constant = flow_expr.getConstant()
        self.movement_vars = [movement_expr.getVar(i) for i in range(movement_expr.size())]
        self.movement_coefficients = [movement_expr.getCoeff(i) for i in range(movement_expr.size())]
        self.movement_constant = movement_expr.getConstant()
        self.extract_callback_metrics = extract_callback_metrics
        self.best_incumbent = None
        self.best_live_incumbent = None
        self.best_flow_candidate = None
        self.best_observed_weighted_objective = None
        self.last_checkpoint = None
        self.best_phase_bound = None
        self.best_phase_bound_runtime = None
        self.phase1_snapshot = None
        self.phase_transition_runtime = None
        self.extension_used = False
        self.extension_announced = False
        self.stop_reason = ""
        self.counterexample = None
        self.first_flow_proof = None
        self.first_flow_proof_incumbent = None
        self.flow_proof_invalidated = False
        self.flow_proof_invalidation_reason = ""
        self.best_flow_proof_bound = None
        self.last_flow_proof_checked_flowtime = None
        self.flow_proof_check_count = 0
        self.error = None

    @staticmethod
    def _finite(value):
        return math.isfinite(value) and abs(value) < GRB.INFINITY

    def _read_candidate(self, model, runtime):
        values = model.cbGetSolution(self.flow_vars)
        flow = _integer_value(self.flow_constant + sum(
            coefficient * value for coefficient, value in zip(self.flow_coefficients, values)))
        reported_objective = model.cbGet(GRB.Callback.MIPSOL_OBJ)
        if (self.best_incumbent is not None and self.best_flow_candidate is not None
                and reported_objective >= self.best_incumbent["scaled_objective"] - CERTIFICATE_GAP_MARGIN
                and flow >= self.best_flow_candidate["flowtime"]):
            return None
        # Round F and M separately, just as the final solution extractor does.
        # R times a tiny F residual can be materially larger than the tolerance
        # for checking the already-scaled MIPSOL objective's integrality.
        movement_values = model.cbGetSolution(self.movement_vars)
        movements = _integer_value(self.movement_constant + sum(
            coefficient * value for coefficient, value in zip(self.movement_coefficients, movement_values)))
        objective = self.flow_weight * flow + movements
        if not self._finite(reported_objective) or abs(reported_objective - objective) > (self.flow_weight + 1) * 1e-4 + CERTIFICATE_GAP_MARGIN:
            raise ValueError("Callback weighted objective disagrees with its flow and movement expressions")
        if flow < self.distance_sum or movements < self.distance_sum:
            raise ValueError("Callback objective components are below the target-distance lower bound")
        candidate = dict(has_solution=True, flowtime=flow, movements=movements,
                         scaled_objective=objective, incumbent_runtime=runtime)
        if self.best_observed_weighted_objective is None or objective < self.best_observed_weighted_objective:
            self.best_observed_weighted_objective = objective
        if self.extract_callback_metrics is not None:
            candidate.update(self.extract_callback_metrics(model))
        if (self.best_live_incumbent is None
                or objective < self.best_live_incumbent["scaled_objective"]):
            self.best_live_incumbent = candidate
        if self.best_flow_candidate is None or flow < self.best_flow_candidate["flowtime"]:
            self.best_flow_candidate = candidate.copy()
        return candidate

    def _snapshot(self, runtime, status):
        checkpoint = self.last_checkpoint or {}
        snapshot = dict(has_solution=False, flowtime=None, movements=None, makespan=None,
                        objective=None, best_bound=None, absolute_gap=None,
                        scaled_objective=None, scaled_best_bound=self.best_phase_bound,
                        scaled_absolute_gap=None, weighted_proven=False, bound_consistent=True,
                        runtime=runtime, status_name=status, work=checkpoint.get("work"),
                        node_count=checkpoint.get("node_count"),
                        bound_checkpoint_runtime=self.best_phase_bound_runtime,
                        statistics_checkpoint_runtime=checkpoint.get("runtime"),
                        weight_scale=self.flow_weight, snapshot_source="CALLBACK_OBSERVATIONS")
        if self.best_incumbent is not None:
            snapshot.update(self.best_incumbent)
            snapshot["objective"] = snapshot["scaled_objective"] / self.flow_weight
        bound = snapshot["scaled_best_bound"]
        if bound is not None:
            snapshot["best_bound"] = bound / self.flow_weight
            if snapshot["has_solution"]:
                gap = abs(snapshot["scaled_objective"] - bound)
                snapshot.update(scaled_absolute_gap=gap, absolute_gap=gap / self.flow_weight,
                                bound_consistent=bound <= snapshot["scaled_objective"] + CERTIFICATE_GAP_MARGIN)
                snapshot["weighted_proven"] = snapshot["bound_consistent"] and gap < 1 - CERTIFICATE_GAP_MARGIN
        return snapshot

    def _proof(self, incumbent, bound, status="TIME_LIMIT"):
        candidate = dict(incumbent or {}, scaled_best_bound=bound, status_name=status,
                         weight_scale=self.flow_weight)
        # A first-phase consistency flag pertains to its old bound. Recheck the
        # supplied bound against the same frozen feasible objective instead.
        candidate.pop("bound_consistent", None)
        proof = gap_flow_certificate(
            candidate, self.context["targets"], self.context["outputs"],
            self.context["cell_count"], self.context["escort_count"],
            self.context["physical_horizon"], self.flow_weight)
        if (candidate.get("has_solution") and status in RELIABLE_FINISHED_STATUSES
                and bound is not None and proof["reason"] in {"", "WEIGHTED_GAP_TOO_LARGE"}):
            objective = self.flow_weight * candidate["flowtime"] + candidate["movements"]
            if self._finite(bound) and bound <= objective + CERTIFICATE_GAP_MARGIN and abs(objective - bound) < 1 - CERTIFICATE_GAP_MARGIN:
                params = safe_weight_parameters(
                    self.context["targets"], self.context["outputs"],
                    self.context["cell_count"], self.context["escort_count"], candidate["flowtime"])
                if self.flow_weight >= params["flow_weight"] and self.context["physical_horizon"] >= params["flow_horizon"]:
                    proof.update(flow_proven=True, proof_source="weighted_optimum",
                                 reason="", flow_lower_bound=candidate["flowtime"])
        return proof

    def _freeze(self, runtime):
        if self.phase1_snapshot is None:
            self.phase1_snapshot = self._snapshot(self.phase_limit, "BUDGET_REACHED")
            self.phase_transition_runtime = runtime
            method = self._cached_flow_proof_method(self.phase1_snapshot, self.phase1_snapshot.get("scaled_best_bound"))
            self.phase1_snapshot.update(flow_proven=bool(method), flow_proof_source=method or "")
            self.extension_used = (
                bool(self.phase1_snapshot.get("has_solution") and not method)
                if self.stop_on_flow_proof else self.extension_limit > 0)

    def _build_flow_criterion(self, flow):
        """Only F changes either mathematical certificate threshold."""
        better_horizon = max(0, flow - 1 - self.distance_sum + self.distance_max)
        flow_horizon = max(0, flow - self.distance_sum + self.distance_max)
        return dict(
            distance_optimal=flow == self.distance_sum,
            gap_horizon_sufficient=self.context["physical_horizon"] >= better_horizon,
            bound_threshold=(self.flow_weight * (flow - 1) + self.load_count * better_horizon
                             + WEIGHTED_BOUND_MARGIN + CERTIFICATE_GAP_MARGIN),
            weighted_optimum_scope=(self.flow_weight >= self.load_count * flow_horizon - self.distance_sum + 1
                                    and self.context["physical_horizon"] >= flow_horizon),
        )

    def _cached_flow_proof_method(self, incumbent, bound, status="TIME_LIMIT"):
        """Use scalar comparisons with a criterion cached by incumbent F.

        The complete public certificate helper revalidates saved evidence at
        finalization. Callback checks never traverse instance coordinates or
        retrieve a solver solution vector.
        """
        if not incumbent or not incumbent.get("has_solution") or status not in RELIABLE_FINISHED_STATUSES:
            return None
        flow, movements = incumbent.get("flowtime"), incumbent.get("movements")
        if (not isinstance(flow, (int, float)) or not isinstance(movements, (int, float))
                or not math.isfinite(flow) or not math.isfinite(movements)
                or flow < self.distance_sum or movements < self.distance_sum
                or int(flow) != flow or int(movements) != movements
                or incumbent.get("weight_scale", self.flow_weight) != self.flow_weight):
            return None
        objective = self.flow_weight * flow + movements
        reported = incumbent.get("scaled_objective")
        if reported is not None and (not self._finite(reported) or abs(reported - objective) > CERTIFICATE_GAP_MARGIN):
            return None
        criterion = self.flow_criteria.get(flow)
        if criterion is None:
            criterion = self._build_flow_criterion(flow)
            self.flow_criteria[flow] = criterion
        bound_valid = bound is not None and self._finite(bound) and bound <= objective + CERTIFICATE_GAP_MARGIN
        if (bound_valid and criterion["weighted_optimum_scope"]
                and abs(objective - bound) < 1 - CERTIFICATE_GAP_MARGIN):
            return "weighted_optimum"
        if criterion["distance_optimal"]:
            return "distance_bound"
        if bound_valid and criterion["gap_horizon_sufficient"] and bound > criterion["bound_threshold"]:
            return "weighted_gap"
        return None

    def _invalidate_first_flow_proof(self, reason):
        self.first_flow_proof = None
        self.first_flow_proof_incumbent = None
        self.flow_proof_invalidated = True
        self.flow_proof_invalidation_reason = reason

    def _observe_flow_proof(self, incumbent, bound, checkpoint, source, *,
                            force=False, status="TIME_LIMIT"):
        """Time the first certificate at an actual supported solver observation.

        A new best flow or a strictly stronger finite bound can change the
        certificate. All other callback observations skip criterion work.
        The strongest observed bound is retained before the first candidate.
        """
        objective = incumbent.get("scaled_objective") if incumbent and incumbent.get("has_solution") else None
        if self.best_observed_weighted_objective is not None:
            objective = (self.best_observed_weighted_objective if objective is None else
                         min(objective, self.best_observed_weighted_objective))
        if self.first_flow_proof is not None:
            if incumbent is not None and incumbent.get("has_solution") and incumbent["flowtime"] < self.first_flow_proof["flowtime"]:
                self._invalidate_first_flow_proof("BETTER_FLOW_AFTER_RECORDED_PROOF")
            elif (objective is not None and self.first_flow_proof["scaled_bound"] is not None
                    and self.first_flow_proof["scaled_bound"] > objective + CERTIFICATE_GAP_MARGIN):
                self._invalidate_first_flow_proof("INCONSISTENT_RECORDED_FLOW_BOUND")
            return
        if self.flow_proof_invalidated:
            return
        runtime, nodes = checkpoint["runtime"], checkpoint.get("node_count")
        # Discard an incompatible pre-candidate bound rather than letting it
        # hide later valid bound improvements. No comparison uses an epsilon
        # for deciding whether a valid bound improved.
        if (objective is not None and self.best_flow_proof_bound is not None
                and self.best_flow_proof_bound > objective + CERTIFICATE_GAP_MARGIN):
            self.best_flow_proof_bound = None
        bound_valid = (bound is not None and self._finite(bound)
                       and (objective is None or bound <= objective + CERTIFICATE_GAP_MARGIN))
        bound_improved = bound_valid and (self.best_flow_proof_bound is None or bound > self.best_flow_proof_bound)
        if bound_improved:
            self.best_flow_proof_bound = bound
        if not incumbent or not incumbent.get("has_solution"):
            return
        flow = incumbent["flowtime"]
        flow_improved = self.last_flow_proof_checked_flowtime is None or flow < self.last_flow_proof_checked_flowtime
        if not force and not bound_improved and not flow_improved:
            return
        self.last_flow_proof_checked_flowtime = flow
        self.flow_proof_check_count += 1
        proof_bound = self.best_flow_proof_bound
        method = self._cached_flow_proof_method(incumbent, proof_bound, status)
        if method:
            self.first_flow_proof_incumbent = dict(incumbent)
            self.first_flow_proof = dict(
                runtime=runtime, node_count=nodes, work=checkpoint.get("work"),
                flowtime=incumbent["flowtime"], scaled_bound=proof_bound,
                source=source, method=method,
            )

    def _stop(self, model, reason):
        self.stop_reason = reason
        model.terminate()

    def _stop_at_proved_flow(self, model):
        # Reuse the recorded event-driven proof. A nonimproving MIPSOL can
        # prove a witness's F without making it the weighted incumbent, so
        # only stop when the solution that will be reported has that same F.
        if (self.stop_at_flow_proof and self.first_flow_proof is not None
                and self.best_live_incumbent is not None
                and self.best_live_incumbent["flowtime"] == self.first_flow_proof["flowtime"]):
            self._stop(model, "FLOW_PROVEN_EARLY")
            return True
        return False

    def __call__(self, model, where):
        if self.error is not None or self.stop_reason or where not in (GRB.Callback.MIP, GRB.Callback.MIPSOL):
            return
        try:
            runtime = model.cbGet(GRB.Callback.RUNTIME)
            prefix = "MIP" if where == GRB.Callback.MIP else "MIPSOL"
            bound = model.cbGet(getattr(GRB.Callback, prefix + "_OBJBND"))
            if not self._finite(bound):
                bound = None
            checkpoint = dict(runtime=runtime, scaled_best_bound=bound,
                              node_count=model.cbGet(getattr(GRB.Callback, prefix + "_NODCNT")),
                              work=model.cbGet(GRB.Callback.WORK))
            before_cutoff = runtime <= self.phase_limit
            if not before_cutoff:
                # Freeze before reading a candidate first observed after the cap.
                self._freeze(runtime)
            candidate = self._read_candidate(model, runtime) if where == GRB.Callback.MIPSOL else None
            self._observe_flow_proof(
                self.best_flow_candidate, bound, checkpoint,
                "CALLBACK_" + prefix)
            if before_cutoff:
                self.last_checkpoint = checkpoint
                if bound is not None and (self.best_phase_bound is None or bound > self.best_phase_bound):
                    self.best_phase_bound = bound
                    self.best_phase_bound_runtime = runtime
                if candidate is not None and (self.best_incumbent is None or
                        candidate["scaled_objective"] < self.best_incumbent["scaled_objective"]):
                    self.best_incumbent = candidate
                self._stop_at_proved_flow(model)
                return
            if self._stop_at_proved_flow(model):
                return
            if not self.phase1_snapshot["has_solution"]:
                if self.stop_on_flow_proof:
                    self._stop(model, "NO_PHASE1_SOLUTION")
                return
            target_flow = self.phase1_snapshot["flowtime"]
            if self.best_flow_candidate is not None and self.best_flow_candidate["flowtime"] < target_flow:
                self.counterexample = self.best_flow_candidate.copy()
                if self.stop_on_flow_proof:
                    self._stop(model, "BETTER_FLOW_FOUND")
                return
            if self.phase1_snapshot["flow_proven"]:
                if self.stop_on_flow_proof:
                    self._observe_flow_proof(
                        self.phase1_snapshot, self.phase1_snapshot.get("scaled_best_bound"),
                        checkpoint, "CALLBACK_" + prefix, force=True)
                    self._stop(model, "FLOW_PROVEN_AT_PHASE_LIMIT")
                return
            if (self.stop_on_flow_proof and self.first_flow_proof is not None
                    and self.first_flow_proof["flowtime"] == target_flow):
                self._stop(model, "FLOW_PROVEN_DURING_EXTENSION" if self.extension_used else "FLOW_PROVEN_AT_PHASE_LIMIT")
                return
            if not self.extension_announced:
                self.extension_announced = True
                purpose = f"to certify fixed flow {target_flow}" if self.stop_on_flow_proof else "for both objectives"
                print(f"Weighted phase ended at {self.phase_limit:g}s; continuing the same search "
                      f"{purpose} (total cap {self.phase_limit + self.extension_limit:g}s).", flush=True)
        except Exception as error:
            self.error = error
            model.terminate()

    def raise_if_failed(self):
        if self.error is not None:
            raise RuntimeError("Continuous weighted-search callback failed") from self.error

    def finalize(self, result):
        """Return metadata while retaining the actual final solver result."""
        self.raise_if_failed()
        runtime = result["runtime"]
        if self.phase1_snapshot is None and runtime <= self.phase_limit:
            self.phase1_snapshot = dict(result, snapshot_source="SOLVE_FINISHED",
                                        bound_checkpoint_runtime=runtime)
        elif self.phase1_snapshot is None:
            self._freeze(None)
        snapshot = self.phase1_snapshot
        final_bound = result.get("scaled_best_bound")
        status = result.get("status_name")
        # The final weighted incumbent can have a larger F than a nonincumbent
        # best-flow witness. Its smaller objective still constrains every
        # valid weighted lower bound, regardless of which witness we certify.
        if result.get("has_solution"):
            final_objective = result.get("scaled_objective")
            flow, movements = result.get("flowtime"), result.get("movements")
            if (isinstance(flow, (int, float)) and isinstance(movements, (int, float))
                    and math.isfinite(flow) and math.isfinite(movements)
                    and flow >= self.distance_sum and movements >= self.distance_sum
                    and int(flow) == flow and int(movements) == movements
                    and result.get("weight_scale", self.flow_weight) == self.flow_weight
                    and final_objective is not None and self._finite(final_objective)
                    and abs(final_objective - (self.flow_weight * flow + movements)) <= CERTIFICATE_GAP_MARGIN
                    and (self.best_observed_weighted_objective is None or final_objective < self.best_observed_weighted_objective)):
                self.best_observed_weighted_objective = final_objective
        final_candidate = result if result.get("has_solution") else None
        if (self.best_flow_candidate is not None and
                (final_candidate is None or self.best_flow_candidate["flowtime"] < final_candidate["flowtime"])):
            final_candidate = self.best_flow_candidate
        if status not in RELIABLE_FINISHED_STATUSES:
            if self.first_flow_proof is not None:
                self._invalidate_first_flow_proof("UNRELIABLE_FINAL_STATUS")
        else:
            self._observe_flow_proof(
                final_candidate, final_bound,
                dict(runtime=runtime, node_count=result.get("node_count"), work=result.get("work")),
                "SOLVE_FINISHED", force=True, status=status)
            if self.first_flow_proof is not None:
                saved = self._proof(self.first_flow_proof_incumbent,
                                    self.first_flow_proof["scaled_bound"], status)
                if not saved["flow_proven"]:
                    self._invalidate_first_flow_proof("RECORDED_FLOW_PROOF_FAILED_VALIDATION")
        initial_proof = self._proof(snapshot, final_bound, status)
        final_proof = self._proof(result, final_bound, status)
        if result.get("has_solution") and snapshot.get("has_solution") and result["flowtime"] < snapshot["flowtime"]:
            self.counterexample = dict(has_solution=True, flowtime=result["flowtime"],
                                      movements=result["movements"], makespan=result.get("makespan"),
                                      scaled_objective=result["scaled_objective"])
        if (self.best_flow_candidate is not None and snapshot.get("has_solution")
                and self.best_flow_candidate["flowtime"] < snapshot["flowtime"]
                and (self.counterexample is None or self.best_flow_candidate["flowtime"] < self.counterexample["flowtime"])):
            self.counterexample = self.best_flow_candidate.copy()
        if self.counterexample is not None:
            initial_proof.update(flow_proven=False, proof_source="", reason="BETTER_FLOW_FOUND")
        phase1_proof = self._proof(snapshot, snapshot.get("scaled_best_bound"), status)
        if self.counterexample is not None:
            phase1_proof.update(flow_proven=False, proof_source="", reason="BETTER_FLOW_FOUND")
        if self.flow_proof_invalidated and self.flow_proof_invalidation_reason != "UNRELIABLE_FINAL_STATUS":
            for proof in (initial_proof, final_proof, phase1_proof):
                proof.update(flow_proven=False, proof_source="", reason=self.flow_proof_invalidation_reason)
        snapshot["weighted_proven"] = bool(snapshot.get("weighted_proven") and status in RELIABLE_FINISHED_STATUSES)
        snapshot.update(flow_proven=phase1_proof["flow_proven"],
                        flow_proof_source=phase1_proof["proof_source"])
        stop_reason = self.stop_reason
        if not stop_reason:
            if not self.stop_on_flow_proof:
                stop_reason = "WEIGHTED_OPTIMUM" if result.get("weighted_proven") else status or "UNKNOWN"
            elif self.counterexample is not None:
                stop_reason = "BETTER_FLOW_FOUND"
            elif result.get("weighted_proven"):
                stop_reason = "WEIGHTED_OPTIMUM"
            elif initial_proof["flow_proven"]:
                stop_reason = "FLOW_PROVEN_AT_STOP"
            elif not snapshot.get("has_solution"):
                stop_reason = "NO_PHASE1_SOLUTION"
            else:
                stop_reason = status or "UNKNOWN"
        first_proof = self.first_flow_proof or {}
        return dict(
            phase1_snapshot=snapshot,
            phase1_snapshot_missing=snapshot["snapshot_source"] != "SOLVE_FINISHED" and self.last_checkpoint is None,
            phase_transition_runtime=self.phase_transition_runtime,
            extension_used=self.extension_used,
            extension_runtime=max(0.0, runtime - self.phase_limit),
            phase_transition_delay=(None if self.phase_transition_runtime is None
                                    else max(0.0, self.phase_transition_runtime - self.phase_limit)),
            stop_reason=stop_reason, optimization_calls=1,
            flow_proven=initial_proof["flow_proven"], flow_proof_source=initial_proof["proof_source"],
            flow_gap_certificate=initial_proof,
            final_flow_proven=final_proof["flow_proven"], final_flow_proof_source=final_proof["proof_source"],
            final_flow_gap_certificate=final_proof,
            counterexample=self.counterexample is not None, flow_counterexample=self.counterexample,
            first_flow_proof_runtime=first_proof.get("runtime"),
            first_flow_proof_node_count=first_proof.get("node_count"),
            first_flow_proof_work=first_proof.get("work"),
            first_flow_proof_flowtime=first_proof.get("flowtime"),
            first_flow_proof_source=first_proof.get("source"),
            first_flow_proof_method=first_proof.get("method"),
            first_flow_proof_scaled_bound=first_proof.get("scaled_bound"),
            flow_proof_invalidated=self.flow_proof_invalidated,
            flow_proof_invalidation_reason=self.flow_proof_invalidation_reason,
        )
