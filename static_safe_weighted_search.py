"""One continuous weighted MIP search with a fixed first-phase flow target.

The first-phase incumbent and bound are the last observations at or before its
runtime cap. A later callback never replaces that reported incumbent. The
search tree, objective, and focus remain unchanged throughout the extension.
"""

import math

from gurobipy import GRB

from static_lexicographic import CERTIFICATE_GAP_MARGIN, _integer_value
from static_weighted_certification import (
    RELIABLE_FINISHED_STATUSES, gap_flow_certificate, safe_weight_parameters,
)


class SafeWeightedContinuation:
    def __init__(self, flow_expr, movement_expr, flow_weight, context,
                 extract_callback_metrics=None):
        self.context = dict(context)
        self.flow_weight = flow_weight
        self.phase_limit = float(context["weighted_time_limit"])
        self.extension_limit = float(context["extension_time_limit"])
        self.flow_vars = [flow_expr.getVar(i) for i in range(flow_expr.size())]
        self.flow_coefficients = [flow_expr.getCoeff(i) for i in range(flow_expr.size())]
        self.flow_constant = flow_expr.getConstant()
        self.movement_vars = [movement_expr.getVar(i) for i in range(movement_expr.size())]
        self.movement_coefficients = [movement_expr.getCoeff(i) for i in range(movement_expr.size())]
        self.movement_constant = movement_expr.getConstant()
        self.extract_callback_metrics = extract_callback_metrics
        self.best_incumbent = None
        self.best_flow_candidate = None
        self.last_checkpoint = None
        self.best_phase_bound = None
        self.best_phase_bound_runtime = None
        self.phase1_snapshot = None
        self.phase_transition_runtime = None
        self.extension_used = False
        self.extension_announced = False
        self.stop_reason = ""
        self.counterexample = None
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
        if flow < 0 or movements < 0:
            raise ValueError("Callback returned invalid integer weighted objective components")
        candidate = dict(has_solution=True, flowtime=flow, movements=movements,
                         scaled_objective=objective, incumbent_runtime=runtime)
        if self.extract_callback_metrics is not None:
            candidate.update(self.extract_callback_metrics(model))
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
        if candidate.get("has_solution") and status in RELIABLE_FINISHED_STATUSES and bound is not None:
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
            proof = self._proof(self.phase1_snapshot, self.phase1_snapshot.get("scaled_best_bound"))
            self.phase1_snapshot.update(flow_proven=proof["flow_proven"],
                                        flow_proof_source=proof["proof_source"])
            self.extension_used = bool(self.phase1_snapshot.get("has_solution") and not proof["flow_proven"])

    def _stop(self, model, reason):
        self.stop_reason = reason
        model.terminate()

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
            if before_cutoff:
                self.last_checkpoint = checkpoint
                if bound is not None and (self.best_phase_bound is None or bound > self.best_phase_bound):
                    self.best_phase_bound = bound
                    self.best_phase_bound_runtime = runtime
                if candidate is not None and (self.best_incumbent is None or
                        candidate["scaled_objective"] < self.best_incumbent["scaled_objective"]):
                    self.best_incumbent = candidate
                return
            if not self.phase1_snapshot["has_solution"]:
                self._stop(model, "NO_PHASE1_SOLUTION")
                return
            target_flow = self.phase1_snapshot["flowtime"]
            if self.best_flow_candidate is not None and self.best_flow_candidate["flowtime"] < target_flow:
                self.counterexample = self.best_flow_candidate.copy()
                self._stop(model, "BETTER_FLOW_FOUND")
                return
            if self.phase1_snapshot["flow_proven"]:
                self._stop(model, "FLOW_PROVEN_AT_PHASE_LIMIT")
                return
            proof = self._proof(self.phase1_snapshot, bound)
            if proof["flow_proven"]:
                self._stop(model, "FLOW_PROVEN_DURING_EXTENSION" if self.extension_used else "FLOW_PROVEN_AT_PHASE_LIMIT")
                return
            if not self.extension_announced:
                self.extension_announced = True
                print(f"Weighted phase ended at {self.phase_limit:g}s; continuing the same search "
                      f"to certify fixed flow {target_flow} (total cap "
                      f"{self.phase_limit + self.extension_limit:g}s).", flush=True)
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
        initial_proof = self._proof(snapshot, final_bound, status)
        final_proof = self._proof(result, final_bound, status)
        if result.get("has_solution") and snapshot.get("has_solution") and result["flowtime"] < snapshot["flowtime"]:
            self.counterexample = dict(has_solution=True, flowtime=result["flowtime"],
                                      movements=result["movements"], makespan=result.get("makespan"),
                                      scaled_objective=result["scaled_objective"])
        if self.counterexample is not None:
            initial_proof.update(flow_proven=False, proof_source="", reason="BETTER_FLOW_FOUND")
        phase1_proof = self._proof(snapshot, snapshot.get("scaled_best_bound"), status)
        snapshot["weighted_proven"] = bool(snapshot.get("weighted_proven") and status in RELIABLE_FINISHED_STATUSES)
        snapshot.update(flow_proven=phase1_proof["flow_proven"],
                        flow_proof_source=phase1_proof["proof_source"])
        stop_reason = self.stop_reason
        if not stop_reason:
            if self.counterexample is not None:
                stop_reason = "BETTER_FLOW_FOUND"
            elif result.get("weighted_proven"):
                stop_reason = "WEIGHTED_OPTIMUM"
            elif initial_proof["flow_proven"]:
                stop_reason = "FLOW_PROVEN_AT_STOP"
            elif not snapshot.get("has_solution"):
                stop_reason = "NO_PHASE1_SOLUTION"
            else:
                stop_reason = status or "UNKNOWN"
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
        )
