"""Integer-scaled weighted optimization and independent flow-time certificates.

These are separate solves with separate budgets. The certificate caller must
provide a horizon containing every schedule that could improve the candidate's
flow time. Neither an objective cutoff nor a fixed-flow constraint is added.
"""

import math
from numbers import Integral
import time

from gurobipy import GRB

from static_lexicographic import CERTIFICATE_GAP_MARGIN, _integer_value


WEIGHT_SCALE = 100
FLOW_CERTIFICATE_GAP = 0.999
WEIGHTED_BOUND_MARGIN = 0.001
UNRELIABLE_STATUSES = {"NUMERIC", "SUBOPTIMAL", "INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED"}
RELIABLE_FINISHED_STATUSES = {
    "OPTIMAL", "ITERATION_LIMIT", "NODE_LIMIT", "TIME_LIMIT", "SOLUTION_LIMIT",
    "INTERRUPTED", "USER_OBJ_LIMIT", "WORK_LIMIT", "MEM_LIMIT",
}


def _weight_scale(value):
    if isinstance(value, bool) or not isinstance(value, Integral) or value <= 0:
        raise ValueError("weight_scale must be a positive integer")
    return int(value)


def validate_objective_mode(config):
    mode = config.objective_mode
    scale = _weight_scale(getattr(config, "weight_scale", WEIGHT_SCALE))
    extension = getattr(config, "flow_proof_extension_time_limit", None)
    if extension is not None:
        if mode != "weighted_integer":
            raise ValueError("A flow-proof extension requires weighted_integer mode")
        if not math.isfinite(extension) or extension < 0:
            raise ValueError("The flow-proof extension time limit must be finite and nonnegative")
        if config.time_limit is None or not math.isfinite(config.time_limit) or config.time_limit <= 0:
            raise ValueError("A flow-proof extension requires a positive finite weighted time limit")
        if not isinstance(getattr(config, "stop_on_flow_proof", True), bool):
            raise ValueError("stop_on_flow_proof must be a boolean")
    if mode not in {"legacy", "weighted_integer", "flow_certificate"}:
        raise ValueError(f"Unknown objective mode: {mode}")
    if mode == "legacy":
        if config.certification_target is not None:
            raise ValueError("certification_target requires flow_certificate mode")
        return
    if config.lp or config.lexicographic:
        raise ValueError("Weighted/certification modes require a non-lexicographic integer model")
    if config.beta != 1 or config.gamma != 1 / scale or getattr(config, "alpha", 0) != 0:
        raise ValueError("Weighted/certification modes require F + M / weight_scale")
    if (getattr(config, "retrieval_mode", "leave") not in {"leave", "continue"}
            or getattr(config, "move_method", "BM") != "BM"):
        raise ValueError("Weighted/certification modes support BM leave or continue retrieval")
    target = config.certification_target
    if mode == "flow_certificate":
        if target is None or not math.isfinite(target) or target < 0 or int(target) != target:
            raise ValueError("Flow certification requires a nonnegative integer target")
    elif target is not None:
        raise ValueError("certification_target requires flow_certificate mode")


def certification_eligibility(weighted, gap_threshold=0.1):
    """Threshold is an absolute gap in original F + 0.01 M units."""
    if not math.isfinite(gap_threshold) or not 0 < gap_threshold < 1:
        raise ValueError("Certification gap threshold must be strictly between zero and one")
    if not weighted.get("has_solution"):
        return False, "NO_WEIGHTED_SOLUTION"
    if weighted.get("status_name") in UNRELIABLE_STATUSES:
        return False, "UNRELIABLE_WEIGHTED_STATUS"
    if weighted.get("bound_consistent") is False:
        return False, "INCONSISTENT_WEIGHTED_BOUND"
    if weighted.get("weighted_proven"):
        return True, ""
    gap = weighted.get("absolute_gap")
    if gap is not None and math.isfinite(gap) and gap < gap_threshold:
        return True, ""
    return False, "WEIGHTED_GAP_TOO_LARGE"


def sufficient_flow_horizon(targets, outputs, candidate_flow):
    """Conservative common EF/LF T covering all schedules with F <= candidate.

    Each target i needs at least d_i movements to reach an output, so its
    arrival time cannot exceed candidate - sum(d_j) + d_i. Setting T to this
    physical arrival horizon is sufficient in both formulations (EF also
    permits the tighter T = max(0, H - 1)).
    """
    distances = [min(abs(x - ox) + abs(y - oy) for ox, oy in outputs) for x, y in targets]
    return max(0, int(candidate_flow) - sum(distances) + max(distances, default=0))


def _nonnegative_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float, Integral)):
        raise ValueError(f"{name} must be a nonnegative integer")
    if not math.isfinite(value) or value < 0 or int(value) != value:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _safe_weight_instance(targets, outputs, cell_count, escort_count):
    cell_count = _nonnegative_integer(cell_count, "cell_count")
    escort_count = _nonnegative_integer(escort_count, "escort_count")
    targets, outputs = tuple(targets), tuple(outputs)
    if cell_count == 0 or escort_count > cell_count:
        raise ValueError("A positive cell count and no more escorts than cells are required")
    load_count = cell_count - escort_count
    if len(targets) > load_count:
        raise ValueError("Target count exceeds the initial number of loads")
    if targets and not outputs:
        raise ValueError("At least one output is required for retrieval")
    for point in targets + outputs:
        if len(point) != 2 or any(isinstance(v, bool) or not isinstance(v, Integral) for v in point):
            raise ValueError("Targets and outputs must have integer cell coordinates")
    distances = [min(abs(x - ox) + abs(y - oy) for ox, oy in outputs) for x, y in targets]
    return load_count, sum(distances), max(distances, default=0)


def safe_weight_parameters(targets, outputs, cell_count, escort_count, feasible_flow):
    """Choose R so an optimum of R F + M is globally lexicographic.

    A feasible flow bound gives H containing all minimum-flow schedules.
    Some minimum-flow schedule has no movements after its last retrieval,
    hence M <= (N - e) H. Every feasible plan needs at least D = sum(d_i)
    movements, so R = (N - e) H - D + 1 makes any worse flow inferior
    to that schedule. The weighted model must also contain this horizon.
    The movement bound does not apply to arbitrary redundant schedules.
    """
    load_count, distance_sum, distance_max = _safe_weight_instance(
        targets, outputs, cell_count, escort_count)
    feasible_flow = _nonnegative_integer(feasible_flow, "feasible_flow")
    if feasible_flow < distance_sum:
        raise ValueError("Feasible flow cannot be below the distance lower bound")
    horizon = max(0, feasible_flow - distance_sum + distance_max)
    movement_bound = load_count * horizon
    return dict(flow_weight=movement_bound - distance_sum + 1,
                movement_bound=movement_bound, movement_lower_bound=distance_sum,
                flow_horizon=horizon)


def gap_flow_certificate(weighted, targets, outputs, cell_count, escort_count,
                         physical_horizon, flow_weight):
    """Certify an incumbent's global flow from the weighted lower bound.

    If better flow existed, a minimum-flow plan with F <= F_inc - 1 would
    finish by H_minus and need at most (N - e) H_minus movements. Its
    weighted objective would be <= R (F_inc - 1) + U_minus, so a strictly
    larger valid lower bound rules it out. A 0.001 margin protects the
    integer threshold. An insufficient weighted horizon cannot prove this.
    """
    flow_weight = _weight_scale(flow_weight)
    physical_horizon = _nonnegative_integer(physical_horizon, "physical_horizon")
    load_count, distance_sum, distance_max = _safe_weight_instance(
        targets, outputs, cell_count, escort_count)
    result = dict(flow_proven=False, proof_source="", reason="",
                  flow_lower_bound=distance_sum, better_flow_horizon=None,
                  better_flow_movement_bound=None, weighted_bound_threshold=None,
                  weighted_gap_threshold=None, scaled_absolute_gap=None,
                  horizon_sufficient=False)

    def reject(reason):
        result["reason"] = reason
        return result

    if not weighted.get("has_solution"):
        return reject("NO_WEIGHTED_SOLUTION")
    if weighted.get("status_name") not in RELIABLE_FINISHED_STATUSES:
        return reject("UNRELIABLE_WEIGHTED_STATUS")
    try:
        flow = _nonnegative_integer(weighted.get("flowtime"), "flowtime")
        movements = _nonnegative_integer(weighted.get("movements"), "movements")
    except ValueError:
        return reject("INVALID_WEIGHTED_INCUMBENT")
    if flow < distance_sum:
        return reject("FLOW_BELOW_DISTANCE_BOUND")
    if movements < distance_sum:
        return reject("MOVEMENTS_BELOW_DISTANCE_BOUND")
    if weighted.get("weight_scale", flow_weight) != flow_weight:
        return reject("INCONSISTENT_WEIGHT_SCALE")
    integer_objective = flow_weight * flow + movements
    reported_objective = weighted.get("scaled_objective")
    if reported_objective is not None and (
            not math.isfinite(reported_objective)
            or abs(reported_objective - integer_objective) > CERTIFICATE_GAP_MARGIN):
        return reject("INCONSISTENT_WEIGHTED_OBJECTIVE")
    if flow == distance_sum:
        result.update(flow_proven=True, proof_source="distance_bound",
                      flow_lower_bound=flow, horizon_sufficient=True)
        return result

    better_flow = flow - 1
    horizon = max(0, better_flow - distance_sum + distance_max)
    movement_bound = load_count * horizon
    result.update(
        better_flow_horizon=horizon, better_flow_movement_bound=movement_bound,
        weighted_bound_threshold=flow_weight * better_flow + movement_bound + WEIGHTED_BOUND_MARGIN,
        weighted_gap_threshold=flow_weight + movements - movement_bound - WEIGHTED_BOUND_MARGIN,
        horizon_sufficient=physical_horizon >= horizon,
    )
    bound = weighted.get("scaled_best_bound")
    if "scaled_best_bound" not in weighted and weighted.get("best_bound") is not None:
        bound = weighted["best_bound"] * flow_weight
    if bound is None or not math.isfinite(bound):
        return reject("NO_FINITE_WEIGHTED_BOUND")
    result["scaled_absolute_gap"] = abs(integer_objective - bound)
    if weighted.get("bound_consistent") is False or bound > integer_objective + CERTIFICATE_GAP_MARGIN:
        return reject("INCONSISTENT_WEIGHTED_BOUND")
    if not result["horizon_sufficient"]:
        return reject("INSUFFICIENT_WEIGHTED_HORIZON")
    if bound > result["weighted_bound_threshold"] + CERTIFICATE_GAP_MARGIN:
        result.update(flow_proven=True, proof_source="weighted_gap", flow_lower_bound=flow)
        return result
    return reject("WEIGHTED_GAP_TOO_LARGE")


def solve_weighted_or_certificate(model, flow_expr, movement_expr, extract_solution, *,
                                  status_name, solve_start, mode, target=None,
                                  weight_scale=WEIGHT_SCALE, flow_proof_context=None,
                                  extract_callback_metrics=None):
    """Optimize once, extract a numerical certificate, and dispose the model."""
    result = dict(has_solution=False, makespan=None, flowtime=None, movements=None,
                  objective=None, best_bound=None, absolute_gap=None,
                  animation_moves=None, weighted_proven=False, flow_proven=False,
                  counterexample=False, bound_consistent=True, user_cut_time=0.0,
                  scaled_objective=None, scaled_best_bound=None, scaled_absolute_gap=None)
    try:
        weight_scale = _weight_scale(weight_scale)
        result["weight_scale"] = weight_scale
        model.Params.MIPGap = 0.0
        model.Params.MIPGapAbs = .999
        model.Params.Cutoff = GRB.INFINITY
        model.Params.BestBdStop = GRB.INFINITY
        model.Params.BestObjStop = -GRB.INFINITY
        if mode == "weighted_integer":
            model.setObjective(weight_scale * flow_expr + movement_expr, GRB.MINIMIZE)
        elif mode == "flow_certificate":
            model.setObjective(flow_expr, GRB.MINIMIZE)
            # The postcheck below remains authoritative, not the solver status.
            model.Params.BestBdStop = target - FLOW_CERTIFICATE_GAP + 2 * CERTIFICATE_GAP_MARGIN
            model.Params.BestObjStop = target - 1 + CERTIFICATE_GAP_MARGIN
        else:
            raise ValueError(f"Unknown objective mode: {mode}")
        continuation = None
        model_build_time = time.perf_counter() - solve_start
        if flow_proof_context is not None:
            if mode != "weighted_integer":
                raise ValueError("A flow-proof extension requires weighted_integer mode")
            from static_safe_weighted_search import SafeWeightedContinuation
            continuation = SafeWeightedContinuation(
                flow_expr, movement_expr, weight_scale, flow_proof_context,
                extract_callback_metrics=extract_callback_metrics)
            model.Params.TimeLimit = (flow_proof_context["weighted_time_limit"]
                                      + flow_proof_context["extension_time_limit"])
            model.optimize(continuation)
            continuation.raise_if_failed()
        else:
            model.optimize()
        result.update(status_name=status_name(model.Status), runtime=model.Runtime,
                      work=model.Work, cpu_time=time.perf_counter() - solve_start,
                      node_count=getattr(model, "NodeCount", None))
        bound = getattr(model, "ObjBound", None)
        valid_bound = bound is not None and math.isfinite(bound)
        reliable = result["status_name"] in RELIABLE_FINISHED_STATUSES
        if model.SolCount:
            result.update(extract_solution(), has_solution=True,
                          flowtime=_integer_value(flow_expr.getValue()),
                          movements=_integer_value(movement_expr.getValue()))
        if mode == "weighted_integer":
            result["scaled_best_bound"] = bound if valid_bound else None
            result["best_bound"] = bound / weight_scale if valid_bound else None
            if result["has_solution"]:
                integer_objective = weight_scale * result["flowtime"] + result["movements"]
                result["scaled_objective"] = integer_objective
                result["objective"] = integer_objective / weight_scale
                if valid_bound:
                    gap = abs(integer_objective - bound)
                    result["scaled_absolute_gap"] = gap
                    result["absolute_gap"] = gap / weight_scale
                    result["bound_consistent"] = bound <= integer_objective + CERTIFICATE_GAP_MARGIN
                    result["weighted_proven"] = (reliable and result["bound_consistent"]
                                                  and gap < 1 - CERTIFICATE_GAP_MARGIN)
        else:
            result["best_bound"] = bound if valid_bound else None
            result["objective"] = result["flowtime"]
            result["counterexample"] = reliable and result["has_solution"] and result["flowtime"] < target
            result["flow_proven"] = bool(
                reliable and valid_bound and not result["counterexample"]
                and target - bound < FLOW_CERTIFICATE_GAP - CERTIFICATE_GAP_MARGIN
                and bound <= target + CERTIFICATE_GAP_MARGIN
            )
            if result["has_solution"] and valid_bound:
                result["absolute_gap"] = abs(result["flowtime"] - bound)
        if continuation is not None:
            result.update(continuation.finalize(result))
            proof_runtime = result.get("first_flow_proof_runtime")
            result["first_flow_proof_cpu_time"] = (
                None if proof_runtime is None else model_build_time + proof_runtime)
            snapshot = result["phase1_snapshot"]
            snapshot["cpu_time_is_estimate"] = snapshot["snapshot_source"] != "SOLVE_FINISHED"
            snapshot["cpu_time"] = (model_build_time + snapshot["runtime"]
                                    if snapshot["cpu_time_is_estimate"] else result["cpu_time"])
        return result
    finally:
        model.dispose()
