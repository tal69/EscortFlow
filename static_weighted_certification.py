"""Integer-scaled weighted optimization and independent flow-time certificates.

These are separate solves with separate budgets. The certificate caller must
provide a horizon containing every schedule that could improve the candidate's
flow time. Neither an objective cutoff nor a fixed-flow constraint is added.
"""

import math
import time

from gurobipy import GRB

from static_lexicographic import CERTIFICATE_GAP_MARGIN, _integer_value


WEIGHT_SCALE = 100
FLOW_CERTIFICATE_GAP = 0.999
UNRELIABLE_STATUSES = {"NUMERIC", "SUBOPTIMAL", "INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED"}


def validate_objective_mode(config):
    mode = config.objective_mode
    if mode not in {"legacy", "weighted_integer", "flow_certificate"}:
        raise ValueError(f"Unknown objective mode: {mode}")
    if mode == "legacy":
        if config.certification_target is not None:
            raise ValueError("certification_target requires flow_certificate mode")
        return
    if config.lp or config.lexicographic:
        raise ValueError("Weighted/certification modes require a non-lexicographic integer model")
    if config.beta != 1 or config.gamma != .01 or getattr(config, "alpha", 0) != 0:
        raise ValueError("Weighted/certification modes require F + 0.01 M")
    if getattr(config, "retrieval_mode", "leave") != "leave" or getattr(config, "move_method", "BM") != "BM":
        raise ValueError("Weighted/certification modes support BM leave retrieval")
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


def solve_weighted_or_certificate(model, flow_expr, movement_expr, extract_solution, *,
                                  status_name, solve_start, mode, target=None):
    """Optimize once, extract a numerical certificate, and dispose the model."""
    result = dict(has_solution=False, makespan=None, flowtime=None, movements=None,
                  objective=None, best_bound=None, absolute_gap=None,
                  animation_moves=None, weighted_proven=False, flow_proven=False,
                  counterexample=False, bound_consistent=True, user_cut_time=0.0)
    try:
        model.Params.MIPGap = 0.0
        model.Params.MIPGapAbs = .999
        model.Params.Cutoff = GRB.INFINITY
        model.Params.BestBdStop = GRB.INFINITY
        model.Params.BestObjStop = -GRB.INFINITY
        if mode == "weighted_integer":
            model.setObjective(WEIGHT_SCALE * flow_expr + movement_expr, GRB.MINIMIZE)
        elif mode == "flow_certificate":
            model.setObjective(flow_expr, GRB.MINIMIZE)
            # The postcheck below remains authoritative, not the solver status.
            model.Params.BestBdStop = target - FLOW_CERTIFICATE_GAP + 2 * CERTIFICATE_GAP_MARGIN
            model.Params.BestObjStop = target - 1 + CERTIFICATE_GAP_MARGIN
        else:
            raise ValueError(f"Unknown objective mode: {mode}")
        model.optimize()
        result.update(status_name=status_name(model.Status), runtime=model.Runtime,
                      work=model.Work, cpu_time=time.perf_counter() - solve_start)
        bound = getattr(model, "ObjBound", None)
        valid_bound = bound is not None and math.isfinite(bound)
        reliable = result["status_name"] not in UNRELIABLE_STATUSES
        if model.SolCount:
            result.update(extract_solution(), has_solution=True,
                          flowtime=_integer_value(flow_expr.getValue()),
                          movements=_integer_value(movement_expr.getValue()))
        if mode == "weighted_integer":
            result["best_bound"] = bound / WEIGHT_SCALE if valid_bound else None
            if result["has_solution"]:
                integer_objective = WEIGHT_SCALE * result["flowtime"] + result["movements"]
                result["objective"] = integer_objective / WEIGHT_SCALE
                if valid_bound:
                    gap = abs(integer_objective - bound)
                    result["absolute_gap"] = gap / WEIGHT_SCALE
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
        return result
    finally:
        model.dispose()
