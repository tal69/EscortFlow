"""Two explicit MIP solves for integer flow time, then integer load movements.

Time limits count Gurobi Runtime, excluding model construction and reporting.
The phase-one incumbent is retained if phase two cannot return a solution.
"""

import math
import time

from gurobipy import GRB


ABSOLUTE_GAP_LIMIT = 1.0
LEX_RESULT_HEADER = (
    "makespan,flowtime,#load movements,Flowtime LB,Movement LB,Wall Clock Time,"
    "Work,User Cut Time,Phase 1 Status,Phase 1 Proven,Phase 1 Flowtime,"
    "Phase 1 Absolute Gap,Phase 1 Solve Time,Phase 1 Work,Phase 2 Status,"
    "Phase 2 Proven,Phase 2 Absolute Gap,Phase 2 Solve Time,Phase 2 Work,"
    "Lexicographic Optimal,Solver Status,Phase 2 Skip Reason,"
    "Requested Phase 1 Time Limit,Total Time Limit,"
    "Effective Phase 1 Time Limit,Phase 2 Time Limit"
)


def _integer_value(value):
    rounded = int(round(value))
    if abs(value - rounded) > 1e-4:
        raise ValueError(f"Lexicographic objective is not integral: {value}")
    return rounded


def _phase_result(model, status_name, time_limit):
    objective = _integer_value(model.ObjVal) if model.SolCount else None
    bound = getattr(model, "ObjBound", None)
    gap = None
    if objective is not None and bound is not None and math.isfinite(bound):
        gap = abs(objective - bound)
    status = status_name(model.Status)
    proven = (
        gap is not None and gap < ABSOLUTE_GAP_LIMIT
        and status not in {"NUMERIC", "SUBOPTIMAL", "INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED"}
    )
    return {
        "status_name": status,
        "objective": objective,
        "best_bound": bound,
        "absolute_gap": gap,
        "proven_optimal": proven,
        "runtime": model.Runtime,
        "work": model.Work,
        "time_limit": time_limit,
    }


def solve_lexicographic(
    model, flowtime_expr, movement_expr, extract_solution, *, status_name,
    solve_start, time_limit=None, phase1_time_limit=None, work_limit=None,
    callback=None,
):
    """Minimize flow time, fix its best integer incumbent, then minimize moves.

    Phase two also runs when phase one hits a limit with an incumbent. Only
    proven optima in BOTH phases certify a lexicographic optimum. The total
    time/work budgets cover both optimize calls. The caller transfers ownership
    of the model to this function, which always disposes it.
    """
    result = {
        "has_solution": False, "makespan": None, "flowtime": None,
        "movements": None, "objective": None, "best_bound": None,
        "animation_moves": None, "user_cut_time": 0.0,
        "lexicographic_optimal": False, "phase2_skip_reason": "",
        "time_limit": time_limit, "phase1_time_limit": phase1_time_limit,
        "phase2": {"status_name": "NOT_RUN", "proven_optimal": False,
                   "runtime": 0.0, "work": 0.0, "time_limit": None},
    }

    def snapshot():
        return dict(
            extract_solution(), has_solution=True,
            flowtime=_integer_value(flowtime_expr.getValue()),
            movements=_integer_value(movement_expr.getValue()),
        )

    def optimize():
        if callback is None:
            model.optimize()
        else:
            model.optimize(callback)

    try:
        for name, value in (("time_limit", time_limit),
                            ("phase1_time_limit", phase1_time_limit),
                            ("work_limit", work_limit)):
            if value is not None and (not math.isfinite(value) or value < 0):
                raise ValueError(f"{name} must be finite and nonnegative")
        caps = [cap for cap in (time_limit, phase1_time_limit) if cap is not None]
        first_limit = min(caps) if caps else None
        # Relative gap termination could otherwise accept a gap of one or more.
        model.Params.MIPGap = 0.0
        model.Params.MIPGapAbs = ABSOLUTE_GAP_LIMIT
        model.Params.TimeLimit = first_limit if first_limit is not None else GRB.INFINITY
        model.Params.WorkLimit = work_limit if work_limit is not None else GRB.INFINITY
        model.setObjective(flowtime_expr, GRB.MINIMIZE)
        optimize()
        phase1 = result["phase1"] = _phase_result(model, status_name, first_limit)
        result["status_name"] = phase1["status_name"]
        result["best_bound"] = phase1["best_bound"]
        result["objective"] = phase1["objective"]
        if model.SolCount:
            result.update(snapshot())

        remaining_time = None if time_limit is None else max(0.0, time_limit - phase1["runtime"])
        remaining_work = None if work_limit is None else max(0.0, work_limit - phase1["work"])
        result["phase2"]["time_limit"] = remaining_time
        if not result["has_solution"]:
            result["phase2_skip_reason"] = "NO_PHASE1_INCUMBENT"
        elif phase1["status_name"] in {"INTERRUPTED", "NUMERIC", "SUBOPTIMAL"}:
            result["phase2_skip_reason"] = "PHASE1_" + phase1["status_name"]
        elif remaining_time is not None and remaining_time <= 0:
            result["phase2_skip_reason"] = "TOTAL_TIME_LIMIT"
        elif remaining_work is not None and remaining_work <= 0:
            result["phase2_skip_reason"] = "TOTAL_WORK_LIMIT"
        else:
            # Save a full feasible start before changing the model. Keep the
            # extracted plan above even if phase two stops before accepting it.
            variables = model.getVars()
            start_values = model.getAttr("X", variables)
            model.addConstr(flowtime_expr == result["flowtime"], name="lex_fixed_flowtime")
            model.setObjective(movement_expr, GRB.MINIMIZE)
            model.NumStart = 1
            model.setAttr("Start", variables, start_values)
            model.Params.TimeLimit = remaining_time if remaining_time is not None else GRB.INFINITY
            model.Params.WorkLimit = remaining_work if remaining_work is not None else GRB.INFINITY
            optimize()
            phase2 = result["phase2"] = _phase_result(model, status_name, remaining_time)
            if model.SolCount:
                candidate = snapshot()
                if candidate["flowtime"] != result["flowtime"]:
                    raise RuntimeError("Phase two changed the fixed flow time")
                if candidate["movements"] <= result["movements"]:
                    result.update(candidate)
            result["objective"] = result["movements"]
            result["best_bound"] = phase2["best_bound"]
            result["lexicographic_optimal"] = phase1["proven_optimal"] and phase2["proven_optimal"]
            if result["lexicographic_optimal"]:
                result["status_name"] = "OPTIMAL"
            elif phase2["proven_optimal"]:
                result["status_name"] = "FLOWTIME_NOT_PROVEN"
            else:
                result["status_name"] = phase2["status_name"]
        if result["phase2_skip_reason"] and result["status_name"] == "OPTIMAL":
            result["status_name"] = "PHASE2_NOT_RUN"
        result["work"] = phase1["work"] + result["phase2"]["work"]
        result["cpu_time"] = time.perf_counter() - solve_start
        return result
    finally:
        model.dispose()


def build_lex_csv_suffix(result):
    """Report each phase's bounds and certificate in its own objective units."""
    p1, p2 = result.get("phase1", {}), result.get("phase2", {})
    values = [
        result.get("makespan"), result.get("flowtime"), result.get("movements"),
        p1.get("best_bound"), p2.get("best_bound"), result.get("cpu_time"),
        result.get("work"), result.get("user_cut_time", 0.0),
        p1.get("status_name", "NOT_RUN"), p1.get("proven_optimal", False),
        p1.get("objective"), p1.get("absolute_gap"), p1.get("runtime"), p1.get("work"),
        p2.get("status_name", "NOT_RUN"), p2.get("proven_optimal", False),
        p2.get("absolute_gap"), p2.get("runtime"), p2.get("work"),
        result.get("lexicographic_optimal", False), result.get("status_name", "ERROR"),
        result.get("phase2_skip_reason", ""), result.get("phase1_time_limit"),
        result.get("time_limit"), p1.get("time_limit"), p2.get("time_limit"),
    ]

    def format_value(value):
        if value is None or isinstance(value, float) and not math.isfinite(value):
            return "-"
        if isinstance(value, bool):
            return str(int(value))
        if isinstance(value, float):
            return f"{value:.10g}"
        return str(value)

    return "," + ",".join(format_value(value) for value in values)
