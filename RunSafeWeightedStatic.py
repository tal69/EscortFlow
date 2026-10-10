#!/usr/bin/env python3
"""Safe integer PBS objective with flow-proof timing in one weighted search.

The first budget's incumbent remains the primary experimental result. If its
flow time is unproved, the same weighted search continues for another budget.
The final search result is reported separately, including any counterexample
with better flow time. Both formulations receive the common greedy start.
"""

import argparse
import csv
from datetime import datetime
import itertools
import math
from pathlib import Path
import platform
import sys
import time

from PBSCom import GeneretaeRandomInstance
import OneStepHeuristic_v2
from static_integrated_lp import LP_FIELDS, attach_lp, embedded_value
from RunWeightedStatic import parse_range, positive_number
from static_weighted_certification import (
    CERTIFICATE_GAP_MARGIN, RELIABLE_FINISHED_STATUSES, gap_flow_certificate,
    safe_weight_parameters,
)

PROTOCOL = "safe_integer_flow_timing_v5"
FLOW_PROOF_CHECK_MODE = "bound_or_flow_change"
SOLUTION_FIELDS = (
    "has_solution", "makespan", "flowtime", "movements", "objective", "best_bound",
    "absolute_gap", "scaled_objective", "scaled_best_bound", "scaled_absolute_gap", "weighted_proven",
)
FIELDNAMES = [
    "host", "date", "protocol", "formulation", "Lx x Ly", "#IOs", "# Escorts", "#Loads",
    "IOs", "Escorts", "Target Loads", "seed", "retrieval_mode", "movement_mode",
    "flow_weight", "movement_integer_weight", "movement_weight", "objective_units",
    "scaled_objective_units", "safe_movement_bound", "safe_movement_lower_bound",
    "safe_flow_horizon", "warmstart", "threads",
    "weighted_horizon", "weighted_physical_horizon", "weighted_global_scope", "greedy_makespan",
    "greedy_flowtime", "greedy_movements", "greedy_objective", "greedy_scaled_objective", "greedy_time",
    "weighted_time_limit", "extension_time_limit", "total_time_limit", "Solver Status",
    "flow_proof_check_mode",
    *SOLUTION_FIELDS, "weighted_bound_consistent", "weighted_runtime", "weighted_cpu_time",
    "weighted_cpu_time_is_estimate", "weighted_work", "phase1_flow_proven", "phase1_flow_proof_source",
    "phase1_incumbent_runtime", "phase1_snapshot_source", "phase1_snapshot_missing",
    "phase1_bound_checkpoint_runtime", "phase1_statistics_checkpoint_runtime",
    "phase_transition_runtime", "phase_transition_delay", "phase1_node_count",
    "extension_used", "extension_runtime", "search_stop_reason", "optimization_calls",
    "final_status", *(f"final_{key}" for key in SOLUTION_FIELDS),
    "final_bound_consistent", "final_runtime", "final_cpu_time", "final_work",
    "final_flow_proven", "final_flow_proof_source", "final_lexicographic_proven",
    "gap_flow_proven", "gap_certificate_reason", "gap_certificate_horizon_sufficient",
    "better_flow_horizon", "better_flow_movement_bound", "weighted_bound_threshold",
    "weighted_gap_threshold", "proof_scaled_best_bound", "proof_scaled_absolute_gap",
    "flow_lower_bound", "flow_proven", "flow_proof_source", "counterexample",
    "counterexample_flowtime", "counterexample_movements", "counterexample_scaled_objective",
    "lexicographic_proven", "first_flow_proof_runtime", "first_flow_proof_cpu_time",
    "first_flow_proof_node_count", "first_flow_proof_flowtime", "first_flow_proof_source",
    "first_flow_proof_scaled_bound", "first_flow_proof_work", "first_flow_proof_method",
    "flow_proof_invalidated", "flow_proof_invalidation_reason", "total_wall_time", "error",
    *LP_FIELDS,
]


def make_solver(args, *, flow_weight):
    retrieval_mode = getattr(args, "retrieval_mode", "leave")
    common = dict(
        Lx=args.Lx, Ly=args.Ly, output_cells=tuple(args.outputs),
        beta=1, gamma=1 / flow_weight, weight_scale=flow_weight,
        time_limit=args.weighted_time_limit, threads=args.threads, mip_focus=0,
        objective_mode="weighted_integer", flow_proof_extension_time_limit=args.extension_time_limit,
        stop_on_flow_proof=True,
    )
    if args.formulation == "escortflow":
        from escort_flow_static_gurobi import StaticGurobiConfig, StaticEscortFlowGurobiSolver
        return StaticEscortFlowGurobiSolver(StaticGurobiConfig(retrieval_mode=retrieval_mode, **common))
    from load_flow_static_gurobi import LoadFlowStaticGurobiConfig, LoadFlowStaticGurobiSolver
    return LoadFlowStaticGurobiSolver(LoadFlowStaticGurobiConfig(
        move_method="BM", alpha=0, retrieval_mode=retrieval_mode, **common))


def _check_weighted_proof(weighted, coefficient, greedy_flow, greedy_objective,
                          movement_bound, movement_lower_bound, scope):
    """Validate native bound data before interpreting a proof globally."""
    def nonnegative_integer(value):
        return (not isinstance(value, bool) and isinstance(value, (int, float))
                and math.isfinite(value) and value >= 0 and int(value) == value)

    if (any(not nonnegative_integer(value) for value in
            (coefficient, greedy_flow, greedy_objective, movement_bound, movement_lower_bound))
            or coefficient == 0 or movement_bound < movement_lower_bound
            or greedy_flow < movement_lower_bound
            or greedy_objective < coefficient * greedy_flow + movement_lower_bound):
        raise ValueError("Invalid safe-weight proof parameters")
    if (not weighted.get("has_solution") or any(
            not nonnegative_integer(weighted.get(key))
            or weighted[key] < movement_lower_bound for key in ("flowtime", "movements"))):
        raise ValueError("Invalid weighted optimality certificate")
    objective = coefficient * weighted["flowtime"] + weighted["movements"]
    bound = weighted.get("scaled_best_bound")
    # BUDGET_REACHED labels a frozen, valid callback checkpoint rather than
    # the status of a second optimization call.
    reliable = RELIABLE_FINISHED_STATUSES | {"BUDGET_REACHED"}
    if (weighted["status_name"] not in reliable or weighted.get("bound_consistent") is False
            or weighted.get("weight_scale") != coefficient
            or weighted.get("scaled_objective") != objective
            or bound is None or not math.isfinite(bound)
            or bound > objective + CERTIFICATE_GAP_MARGIN
            or abs(objective - bound) >= 1 - CERTIFICATE_GAP_MARGIN):
        raise ValueError("Invalid weighted optimality certificate")
    if (not scope or coefficient <= movement_bound - movement_lower_bound
            or weighted["flowtime"] > greedy_flow
            or objective > greedy_objective):
        raise ValueError("Proven weighted solution contradicts the safe-weight scope or known greedy solution")


def run_instance(args, seed, escort_count, solver_factory=make_solver):
    start = time.perf_counter()
    retrieval_mode = getattr(args, "retrieval_mode", "leave")
    locations = sorted(itertools.product(range(args.Lx), range(args.Ly)))
    targets, escorts = GeneretaeRandomInstance(seed, locations, escort_count, args.loads)
    row = dict.fromkeys(FIELDNAMES, "")
    row.update({
        "host": platform.node(), "date": datetime.now().astimezone().isoformat(), "protocol": PROTOCOL,
        "formulation": args.formulation, "Lx x Ly": f"{args.Lx}x{args.Ly}",
        "#IOs": len(args.outputs), "# Escorts": escort_count, "#Loads": args.loads,
        "IOs": repr(args.outputs), "Escorts": repr(escorts), "Target Loads": repr(targets),
        "seed": seed, "retrieval_mode": retrieval_mode, "movement_mode": "BM", "warmstart": 1,
        "threads": args.threads, "movement_integer_weight": 1,
        "objective_units": "F+M/R", "scaled_objective_units": "R*F+M",
        "weighted_time_limit": args.weighted_time_limit, "extension_time_limit": args.extension_time_limit,
        "total_time_limit": args.weighted_time_limit + args.extension_time_limit,
        "flow_proof_check_mode": FLOW_PROOF_CHECK_MODE,
        "extension_used": 0, "extension_runtime": 0, "optimization_calls": 0,
        "gap_flow_proven": 0, "flow_proven": 0, "counterexample": 0, "lexicographic_proven": 0,
        "final_flow_proven": 0, "final_lexicographic_proven": 0,
        "lp_requested": 0, "lp_status": "NOT_RUN",
    })
    print(f"[{row['date']}] seed={seed} escorts={escort_count}: greedy and safe weighted search starting", flush=True)
    try:
        greedy_start = time.perf_counter()
        trace = OneStepHeuristic_v2.SolveGreedy(
            args.Lx, args.Ly, set(args.outputs), set(targets), set(escorts),
            verbal=False, max_steps=max(1, 4 * args.loads * (args.Lx + args.Ly - 2) + 1),
            retrieval_mode=retrieval_mode, return_trace=True,
        )
        makespan, flow, movements, _, escort_history, target_history = trace
        parameters = safe_weight_parameters(targets, args.outputs, cell_count=args.Lx * args.Ly,
                                            escort_count=escort_count, feasible_flow=flow)
        coefficient = parameters["flow_weight"]
        movement_lower_bound = parameters["movement_lower_bound"]
        if (isinstance(movements, bool) or not isinstance(movements, (int, float))
                or not math.isfinite(movements)
                or int(movements) != movements or movements < movement_lower_bound):
            raise ValueError("Greedy movements violate the distance lower bound")
        safe_horizon = parameters["flow_horizon"]
        # Compare the formulations over the same latest retrieval time. EF
        # arrivals at movement index T occur at T+1; LF retrievals use index T.
        # Preserve the complete greedy trace, including its output-service tail.
        physical_horizon = max(safe_horizon, makespan + 1)
        horizon = physical_horizon - 1 if args.formulation == "escortflow" else physical_horizon
        scope = physical_horizon >= safe_horizon
        row.update(
            flow_weight=coefficient, movement_weight=1 / coefficient,
            safe_movement_bound=parameters["movement_bound"], safe_flow_horizon=safe_horizon,
            safe_movement_lower_bound=movement_lower_bound,
            greedy_makespan=makespan, greedy_flowtime=flow, greedy_movements=movements,
            greedy_objective=flow + movements / coefficient, greedy_scaled_objective=coefficient * flow + movements,
            greedy_time=time.perf_counter() - greedy_start, weighted_horizon=horizon,
            weighted_physical_horizon=physical_horizon, weighted_global_scope=int(scope),
        )
        print(f"seed={seed} escorts={escort_count}: R={coefficient}, T={horizon}, "
              f"weighted budget={args.weighted_time_limit:g}s, conditional same-tree extension="
              f"{args.extension_time_limit:g}s, warm start=greedy", flush=True)
        solver = solver_factory(args, flow_weight=coefficient)
        try:
            warmstart = solver.build_warmstart_from_trace(targets, escorts, horizon, target_history, escort_history)
            search = solver.solve(targets, escorts, horizon, warmstart=warmstart)
        finally:
            solver.close()
        final = search
        snapshot = search["phase1_snapshot"]
        weighted = dict(snapshot)
        row["Solver Status"] = weighted["status_name"]
        row["final_status"] = final["status_name"]
        for source, prefix in ((weighted, ""), (final, "final_")):
            for key in SOLUTION_FIELDS:
                value = source.get(key)
                row[f"{prefix}{key}"] = int(value) if isinstance(value, bool) else value
        row.update(weighted_bound_consistent=int(bool(weighted.get("bound_consistent", True))),
                   final_bound_consistent=int(bool(final.get("bound_consistent", True))),
                   weighted_runtime=weighted.get("runtime"), weighted_work=weighted.get("work"),
                   weighted_cpu_time=weighted.get("cpu_time"),
                   weighted_cpu_time_is_estimate=int(bool(weighted.get("cpu_time_is_estimate"))),
                   phase1_flow_proven=int(bool(weighted.get("flow_proven"))),
                   phase1_flow_proof_source=weighted.get("flow_proof_source", ""),
                   phase1_incumbent_runtime=weighted.get("incumbent_runtime"),
                   phase1_snapshot_source=weighted.get("snapshot_source", ""),
                   phase1_snapshot_missing=int(bool(search.get("phase1_snapshot_missing"))),
                   phase1_bound_checkpoint_runtime=snapshot.get("bound_checkpoint_runtime"),
                   phase1_statistics_checkpoint_runtime=snapshot.get("statistics_checkpoint_runtime"),
                   phase_transition_runtime=search.get("phase_transition_runtime"),
                   phase_transition_delay=search.get("phase_transition_delay"),
                   phase1_node_count=snapshot.get("node_count"),
                   extension_used=int(bool(search.get("extension_used"))),
                   extension_runtime=search.get("extension_runtime", 0),
                   search_stop_reason=search.get("stop_reason", ""),
                   optimization_calls=search.get("optimization_calls", 0))
        for key in ("runtime", "cpu_time", "work"):
            row[f"final_{key}"] = final.get(key)
        for key in ("runtime", "cpu_time", "node_count", "flowtime", "source", "scaled_bound", "work", "method"):
            row[f"first_flow_proof_{key}"] = final.get(f"first_flow_proof_{key}")
        row["flow_proof_invalidated"] = int(bool(final.get("flow_proof_invalidated")))
        row["flow_proof_invalidation_reason"] = final.get("flow_proof_invalidation_reason", "")
        if row["optimization_calls"] != 1:
            raise ValueError("Safe weighted continuation must use exactly one optimization call")

        # Use the final bound to assess the frozen incumbent. This changes
        # neither its saved first-budget objective nor its saved first-budget gap.
        proof_input = dict(weighted, status_name=final["status_name"],
                           scaled_best_bound=final.get("scaled_best_bound"),
                           best_bound=final.get("best_bound"),
                           bound_consistent=final.get("bound_consistent", True))
        certificate = gap_flow_certificate(proof_input, targets, args.outputs, args.Lx * args.Ly,
                                           escort_count, physical_horizon, coefficient)
        proof_bound = final.get("scaled_best_bound")
        proof_gap = None
        if weighted.get("has_solution") and proof_bound is not None and math.isfinite(proof_bound):
            proof_gap = abs(coefficient * weighted["flowtime"] + weighted["movements"] - proof_bound)
        witness = search.get("flow_counterexample")
        if (witness is None and final["status_name"] in RELIABLE_FINISHED_STATUSES
                and weighted.get("has_solution") and final.get("has_solution")
                and final["flowtime"] < weighted["flowtime"]):
            witness = final
        counterexample = bool(witness is not None and weighted.get("has_solution")
                              and witness.get("has_solution")
                              and witness["flowtime"] < weighted["flowtime"])
        if counterexample:
            for key in ("flowtime", "movements", "scaled_objective"):
                row[f"counterexample_{key}"] = witness.get(key)
        row.update(counterexample=int(counterexample),
                   gap_flow_proven=int(bool(certificate["flow_proven"]) and not counterexample),
                   gap_certificate_reason=certificate["reason"],
                   gap_certificate_horizon_sufficient=int(certificate["horizon_sufficient"]),
                   proof_scaled_best_bound=final.get("scaled_best_bound"),
                   proof_scaled_absolute_gap=proof_gap)
        for key in ("better_flow_horizon", "better_flow_movement_bound", "weighted_bound_threshold",
                    "weighted_gap_threshold", "flow_lower_bound"):
            row[key] = certificate.get(key)
        if weighted.get("weighted_proven"):
            _check_weighted_proof(weighted, coefficient, flow, row["greedy_scaled_objective"],
                                  parameters["movement_bound"], movement_lower_bound, scope)
            if counterexample:
                raise ValueError("Final solution contradicts the frozen weighted optimality certificate")
            row.update(flow_proven=1, lexicographic_proven=1, flow_proof_source="safe_weighted_optimum",
                       flow_lower_bound=weighted["flowtime"])
        elif certificate["flow_proven"] and not counterexample:
            row.update(flow_proven=1, flow_proof_source=certificate["proof_source"])
        # A later bound can also prove the frozen movement count optimal. Its
        # first-budget weighted_proven flag and gap still remain unchanged.
        if (not counterexample and row["flow_proven"]
                and proof_gap is not None and proof_gap < 1 - CERTIFICATE_GAP_MARGIN):
            _check_weighted_proof(proof_input, coefficient, flow, row["greedy_scaled_objective"],
                                  parameters["movement_bound"], movement_lower_bound, scope)
            row["lexicographic_proven"] = 1
        final_certificate = gap_flow_certificate(final, targets, args.outputs, args.Lx * args.Ly,
                                                 escort_count, physical_horizon, coefficient)
        row.update(final_flow_proven=int(bool(final_certificate["flow_proven"])),
                   final_flow_proof_source=final_certificate["proof_source"])
        if final.get("weighted_proven"):
            known_flow = row.get("first_flow_proof_flowtime")
            if ((witness is not None and witness.get("flowtime", final["flowtime"]) < final["flowtime"])
                    or (known_flow is not None and known_flow < final["flowtime"])):
                raise ValueError("Final weighted optimality contradicts a known feasible smaller flow")
            _check_weighted_proof(final, coefficient, flow, row["greedy_scaled_objective"],
                                  parameters["movement_bound"], movement_lower_bound, scope)
            row.update(final_flow_proven=1, final_lexicographic_proven=1,
                       final_flow_proof_source="safe_weighted_optimum")
        if (row["flow_proof_invalidated"]
                and row["flow_proof_invalidation_reason"] != "UNRELIABLE_FINAL_STATUS"):
            raise ValueError(f"Flow proof invalidated: {row['flow_proof_invalidation_reason']}")
        print(f"seed={seed} escorts={escort_count}: first-budget F={weighted.get('flowtime')} "
              f"M={weighted.get('movements')} integer gap={weighted.get('scaled_absolute_gap')}; "
              f"final F={final.get('flowtime')} M={final.get('movements')} "
              f"flow first proved at {row['first_flow_proof_runtime']} solver seconds; "
              f"total solve={row['final_runtime']}s, stop={row['search_stop_reason']}", flush=True)
    except Exception as exc:
        row["error"] = f"{type(exc).__name__}: {exc}"
        for key in ("weighted_proven", "final_weighted_proven", "phase1_flow_proven",
                    "flow_proven", "gap_flow_proven", "lexicographic_proven",
                    "final_flow_proven", "final_lexicographic_proven"):
            row[key] = 0
        for key in ("runtime", "cpu_time", "node_count", "flowtime", "source", "scaled_bound", "work", "method"):
            row[f"first_flow_proof_{key}"] = None
        if not row["Solver Status"]:
            row["Solver Status"] = "ERROR"
        print(f"seed={seed} escorts={escort_count}: ERROR: {row['error']}", flush=True)
    row["total_wall_time"] = time.perf_counter() - start
    print(f"seed={seed} escorts={escort_count}: complete, original flow_proven={row['flow_proven']} "
          f"({row['flow_proof_source']}), counterexample={row['counterexample']}, "
          f"wall={row['total_wall_time']:.3f}s", flush=True)
    return row


def _validate_merge_rows(rows, path, weighted_limit, extension_limit):
    """Reject incompatible protocols and unsafe coefficients on either side."""
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed CSV row: {path}")
        if row["retrieval_mode"] not in {"leave", "continue"}:
            raise ValueError(f"Unsupported retrieval mode: {path}")
        if row["error"] or row["Solver Status"] in {"", "ERROR"}:
            raise ValueError(f"Solver error at seed {row['seed']}, escorts {row['# Escorts']}: {path}")
        if row["protocol"] != PROTOCOL:
            raise ValueError(f"Incompatible CSV protocol: {path}")
        for key, expected_value in (("weighted_time_limit", weighted_limit),
                                    ("extension_time_limit", extension_limit),
                                    ("total_time_limit", float(weighted_limit) + float(extension_limit))):
            if float(row[key]) != float(expected_value):
                raise ValueError(f"Unexpected {key}: {path}")
        coefficient = int(row["flow_weight"])
        movement_upper = int(row["safe_movement_bound"])
        movement_lower = int(row["safe_movement_lower_bound"])
        if (row["warmstart"] != "1"
                or row["movement_integer_weight"] != "1"
                or movement_lower < 0 or movement_upper < movement_lower
                or coefficient != movement_upper - movement_lower + 1
                or row["weighted_global_scope"] != "1" or row["optimization_calls"] != "1"
                or row["flow_proof_check_mode"] != FLOW_PROOF_CHECK_MODE
                or float(row["movement_weight"]) != 1 / coefficient):
            raise ValueError(f"Unexpected warm-start, continuation, or safe objective setting: {path}")
        component_fields = ["greedy_flowtime", "greedy_movements"]
        for prefix in ("", "final_"):
            if row[prefix + "has_solution"] == "1":
                component_fields.extend((prefix + "flowtime", prefix + "movements"))
        for key in component_fields:
            value = float(row[key])
            if not math.isfinite(value) or value < movement_lower or int(value) != value:
                raise ValueError(f"Invalid {key} relative to distance lower bound: {path}")
        embedded_value(row)


def merge_batch(source, destination, seeds, escorts, weighted_limit, extension_limit):
    """Validate complete instance coverage and continuation protocol before append."""
    source, destination = Path(source), Path(destination)
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDNAMES:
            raise ValueError(f"Unexpected CSV schema: {source}")
        rows = list(reader)
    _validate_merge_rows(rows, source, weighted_limit, extension_limit)
    modes = {row["retrieval_mode"] for row in rows}
    if len(modes) != 1:
        raise ValueError(f"Cannot merge mixed retrieval modes: {source}")
    expected = {(seed, count) for seed in parse_range(seeds, minimum=0)
                for count in parse_range(escorts, minimum=1)}
    actual = [(int(row["seed"]), int(row["# Escorts"])) for row in rows]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError(f"Missing or duplicate instance results: {source}")
    exists = destination.exists()
    if exists:
        with destination.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != FIELDNAMES:
                raise ValueError(f"Cannot merge incompatible CSV schemas: {destination}")
            existing = list(reader)
        _validate_merge_rows(existing, destination, weighted_limit, extension_limit)
        if existing and {row["retrieval_mode"] for row in existing} != modes:
            raise ValueError(f"Cannot merge mixed retrieval modes: {destination}")
        key_fields = ("formulation", "Lx x Ly", "IOs", "# Escorts", "#Loads", "seed")
        keys = [tuple(row[key] for key in key_fields) for row in existing + rows]
        if len(set(keys)) != len(keys):
            raise ValueError(f"Cannot merge duplicate instance results: {destination}")
    with destination.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if not exists:
            writer.writeheader()
        writer.writerows(rows)
    print(f"Validated and saved {len(rows)} instance results to {destination}", flush=True)


def nonnegative_number(value):
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise argparse.ArgumentTypeError("Expected a finite nonnegative number")
    return result


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formulation", choices=("escortflow", "loadflow"), required=True)
    parser.add_argument("-x", dest="Lx", type=int, required=True)
    parser.add_argument("-y", dest="Ly", type=int, required=True)
    parser.add_argument("-O", "--outputs", type=int, nargs="+", required=True)
    parser.add_argument("-e", "--escorts", default="3-8")
    parser.add_argument("-l", "--loads", type=int, default=1)
    parser.add_argument("-r", "--seeds", default="1-100")
    parser.add_argument("-m", "--retrieval-mode", choices=("leave", "continue"), default="leave")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--weighted-time-limit", type=positive_number, default=300,
                        help="initial weighted search budget in solver seconds (default: 300)")
    parser.add_argument("--extension-time-limit", "--certification-time-limit", dest="extension_time_limit",
                        type=nonnegative_number, default=300,
                        help="additional search budget only if first-phase flow is unproved (default: 300)")
    parser.add_argument("-f", "--output", type=Path, required=True)
    parser.add_argument('--lp', action='store_true', help='Generate the same instances and solve their continuous relaxations without integer CSV inputs')
    parser.add_argument('--with-lp', action='store_true', help='After each integer instance, save its continuous LP bound in the same CSV row')
    parser.add_argument('--lp-threads', type=int, default=1, help='Threads for each integrated LP solve (default 1)')
    parser.add_argument('--lp-time-limit', type=positive_number, default=300, help='Integrated LP budget, separate from integer budgets')
    parser.add_argument('--lp-protocol', choices=['v4', 'v5'], default='v5', help='LP coefficient/horizon conventions: archived v4 or current v5')
    parser.add_argument('--lp-workers', type=int, default=1)
    parser.add_argument('--lp-retry-time-limit', type=positive_number, default=600)
    parser.add_argument('--resume', action='store_true', help='Resume direct --lp results')
    args = parser.parse_args(argv)
    try:
        if args.Lx <= 0 or args.Ly <= 0 or args.loads <= 0 or args.threads <= 0:
            raise ValueError("Grid dimensions, loads, and threads must be positive")
        if len(args.outputs) % 2:
            raise ValueError("Outputs must be x y coordinate pairs")
        args.outputs = sorted(zip(args.outputs[::2], args.outputs[1::2]))
        if len(set(args.outputs)) != len(args.outputs) or any(
                not (0 <= x < args.Lx and 0 <= y < args.Ly) for x, y in args.outputs):
            raise ValueError("Outputs must be distinct cells inside the grid")
        args.seed_values = parse_range(args.seeds, minimum=0)
        args.escort_values = parse_range(args.escorts, minimum=1)
        if max(args.escort_values) + args.loads > args.Lx * args.Ly:
            raise ValueError("Targets and escorts exceed the number of grid cells")
        if args.resume and not args.lp:
            raise ValueError('--resume requires --lp')
        if args.lp_protocol != 'v5' and not args.lp:
            raise ValueError('--lp-protocol v4 requires --lp; integer runs use v5')
        if args.with_lp and args.lp:
            raise ValueError('Choose --with-lp for integer plus LP, or --lp for LP only')
        if args.lp_threads <= 0:
            raise ValueError('LP threads must be positive')
        if args.output.exists() and not (args.lp and args.resume):
            raise ValueError(f"Output file already exists: {args.output}")
    except ValueError as exc:
        parser.error(str(exc))
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.lp:
        from static_generated_lp import run_standard
        try:
            return run_standard(args, args.formulation)
        except (ValueError, OSError) as exc:
            print(f'Direct LP error: {exc}', file=sys.stderr)
            return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lp_incomplete = False
    with args.output.open("x+", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        handle.flush()
        for escort_count in args.escort_values:
            for seed in args.seed_values:
                row = run_instance(args, seed, escort_count)
                if args.with_lp and not row['error']:
                    row.update(lp_requested=1, lp_status='PENDING', lp_threads=args.lp_threads,
                               lp_time_limit=args.lp_time_limit, lp_retry_time_limit=args.lp_retry_time_limit)
                position = handle.tell()
                writer.writerow(row)
                handle.flush()
                if row["error"]:
                    return 1
                if args.with_lp:
                    # The integer checkpoint is already durable if the LP is
                    # interrupted. Only this last row is replaced after solving.
                    attach_lp(row, args, Path(__file__).resolve().parent)
                    handle.seek(position)
                    writer.writerow(row)
                    handle.flush()
                    handle.truncate()
                    handle.flush()
                    lp_incomplete |= row['lp_status'] != 'OPTIMAL'
    # Distinguish saved, valid integer rows with missing LPs from integer errors.
    return 3 if lp_incomplete else 0


if __name__ == "__main__":
    raise SystemExit(main())
