#!/usr/bin/env python3
"""Weighted static PBS experiments, with a separate flow-time certificate.

Both weighted formulations receive the same deterministic greedy physical plan.
Certification starts from the saved weighted incumbent. The weighted candidate
is retained even if certification finds a better flow time. See
weighted_flow_certification_notes.md for proof scope.
"""

import argparse
import csv
from datetime import datetime
import itertools
import json
import math
from pathlib import Path
import platform
import time

from PBSCom import GeneretaeRandomInstance, str2range
import OneStepHeuristic_v2
from static_weighted_certification import certification_eligibility, sufficient_flow_horizon


FIELDNAMES = [
    "host", "date", "formulation", "Lx x Ly", "#IOs", "# Escorts", "#Loads",
    "IOs", "Escorts", "Target Loads", "seed", "retrieval_mode", "movement_mode",
    "movement_weight", "warmstart", "threads", "weighted_horizon",
    "certification_horizon", "global_flow_horizon_requirement",
    "greedy_makespan", "greedy_flowtime", "greedy_movements",
    "greedy_objective", "greedy_time", "weighted_time_limit",
    "certification_time_limit", "certification_gap_threshold", "Solver Status",
    "has_solution", "makespan", "flowtime", "movements", "objective",
    "best_bound", "absolute_gap", "scaled_objective", "scaled_best_bound",
    "scaled_absolute_gap", "weighted_bound_consistent", "weighted_proven", "weighted_cpu_time",
    "weighted_runtime", "weighted_work", "certification_eligible",
    "certification_warmstart_source", "certification_start_removed_post_retrieval_movements",
    "certification_skip_reason", "certification_status", "certification_cpu_time",
    "certification_runtime", "certification_work", "flow_lower_bound",
    "certification_best_flowtime", "certification_best_movements",
    "flow_proven", "counterexample", "counterexample_file", "lexicographic_proven",
    "lexicographic_proven_within_weighted_horizon", "total_wall_time", "error",
]


def positive_number(value):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("must be a finite positive number")
    return value


def certification_threshold(value):
    value = positive_number(value)
    if value >= 1:
        raise argparse.ArgumentTypeError("must be strictly between zero and one")
    return value


def parse_range(value, *, minimum):
    values = list(str2range(value))
    if not values or len(values) != len(set(values)) or min(values) < minimum:
        raise ValueError(f"Range must contain distinct integers >= {minimum}: {value}")
    return values


def make_solver(args, *, objective_mode, time_limit, certification_target=None):
    common = dict(
        Lx=args.Lx, Ly=args.Ly, output_cells=tuple(args.outputs),
        beta=1, gamma=0.01, time_limit=time_limit, threads=args.threads,
        mip_focus=0 if objective_mode == "weighted_integer" else 3,
        objective_mode=objective_mode, certification_target=certification_target,
    )
    if args.formulation == "escortflow":
        from escort_flow_static_gurobi import StaticGurobiConfig, StaticEscortFlowGurobiSolver
        return StaticEscortFlowGurobiSolver(StaticGurobiConfig(retrieval_mode="leave", **common))
    from load_flow_static_gurobi import LoadFlowStaticGurobiConfig, LoadFlowStaticGurobiSolver
    return LoadFlowStaticGurobiSolver(LoadFlowStaticGurobiConfig(move_method="BM", alpha=0, **common))


def run_instance(args, seed, escort_count, weighted_solver, solver_factory=make_solver):
    start = time.perf_counter()
    locations = sorted(itertools.product(range(args.Lx), range(args.Ly)))
    targets, escorts = GeneretaeRandomInstance(seed, locations, escort_count, args.loads)
    row = dict.fromkeys(FIELDNAMES, "")
    row.update({
        "host": platform.node(), "date": datetime.now().astimezone().isoformat(),
        "formulation": args.formulation, "Lx x Ly": f"{args.Lx}x{args.Ly}",
        "#IOs": len(args.outputs), "# Escorts": escort_count, "#Loads": args.loads,
        "IOs": repr(args.outputs), "Escorts": repr(escorts), "Target Loads": repr(targets),
        "seed": seed, "retrieval_mode": "leave", "movement_mode": "BM",
        "movement_weight": 0.01, "warmstart": 1, "threads": args.threads,
        "weighted_time_limit": args.weighted_time_limit,
        "certification_time_limit": args.certification_time_limit,
        "certification_gap_threshold": args.certification_gap_threshold,
        "certification_status": "NOT_RUN", "certification_cpu_time": 0,
        "certification_runtime": 0, "certification_work": 0, "flow_proven": 0,
        "counterexample": 0, "lexicographic_proven": 0,
        "lexicographic_proven_within_weighted_horizon": 0,
    })
    print(f"[{row['date']}] seed={seed} escorts={escort_count}: greedy and weighted solve starting", flush=True)
    try:
        greedy_start = time.perf_counter()
        trace = OneStepHeuristic_v2.SolveGreedy(
            args.Lx, args.Ly, set(args.outputs), set(targets), set(escorts),
            verbal=False, max_steps=max(1, (args.Lx + args.Ly) * args.loads * 20 // escort_count),
            retrieval_mode="leave", return_trace=True,
        )
        makespan, flow, movements, _, escort_history, target_history = trace
        horizon = max(0, makespan - 1) if args.formulation == "escortflow" else makespan + 1
        row.update(greedy_makespan=makespan, greedy_flowtime=flow,
                   greedy_movements=movements, greedy_objective=flow + 0.01 * movements,
                   greedy_time=time.perf_counter() - greedy_start, weighted_horizon=horizon)
        warmstart = weighted_solver.build_warmstart_from_trace(
            targets, escorts, horizon, target_history, escort_history)
        weighted = weighted_solver.solve(targets, escorts, horizon, warmstart=warmstart)
        row["Solver Status"] = weighted["status_name"]
        for key in ("has_solution", "makespan", "flowtime", "movements", "objective",
                    "best_bound", "absolute_gap", "weighted_proven"):
            value = weighted.get(key)
            row[key] = int(value) if isinstance(value, bool) else value
        for key in ("cpu_time", "runtime", "work"):
            row[f"weighted_{key}"] = weighted.get(key)
        row["weighted_bound_consistent"] = int(bool(weighted.get("bound_consistent", True)))
        for key in ("objective", "best_bound", "absolute_gap"):
            value = weighted.get(key)
            row[f"scaled_{key}"] = None if value is None else 100 * value
        if weighted.get("has_solution"):
            row["scaled_objective"] = 100 * weighted["flowtime"] + weighted["movements"]

        eligible, reason = certification_eligibility(weighted, args.certification_gap_threshold)
        row.update(certification_eligible=int(eligible), certification_skip_reason=reason)
        if eligible:
            requirement = sufficient_flow_horizon(targets, args.outputs, weighted["flowtime"])
            certificate_horizon = requirement
            row.update(global_flow_horizon_requirement=requirement,
                       certification_horizon=certificate_horizon)
            print(f"seed={seed} escorts={escort_count}: weighted {weighted['status_name']}, "
                  f"F={weighted['flowtime']} M={weighted['movements']} gap={weighted.get('absolute_gap')}; "
                  f"flow certificate starting (T={certificate_horizon}, "
                  f"separate limit={args.certification_time_limit:g}s, warm start=weighted solution)", flush=True)
            certificate_solver = solver_factory(
                args, objective_mode="flow_certificate", time_limit=args.certification_time_limit,
                certification_target=int(weighted["flowtime"]),
            )
            try:
                certificate_start = certificate_solver.build_warmstart_from_solution(
                    weighted, certificate_horizon)
                row["certification_warmstart_source"] = "weighted_solution"
                row["certification_start_removed_post_retrieval_movements"] = certificate_start.get(
                    "removed_post_retrieval_movements", 0)
                certificate = certificate_solver.solve(
                    targets, escorts, certificate_horizon, warmstart=certificate_start)
            finally:
                certificate_solver.close()
            row.update(certification_status=certificate["status_name"],
                       flow_lower_bound=certificate.get("best_bound"),
                       certification_best_flowtime=certificate.get("flowtime"),
                       certification_best_movements=certificate.get("movements"),
                       flow_proven=int(bool(certificate.get("flow_proven"))),
                       counterexample=int(bool(certificate.get("counterexample"))))
            for key in ("cpu_time", "runtime", "work"):
                row[f"certification_{key}"] = certificate.get(key)
            within_horizon = bool(weighted.get("weighted_proven") and certificate.get("flow_proven"))
            row["lexicographic_proven_within_weighted_horizon"] = int(within_horizon)
            physical_horizon = horizon + 1 if args.formulation == "escortflow" else horizon
            row["lexicographic_proven"] = int(within_horizon and physical_horizon >= requirement)
            if certificate.get("counterexample"):
                witness_path = args.output.with_name(
                    f"{args.output.stem}_e{escort_count}_seed{seed}_counterexample.json")
                witness = dict(formulation=args.formulation, seed=seed, targets=targets,
                               escorts=escorts, outputs=args.outputs,
                               weighted_horizon=horizon, certification_horizon=certificate_horizon,
                               weighted_flowtime=weighted["flowtime"],
                               weighted_movements=weighted["movements"],
                               certificate_flowtime=certificate.get("flowtime"),
                               certificate_movements=certificate.get("movements"),
                               weighted_moves=weighted.get("animation_moves"),
                               certificate_moves=certificate.get("animation_moves"))
                with witness_path.open("x") as handle:
                    json.dump(witness, handle, indent=2, default=int)
                    handle.write("\n")
                row["counterexample_file"] = str(witness_path.resolve())
        else:
            print(f"seed={seed} escorts={escort_count}: certification skipped ({reason})", flush=True)
    except Exception as exc:
        row["error"] = f"{type(exc).__name__}: {exc}"
        if not row["Solver Status"]:
            row["Solver Status"] = "ERROR"
        row["certification_status"] = "ERROR"
        print(f"seed={seed} escorts={escort_count}: ERROR: {row['error']}", flush=True)
    row["total_wall_time"] = time.perf_counter() - start
    print(f"seed={seed} escorts={escort_count}: complete, weighted={row['Solver Status']}, "
          f"flow_proven={row['flow_proven']}, counterexample={row['counterexample']}, "
          f"wall={row['total_wall_time']:.3f}s", flush=True)
    return row


def merge_batch(source, destination, seeds, escorts, weighted_limit, certification_limit, threshold):
    """Validate complete coverage and protocol before appending a layout batch."""
    source, destination = Path(source), Path(destination)
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDNAMES:
            raise ValueError(f"Unexpected CSV schema: {source}")
        rows = list(reader)
    expected = {(seed, count) for seed in parse_range(seeds, minimum=0)
                for count in parse_range(escorts, minimum=1)}
    actual = [(int(row["seed"]), int(row["# Escorts"])) for row in rows]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError(f"Missing or duplicate instance results: {source}")
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed CSV row: {source}")
        if row["error"] or row["Solver Status"] in {"", "ERROR"}:
            raise ValueError(f"Solver error at seed {row['seed']}, escorts {row['# Escorts']}: {source}")
        for key, expected_value in (("weighted_time_limit", weighted_limit),
                                    ("certification_time_limit", certification_limit),
                                    ("certification_gap_threshold", threshold)):
            if float(row[key]) != float(expected_value):
                raise ValueError(f"Unexpected {key}: {source}")
        if row["warmstart"] != "1" or float(row["movement_weight"]) != 0.01:
            raise ValueError(f"Unexpected warm-start or objective setting: {source}")
    exists = destination.exists()
    if exists:
        with destination.open(newline="") as handle:
            if next(csv.reader(handle), None) != FIELDNAMES:
                raise ValueError(f"Cannot merge incompatible CSV schemas: {destination}")
    with destination.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if not exists:
            writer.writeheader()
        writer.writerows(rows)
    print(f"Validated and saved {len(rows)} instance results to {destination}", flush=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formulation", choices=("escortflow", "loadflow"), required=True)
    parser.add_argument("-x", dest="Lx", type=int, required=True)
    parser.add_argument("-y", dest="Ly", type=int, required=True)
    parser.add_argument("-O", "--outputs", type=int, nargs="+", required=True)
    parser.add_argument("-e", "--escorts", default="3-8")
    parser.add_argument("-l", "--loads", type=int, default=1)
    parser.add_argument("-r", "--seeds", default="1-100")
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--weighted-time-limit", type=positive_number, default=300)
    parser.add_argument("--certification-time-limit", type=positive_number, default=300)
    parser.add_argument("--certification-gap-threshold", type=certification_threshold, default=0.1,
                        help="maximum absolute weighted gap eligible for certification, in F+0.01M units (strict; default 0.1)")
    parser.add_argument("-f", "--output", type=Path, required=True)
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
        if args.output.exists():
            raise ValueError(f"Output file already exists: {args.output}")
    except ValueError as exc:
        parser.error(str(exc))
    return args


def main(argv=None):
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    weighted_solver = make_solver(args, objective_mode="weighted_integer", time_limit=args.weighted_time_limit)
    try:
        with args.output.open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
            writer.writeheader()
            handle.flush()
            for escort_count in args.escort_values:
                for seed in args.seed_values:
                    row = run_instance(args, seed, escort_count, weighted_solver)
                    writer.writerow(row)
                    handle.flush()
                    if row["error"]:
                        return 1
    finally:
        weighted_solver.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
