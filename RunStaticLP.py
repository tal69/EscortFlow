#!/usr/bin/env python3
"""Replay continuous PBS relaxations using recorded instances, weights, and horizons.

Accepts RunSafeWeightedStatic CSVs, including archived v4 and current v5 runs.
Only OPTIMAL LP values enter the output. Integer proof callbacks, warm starts,
greedy recomputation, fixed-flow constraints, and objective cutoffs are absent.
"""

import argparse
import ast
import concurrent.futures
import csv
import glob
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

HERE = Path(__file__).resolve().parent
MODEL_FILES = (
    "escort_flow_static_gurobi.py", "load_flow_static_gurobi.py",
    "static_lexicographic.py", "static_weighted_certification.py",
    "static_safe_weighted_search.py",
)
FIELDS = (
    "layout", "escorts", "loads", "seed", "method", "flow_weight", "horizon",
    "lp_objective", "lp_flow", "lp_movements", "status", "elapsed_seconds",
    "physical_horizon", "retrieval_mode", "problem_sha256", "source_sha256",
    "solver_version", "algorithm",
)


class LPNotOptimalError(RuntimeError):
    def __init__(self, message, status='ERROR'):
        super().__init__(message, status)
        self.status = status

    def __str__(self):
        return self.args[0]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def cells(text, name, lx, ly):
    value = ast.literal_eval(text)
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} must be a list of coordinate pairs")
    result = []
    for cell in value:
        if (not isinstance(cell, (list, tuple)) or len(cell) != 2
                or any(type(v) is not int for v in cell)
                or not (0 <= cell[0] < lx and 0 <= cell[1] < ly)):
            raise ValueError(f"Invalid {name} coordinate: {cell}")
        result.append(tuple(cell))
    if len(set(result)) != len(result):
        raise ValueError(f"Duplicate {name} coordinates")
    return sorted(result)


def parse_instance(row):
    """Freeze the model-defining fields; never derive R or T from a new heuristic."""
    lx, ly = map(int, row["Lx x Ly"].split("x"))
    method = row["formulation"]
    mode = row.get("retrieval_mode", "leave")
    coefficient, horizon = int(row["flow_weight"]), int(row["weighted_horizon"])
    if lx <= 0 or ly <= 0 or coefficient <= 0 or horizon < 0:
        raise ValueError("Positive dimensions and R, and a nonnegative horizon, are required")
    if method not in {"loadflow", "escortflow"} or mode not in {"leave", "continue"}:
        raise ValueError("Expected loadflow/escortflow and leave/continue retrieval")
    if row.get("movement_mode", row.get("Moves", "BM")).strip() != "BM":
        raise ValueError("The replay script supports the paper's BM movement regime")
    physical = horizon + (method == "escortflow")
    if row.get("weighted_physical_horizon") and int(row["weighted_physical_horizon"]) != physical:
        raise ValueError("Recorded physical horizon disagrees with the formulation's time indexing")
    if row.get("movement_weight") and not math.isclose(
            float(row["movement_weight"]), 1 / coefficient, rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError("Recorded movement weight is inconsistent with R")
    outputs = cells(row["IOs"], "outputs", lx, ly)
    targets = cells(row["Target Loads"], "targets", lx, ly)
    escorts = cells(row["Escorts"], "escorts", lx, ly)
    if not outputs or not escorts or set(targets) & set(escorts):
        raise ValueError("Outputs and escorts are required; targets cannot occupy escort cells")
    if len(targets) != int(row["#Loads"]) or len(escorts) != int(row["# Escorts"]):
        raise ValueError("Recorded target/escort counts disagree with the coordinates")
    result = dict(layout=f"{lx}x{ly}", Lx=lx, Ly=ly, method=method, seed=int(row["seed"]),
                  outputs=outputs, targets=targets, escort_cells=escorts, escorts=len(escorts), loads=len(targets),
                  flow_weight=coefficient, horizon=horizon, physical_horizon=physical,
                  retrieval_mode=mode)
    result["problem_sha256"] = fingerprint(result)
    return result


def parse_range(text):
    """Comma-separated integers and inclusive ranges, e.g., 1,3-8."""
    result = set()
    for part in text.split(","):
        bounds = part.split("-")
        if len(bounds) == 1:
            start = stop = int(bounds[0])
        elif len(bounds) == 2:
            start, stop = map(int, bounds)
        else:
            raise ValueError("Use comma-separated integers or inclusive ranges")
        if start < 0 or stop < start:
            raise ValueError("Ranges must be nonnegative and increasing")
        result.update(range(start, stop + 1))
    return result


def read_instances(paths, seeds=None, escorts=None, formulation=None):
    instances, labels = [], set()
    for path in paths:
        with path.open(newline="") as handle:
            for line, row in enumerate(csv.DictReader(handle), 2):
                try:
                    problem = parse_instance(row)
                    if ((seeds is not None and problem["seed"] not in seeds)
                            or (escorts is not None and problem["escorts"] not in escorts)
                            or (formulation and problem["method"] != formulation)):
                        continue
                    # A duplicate label must not hide a different R, horizon, or initial state.
                    label = (problem["method"], problem["layout"], problem["escorts"],
                             problem["loads"], problem["seed"], problem["retrieval_mode"])
                    if label in labels:
                        raise ValueError(f"Duplicate instance label {label}; run campaigns separately")
                    labels.add(label)
                    instances.append(problem)
                except (ValueError, KeyError, SyntaxError, TypeError) as exc:
                    raise ValueError(f"{path.name}:{line}: {exc}") from exc
    if not instances:
        raise ValueError("No instances match the inputs and filters")
    return sorted(instances, key=lambda p: (p["seed"], p["loads"], p["escorts"], p["layout"],
                                          p["method"], p["retrieval_mode"]))


def validate_value(row):
    if row["status"] != "OPTIMAL":
        raise ValueError("Only OPTIMAL LP results can be used or resumed")
    objective, flow, movements = (float(row[k]) for k in ("lp_objective", "lp_flow", "lp_movements"))
    if any(not math.isfinite(v) or v < -1e-6 for v in (objective, flow, movements)):
        raise ValueError("LP objective components must be finite and nonnegative")
    if not math.isclose(objective, flow + movements / int(row["flow_weight"]),
                        rel_tol=1e-8, abs_tol=1e-6):
        raise ValueError("LP objective disagrees with F + M/R")


def initialize_worker(source_dir, log_dir):
    # Each process gets a separate solver log and exactly the selected source tree.
    sys.path.insert(0, source_dir)
    for name in MODEL_FILES:
        sys.modules.pop(Path(name).stem, None)
    importlib.invalidate_caches()
    log = os.open(str(Path(log_dir) / f"worker_{os.getpid()}.log"),
                  os.O_CREAT | os.O_WRONLY | os.O_APPEND, 0o600)
    os.dup2(log, 1)
    os.dup2(log, 2)
    os.close(log)


def solve_instance(problem, settings):
    gp = importlib.import_module("gurobipy")
    if problem["method"] == "loadflow":
        module = importlib.import_module("load_flow_static_gurobi")
        constructor, config_class = module.LoadFlowStaticGurobiSolver, module.LoadFlowStaticGurobiConfig
        extra = dict(move_method="BM", alpha=0)
    else:
        module = importlib.import_module("escort_flow_static_gurobi")
        constructor, config_class = module.StaticEscortFlowGurobiSolver, module.StaticGurobiConfig
        extra = {}
    config = dict(Lx=problem["Lx"], Ly=problem["Ly"], output_cells=tuple(problem["outputs"]),
                  retrieval_mode=problem["retrieval_mode"], beta=1, gamma=1 / problem["flow_weight"],
                  lp=True, threads=settings["threads"], **extra)
    start = time.perf_counter()
    algorithm = "automatic"
    for attempt, limit in enumerate((settings["time_limit"], settings["retry_time_limit"])):
        solver = constructor(config_class(time_limit=limit, **config))
        try:
            if attempt:
                algorithm = "barrier_retry"
                solver.env.setParam("Method", 2)
                solver.env.setParam("Crossover", 0)
            result = solver.solve(problem["targets"], problem["escort_cells"], problem["horizon"])
        finally:
            solver.close()
        if result["status_name"] != "TIME_LIMIT" or not settings["retry_time_limit"] or attempt:
            break
    if result["status_name"] != "OPTIMAL" or not result["has_solution"]:
        raise LPNotOptimalError(f"{problem['method']} {problem['layout']} e={problem['escorts']} "
                                f"loads={problem['loads']} seed={problem['seed']}: LP status {result['status_name']}",
                                result['status_name'])
    row = {key: problem[key] for key in ("layout", "escorts", "loads", "seed", "flow_weight", "method",
                                       "horizon", "physical_horizon", "retrieval_mode", "problem_sha256")}
    row.update(lp_objective=result["objective"], lp_flow=result["flowtime"],
               lp_movements=result["movements"], status=result["status_name"],
               elapsed_seconds=time.perf_counter() - start, source_sha256=settings["source_sha256"],
               solver_version=".".join(map(str, gp.gurobi.version())), algorithm=algorithm)
    validate_value(row)
    return row


def reference_key(row):
    # Archived one-target reference CSVs predate explicit loads/mode columns.
    return (row["layout"], int(row["escorts"]), int(row.get("loads", 1)),
            int(row["seed"]), row["method"], row.get("retrieval_mode", "leave"))


def read_reference(path):
    reference = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            validate_value(row)
            key = reference_key(row)
            if key in reference:
                raise ValueError(f"Duplicate reference key {key}")
            reference[key] = row
    return reference


def check_reference(row, reference):
    expected = reference[reference_key(row)]
    if (int(row["flow_weight"]), int(row["horizon"])) != (
            int(expected["flow_weight"]), int(expected["horizon"])):
        raise ValueError("Reference R or horizon differs from the replay")
    # Fractional F and M can vary across alternative optimal solutions. Compare Z only.
    if not math.isclose(float(row["lp_objective"]), float(expected["lp_objective"]),
                        rel_tol=1e-8, abs_tol=1e-6):
        raise ValueError(f"LP objective differs from reference for {reference_key(row)}")


def positive_number(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("Expected a finite positive number")
    return number


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True, help="recorded integer CSVs or quoted glob patterns")
    parser.add_argument("-f", "--output", type=Path, required=True, help="fresh LP result CSV")
    parser.add_argument("--source-dir", type=Path, default=HERE, help="complete model snapshot (default: this code directory)")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1, help="solver threads per worker (default: 1)")
    parser.add_argument("--time-limit", type=positive_number, default=300)
    parser.add_argument("--retry-time-limit", type=positive_number, default=600,
                        help="barrier retry budget after TIME_LIMIT (default: 600)")
    parser.add_argument("--no-retry", action="store_true")
    parser.add_argument("--seeds", help="optional filter, e.g., 1-100 or 1,44")
    parser.add_argument("--escorts", help="optional filter, e.g., 3-8")
    parser.add_argument("--formulation", choices=("loadflow", "escortflow"))
    parser.add_argument("--resume", action="store_true", help="verify the manifest and skip saved OPTIMAL results")
    parser.add_argument("--extend", action="store_true",
                        help="with --resume, allow additional instances; all previous instances and model sources must remain unchanged")
    parser.add_argument("--check-against", type=Path, help="verify LP objectives against a reference LP CSV")
    parser.add_argument('--fill-input-lp', action='store_true',
                        help='fill missing LP columns in new integrated integer CSVs after retry; historical CSVs stay read-only')
    args = parser.parse_args(argv)
    try:
        if args.workers <= 0 or args.threads <= 0:
            raise ValueError("Workers and threads must be positive")
        if args.extend and not args.resume:
            raise ValueError("--extend requires --resume")
        paths = []
        for pattern in args.input:
            matches = sorted(glob.glob(pattern))
            if not matches:
                raise ValueError(f"Input does not exist or match any files: {pattern}")
            paths.extend(Path(p).resolve() for p in matches)
        if len(paths) != len(set(paths)):
            raise ValueError("An input CSV was supplied more than once")
        args.paths = paths
        args.source_dir = args.source_dir.resolve()
        args.seed_values = parse_range(args.seeds) if args.seeds else None
        args.escort_values = parse_range(args.escorts) if args.escorts else None
        for name in MODEL_FILES:
            if not (args.source_dir / name).is_file():
                raise ValueError(f"Incomplete model source tree: missing {name}")
        if args.output.resolve() in paths or (args.check_against is not None
                                             and args.output.resolve() == args.check_against.resolve()):
            raise ValueError("The output must not overwrite an input or reference")
        if args.output.exists() and not args.resume:
            raise ValueError("Output already exists; use a fresh file or --resume")
        if args.no_retry:
            args.retry_time_limit = None
    except ValueError as exc:
        parser.error(str(exc))
    return args


def run(args):
    instances = read_instances(args.paths, args.seed_values, args.escort_values, args.formulation)
    sources = {name: sha256(args.source_dir / name) for name in MODEL_FILES}
    source_hash = fingerprint(sources)
    identity = dict(objective="F + M/R", model_sources=sources,
                    problem_sha256=[p["problem_sha256"] for p in instances])
    manifest_path = Path(str(args.output) + ".manifest.json")
    reference = read_reference(args.check_against) if args.check_against else None
    if reference is not None:
        for problem in instances:
            expected = reference.get(reference_key(problem))
            if (expected is None or int(expected["flow_weight"]) != problem["flow_weight"]
                    or int(expected["horizon"]) != problem["horizon"]):
                raise ValueError(f"Missing or incompatible reference: {reference_key(problem)}")
    completed = {}
    problems = {p["problem_sha256"]: p for p in instances}
    if args.output.exists():
        if not manifest_path.exists():
            raise ValueError("Resume refused: model sources or selected instances changed, or manifest missing")
        manifest = json.loads(manifest_path.read_text())
        previous = manifest["identity"]
        if args.extend:
            if (previous["objective"] != identity["objective"]
                    or previous["model_sources"] != identity["model_sources"]
                    or not set(previous["problem_sha256"]) <= set(identity["problem_sha256"])):
                raise ValueError("Resume refused: previous instances were removed/changed or model sources changed")
        elif previous != identity:
            raise ValueError("Resume refused: model sources or selected instances changed, or manifest missing")
        with args.output.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != list(FIELDS):
                raise ValueError("Resume refused: incompatible LP CSV schema")
            for row in reader:
                validate_value(row)
                key = row["problem_sha256"]
                if key not in problems or key in completed or row["source_sha256"] != source_hash:
                    raise ValueError("Resume refused: unknown, duplicate, or changed-source result")
                expected = problems[key]
                if (reference_key(row) != reference_key(expected)
                        or int(row["horizon"]) != expected["horizon"]
                        or int(row["flow_weight"]) != expected["flow_weight"]
                        or int(row["physical_horizon"]) != expected["physical_horizon"]
                        or row["retrieval_mode"] != expected["retrieval_mode"]):
                    raise ValueError("Resume refused: saved metadata disagrees with the instance")
                if reference is not None:
                    check_reference(row, reference)
                completed[key] = row
        if args.extend:
            manifest.update(identity=identity, inputs=[dict(file=p.name, sha256=sha256(p)) for p in args.paths])
    else:
        if manifest_path.exists():
            raise ValueError("Manifest already exists without its output; choose a fresh output")
        manifest = dict(identity=identity, inputs=[dict(file=p.name, sha256=sha256(p)) for p in args.paths],
                        python=platform.python_version(), platform=platform.platform(), sessions=[])
    # Reuse the optimal relaxations already saved in each integer CSV row.
    # The replay remains a recovery tool for older campaigns and failed LPs.
    from static_integrated_lp import read_embedded
    imported = []
    for key, row in read_embedded(args.paths, problems, source_hash).items():
        if reference is not None:
            check_reference(row, reference)
        if key in completed:
            if not math.isclose(float(row['lp_objective']), float(completed[key]['lp_objective']),
                                abs_tol=1e-6, rel_tol=1e-8):
                raise ValueError('Embedded and cached LP objectives disagree')
        else:
            completed[key] = row
            imported.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    settings = dict(threads=args.threads, time_limit=args.time_limit, retry_time_limit=args.retry_time_limit,
                    source_sha256=source_hash)
    session = dict(started_unix=time.time(), workers=args.workers,
                   python=platform.python_version(), platform=platform.platform(), **settings)
    manifest["sessions"].append(session)

    def save_progress():
        manifest.update(total=len(instances), completed=len(completed), updated_unix=time.time())
        temporary_manifest = Path(str(manifest_path) + ".tmp")
        temporary_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
        temporary_manifest.replace(manifest_path)

    save_progress()
    if imported:
        with args.output.open('a' if args.output.exists() else 'x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            if handle.tell() == 0:
                writer.writeheader()
            writer.writerows(imported)
            handle.flush()
        print(f'Reused {len(imported)} optimal LP bounds from integer CSV rows.', flush=True)
    pending = [p for p in instances if p["problem_sha256"] not in completed]
    print(f"LP replay: {len(completed)}/{len(instances)} saved; {len(pending)} to solve.", flush=True)
    failures = []
    if pending:
        logs = Path(str(args.output) + ".logs")
        logs.mkdir(exist_ok=True)
        with args.output.open("a" if args.output.exists() else "x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            if handle.tell() == 0:
                writer.writeheader()
                handle.flush()
            with concurrent.futures.ProcessPoolExecutor(
                    max_workers=args.workers, initializer=initialize_worker,
                    initargs=(str(args.source_dir), str(logs))) as pool:
                iterator = iter(pending)
                futures = {pool.submit(solve_instance, p, settings): p
                           for p in [next(iterator, None) for _ in range(2 * args.workers)] if p is not None}
                while futures:
                    done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                    for future in done:
                        problem = futures.pop(future)
                        try:
                            row = future.result()
                            if reference is not None:
                                check_reference(row, reference)
                        except Exception as exc:
                            failures.append(dict(problem_sha256=problem["problem_sha256"], error=str(exc)))
                            print(str(exc), file=sys.stderr, flush=True)
                        else:
                            writer.writerow(row)
                            handle.flush()
                            completed[row["problem_sha256"]] = row
                        following = next(iterator, None)
                        if following is not None:
                            futures[pool.submit(solve_instance, following, settings)] = following
                    save_progress()
                    print(f"LP replay: {len(completed)}/{len(instances)} OPTIMAL.", flush=True)
    if getattr(args, 'fill_input_lp', False):
        from static_integrated_lp import fill_inputs
        filled = fill_inputs(args.paths, completed, source_hash, settings)
        manifest['inputs'] = [dict(file=p.name, sha256=sha256(p)) for p in args.paths]
        print(f'Filled {filled} missing LP values in integrated integer CSVs.', flush=True)
    session.update(finished_unix=time.time(), failures=failures)
    save_progress()
    if failures:
        Path(str(args.output) + ".failures.json").write_text(json.dumps(failures, indent=2) + "\n")
        print("Incomplete replay. OPTIMAL values are saved; use --resume to retry missing cases.", file=sys.stderr)
        return 1
    print(f"Complete: {len(completed)} OPTIMAL LP values" + (", matching the reference." if reference else "."))
    return 0


def main(argv=None):
    args = parse_args(argv)
    try:
        return run(args)
    except (ValueError, KeyError, OSError) as exc:
        print(f"LP replay error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
