#!/usr/bin/env python3
"""Reproduce the revised formulation paper's numerical tables in one command.

Examples: python3 ReproducePaper.py; python3 ReproducePaper.py 1-10.
Runs both formulations, leave/continue retrieval, integer searches and LPs,
then writes CSV metrics and portrait LaTeX table fragments. See README.md.
"""

import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from RunStaticCampaign import LAYOUTS, PREFLIGHT, automatic_threads, finite_number
from PaperTables import build_tables

SOURCE_FILES = (
    "ReproducePaper.py", "PaperTables.py", "RunStaticCampaign.py",
    "RunSafeWeightedStatic.py", "RunWeightedStatic.py", "RunStaticLP.py",
    "PBSCom.py", "OneStepHeuristic_v2.py", "static_generated_lp.py", "static_integrated_lp.py", "static_lexicographic.py",
    "static_weighted_certification.py", "static_safe_weighted_search.py",
    "escort_flow_static_gurobi.py", "load_flow_static_gurobi.py",
)
METHODS = ("loadflow", "escortflow")


def paper_configurations(layouts):
    """Union of the revised leave benchmarks and the new continue campaign."""
    for lx, ly, occupancy_escorts, outputs in LAYOUTS:
        if f"{lx}x{ly}" not in layouts:
            continue
        for mode in ("leave", "continue"):
            for loads in ((1, 2, 4, 6) if mode == "leave" else (2, 4, 6)):
                escorts = list(range(3, 9)) if loads == 1 else [8, 12, 16]
                if loads == 6 and mode == "leave":
                    escorts.append(20)
                for count in sorted(set(escorts)):
                    yield dict(lx=lx, ly=ly, escorts=count, loads=loads, outputs=list(outputs), mode=mode)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def integer_command(python, source, destination, config, seeds, settings, method):
    return [python, "-u", str(source / "RunSafeWeightedStatic.py"),
            "--formulation", method, "-x", str(config["lx"]), "-y", str(config["ly"]),
            "-O", *map(str, config["outputs"]), "-e", str(config["escorts"]),
            "-l", str(config["loads"]), "--retrieval-mode", config["mode"],
            "--seeds", seeds, "--threads", str(settings["threads"]),
            "--weighted-time-limit", str(settings["weighted_time_limit"]),
            "--extension-time-limit", str(settings["extension_time_limit"]),
            "-f", str(destination), '--with-lp', '--lp-threads', str(settings.get('lp_threads', 1)),
            '--lp-time-limit', str(settings.get('lp_time_limit', 300)),
            '--lp-retry-time-limit', str(settings.get('lp_retry_time_limit', 600))]


def stem(config, method):
    return f"{method}_{config['lx']}x{config['ly']}_l{config['loads']}_e{config['escorts']}"


def lp_command(python, source, root, mode, settings):
    return [python, "-u", str(source / "RunStaticLP.py"), "--input",
            *(str(root / mode / f"{method}.csv") for method in METHODS),
            "--source-dir", str(source), "--workers", "1", "--threads", str(settings["lp_threads"]),
            "--time-limit", str(settings["lp_time_limit"]),
            "--retry-time-limit", str(settings["lp_retry_time_limit"]),
            "-f", str(root / mode / "lp_results.csv"), "--resume", '--fill-input-lp']


def write_json(path, value):
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def run_integer_batches(python, root, manifest):
    # Import only after selecting the frozen tree, including for resume.
    from RunSafeWeightedStatic import FIELDNAMES, _validate_merge_rows
    from RunStaticLP import parse_instance
    source = root / "source"
    settings, seeds = manifest["settings"], manifest["seed_values"]

    def read_checked(path, config, method):
        if not path.exists():
            return []
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != FIELDNAMES:
                raise ValueError(f"Incompatible or truncated CSV header: {path}")
            rows = list(reader)
        _validate_merge_rows(rows, path, settings["weighted_time_limit"], settings["extension_time_limit"])
        actual = []
        for row in rows:
            problem = parse_instance(row)
            if (problem["layout"] != f"{config['lx']}x{config['ly']}"
                    or problem["method"] != method or problem["retrieval_mode"] != config["mode"]
                    or problem["loads"] != config["loads"] or problem["escorts"] != config["escorts"]
                    or problem["outputs"] != sorted(zip(config["outputs"][::2], config["outputs"][1::2]))
                    or int(row["threads"]) != settings["threads"]
                    or row["has_solution"] != "1" or row["final_has_solution"] != "1"):
                raise ValueError(f"Saved configuration or settings differ: {path}")
            actual.append(problem["seed"])
        if len(actual) != len(set(actual)) or not set(actual) <= set(seeds):
            raise ValueError(f"Duplicate or unexpected saved seeds: {path}")
        return rows

    def save_completed(path, rows):
        temporary = path.with_suffix(".csv.tmp")
        with temporary.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda row: int(row["seed"])))
        temporary.replace(path)

    for config in manifest["configurations"]:
        for method in METHODS:
            path = root / config["mode"] / "parts" / (stem(config, method) + ".csv")
            staging = path.with_suffix(".pending.csv")
            saved = read_checked(path, config, method)
            if staging.exists():
                keyed = {int(row["seed"]): row for row in saved}
                for row in read_checked(staging, config, method):
                    seed = int(row["seed"])
                    if seed in keyed and keyed[seed] != row:
                        raise ValueError("Conflicting saved rows for seed " + str(seed))
                    keyed[seed] = row
                saved = list(keyed.values())
                save_completed(path, saved)
                staging.unlink()
            pending = sorted(set(seeds) - {int(row["seed"]) for row in saved})
            if not pending:
                print(f"Reusing complete batch: {config['mode']}/{path.name}", flush=True)
                continue
            # The standard runner intentionally requires a fresh file. Keep
            # missing seeds separately until their rows have been validated.
            command = integer_command(python, source, staging, config, ",".join(map(str, pending)), settings, method)
            log = root / config["mode"] / "logs" / (stem(config, method) + ".log")
            print(f"[{datetime.now().astimezone().isoformat()}] {config['mode']} {path.stem}: "
                  f"{len(pending)} seeds to run (log: {log})", flush=True)
            with log.open("a") as handle:
                handle.write("\n" + shlex.join(command) + "\n")
                handle.flush()
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
                if result.returncode not in (0, 3):
                    result.check_returncode()
            new_rows = read_checked(staging, config, method)
            if {int(row["seed"]) for row in new_rows} != set(pending):
                raise ValueError("Missing results in completed batch: " + str(staging))
            save_completed(path, saved + new_rows)
            staging.unlink()
    # Rebuild merged inputs deterministically so resume cannot append duplicates.
    from RunSafeWeightedStatic import merge_batch
    for mode in ("leave", "continue"):
        for method in METHODS:
            destination = root / mode / f"{method}.csv"
            temporary = destination.with_suffix(".csv.tmp")
            temporary.unlink(missing_ok=True)
            for config in manifest["configurations"]:
                if config["mode"] == mode:
                    merge_batch(root / mode / "parts" / (stem(config, method) + ".csv"), temporary,
                                manifest["seeds"], str(config["escorts"]),
                                settings["weighted_time_limit"], settings["extension_time_limit"])
            temporary.replace(destination)


def parse_args(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("seed_range", nargs="?", help="optional inclusive seed range, e.g. 1-10")
    parser.add_argument("--seeds", "-r", help="same as the positional seed range; default 1-100")
    parser.add_argument("--python", default=os.environ.get("PYTHON") or sys.executable,
                        help="Python 3.10+ with NumPy, gurobipy, and a working license")
    parser.add_argument("--threads", type=int, default=automatic_threads(), help="threads per integer solve")
    parser.add_argument("--lp-threads", type=int, default=1, help="threads per LP; LPs run sequentially")
    parser.add_argument("--layouts", nargs="+", choices=[f"{x}x{y}" for x, y, _, _ in LAYOUTS],
                        default=[f"{x}x{y}" for x, y, _, _ in LAYOUTS])
    parser.add_argument("--weighted-time-limit", type=finite_number, default=300)
    parser.add_argument("--extension-time-limit", type=finite_number, default=300)
    parser.add_argument("--lp-time-limit", type=finite_number, default=300)
    parser.add_argument("--lp-retry-time-limit", type=finite_number, default=600)
    parser.add_argument("--reference-input", nargs="+", type=Path,
                        help="optional previous integer CSVs to include in the common best-known solutions")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--output-dir", type=Path, help="fresh output directory")
    group.add_argument("--resume", type=Path, help="resume a results directory using its saved settings and sources")
    group.add_argument("--tables-only", type=Path, help="rebuild tables from a completed results directory; no solves")
    parser.add_argument("--dry-run", action="store_true", help="print plan and commands; no files, imports of Gurobi, or solves")
    args = parser.parse_args(argv)
    if args.seed_range and args.seeds:
        parser.error("Supply the seed range once, either positionally or using --seeds")
    from PBSCom import str2range
    try:
        args.seeds = args.seed_range or args.seeds or "1-100"
        args.seed_values = sorted(str2range(args.seeds))
        if (not args.seed_values or min(args.seed_values) < 0
                or len(args.seed_values) != len(set(args.seed_values))):
            raise ValueError("Seeds must be distinct nonnegative integers")
        if len(set(args.layouts)) != len(args.layouts):
            raise ValueError("Select each layout only once")
        if min(args.threads, args.lp_threads, args.weighted_time_limit,
               args.lp_time_limit, args.lp_retry_time_limit) <= 0:
            raise ValueError("Threads and initial/LP time limits must be positive")
    except ValueError as exc:
        parser.error(str(exc))
    if (args.resume or args.tables_only) and (args.seed_range or args.seeds != "1-100" or args.reference_input):
        parser.error("Resume/table rebuild uses saved seeds and reference inputs")
    return args


def main(argv=None):
    args = parse_args(argv)
    here = Path(__file__).resolve().parent
    existing = args.resume or args.tables_only
    root = (existing or args.output_dir or here / (
        "results_paper_" + datetime.now().astimezone().strftime("%Y%m%d_%H%M%S") + f"_{os.getpid()}")).resolve()
    source = root / "source"
    if existing:
        manifest = json.loads((root / "manifest.json").read_text())
        for name, expected in manifest["source_sha256"].items():
            if sha256(source / name) != expected:
                raise ValueError("Frozen source changed: " + name)
        for item in manifest["reference_inputs"]:
            if sha256(root / item["file"]) != item["sha256"]:
                raise ValueError("Saved reference input changed: " + item["file"])
        # Execute the saved launcher, even if the user's checkout changed later.
        if here != source:
            option = "--tables-only" if args.tables_only else "--resume"
            subprocess.run([args.python, "-u", str(source / "ReproducePaper.py"), option, str(root),
                            "--python", args.python, *(["--dry-run"] if args.dry_run else [])], check=True)
            return 0
    else:
        settings = {name: getattr(args, name) for name in (
            "threads", "lp_threads", "weighted_time_limit", "extension_time_limit",
            "lp_time_limit", "lp_retry_time_limit")}
        manifest = dict(format_version=1, seeds=args.seeds, seed_values=args.seed_values,
                        settings=settings, configurations=list(paper_configurations(args.layouts)), reference_inputs=[])
    count = len(manifest["configurations"]) * len(manifest["seed_values"])
    print(f"Paper reproduction: {count} matched instances; {2 * count} integer runs and {2 * count} LP solves.", flush=True)
    print("Results: " + str(root), flush=True)
    if args.dry_run:
        for config in manifest["configurations"]:
            for method in METHODS:
                print(shlex.join(integer_command(args.python, source, root / config["mode"] / "parts" /
                    (stem(config, method) + ".csv"), config, manifest["seeds"], manifest["settings"], method)))
        for mode in ("leave", "continue"):
            print(shlex.join(lp_command(args.python, source, root, mode, manifest["settings"])))
        return 0
    if not existing:
        if root.exists():
            raise ValueError("Output directory already exists; use --resume or a fresh directory")
        for mode in ("leave", "continue"):
            for name in ("parts", "logs"):
                (root / mode / name).mkdir(parents=True, exist_ok=True)
        source.mkdir()
        for name in SOURCE_FILES:
            shutil.copy2(here / name, source / name)
        for name in ("README.md", "requirements.txt", "LP_REPRODUCIBILITY.md"):
            shutil.copy2(here / name, source / name)
        manifest.update(started=datetime.now().astimezone().isoformat(),
                        source_sha256={name: sha256(source / name) for name in SOURCE_FILES}, status="prepared")
        for index, path in enumerate(args.reference_input or []):
            destination = root / "reference_inputs" / f"{index}_{path.name}"
            destination.parent.mkdir(exist_ok=True)
            shutil.copy2(path, destination)
            manifest["reference_inputs"].append(dict(file=str(destination.relative_to(root)), sha256=sha256(destination)))
        write_json(root / "manifest.json", manifest)
        # All solving and reporting use the snapshot made before the first job.
        subprocess.run([args.python, "-u", str(source / "ReproducePaper.py"), "--resume", str(root),
                        "--python", args.python], check=True)
        return 0
    if not args.tables_only:
        preflight = subprocess.run([args.python, "-u", "-c", PREFLIGHT, str(source), manifest["seeds"]],
                                   capture_output=True, text=True)
        (root / "preflight.log").write_text(preflight.stdout + preflight.stderr)
        if preflight.returncode:
            raise RuntimeError("Preflight failed; inspect " + str(root / "preflight.log"))
        environment = json.loads(preflight.stdout)
        manifest.setdefault("sessions", []).append(dict(started=datetime.now().astimezone().isoformat(), **environment))
        manifest["status"] = "integer_running"
        write_json(root / "manifest.json", manifest)
        run_integer_batches(args.python, root, manifest)
        manifest["status"] = "lp_collecting"
        write_json(root / "manifest.json", manifest)
        for mode in ("leave", "continue"):
            command = lp_command(args.python, source, root, mode, manifest["settings"])
            print(f"Collecting {mode} LP bounds; retrying only missing cases.", flush=True)
            with (root / mode / "logs" / "lp.log").open("a") as handle:
                subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        from static_integrated_lp import fill_inputs
        from RunStaticLP import read_reference, fingerprint, MODEL_FILES
        source_hash = fingerprint({name: sha256(source/name) for name in MODEL_FILES})
        for mode in ('leave', 'continue'):
            completed = {r['problem_sha256']: r for r in read_reference(root/mode/'lp_results.csv').values()}
            fill_inputs((root/mode/'parts').glob('*.csv'), completed, source_hash,
                        dict(threads=manifest['settings']['lp_threads'], time_limit=manifest['settings']['lp_time_limit'],
                             retry_time_limit=manifest['settings']['lp_retry_time_limit']))
    build_tables(root, manifest)
    manifest.update(status="complete", finished=datetime.now().astimezone().isoformat())
    write_json(root / "manifest.json", manifest)
    print("Complete. CSV summaries and LaTeX fragments: " + str(root / "tables"), flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        print("Completed results are preserved. Resume with --resume RESULTS_DIRECTORY.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("Interrupted. Completed rows are preserved; use --resume RESULTS_DIRECTORY.", file=sys.stderr)
        sys.exit(130)
