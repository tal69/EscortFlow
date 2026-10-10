#!/usr/bin/env python3
"""Shared launcher for matched static SBM experiments on Linux or macOS.

Run with a Python environment containing NumPy, gurobipy, and a Gurobi license,
or select that environment with --python. No experiment is launched by --dry-run.
"""

import argparse
import ast
import csv
from datetime import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys


LAYOUTS = (
    (13, 7, 27, (6, 0)),
    (10, 10, 30, (0, 0)),
    (16, 10, 48, (4, 0, 11, 0)),
    (27, 10, 81, (4, 0, 13, 0, 22, 0)),
)
FORMULATIONS = ("escortflow", "loadflow")
CAMPAIGNS = {
    "occupancy70": dict(prefix="occupancy70", directory="70percent", mode="leave", loads=(4,),
                        description="Four-target leave-mode cases at approximately 70% occupancy."),
    "continue": dict(prefix="continue", directory="continue", mode="continue", loads=(2, 4, 6),
                     description="Continue-mode cases with 2, 4 and 6 targets, at every Table 2(b) escort count including approximately 70% occupancy."),
    "target_counts": dict(prefix="table2b_targets", directory="table2b_targets", mode="leave", loads=(2, 6),
                         escorts_by_loads={2: (8, 12, 16), 6: (8, 12, 16, 20)}, lp=True,
                         description="Tables 2 and 3 leave-mode cases: 2 targets with 8/12/16 escorts, "
                                     "6 targets with 8/12/16/20 escorts, with an LP bound in each instance row."),
}


def configurations(campaign, layouts):
    for lx, ly, occupancy_escorts, outputs in LAYOUTS:
        if "{}x{}".format(lx, ly) not in layouts:
            continue
        for loads in CAMPAIGNS[campaign]["loads"]:
            escort_counts = CAMPAIGNS[campaign].get("escorts_by_loads", {}).get(loads)
            if escort_counts is None:
                escort_counts = ((occupancy_escorts,) if campaign == "occupancy70"
                                 else (8, 12, 16, occupancy_escorts))
            for escorts in escort_counts:
                yield dict(lx=lx, ly=ly, escorts=escorts, loads=loads, outputs=outputs,
                           occupancy=(lx * ly - escorts) / (lx * ly))

PREFLIGHT = r'''
import json, os, platform, subprocess, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
print("Preflight interpreter: {} (Python {})".format(sys.executable, sys.version.split()[0]),
      file=sys.stderr, flush=True)
if sys.version_info < (3, 10):
    raise SystemExit("The solver environment requires Python 3.10 or newer. "
                     "Select the intended interpreter with --python /path/to/python.")
import numpy
import gurobipy as gp
from RunSafeWeightedStatic import PROTOCOL, parse_range
seeds = parse_range(sys.argv[2], minimum=0)
hardware = {}
if platform.system() == "Darwin":
    for key in ("hw.model", "machdep.cpu.brand_string", "hw.memsize",
                "hw.physicalcpu", "hw.logicalcpu", "hw.perflevel0.physicalcpu",
                "hw.perflevel1.physicalcpu"):
        result = subprocess.run(["/usr/sbin/sysctl", "-n", key],
                                capture_output=True, text=True)
        if result.returncode == 0:
            hardware[key] = result.stdout.strip()
    memory = int(hardware["hw.memsize"]) / 2**30 if "hw.memsize" in hardware else None
else:
    try:
        cpuinfo = Path("/proc/cpuinfo").read_text()
        processors = [dict((key.strip(), value.strip()) for key, value in
                           (line.split(":", 1) for line in block.splitlines() if ":" in line))
                      for block in cpuinfo.strip().split("\n\n")]
        models = sorted({p["model name"] for p in processors if "model name" in p})
        if models:
            hardware["cpu_models"] = models
        cores = {(p["physical id"], p["core id"]) for p in processors
                 if "physical id" in p and "core id" in p}
        if cores:
            hardware["physical_cores"] = len(cores)
    except OSError:
        pass
    try:
        memory = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 2**30
    except (ValueError, OSError):
        memory = None
with gp.Env(empty=True) as env:
    env.setParam("OutputFlag", 0)
    env.start()
    with gp.Model(env=env) as model:
        model.addVar(lb=0, obj=1)
        model.optimize()
        if model.Status != gp.GRB.OPTIMAL:
            raise SystemExit("Gurobi license/solver preflight failed")
print(json.dumps(dict(host=platform.node(), platform=platform.platform(),
    processor=platform.processor(), logical_cpus=os.cpu_count(), hardware=hardware,
    physical_memory_gib=memory, python=sys.version, executable=sys.executable,
    numpy=numpy.__version__, gurobipy=getattr(gp, "__version__", None),
    gurobi=".".join(map(str, gp.gurobi.version())),
    protocol=PROTOCOL, seed_values=seeds)))
'''

MERGE = r'''
import ast, csv, sys
sys.path.insert(0, sys.argv[1])
from RunSafeWeightedStatic import merge_batch
source, destination, seeds, escorts, cutoff, extension = sys.argv[2:8]
formulation, grid, outputs, loads, retrieval_mode = sys.argv[8:13]
with open(source, newline="") as handle:
    for row in csv.DictReader(handle):
        if (row["formulation"] != formulation or row["Lx x Ly"] != grid
                or row["#Loads"] != loads or row["# Escorts"] != escorts
                or sorted(ast.literal_eval(row["IOs"])) != sorted(ast.literal_eval(outputs))
                or row["retrieval_mode"] != retrieval_mode or row["movement_mode"] != "BM"):
            raise SystemExit("Unexpected experiment configuration in " + source)
merge_batch(source, destination, seeds, escorts, cutoff, extension)
'''

FILL_PARTS = r'''
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from static_integrated_lp import fill_inputs, source_fingerprint
from RunStaticLP import read_reference
completed = {r['problem_sha256']: r for r in read_reference(Path(sys.argv[2])).values()}
fill_inputs(sys.argv[6:], completed, source_fingerprint(sys.argv[1]),
            dict(threads=int(sys.argv[3]), time_limit=float(sys.argv[4]), retry_time_limit=float(sys.argv[5])))
'''


def automatic_threads():
    """Use Apple performance cores, or the Linux campaign's 16-thread setting."""
    if platform.system() == "Darwin":
        try:
            result = subprocess.run(
                ["/usr/sbin/sysctl", "-n", "hw.perflevel0.physicalcpu"],
                capture_output=True, text=True, check=True,
            )
            count = int(result.stdout.strip())
            if count > 0:
                return count
        except (OSError, ValueError, subprocess.CalledProcessError):
            pass
    return min(16, os.cpu_count() or 16)


def finite_number(value):
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("Expected a finite nonnegative number")
    return number


def commands(args, source_dir, result_dir):
    settings = CAMPAIGNS[args.campaign]
    for config in configurations(args.campaign, args.layouts):
        lx, ly, escorts, outputs, loads = (config[k] for k in ("lx", "ly", "escorts", "outputs", "loads"))
        for formulation in FORMULATIONS:
            stem = "{}_{}_{}x{}".format(settings["prefix"], formulation, lx, ly)
            if args.campaign != "occupancy70":
                stem += "_l{}_e{}".format(loads, escorts)
            command = [args.python, "-u", str(source_dir / "RunSafeWeightedStatic.py"),
                       "--formulation", formulation, "-x", str(lx), "-y", str(ly),
                       "-O", *map(str, outputs), "-e", str(escorts), "-l", str(loads),
                       "--retrieval-mode", settings["mode"],
                       "-r", args.seeds, "--threads", str(args.threads),
                       "--weighted-time-limit", str(args.weighted_time_limit),
                       "--extension-time-limit", str(args.extension_time_limit),
                       "-f", str(result_dir / "parts" / (stem + ".csv"))]
            if args.lp:
                command += ['--with-lp', '--lp-threads', str(args.lp_threads),
                            '--lp-time-limit', str(args.lp_time_limit),
                            '--lp-retry-time-limit', str(args.lp_retry_time_limit)]
            yield dict(config, formulation=formulation, stem=stem, command=command)


def lp_command(args, source_dir, result_dir):
    """Replay both merged integer CSVs with their recorded R and physical horizon."""
    prefix = CAMPAIGNS[args.campaign]["prefix"]
    return [args.python, "-u", str(source_dir / "RunStaticLP.py"),
            "--input", *(str(result_dir / (prefix + "_" + method + ".csv"))
                         for method in FORMULATIONS),
            "--source-dir", str(source_dir), "--workers", str(args.lp_workers),
            "--threads", str(args.lp_threads), "--time-limit", str(args.lp_time_limit),
            "--retry-time-limit", str(args.lp_retry_time_limit),
            "-f", str(result_dir / "lp_results.csv"), '--fill-input-lp']


def validate_pairs(result_dir, seed_values, configs, prefix):
    """Require identical instances and settings across the two formulations."""
    shared = ("protocol", "threads", "warmstart", "retrieval_mode", "movement_mode",
              "flow_weight", "movement_integer_weight", "safe_movement_bound",
              "safe_movement_lower_bound", "safe_flow_horizon",
              "weighted_physical_horizon", "weighted_global_scope",
              "greedy_flowtime", "greedy_movements", "greedy_makespan",
              "weighted_time_limit", "extension_time_limit", "flow_proof_check_mode")
    tables = {}
    for formulation in FORMULATIONS:
        with (result_dir / (prefix + "_" + formulation + ".csv")).open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        keyed = {}
        for row in rows:
            key = (row["Lx x Ly"], row["# Escorts"], row["#Loads"], row["seed"],
                   *(tuple(sorted(ast.literal_eval(row[field])))
                     for field in ("IOs", "Escorts", "Target Loads")))
            if key in keyed:
                raise ValueError("Duplicate paired instance: " + repr(key[:4]))
            keyed[key] = row
        expected = {(str(c["lx"]) + "x" + str(c["ly"]), str(c["escorts"]), str(c["loads"]), str(seed))
                    for c in configs for seed in seed_values}
        actual = [key[:4] for key in keyed]
        if len(actual) != len(expected) or set(actual) != expected:
            raise ValueError("Missing or unexpected configurations for " + formulation)
        tables[formulation] = keyed
    ef, lf = (tables[name] for name in FORMULATIONS)
    if ef.keys() != lf.keys():
        raise ValueError("The two formulations have different initial instances")
    for key in ef:
        for field in shared:
            if ef[key][field] != lf[key][field]:
                raise ValueError("Different paired {} at {}".format(field, key[:4]))
    report = dict(matched_instances=len(ef), solver_runs=2 * len(ef),
                  shared_fields_checked=list(shared), passed=True)
    (result_dir / "pairing.json").write_text(json.dumps(report, indent=2) + "\n")
    print("Verified {} matched instances and common solver settings.".format(len(ef)), flush=True)


def main(campaign, argv=None):
    settings = CAMPAIGNS[campaign]
    parser = argparse.ArgumentParser(description=settings["description"] + " " + __doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--python", default=os.environ.get("PYTHON") or sys.executable,
                        help="solver Python executable; defaults to PYTHON when set, "
                             "otherwise the interpreter running this launcher")
    parser.add_argument("--threads", type=int,
                        default=os.environ.get("NUM_THREADS") or automatic_threads(),
                        help="threads per solve; auto default uses Mac performance cores or up to 16 on Linux")
    parser.add_argument("--seeds", default="1-100", help="initial-state seed range")
    parser.add_argument("--layouts", nargs="+", choices=["{}x{}".format(x, y) for x, y, _, _ in LAYOUTS],
                        default=["{}x{}".format(x, y) for x, y, _, _ in LAYOUTS],
                        help="layouts to include; default is all four")
    parser.add_argument("--weighted-time-limit", type=finite_number, default=300,
                        help="main comparison cutoff in solver seconds")
    parser.add_argument("--extension-time-limit", type=finite_number, default=300,
                        help="extra solver seconds only when cutoff flow is unproved")
    parser.add_argument("--lp", action=argparse.BooleanOptionalAction,
                        default=settings.get("lp", False),
                        help="save the matched LP bound after each integer instance (default on for 2/6 targets)")
    parser.add_argument("--lp-workers", type=int, default=1,
                        help="parallel LP recovery processes for missing bounds; integrated LPs run sequentially")
    parser.add_argument("--lp-threads", type=int, default=1, help="solver threads per LP process")
    parser.add_argument("--lp-time-limit", type=finite_number, default=300,
                        help="initial LP solver seconds per instance")
    parser.add_argument("--lp-retry-time-limit", type=finite_number, default=600,
                        help="barrier retry seconds after an LP time limit")
    parser.add_argument("--output-dir", type=Path, help="new results directory")
    parser.add_argument("--dry-run", action="store_true",
                        help="print every batch command without creating files or solving")
    args = parser.parse_args(argv)
    args.campaign = campaign
    if args.threads <= 0 or args.weighted_time_limit <= 0:
        parser.error("Threads and the main time limit must be positive")
    if (args.lp_workers <= 0 or args.lp_threads <= 0
            or args.lp_time_limit <= 0 or args.lp_retry_time_limit <= 0):
        parser.error("LP workers, threads, and time limits must be positive")
    if len(set(args.layouts)) != len(args.layouts):
        parser.error("Each layout may be selected only once")
    source_dir = Path(__file__).resolve().parent
    stamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    result_dir = (args.output_dir or source_dir / (
        "results_{}_{}_{}".format(settings["directory"], stamp, os.getpid()))).resolve()
    if args.dry_run:
        for batch in commands(args, source_dir, result_dir):
            print(shlex.join(batch["command"]))
        if args.lp:
            print(shlex.join(lp_command(args, source_dir, result_dir)))
        return 0
    python_path = shutil.which(args.python)
    if python_path is None:
        parser.error("Python executable not found: " + args.python)
    args.python = python_path
    if result_dir.exists():
        parser.error("Output directory already exists: " + str(result_dir))
    print("Launcher interpreter: {} (Python {})".format(
        sys.executable, platform.python_version()), flush=True)
    print("Selected solver interpreter: {}".format(args.python), flush=True)
    result_dir.mkdir(parents=True, exist_ok=False)
    for name in ("parts", "logs", "source"):
        (result_dir / name).mkdir()
    # Freeze local Python modules so edits during the campaign cannot change a
    # later batch's formulation, coefficient, warm start, or reporting protocol.
    runtime_dir = result_dir / "source"
    for path in source_dir.glob("*.py"):
        if path.is_file():
            shutil.copy2(path, runtime_dir / path.name)
    preflight = subprocess.run([args.python, "-u", "-c", PREFLIGHT,
                                str(runtime_dir), args.seeds], capture_output=True, text=True)
    (result_dir / "preflight.log").write_text(preflight.stdout + preflight.stderr)
    if preflight.returncode:
        print(preflight.stdout + preflight.stderr, file=sys.stderr)
        raise RuntimeError("Preflight failed; inspect " + str(result_dir / "preflight.log"))
    metadata = json.loads(preflight.stdout)
    configs = list(configurations(campaign, args.layouts))
    metadata.update(started=datetime.now().astimezone().isoformat(), campaign=campaign,
                    threads=args.threads, weighted_time_limit=args.weighted_time_limit,
                    extension_time_limit=args.extension_time_limit,
                    objective="R*F+M with the sufficient coefficient recorded in each row",
                    comparison="Saved incumbent and bound by the initial cutoff; extensions separate",
                    retrieval_mode=settings["mode"], movement_mode="BM", target_loads=list(settings["loads"]),
                    warm_start="Common complete greedy plan in both formulations",
                    solver_runs=2 * len(configs) * len(metadata["seed_values"]),
                    lp=dict(enabled=args.lp, schedule='after_each_integer_instance', workers=1,
                            recovery_workers=args.lp_workers, threads=args.lp_threads,
                            time_limit=args.lp_time_limit, retry_time_limit=args.lp_retry_time_limit,
                            objective="F+M/R, using each recorded R and horizon",
                            solver_runs=2 * len(configs) * len(metadata["seed_values"]) if args.lp else 0),
                    configurations=[dict(c, outputs=list(zip(c["outputs"][::2], c["outputs"][1::2])))
                                    for c in configs])
    for option in (["rev-parse", "HEAD"], ["status", "--short"]):
        try:
            result = subprocess.run(["git", "-C", str(source_dir), *option],
                                    capture_output=True, text=True, timeout=10)
            metadata["git_" + option[0]] = result.stdout.strip() or result.stderr.strip()
        except (OSError, subprocess.TimeoutExpired) as exc:
            metadata["git_" + option[0]] = "Unavailable: " + str(exc)
    metadata["source_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted(runtime_dir.glob("*.py"))}
    (result_dir / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("Host: {}; protocol: {}; threads: {}".format(
        metadata["host"], metadata["protocol"], args.threads), flush=True)
    print("Python: {}; Gurobi: {}; memory GiB: {}".format(
        metadata["executable"], metadata["gurobi"], metadata["physical_memory_gib"]), flush=True)
    print("{} {} solver runs, sequential batches".format(
        settings["description"], metadata["solver_runs"]), flush=True)
    print("First-phase cap: {}s; conditional extension: {}s".format(
        args.weighted_time_limit, args.extension_time_limit), flush=True)
    print("Results: " + str(result_dir), flush=True)
    batches = list(commands(args, runtime_dir, result_dir))
    (result_dir / "commands.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" +
        "\n".join(shlex.join(batch["command"]) for batch in batches) + "\n")
    if args.lp:
        relaxation_command = lp_command(args, runtime_dir, result_dir)
        fill_command = [args.python, '-u', '-c', FILL_PARTS, str(runtime_dir), str(result_dir/'lp_results.csv'),
                        str(args.lp_threads), str(args.lp_time_limit), str(args.lp_retry_time_limit),
                        *(str(result_dir/'parts'/(batch['stem']+'.csv')) for batch in batches)]
        # This separate script also retries only missing LP values after a failure
        # or interruption, using the frozen solver modules and input fingerprints.
        (result_dir / "lp_commands.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" +
            shlex.join(relaxation_command + ["--resume"]) + "\n" + shlex.join(fill_command) + "\n")
    for batch in batches:
        lx, ly, escorts, outputs, loads, formulation, stem, command = (
            batch[k] for k in ("lx", "ly", "escorts", "outputs", "loads", "formulation", "stem", "command"))
        log = result_dir / "logs" / (stem + ".log")
        print("[{}] {} {}x{}, {} targets, {} escorts, {}: starting (log: {})".format(
            datetime.now().astimezone().isoformat(), formulation, lx, ly, loads, escorts, settings["mode"], log), flush=True)
        with log.open("w") as handle:
            result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode not in (0, 3):
            result.check_returncode()
        if result.returncode == 3:
            print('Integer rows saved; missing LP bounds will be retried after the campaign.', flush=True)
        batch_csv = result_dir / "parts" / (stem + ".csv")
        merged_csv = result_dir / (settings["prefix"] + "_" + formulation + ".csv")
        subprocess.run([args.python, "-u", "-c", MERGE, str(runtime_dir), str(batch_csv),
                        str(merged_csv), args.seeds, str(escorts), str(args.weighted_time_limit),
                        str(args.extension_time_limit), formulation, "{}x{}".format(lx, ly),
                        repr(list(zip(outputs[::2], outputs[1::2]))), str(loads), settings["mode"]], check=True)
        print("[{}] {} {}x{}: completed".format(
            datetime.now().astimezone().isoformat(), formulation, lx, ly), flush=True)
    validate_pairs(result_dir, metadata["seed_values"], configs, settings["prefix"])
    if args.lp:
        log = result_dir / "logs" / (settings["prefix"] + "_lp.log")
        print("Collecting saved LP bounds; retrying only missing cases (log: {})".format(log), flush=True)
        with log.open("w") as handle:
            subprocess.run(relaxation_command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        subprocess.run(fill_command, check=True)
        print("LP relaxations complete. Saved OPTIMAL values: " + str(result_dir / "lp_results.csv"), flush=True)
    print("Completed. Results: " + str(result_dir), flush=True)
    return 0


def run_campaign(campaign):
    try:
        sys.exit(main(campaign))
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("Interrupted. Already written CSV rows and logs are preserved.", file=sys.stderr)
        sys.exit(130)
