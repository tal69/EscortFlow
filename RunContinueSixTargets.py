#!/usr/bin/env python3
"""Six-target continue campaign: 12/16/20 escorts; LF 300+0, EF 300+300 seconds."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
LAYOUTS = {"13x7": (6, 0), "10x10": (0, 0), "16x10": (4, 0, 11, 0), "27x10": (4, 0, 13, 0, 22, 0)}
ESCORTS = (12, 16, 20)
EXTENSIONS = {"loadflow": 0, "escortflow": 300}
PROTOCOL = "safe_integer_flow_timing_v6"
PREFLIGHT = """
import json, platform, sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
import numpy, gurobipy as gp
from RunSafeWeightedStatic import PROTOCOL
from static_integrated_lp import source_fingerprint
assert PROTOCOL == 'safe_integer_flow_timing_v6'
with gp.Env(empty=True) as env:
    env.setParam('OutputFlag',0)
    env.start()
    with gp.Model(env=env) as model:
        model.addVar(lb=0,obj=1)
        model.optimize()
        assert model.Status == gp.GRB.OPTIMAL
print(json.dumps(dict(python=sys.executable,python_version=sys.version,
    numpy_version=numpy.__version__,gurobi_version=gp.gurobi.version(),
    platform=platform.platform(),host=platform.node(),
    model_source_sha256=source_fingerprint(Path(sys.argv[1])))))
"""


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def seeds(text):
    result = []
    for token in text.split(","):
        if not re.fullmatch(r"[0-9]+(?:-[0-9]+)?", token):
            raise ValueError("Seeds must be comma-separated integers or ranges within 1-100")
        values = list(map(int, token.split("-")))
        first, last = values[0], values[-1]
        if not 1 <= first <= last <= 100:
            raise ValueError("Seeds must be within 1-100, in increasing ranges")
        result.extend(range(first, last + 1))
    if not result or len(result) != len(set(result)):
        raise ValueError("Empty or duplicate seed choices")
    return sorted(result)


def commands(output, layouts, seed_text):
    for layout in layouts:
        lx, ly = layout.split("x")
        for escort_count in ESCORTS:
            for method, extension in EXTENSIONS.items():
                label = f"{method}_{layout}_t6_e{escort_count}"
                yield label, [sys.executable, "-u", str(output / "source/SolveStatic.py"),
                    "--formulation", method, "-x", lx, "-y", ly, "-O",
                    *map(str, LAYOUTS[layout]), "-l", "6", "-e", str(escort_count),
                    "-r", seed_text, "-m", "continue", "--threads", "16",
                    "--weighted-time-limit", "300", "--extension-time-limit", str(extension),
                    "--with-lp", "--lp-threads", "1", "--lp-time-limit", "300",
                    "--lp-retry-time-limit", "600", "-f", str(output / "parts" / (label + ".csv"))]


def validate_source(source):
    env = json.loads((source.parent / "environment.json").read_text())
    if env["protocol"] != PROTOCOL or env["retrieval_mode"] != "continue":
        raise ValueError("Source must be an audited v6 continue-mode snapshot")
    for name, expected in env["source_sha256"].items():
        path = source / name
        if Path(name).name != name or path.is_symlink() or not path.is_file() or sha(path) != expected:
            raise ValueError("Frozen source integrity check failed: " + name)
    for name in ("SolveStatic.py", "RunSafeWeightedStatic.py",
                 "escort_flow_static_gurobi.py", "load_flow_static_gurobi.py"):
        if name not in env["source_sha256"]:
            raise ValueError("Missing solver source: " + name)
    return env


def freeze_source(source, output, env):
    destination = output / "source"
    destination.mkdir()
    for name, expected in env["source_sha256"].items():
        shutil.copy2(source / name, destination / name)
        if sha(destination / name) != expected:
            raise ValueError("Source changed while copying: " + name)
    for name in ("RunContinueSixTargets.py", "RunContinueSixTargets.sh", "RunContinueSixTargets.md"):
        shutil.copy2(HERE / name, destination / name)
    return destination


def active_study_processes(project):
    """Inspect exact study script arguments and working directories, without signals."""
    if not Path("/proc").exists():
        return []
    active = []
    driver_names = {"RunContinue.sh", "RunContinueSixTargets.sh", "RunContinueSixTargets.py"}
    worker_names = {"SolveStatic.py", "RunSafeWeightedStatic.py", "RunStaticLP.py"}
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) in {os.getpid(), os.getppid()}:
            continue
        try:
            args = [x.decode() for x in (proc / "cmdline").read_bytes().split(b"\0") if x]
            cwd = (proc / "cwd").resolve()
            match = False
            for arg in args:
                name = Path(arg).name
                if name in driver_names | worker_names:
                    path = (cwd / arg).resolve()
                    if path.is_relative_to(project):
                        match = True
            if match:
                active.append(dict(pid=int(proc.name), cwd=str(cwd), argv=args))
        except (OSError, UnicodeDecodeError):
            continue
    return active


def check_batch(path, expected_seeds, method):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError("Malformed/truncated output: " + str(path))
    if len(rows) != len(expected_seeds) or {int(r["seed"]) for r in rows} != set(expected_seeds):
        raise ValueError("Incomplete or duplicate seed coverage: " + str(path))
    for row in rows:
        if (row["error"] or row["has_solution"] != "1" or row["final_has_solution"] != "1"
                or row["protocol"] != PROTOCOL or row["retrieval_mode"] != "continue"
                or row["movement_mode"] != "BM" or row["#Loads"] != "6"
                or int(row["# Escorts"]) not in ESCORTS or row["formulation"] != method
                or float(row["weighted_time_limit"]) != 300
                or float(row["extension_time_limit"]) != EXTENSIONS[method]
                or row["phase1_snapshot_missing"] != "0"):
            raise ValueError("Failed or incompatible integer record: " + str(path))
    return sum(r["lp_status"] != "OPTIMAL" or bool(r["lp_error"]) for r in rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path, default=HERE / "results_continue_v6_six_targets",
                        help="new output directory; existing directories are refused")
    parser.add_argument("--dry-run", action="store_true", help="print commands without importing Gurobi or writing files")
    parser.add_argument("--source-dir", type=Path, default=HERE / "results_continue_v6/source",
                        help="audited frozen v6 source, with environment.json in its parent")
    parser.add_argument("--seeds", default="1-100", help="production: 1-100; subset for a separate pilot or interrupted batches")
    parser.add_argument("--layouts", nargs="+", choices=tuple(LAYOUTS), default=list(LAYOUTS))
    args = parser.parse_args(argv)
    selected_seeds = seeds(args.seeds)
    if len(args.layouts) != len(set(args.layouts)):
        parser.error("Repeated layouts")
    output, source = args.output.resolve(), args.source_dir.resolve()
    plan = list(commands(output, args.layouts, args.seeds))
    if args.dry_run:
        for _, command in plan:
            print(shlex.join(command))
        return 0
    if output.exists():
        parser.error("Output directory already exists: " + str(output))
    project = HERE.parent.parent if HERE.name == "source" else HERE
    active = active_study_processes(project)
    if active:
        parser.error("Existing study processes are still active; this launcher does not stop them:\n"
                     + json.dumps(active, indent=2))
    source_env = validate_source(source)
    output.mkdir(parents=True, exist_ok=False)
    state = dict(started_utc=now(), status="preflight", completed_batches=[], pending_lp_records=0)

    def save_state():
        state["updated_utc"] = now()
        (output / "run_state.json").write_text(json.dumps(state, indent=2) + "\n")

    save_state()
    try:
        runtime = freeze_source(source, output, source_env)
        (output / "parts").mkdir()
        (output / "logs").mkdir()
        check = subprocess.run([sys.executable, "-c", PREFLIGHT, str(runtime)], capture_output=True, text=True)
        (output / "preflight.log").write_text(check.stdout + check.stderr)
        check.check_returncode()
        software = json.loads(check.stdout.splitlines()[-1])
        metadata = dict(**software, started_utc=state["started_utc"], protocol=PROTOCOL,
            campaign="continue_six_targets_12_16_20_lf300_ef600",
            source_snapshot=str(source), original_source_sha256=source_env["source_sha256"],
            seed_values=selected_seeds, layouts=args.layouts, targets=[6], escorts=list(ESCORTS),
            retrieval_mode="continue", movement_mode="BM", threads=16, initial_time_limit=300,
            extension_time_limit_by_method=EXTENSIONS, with_lp=True, lp_threads=1,
            lp_time_limit=300, lp_retry_time_limit=600, configurations=len(plan)//2,
            matched_instances=len(plan)//2*len(selected_seeds), integer_runs=len(plan)*len(selected_seeds),
            method_order=list(EXTENSIONS), logical_cpus=os.cpu_count(),
            source_sha256={p.name:sha(p) for p in sorted(runtime.iterdir()) if p.is_file()})
        if Path("/proc/cpuinfo").exists():
            metadata.update(cpu_info=Path("/proc/cpuinfo").read_text(),
                            memory_info=Path("/proc/meminfo").read_text())
        (output / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n")
        (output / "commands.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n"
                                            + "\n".join(shlex.join(c) for _, c in plan) + "\n")
        state["status"] = "running"
        save_state()
        for label, command in plan:
            state["current_batch"] = label
            save_state()
            log = output / "logs" / (label + ".log")
            print(f"Starting {label}; log: {log}", flush=True)
            with log.open("w") as handle:
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
            if result.returncode not in (0, 3):
                raise RuntimeError(f"{label} exited {result.returncode}; inspect {log}")
            method = command[command.index("--formulation") + 1]
            pending = check_batch(Path(command[-1]), selected_seeds, method)
            state["pending_lp_records"] += pending
            state["completed_batches"].append(dict(batch=label, exit_code=result.returncode,
                                                    pending_lp_records=pending, completed_utc=now()))
            save_state()
        state.pop("current_batch", None)
        state["status"] = "integer_batches_complete_lp_pending" if state["pending_lp_records"] else "complete"
        state["finished_utc"] = now()
        save_state()
        print(f"{state['status']}: {output}", flush=True)
        return 0
    except BaseException as exc:
        state["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        state["error"] = str(exc) or type(exc).__name__
        save_state()
        raise


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
