#!/usr/bin/env bash
# Table 2(a,b), weighted objective with greedy starts and separate flow proofs.
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON:-python3}
threads=${NUM_THREADS:-12}
seeds=1-100
part=both
weighted_limit=300
certification_limit=300
certification_gap=0.1
dry_run=0
result_dir="${script_dir}/results_table2_weighted_$(date +%Y%m%d_%H%M%S)_$$"

usage() {
    cat <<'USAGE'
Usage: bash RunTable2Weighted.sh [options]

Run weighted escort-flow and load-flow Table 2(a,b) experiments with the same
greedy warm starts. Each weighted solve gets 300 solver seconds. Eligible
solutions get a separate flow-time certification solve, capped at 300 seconds
and warm-started from the weighted incumbent.
Certification runs if weighted optimality is proven or the absolute weighted
gap is strictly below 0.1 in the original F + 0.01 M units.

Options:
  --python PATH                   Python with numpy, gurobipy, and Gurobi license
  --threads N                     Threads per solve (default: $NUM_THREADS or 12)
  --seeds RANGE                   Seed range (default: 1-100)
  --part a|b|both                  Table part (default: both)
  --weighted-time-limit SECONDS   Weighted solver budget (default: 300)
  --certification-time-limit SEC  Separate certification budget (default: 300)
  --certification-gap-threshold G Absolute weighted gap gate, 0 < G < 1 (default: 0.1)
  --output-dir DIR                New results directory; existing paths rejected
  --dry-run                       Print commands without solving or creating files
  -h, --help                      Show this help

Examples:
  bash RunTable2Weighted.sh --python /opt/conda/envs/pbs/bin/python --threads 12
  nohup bash RunTable2Weighted.sh > table2_weighted.log 2>&1 &
  bash RunTable2Weighted.sh --dry-run

All batches run sequentially: four layouts, both formulations, leave retrieval,
simultaneous block movements. Part (a): one target and 3-8 escorts. Part (b):
four targets and 8/12/16 escorts. Each per-instance log reports stage transitions.
CSV weighted results remain unchanged if certification finds a better flow time.
USAGE
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

while (($#)); do
    case $1 in
        --python|--threads|--seeds|--part|--output-dir|--weighted-time-limit|--certification-time-limit|--certification-gap-threshold)
            (($# >= 2)) || die "Missing value for $1"
            case $1 in
                --python) python_bin=$2 ;;
                --threads) threads=$2 ;;
                --seeds) seeds=$2 ;;
                --part) part=$2 ;;
                --output-dir) result_dir=$2 ;;
                --weighted-time-limit) weighted_limit=$2 ;;
                --certification-time-limit) certification_limit=$2 ;;
                --certification-gap-threshold) certification_gap=$2 ;;
            esac
            shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done
[[ $part == a || $part == b || $part == both ]] || die "--part must be a, b, or both"
[[ $threads =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
[[ -n $seeds ]] || die "--seeds must not be empty"

if (( ! dry_run )); then
    command -v "$python_bin" >/dev/null 2>&1 || die "Python executable not found: $python_bin"
    [[ ! -e $result_dir ]] || die "Output directory already exists: $result_dir"
    mkdir -p -- "$result_dir/parts" "$result_dir/logs"
    result_dir=$(cd -- "$result_dir" && pwd)
    if ! "$python_bin" - "$script_dir" "$threads" "$seeds" "$part" "$weighted_limit" "$certification_limit" "$certification_gap" <<'PY' > "$result_dir/environment.txt" 2>&1
import datetime
import hashlib
import os
from pathlib import Path
import platform
import subprocess
import sys

print("Started:", datetime.datetime.now().astimezone().isoformat(), flush=True)
print("Host:", platform.node(), flush=True)
print("Platform:", platform.platform(), flush=True)
print("Processor:", platform.processor(), flush=True)
print("Logical CPUs:", os.cpu_count(), flush=True)
try:
    print("Physical memory GiB:", os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 2**30, flush=True)
except (ValueError, OSError):
    print("Physical memory GiB: unavailable", flush=True)
print("Python:", sys.version.replace("\n", " "), flush=True)
print("Executable:", sys.executable, flush=True)
if sys.version_info < (3, 10):
    raise SystemExit("Python 3.10 or newer is required; 3.11 is recommended")
import numpy
import gurobipy as gp
sys.path.insert(0, sys.argv[1])
from RunWeightedStatic import parse_range, positive_number, certification_threshold
seed_values = parse_range(sys.argv[3], minimum=0)
weighted = positive_number(sys.argv[5])
certificate = positive_number(sys.argv[6])
threshold = certification_threshold(sys.argv[7])
print("NumPy:", numpy.__version__)
print("Gurobi:", ".".join(map(str, gp.gurobi.version())))
print("Threads:", sys.argv[2])
print("Seeds:", sys.argv[3], "count:", len(seed_values))
print("Table part:", sys.argv[4])
print(f"Weighted cap: {weighted:g} seconds; independent certification cap: {certificate:g} seconds")
print(f"Certification: weighted proven or absolute weighted gap < {threshold:g} (F+0.01M units)")
print("Objective: F+0.01M, implemented as integer 100F+M; MIPGap=0; MIPGapAbs=0.999 in scaled units")
print("Mode: leave; movement: BM; sequential jobs")
print("Warm starts: common complete greedy plan for weighted solves; saved weighted incumbent for certification")
print("MIPFocus: weighted 0; certification 3. Timings exclude the other solve's budget.")
source = Path(sys.argv[1])
for arguments in (["rev-parse", "HEAD"], ["status", "--short"]):
    try:
        result = subprocess.run(["git", "-C", str(source), *arguments], capture_output=True, text=True, check=False)
        print("Git", " ".join(arguments) + ":", result.stdout.strip() or result.stderr.strip())
    except OSError as exc:
        print("Git metadata unavailable:", exc)
for name in ("RunTable2Weighted.sh", "RunWeightedStatic.py", "static_weighted_certification.py",
             "static_lexicographic.py", "escort_flow_static_gurobi.py", "load_flow_static_gurobi.py",
             "OneStepHeuristic_v2.py", "PBSCom.py"):
    print("SHA256", name, hashlib.sha256((source / name).read_bytes()).hexdigest())
with gp.Env(empty=True) as env:
    env.setParam("OutputFlag", 0)
    env.start()
    with gp.Model(env=env) as model:
        model.addVar(lb=0, obj=1)
        model.optimize()
        if model.Status != gp.GRB.OPTIMAL:
            raise SystemExit("Gurobi license/solver check failed")
print("Preflight passed.")
PY
    then
        cat "$result_dir/environment.txt" >&2
        die "Preflight failed; inspect $result_dir/environment.txt"
    fi
    cat "$result_dir/environment.txt"
    printf '#!/usr/bin/env bash\nset -euo pipefail\n' > "$result_dir/commands.sh"
fi

run_layout() {
    local table_part=$1 formulation=$2 lx=$3 ly=$4
    shift 4
    local loads escorts batch_csv merged_csv log_file
    if [[ $table_part == a ]]; then
        loads=1; escorts=3-8
    else
        loads=4; escorts=8-16-4
    fi
    batch_csv="$result_dir/parts/table2${table_part}_${formulation}_${lx}x${ly}.csv"
    merged_csv="$result_dir/table2${table_part}_${formulation}_weighted.csv"
    log_file="$result_dir/logs/table2${table_part}_${formulation}_${lx}x${ly}.log"
    local command=("$python_bin" -u "$script_dir/RunWeightedStatic.py"
        --formulation "$formulation" -x "$lx" -y "$ly" -O "$@"
        -e "$escorts" -l "$loads" -r "$seeds" --threads "$threads"
        --weighted-time-limit "$weighted_limit" --certification-time-limit "$certification_limit"
        --certification-gap-threshold "$certification_gap" -f "$batch_csv")
    if (( dry_run )); then
        printf '%q ' "${command[@]}"
        printf '\n'
        return
    fi
    printf '[%s] Table 2(%s), %s, %sx%s: starting (log: %s)\n' \
        "$(date '+%F %T')" "$table_part" "$formulation" "$lx" "$ly" "$log_file"
    printf '%q ' "${command[@]}" >> "$result_dir/commands.sh"
    printf '\n' >> "$result_dir/commands.sh"
    if ! "${command[@]}" > "$log_file" 2>&1; then
        die "Solver process failed; inspect $log_file"
    fi
    "$python_bin" - "$script_dir" "$batch_csv" "$merged_csv" "$seeds" "$escorts" "$weighted_limit" "$certification_limit" "$certification_gap" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from RunWeightedStatic import merge_batch
merge_batch(*sys.argv[2:])
PY
}

for table_part in a b; do
    [[ $part == both || $part == "$table_part" ]] || continue
    for formulation in escortflow loadflow; do
        run_layout "$table_part" "$formulation" 13 7 6 0
        run_layout "$table_part" "$formulation" 10 10 0 0
        run_layout "$table_part" "$formulation" 16 10 4 0 11 0
        run_layout "$table_part" "$formulation" 27 10 4 0 13 0 22 0
    done
done
if (( ! dry_run )); then
    printf '[%s] Completed Table 2 weighted runs. Results: %s\n' "$(date '+%F %T')" "$result_dir"
fi
