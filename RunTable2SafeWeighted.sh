#!/usr/bin/env bash
# Table 2(a,b), sufficient integer weights and first flow-proof timing.
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON:-python3}
threads=${NUM_THREADS:-16}
seeds=1-100
part=both
weighted_limit=300
extension_limit=300
stop_at_flow_proof=0
dry_run=0
result_dir="${script_dir}/results_table2_safe_weighted_$(date +%Y%m%d_%H%M%S)_$$"

usage() {
    cat <<'USAGE'
Usage: bash RunTable2SafeWeighted.sh [options]

Run escort-flow and load-flow Table 2(a,b) with the same greedy warm starts.
Each instance uses integer R*F+M, with R=(N-e)*H-D+1, where D=sum(d_i) is the
nearest-output distance lower bound and H is a sufficient horizon derived
from its greedy feasible flow time. The first phase attempts both
objectives for up to 300 solver seconds. A gap below 1 proves global
lexicographic optimality. Check flow optimality when the lower bound improves
or a smaller candidate flow is found. Record the first observed proof time.
If flow remains unproved at 300 seconds, extend the same search for up to another
300 seconds, stopping on its flow proof or a lower-flow counterexample.
Preserve the 300-second result and report the final result separately.
One optimize call preserves the search tree, incumbent, cuts, and objective.

Options:
  --python PATH                   Python with numpy, gurobipy, and Gurobi license
  --threads N                     Threads per solve (default: $NUM_THREADS or 16)
  --seeds RANGE                   Seed range (default: 1-100)
  --part a|b|both                  Table part (default: both)
  --weighted-time-limit SECONDS   Initial reporting cutoff (default: 300)
  --extension-time-limit SECONDS  Conditional extra search budget (default: 300)
  --certification-time-limit SEC  Alias for --extension-time-limit
  --stop-at-flow-proof            Stop as soon as flow time is proved; retain all KPIs
  --output-dir DIR                New results directory; existing paths rejected
  --dry-run                       Print commands without solving or creating files
  -h, --help                      Show this help

Examples:
  bash RunTable2SafeWeighted.sh --python /opt/conda/envs/pbs/bin/python --threads 16
  bash RunTable2SafeWeighted.sh --dry-run

All batches run sequentially: four layouts, both formulations, leave retrieval,
simultaneous block movements. Part (a): one target and 3-8 escorts. Part (b):
four targets and 8/12/16 escorts. CSVs report the first flow-proof time and total
solve time. The cutoff and final solutions are saved separately.
USAGE
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

while (($#)); do
    case $1 in
        --python|--threads|--seeds|--part|--output-dir|--weighted-time-limit|--extension-time-limit|--certification-time-limit)
            (($# >= 2)) || die "Missing value for $1"
            case $1 in
                --python) python_bin=$2 ;;
                --threads) threads=$2 ;;
                --seeds) seeds=$2 ;;
                --part) part=$2 ;;
                --output-dir) result_dir=$2 ;;
                --weighted-time-limit) weighted_limit=$2 ;;
                --extension-time-limit|--certification-time-limit) extension_limit=$2 ;;
            esac
            shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        --stop-at-flow-proof) stop_at_flow_proof=1; shift ;;
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
    if ! "$python_bin" - "$script_dir" "$threads" "$seeds" "$part" "$weighted_limit" "$extension_limit" "$stop_at_flow_proof" <<'PY' > "$result_dir/environment.txt" 2>&1
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
from RunSafeWeightedStatic import PROTOCOL, parse_range, positive_number, nonnegative_number
seed_values = parse_range(sys.argv[3], minimum=0)
weighted = positive_number(sys.argv[5])
extension = nonnegative_number(sys.argv[6])
stop_at_flow_proof = bool(int(sys.argv[7]))
print("NumPy:", numpy.__version__)
print("Gurobi:", ".".join(map(str, gp.gurobi.version())))
print("Threads:", sys.argv[2])
print("Seeds:", sys.argv[3], "count:", len(seed_values))
print("Table part:", sys.argv[4])
print(f"First-phase cutoff: {weighted:g} seconds; conditional extension: {extension:g} seconds; total cap: {weighted+extension:g} seconds")
print("Protocol:", PROTOCOL)
print("Objective: integer R*F+M; R=(N-e)*H-D+1 per instance; D=sum(d_i); MIPGap=0; MIPGapAbs=0.999")
print("Movement bounds: lower D=sum(d_i); upper U=(N-e)*H; coefficient R=U-D+1")
print("Horizon: H=greedy_F-sum(d_i)+max(d_i), covering a global lexicographic optimum")
print("Flow proof checks: improved lower bound or smaller candidate flow; mandatory cutoff and final checks")
print("Cached flow criterion; unchanged bounds and movement-only improvements skip proof comparisons")
print("stop_at_flow_proof:", stop_at_flow_proof)
print("Record the first observed flow proof; retain flow, movements, bound, gap and timing KPIs")
print("Extend only if the cutoff flow is unproved; stop on its proof, a lower-flow counterexample, or the total cap")
print("Report first flow-proof time and total runtime, plus cutoff and final solutions separately")
print("Mode: leave; movement: BM; sequential jobs")
print("Warm start: common complete greedy plan; one continuous solve preserves all search state")
print("MIPFocus: 0 throughout. " + ("Flow proof ends the search immediately at its supported callback."
      if stop_at_flow_proof else "Flow proof ends the extension, but does not end the first phase early."))
source = Path(sys.argv[1])
for arguments in (["rev-parse", "HEAD"], ["status", "--short"]):
    try:
        result = subprocess.run(["git", "-C", str(source), *arguments], capture_output=True, text=True, check=False)
        print("Git", " ".join(arguments) + ":", result.stdout.strip() or result.stderr.strip())
    except OSError as exc:
        print("Git metadata unavailable:", exc)
for name in ("RunTable2SafeWeighted.sh", "RunSafeWeightedStatic.py", "RunWeightedStatic.py", "static_weighted_certification.py",
             "static_safe_weighted_search.py", "static_lexicographic.py", "escort_flow_static_gurobi.py", "load_flow_static_gurobi.py",
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
    merged_csv="$result_dir/table2${table_part}_${formulation}_safe_weighted.csv"
    log_file="$result_dir/logs/table2${table_part}_${formulation}_${lx}x${ly}.log"
    local command=("$python_bin" -u "$script_dir/RunSafeWeightedStatic.py"
        --formulation "$formulation" -x "$lx" -y "$ly" -O "$@"
        -e "$escorts" -l "$loads" -r "$seeds" --threads "$threads"
        --weighted-time-limit "$weighted_limit" --extension-time-limit "$extension_limit"
        -f "$batch_csv")
    if (( stop_at_flow_proof )); then
        command+=(--stop-at-flow-proof)
    fi
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
    "$python_bin" - "$script_dir" "$batch_csv" "$merged_csv" "$seeds" "$escorts" "$weighted_limit" "$extension_limit" "$stop_at_flow_proof" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from RunSafeWeightedStatic import merge_batch
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
    printf '[%s] Completed Table 2 safe weighted runs. Results: %s\n' "$(date '+%F %T')" "$result_dir"
fi
