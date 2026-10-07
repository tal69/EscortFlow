#!/usr/bin/env bash
# Current manuscript Table 2(a,b): four layouts, 3,600 instances, both formulations.
# Two-phase integer solves only; the historical weighted LP columns are separate.
# Jobs run sequentially so they do not compete for solver threads.
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON:-python3}
threads=${NUM_THREADS:-12}
seeds=1-100
part=both
dry_run=0
result_dir="${script_dir}/results_table2_lex_$(date +%Y%m%d_%H%M%S)_$$"

usage() {
    cat <<'USAGE'
Usage: bash RunTable2Lex.sh [options]

Run Table 2(a) and (b) with the lexicographic escort-flow and load-flow models.
Each solve has a 270-second phase-one cap and a 300-second total solver budget.
The phase-two cap is 300 minus the actual phase-one solver runtime.

Options:
  --python PATH      Python executable with numpy, gurobipy, and a Gurobi license
                     (default: $PYTHON, otherwise python3; Python 3.11 recommended)
  --threads N        Solver threads per instance (default: $NUM_THREADS or 12)
  --seeds RANGE      Instance seeds in the runners' range syntax (default: 1-100)
  --part a|b|both     Table part to run (default: both)
  --output-dir DIR   New results directory; existing directories are rejected
  --dry-run          Print all solver commands without creating files or solving
  -h, --help         Show this help

Examples:
  bash RunTable2Lex.sh
  bash RunTable2Lex.sh --python /opt/conda/envs/pbs/bin/python --threads 12
  bash RunTable2Lex.sh --dry-run
  nohup bash RunTable2Lex.sh > table2_lex.log 2>&1 &

Defaults match the current manuscript: 13x7, 10x10, 16x10, and 27x10 grids;
leave retrieval with simultaneous block movements; 100 seeds per table row.
Part (a): one target, 3-8 escorts. Part (b): four targets, 8/12/16 escorts.
The existing per-formulation heuristic horizon selection is preserved.
USAGE
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

while (($#)); do
    case $1 in
        --python|--threads|--seeds|--part|--output-dir)
            (($# >= 2)) || die "Missing value for $1"
            case $1 in
                --python) python_bin=$2 ;;
                --threads) threads=$2 ;;
                --seeds) seeds=$2 ;;
                --part) part=$2 ;;
                --output-dir) result_dir=$2 ;;
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

    # Check the same Python interpreter and license used by every experiment.
    if ! "$python_bin" - "$script_dir" "$threads" "$seeds" "$part" <<'PY' > "$result_dir/environment.txt" 2>&1
import datetime
import platform
import sys

print("Started:", datetime.datetime.now().astimezone().isoformat(), flush=True)
print("Host:", platform.node(), flush=True)
print("Platform:", platform.platform(), flush=True)
print("Python:", sys.version.replace("\n", " "), flush=True)
print("Executable:", sys.executable, flush=True)
if sys.version_info < (3, 10):
    raise SystemExit("Python 3.10 or newer is required; 3.11 is recommended")
import numpy
import gurobipy as gp
sys.path.insert(0, sys.argv[1])
from PBSCom import str2range
seed_values = list(str2range(sys.argv[3]))
if not seed_values or len(set(seed_values)) != len(seed_values) or min(seed_values) < 0:
    raise SystemExit("Seeds must form a nonempty range of distinct nonnegative integers")
print("NumPy:", numpy.__version__)
print("Gurobi:", ".".join(map(str, gp.gurobi.version())))
print("Threads:", sys.argv[2])
print("Seeds:", sys.argv[3], "count:", len(seed_values))
print("Table part:", sys.argv[4])
print("Phase 1 cap: 270 seconds; total cap: 300 seconds")
print("Mode: leave; movement: BM; automatic horizons; sequential jobs")
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
fi

run_layout() {
    local table_part=$1 formulation=$2 lx=$3 ly=$4
    shift 4
    local loads escorts runner batch_csv merged_csv log_file
    if [[ $table_part == a ]]; then
        loads=1; escorts=3-8
    else
        loads=4; escorts=8-16-4
    fi
    if [[ $formulation == escortflow ]]; then
        runner=EscortFlowStaticLex.py
    else
        runner=LoadFlowStaticLex.py
    fi
    batch_csv="$result_dir/parts/table2${table_part}_${formulation}_${lx}x${ly}.csv"
    merged_csv="$result_dir/table2${table_part}_${formulation}_lex.csv"
    log_file="$result_dir/logs/table2${table_part}_${formulation}_${lx}x${ly}.log"
    local command=("$python_bin" -u "$script_dir/$runner"
        -x "$lx" -y "$ly" -O "$@" -e "$escorts" -l "$loads" -r "$seeds"
        -m leave --gurobi --phase1_time_limit 270 --time_limit 300
        --num_threads "$threads" --mip_emphasis balanced -f "$batch_csv")
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

    # The runners catch some per-instance exceptions and still exit zero. Check
    # row coverage and status before merging; retain raw batches and their logs.
    "$python_bin" - "$script_dir" "$batch_csv" "$merged_csv" "$seeds" "$escorts" <<'PY'
import csv
from pathlib import Path
import sys

sys.path.insert(0, sys.argv[1])
from PBSCom import str2range
source, destination = map(Path, sys.argv[2:4])
expected = {(seed, escorts) for seed in str2range(sys.argv[4]) for escorts in str2range(sys.argv[5])}
with source.open(newline="") as f:
    records = [[cell.strip() for cell in row] for row in csv.reader(f) if any(cell.strip() for cell in row)]
if not records:
    raise SystemExit(f"Empty result file: {source}")
header = records[0]
rows = [row for row in records[1:] if row != header]
if any(len(row) != len(header) for row in rows):
    raise SystemExit(f"Malformed CSV row: {source}")
items = [dict(zip(header, row)) for row in rows]
actual = [(int(row["seed"]), int(row["# Escorts"])) for row in items]
if len(actual) != len(expected) or set(actual) != expected:
    raise SystemExit(f"Missing or duplicate instance results: {source}")
for row in items:
    if row["Solver Status"] in {"", "-", "ERROR"}:
        raise SystemExit(f"Solver error at seed {row['seed']}, escorts {row['# Escorts']}: {source}")
    if float(row["Requested Phase 1 Time Limit"]) != 270 or float(row["Total Time Limit"]) != 300:
        raise SystemExit(f"Unexpected solver time limits: {source}")
exists = destination.exists()
if exists:
    with destination.open(newline="") as f:
        if next(csv.reader(f)) != header:
            raise SystemExit(f"Cannot merge incompatible CSV schemas: {destination}")
with destination.open("a", newline="") as f:
    writer = csv.writer(f)
    if not exists:
        writer.writerow(header)
    writer.writerows(rows)
print(f"Validated and saved {len(rows)} instance results to {destination}", flush=True)
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
    printf '[%s] Completed Table 2 lexicographic runs. Results: %s\n' "$(date '+%F %T')" "$result_dir"
fi
