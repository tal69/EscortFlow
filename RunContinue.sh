#!/usr/bin/env bash
# Continue-mode SBM study. Each command directly selects the LF or EF solver.
# Change the defaults here, or set these variables before running the script.
# 36 configurations x 100 seeds = 3,600 pairs and 7,200 integer solves.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON:-python3}
seeds=${SEEDS:-1-100}
threads=${THREADS:-16}
layouts=${LAYOUTS:-"13x7 10x10 16x10 27x10"}
target_counts=${TARGET_COUNTS:-"4 2 6"}  # Four targets first, then two and six.
escort_counts=${ESCORT_COUNTS:-"8 12 16"}
time_limit=${TIME_LIMIT:-300}
extension_limit=${EXTENSION_LIMIT:-300}
with_lp=${WITH_LP:-1}
lp_threads=${LP_THREADS:-1}
lp_time_limit=${LP_TIME_LIMIT:-300}
lp_retry_limit=${LP_RETRY_LIMIT:-600}

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    cat <<'HELP'
Usage: bash RunContinue.sh [--dry-run] [NEW_OUTPUT_DIRECTORY]

Defaults: all four layouts; 4, 2, 6 targets; 8, 12, 16 escorts; seeds 1-100;
16 MIP threads; initial 300 seconds; conditional 300-second extension; LPs on.
Variables: PYTHON, SEEDS, THREADS, LAYOUTS, TARGET_COUNTS, ESCORT_COUNTS,
TIME_LIMIT, EXTENSION_LIMIT, WITH_LP, LP_THREADS, LP_TIME_LIMIT, LP_RETRY_LIMIT.

Examples:
  bash RunContinue.sh --dry-run
  PYTHON=/path/to/python bash RunContinue.sh results_continue_v6
  SEEDS=1 TARGET_COUNTS=4 LAYOUTS=13x7 ESCORT_COUNTS=16 THREADS=1 \
    bash RunContinue.sh pilot_continue

SolveStatic.py is the common entry point for both formulations. All solves are
sequential. Results and logs are separate for each configuration
and formulation. Existing output directories are refused. Source is frozen in
OUTPUT/source before solving. See RunContinue.md for design and interrupted runs.
HELP
    exit 0
fi
dry_run=0
if [[ ${1:-} == --dry-run ]]; then dry_run=1; shift; fi
(( $# <= 1 )) || { printf 'Expected at most one output directory.\n' >&2; exit 2; }
result_dir=${1:-"$script_dir/results_continue_v6_$(date -u +%Y%m%d_%H%M%S)_$$"}
[[ $result_dir == /* ]] || result_dir="$(pwd)/$result_dir"
[[ $threads =~ ^[1-9][0-9]*$ && $lp_threads =~ ^[1-9][0-9]*$ && $with_lp =~ ^[01]$ ]]
for layout in $layouts; do
    case $layout in 13x7|10x10|16x10|27x10) ;; *) printf 'Unknown layout: %s\n' "$layout" >&2; exit 2;; esac
done
for targets in $target_counts; do
    case $targets in 2|4|6) ;; *) printf 'Targets must be 2, 4, or 6.\n' >&2; exit 2;; esac
done
for escorts in $escort_counts; do
    case $escorts in 8|12|16) ;; *) printf 'Escorts must be 8, 12, or 16.\n' >&2; exit 2;; esac
done
runtime_dir="$result_dir/source"
if (( ! dry_run )); then
    [[ ! -e $result_dir ]] || { printf 'Output directory already exists: %s\n' "$result_dir" >&2; exit 2; }
    mkdir -p "$runtime_dir" "$result_dir/parts" "$result_dir/logs"
    cp "$script_dir"/*.py "$runtime_dir/"
    cp "$script_dir/RunContinue.sh" "$runtime_dir/"
    cp "$script_dir/README.md" "$script_dir/RunContinue.md" "$runtime_dir/"
    "$python_bin" - "$runtime_dir" "$seeds" "$threads" "$layouts" "$target_counts" \
        "$escort_counts" "$time_limit" "$extension_limit" "$with_lp" "$lp_threads" \
        "$lp_time_limit" "$lp_retry_limit" <<'PY' > "$result_dir/preflight.log"
import hashlib, json, math, os, platform, sys
from datetime import datetime, timezone
from pathlib import Path
(_, source_dir, seed_text, mip_threads, layout_text, target_text, escort_text,
 initial_text, extension_text, lp_enabled, lp_thread_text, lp_time_text, lp_retry_text)=sys.argv
source=Path(source_dir);sys.path.insert(0,str(source))
import numpy, gurobipy as gp
from RunSafeWeightedStatic import PROTOCOL, parse_range
seed_values=parse_range(seed_text, minimum=0)
initial,extension,lp_limit,lp_retry=map(float,(initial_text,extension_text,lp_time_text,lp_retry_text))
assert all(math.isfinite(v) and v>0 for v in (initial,lp_limit,lp_retry)) and math.isfinite(extension) and extension>=0
for values in (layout_text.split(),target_text.split(),escort_text.split()):
    assert values and len(set(values))==len(values), 'Empty or repeated configuration choice'
with gp.Env(empty=True) as env:
    env.setParam('OutputFlag',0);env.start()
    with gp.Model(env=env) as model:
        model.addVar(lb=0,obj=1);model.optimize()
        assert model.Status==gp.GRB.OPTIMAL, 'Gurobi license check failed'
configs=len(layout_text.split())*len(target_text.split())*len(escort_text.split())
metadata=dict(started_utc=datetime.now(timezone.utc).isoformat(),protocol=PROTOCOL,
    python=sys.executable,python_version=sys.version,numpy_version=numpy.__version__,
    gurobi_version=gp.gurobi.version(),platform=platform.platform(),host=platform.node(),
    processor=platform.processor(),logical_cpus=os.cpu_count(),seed_values=seed_values,
    layouts=layout_text.split(),targets=list(map(int,target_text.split())),
    escorts=list(map(int,escort_text.split())),retrieval_mode='continue',movement_mode='BM',
    threads=int(mip_threads),initial_time_limit=initial,extension_time_limit=extension,
    with_lp=bool(int(lp_enabled)),lp_threads=int(lp_thread_text),lp_time_limit=lp_limit,
    lp_retry_time_limit=lp_retry,configurations=configs,matched_instances=configs*len(seed_values),
    integer_runs=2*configs*len(seed_values),method_order=['loadflow','escortflow'],
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.iterdir()) if p.is_file()})
if sys.platform.startswith('linux'):
    metadata['cpu_info']=Path('/proc/cpuinfo').read_text()
    metadata['memory_info']=Path('/proc/meminfo').read_text()
(source.parent/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
print('Preflight passed: '+str(configs*len(seed_values))+' matched instances.')
PY
    printf '#!/usr/bin/env bash\nset -euo pipefail\n' > "$result_dir/commands.sh"
fi

# Print and save the exact solver command, then run it with a separate log.
run_solver() {
    local label=$1 status
    shift
    if (( dry_run )); then printf '%q ' "$@"; printf '\n'; return; fi
    printf 'Starting %s; log: %s/logs/%s.log\n' "$label" "$result_dir" "$label"
    printf '%q ' "$@" >> "$result_dir/commands.sh"; printf '\n' >> "$result_dir/commands.sh"
    if "$@" > "$result_dir/logs/$label.log" 2>&1; then status=0; else status=$?; fi
    case $status in
        0) ;;
        3) printf '%s: integer results saved; some LP bounds remain pending.\n' "$label" >&2;;
        *) printf '%s failed (exit %s); inspect its log. Saved rows are retained.\n' "$label" "$status" >&2; return "$status";;
    esac
}

run_layout() {
    local lx=$1 ly=$2 escorts suffix
    shift 2  # Remaining arguments are the output-cell coordinate pairs.
    case " $layouts " in *" ${lx}x${ly} "*) ;; *) return 0;; esac
    for escorts in $escort_counts; do
        suffix="${lx}x${ly}_t${targets}_e${escorts}"
        local common=(-x "$lx" -y "$ly" -O "$@" -l "$targets" -e "$escorts" -r "$seeds"
            -m continue --threads "$threads" --weighted-time-limit "$time_limit" --extension-time-limit "$extension_limit")
        if (( with_lp )); then
            common+=(--with-lp --lp-threads "$lp_threads" --lp-time-limit "$lp_time_limit"
                     --lp-retry-time-limit "$lp_retry_limit")
        fi
        # Two direct solver calls; the Bash loops define the campaign.
        run_solver "loadflow_$suffix" "$python_bin" -u "$runtime_dir/SolveStatic.py" \
            --formulation loadflow "${common[@]}" -f "$result_dir/parts/loadflow_$suffix.csv"
        run_solver "escortflow_$suffix" "$python_bin" -u "$runtime_dir/SolveStatic.py" \
            --formulation escortflow "${common[@]}" -f "$result_dir/parts/escortflow_$suffix.csv"
    done
}

for targets in $target_counts; do
    run_layout 13 7 6 0
    run_layout 10 10 0 0
    run_layout 16 10 4 0 11 0
    run_layout 27 10 4 0 13 0 22 0
done
if (( ! dry_run )); then printf 'Integer batches completed. Results: %s\n' "$result_dir"; fi
