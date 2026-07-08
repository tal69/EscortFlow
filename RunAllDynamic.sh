#!/usr/bin/env bash
# RunAllDynamic.sh - rerun all dynamic-paper experiments in one go.
#
# Runs, in decreasing order of sensitivity to the 2026-07-08 greedy-heuristic
# fix (see the README Changelog):
#   1. TestHybridRatio.sh   -> HybridRatioTest.csv        (12 runs)
#   2. TestAtten.sh         -> AttentionTest.csv          (12 runs)
#   3. FullFactor9x5.sh     -> FullFactorial9x5-o0.2.csv  (64 runs)
#   4. FullFactor13x7.sh    -> FullFactorial13x7.csv      (64 runs)
#   5. TestIntegrated.sh    -> Modular_vs_integrated.csv  ( 8 runs)
#
# The full battery is long (many hours); run it inside tmux/screen or with:
#   nohup bash RunAllDynamic.sh > runall.log 2>&1 &
#
# Existing CSV outputs are moved aside first (the simulators append!). All
# outputs, per-family logs, and an environment snapshot are collected under
# results_<timestamp>/.
#
# Note: the greedy-vs-optimum estimate (EscortFlowStatic.py --greedy) is not
# part of this battery; run it separately if the appendix numbers are to be
# updated.

set -u
cd "$(dirname "$0")"

STAMP=$(date +%Y%m%d_%H%M)
RESDIR="results_${STAMP}"
mkdir -p "$RESDIR"

# The family scripts invoke `python`; provide a shim if only python3 exists.
if ! command -v python >/dev/null 2>&1; then
  mkdir -p .shim
  ln -sf "$(command -v python3)" .shim/python
  export PATH="$PWD/.shim:$PATH"
fi

# --- environment snapshot (for the reproducibility appendix) ---
{
  echo "date: $(date)"
  echo "host: $(hostname)"
  echo "git commit: $(git rev-parse --short HEAD 2>/dev/null || echo n/a)"
  python3 --version
  python3 - <<'PY'
import numpy
print("numpy:", numpy.__version__)
try:
    import gurobipy
    print("gurobi:", ".".join(map(str, gurobipy.gurobi.version())))
except Exception as exc:  # noqa: BLE001
    print("gurobi: NOT AVAILABLE -", exc)
PY
} | tee "$RESDIR/environment.txt"

# --- preflight ---
python3 - <<'PY' || { echo "ABORT: gurobipy missing or unlicensed"; exit 1; }
import gurobipy as gp
gp.Model()
PY
echo "Preflight: running greedy regression test..."
python3 test_onestep_heuristic.py > "$RESDIR/test_onestep_heuristic.log" 2>&1 \
  || { echo "ABORT: regression test failed (see $RESDIR/test_onestep_heuristic.log)"; exit 1; }
echo "Preflight OK."

run_family () {
  local script=$1; shift
  local status dt t0 c
  # move aside stale outputs (the simulators append)
  for c in "$@"; do
    [ -f "$c" ] && mv "$c" "$RESDIR/stale_${c}"
  done
  echo "=== $(date '+%F %T')  starting $script"
  t0=$SECONDS
  if bash "$script" > "$RESDIR/${script%.sh}.log" 2>&1; then
    status=OK
  else
    status=FAILED
  fi
  dt=$(( SECONDS - t0 ))
  for c in "$@"; do
    [ -f "$c" ] && mv "$c" "$RESDIR/$c"
  done
  printf "%s  %-22s %s  (%d min)\n" "$(date '+%F %T')" "$script" "$status" $((dt/60)) \
    | tee -a "$RESDIR/summary.txt"
}

run_family TestHybridRatio.sh HybridRatioTest.csv
run_family TestAtten.sh       AttentionTest.csv
run_family FullFactor9x5.sh   FullFactorial9x5-o0.2.csv
run_family FullFactor13x7.sh  FullFactorial13x7.csv
run_family TestIntegrated.sh  Modular_vs_integrated.csv

echo "=== all done; outputs in $RESDIR/"
cat "$RESDIR/summary.txt"
