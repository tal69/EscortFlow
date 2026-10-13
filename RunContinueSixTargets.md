# Revised six-target continue campaign

The revised October 13, 2026 plan uses all four layouts, six target loads,
12/16/20 escorts, and seeds 1-100. There are 12 configurations, 1,200 matched
instances, and 2,400 integer runs. LF runs for an initial 300 solver seconds
without an extension. EF uses the same initial 300-second cutoff and a
conditional 300-second extension when FT remains unproved. Both methods keep
searching for MV optimality throughout the initial phase.

The launch script freezes the previously audited code from
`results_continue_v6/source`, verifies every recorded source hash, and preserves
the original campaign. No solver model is modified. Initial states, outputs,
greedy starts, coefficients, physical horizons, and initial 16-thread budgets
match between LF and EF. Each run also requests its matching LP bound, using
one LP thread, 300 seconds, and a 600-second barrier retry when needed.

The existing six-target/eight-escort run must be stopped separately by Tal
before relaunch. The launcher detects active study processes and refuses a
competing run; it never stops or restarts a process.

On bukhan1, from the project directory:

```bash
cd ~/EscortFlow
PYTHON=/home/admina/miniconda3/bin/python3.12 bash RunContinueSixTargets.sh --dry-run results_continue_v6_six_targets
PYTHON=/home/admina/miniconda3/bin/python3.12 bash RunContinueSixTargets.sh results_continue_v6_six_targets
```

Use the launch command in your existing tmux/terminal environment. The new
output directory must not exist. It contains frozen source, `environment.json`,
exact `commands.sh`, per-configuration CSVs and logs, and `run_state.json`
recording running/completed/failed status. Integer batches with pending LPs are
distinguished from fully complete output. Ordinary solver time limits are valid
outcomes. A failed integer record or incomplete seed batch stops the driver.

The default output path is `results_continue_v6_six_targets` beside the script.
Optional `--layouts` and `--seeds` support separate pilots or interrupted-batch
recovery into a fresh directory; they do not resume or overwrite an existing
campaign. `--source-dir` can select another verified frozen v6 source snapshot.
Do not mix six-target/eight-escort rows from the superseded plan into this study.

Paper comparisons use only each method's own initial-cutoff results. Shared
best-known solutions and LP gaps can use LF's available final result and EF's
conditional extension. Publish the different extension policies explicitly.

The completed two-/four-target audit found zero LF-exclusive FT certificates
across 2,400 matched instances, including after extensions. LF extensions added
28 LF certificates already supplied by EF, but improved the shared feasible
pair in four cases. Thus removing them saves time while potentially losing some
best-known solution quality. This evidence does not guarantee the same pattern
for the six-target campaign.
