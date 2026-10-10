# Continue-mode experiment design

The Bash driver `RunContinue.sh` calls `SolveStatic.py` twice per configuration,
once with `--formulation loadflow` and once with `--formulation escortflow`.
Its loops define the campaign directly. `SolveStatic.py` runs the selected
Gurobi model in the current process using the shared corrected protocol.
There is no Python campaign launcher. The README explains individual solves.

## Scope and order

The study addresses reviewer comments R2.1 and R3.3 about multi-target continue
retrieval, and supplies two- and six-target comparisons relevant to R2.4.
Continue retrieval converts a target into an ordinary blocking load at its
first output arrival. That load may subsequently move. No new escort is created
and there is no additional output-service step. Targets initially at outputs
are served at time zero. All experiments use simultaneous block movement (SBM).

| Layout | Output coordinates (x, y) |
| --- | --- |
| 13 x 7 | (6, 0) |
| 10 x 10 | (0, 0) |
| 16 x 10 | (4, 0), (11, 0) |
| 27 x 10 | (4, 0), (13, 0), (22, 0) |

Each layout has 2, 4, and 6 targets, 8, 12, and 16 escorts, and seeds 1-100.
This is 36 configurations, 3,600 matched instances, and 7,200 integer runs.
The default also requests 7,200 LP bounds, with retries if necessary.
Four-target configurations run first, starting with 13 x 7, followed by two
and six targets. Within each configuration LF runs before EF. Both receive
identical generated coordinates, outputs, greedy plans, objective coefficients,
physical horizons, and solver budgets. Seeds identify paired initial states
within a configuration. Escort locations can differ across target counts.
Initial occupancy is `(N-e)/N`, so equal escort counts across layouts do not
mean equal density. There is no added 70%-occupancy experiment.

## Objective, horizon, and budgets

Both methods use `safe_integer_flow_timing_v6`. The feasible common greedy plan
gives flow time `F_g` and makespan `C_g`. If target `i` has nearest-output
Manhattan distance `d_i`, `D=sum(d_i)`, `N` is the cell count, and `e` is the
initial escort count, the automatic parameters are:

```text
K   = N-e
H_g = F_g-D+max(d_i)
R   = K*H_g-D+1
H   = max(H_g, C_g+1)
objective = R*F+M
EF last movement index = H-1
LF last layer index    = H
```

`R` is a sufficient integer weight for prioritizing flow time over movements;
the reported objective is also available as `F+M/R`. The physical horizon covers
every plan with flow time at most the greedy bound and retains the full warm
start. The different backend indices represent the same physical horizon.
No user-supplied movement coefficient or heuristic horizon multiplier is needed.

Each integer solve gets 300 solver seconds initially, with 16 threads. The
solver records the first observed flow-time proof but keeps searching for
movement optimality during that initial budget. Only when flow time is still
unproved does it continue the same search for up to 300 more seconds, stopping
the extension on flow-time proof. Initial-cutoff results and final results have
separate CSV fields. The campaign does not enable `--stop-at-flow-proof`.

After each integer run, its continuous relaxation uses exactly its recorded
coordinates, coefficient, and model horizon. The LP uses one thread, a separate
300-second budget, and a 600-second barrier retry following a time limit.
LP time is excluded from the integer timing metrics. An unfinished LP retains
the integer row, leaves its bound blank, and does not stop later integer batches.
The row records original model variable, constraint, and nonzero counts before
presolve for subsequent model-size comparisons.

## Running the study

Run from the directory containing the scripts (`Code` in the research project,
or the root of a standalone clone) on the numerical-study Linux machine, using
Python with NumPy and a Gurobi license large enough for the models. Use the same machine,
Gurobi version, and thread settings for both methods. Run one campaign at a time
to avoid competing solver workloads. The current v4 process is not stopped by
this script, and its archived results remain provisional for the revised paper.

Preview all 72 batch commands without writing files or starting Gurobi:

```bash
bash RunContinue.sh --dry-run
```

First run a small pilot at the intended production budgets:

```bash
PYTHON=/path/to/licensed/python SEEDS=1-5 LAYOUTS=13x7 TARGET_COUNTS=4 \
  bash RunContinue.sh pilot_continue_v6
```

This pilot contains 15 pairs, 30 integer runs, and their LP bounds. Then launch
the complete study in an existing tmux session:

```bash
PYTHON=/path/to/licensed/python bash RunContinue.sh results_continue_v6
```

If `python3` already belongs to the licensed environment, omit `PYTHON`.
Output directories must be new. All solves run sequentially. Parameters can be
changed explicitly using the variables listed by `bash RunContinue.sh --help`.
For example, `WITH_LP=0` disables LPs, and `THREADS=1` selects one integer thread.
Keep the production defaults when comparing with the 300-second paper protocol.

## Results, interruption, and analysis

The output contains:

- `source/`: frozen Python sources and this Bash script.
- `environment.json`: software, host, CPU/memory information on Linux, seeds,
  scope, budgets, and source hashes.
- `commands.sh`: exact shell-escaped commands issued by the driver.
- `parts/`: one integer CSV per formulation and configuration, including LP fields.
- `logs/`: corresponding logs; `preflight.log` records the license check.

`commands.sh` is an execution record, not an automatic resume script. The driver
refuses an existing output directory and never overwrites completed integer
results. On interruption, already flushed rows are retained. Inspect both
formulations' seed coverage, then run the remaining seeds into a new directory
using the original frozen `source/RunContinue.sh`. If one method already
completed a partially finished configuration, call the frozen `SolveStatic.py`
directly for the other method's missing seeds. Check coordinates, coefficient,
physical horizon, source hashes, protocol, budgets, and duplicate instance keys
before combining outputs. Do not combine v4/v5 rows with corrected v6 results.

To collect saved LP bounds and retry only missing/failed ones using frozen code:

```bash
RESULTS=/absolute/path/to/results_continue_v6
python3 "$RESULTS/source/RunStaticLP.py" \
  --input "$RESULTS/parts/*.csv" --source-dir "$RESULTS/source" \
  --workers 1 --threads 1 --time-limit 300 --retry-time-limit 600 \
  --resume --fill-input-lp -f "$RESULTS/lp_results.csv"
```

Use the same licensed Python environment as the campaign. This command reuses
valid embedded optimal LP bounds after checking their input/source fingerprints.
It does not rerun integer searches.

For the method comparison, use the initial 300-second columns and compare both
flow-time proof and full-objective proof. Report counts, medians, quartiles, and
timeouts per configuration; unresolved times are censored observations. Treat
final extension results separately. For solution-quality tables, use the best
known feasible flow time and movements across the paired final solutions and
mark any unproved values. LP gaps use a common integer reference objective and
each formulation's own relaxation bound. A paired runtime plot can start with
the 300 four-target 13 x 7 instances, on linear axes with the equality line and
winner colors. A runtime pattern alone does not establish the cause of hardness.

The driver and a short local pilot are validated. Full Linux execution,
aggregation, insertion of tables/figures, and revision of numerical claims
remain pending in `tasks.md`.
