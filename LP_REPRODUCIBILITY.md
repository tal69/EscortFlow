# Reproducing the static LP relaxation values

`RunStaticLP.py` is the separate reproducibility runner for the LP lower bounds
underlying the model-specific gaps now reported in Table 2. The archived package
retains its original `table3_lp` name. It runs either formulation with all integer variables made
continuous and minimizes

```text
Z_LP = FT_LP + MV_LP / R.
```

It reads the exact initial coordinates, output cells, coefficient `R`, retrieval
mode, and horizon from each recorded integer experiment row. It does not
regenerate the greedy plan or recalculate the coefficient. The archived one-target calculation used the
archived v4 coefficients; the current v5 coefficient rule is different. The
archived source snapshot and inputs must therefore accompany the reproduction.
The replay also accepts new v5 records without converting their coefficients.

Paths and standard-runner commands in this guide refer to the full `Code/`
directory. For the standalone ZIP, follow its packaged `README.md`.

The physical horizon is `T+1` for escort flow and `T` for load flow. The frozen
single-target inputs match these physical horizons across formulations. Supplying the
same numerical `T` to both formulations would change the comparison.

## All eight four-target LP parts from seeds on the Mac Studio

The four-target integer campaign on Linux can continue while its LP relaxations
run separately on the Mac Studio. From `Code/`, in your licensed Conda environment:

```bash
python -u RunFourTargetLP.py
```

No arguments selects the four layouts (13x7, 10x10, 16x10, 27x10), both
formulations, four targets, 8/12/16 escorts, seeds 1-100, leave mode: eight parts
with 300 LPs each. It calls the **same `static_generated_lp.generated_row`
function used by the standard runners' `--lp` path**, including the original
seed generator and greedy heuristic. It fixes the **v4** coefficient and each
formulation's archived horizon to match the currently running integer campaign.
The future two-/six-target campaign remains v5.

Every available selected row in `Experiment Oct2026/table2b_*.csv` is checked
for identical coordinates, R and horizon. Those files may be missing or partial;
they do not select the LPs. `--input-dir` changes the verification folder. The
integer files are read-only.

Results and immutable source/instance snapshots go into `results_four_target_lp/`.
The combined `lp_results.csv` uses the standard replay schema, so it retains
fractional FT/MV, objective, status, timing, R, both horizon conventions and
source/problem fingerprints. The eight part CSVs, manifests and worker logs
are in `parts/`; `coverage.json` reports complete and missing LPs. Only optimal
LPs are recorded. Repeating the same command resumes with the original frozen
sources and retries missing cases. A changed selection needs a new output folder.

```bash
python RunFourTargetLP.py --dry-run
python RunFourTargetLP.py --check-environment
python -u RunFourTargetLP.py --layouts 13x7 --escorts 16 --seeds 1 \
    --output-dir results_four_target_lp_pilot
```

Use `--seeds 1-100` or another inclusive range, `--workers 1`, `--threads 16`,
`--time-limit 300`, and `--retry-time-limit 600` as needed. Defaults use one LP
at a time and the Mac's performance-core count for threads. Seed generation
requires NumPy in addition to Python 3.10+, Gurobi and a full license. Prefer
`python` on the Mac Studio where `python3` resolves to Apple's older interpreter.
This runner computes LP bounds independently; complete relative-gap means are
aggregated after the missing integer results arrive.

## Four-target Table 2(b) from available CSV rows on the Mac Studio

From `Code/`, set up the Mac's interpreter once, then start the LP campaign:

```bash
bash RunTable2bLP.sh --setup python3.13
bash RunTable2bLP.sh
```

The setup command creates `~/.venvs/escortflow-table2b`, installs
`gurobipy==13.0.3`, validates the full license, and saves the interpreter outside
Dropbox in the Mac's user configuration. It does not start the campaign. New
terminals and tmux windows need no environment activation. To reuse an existing
licensed environment, use `--set-python /absolute/path/to/python` instead.
The runner requires Python 3.10+; a different installed Python version does not
change the interpreter selected by a `python3` command. Use
`bash RunTable2bLP.sh --check-environment` to see the exact path and version and
check Gurobi independently of the experiment inputs. Installing the pip package
does not install a full Gurobi license.

This dedicated wrapper reads the four-target leave-mode Table 2(b) CSVs and
escort counts 8/12/16. It defaults to one LP at a time and the Mac's performance-
core thread count. Results and frozen model/input snapshots are saved in
`results_table2b_lp/`. Repeat the same command to resume after an interruption,
retry failed LPs, or process newly copied integer rows. The wrapper uses
`RunStaticLP.py --resume --extend` after validating the saved source tree.
To generate LPs for seeds whose integer CSV rows are still missing, use the
direct parameter-based runners described below.

Archived four-target v4 physical horizons can differ by one period across
formulations. Each LP uses its own integer run's recorded horizon and coefficient.
The weighted LP objective is `FT+MV/R`. The gap summary uses the common
best-known lexicographic FT/MV pair, including available extension solutions
from both models. Table means require all 100 seeds by default, so unfinished
method groups remain blank.

`lp_results.csv` uses the general replay format. The wrapper also saves
`lp_gaps.csv`, `table2b_lp_summary.csv`, `table2b_lp_columns.tex`, and
`coverage.json`. See the Mac Studio section in [README.md](README.md) for pilot
commands, settings, and output definitions.

## Self-contained package

`reproducibility/table3_lp.zip` contains the runner, all eight frozen input CSVs
(4,800 model solves, 2,400 matched instances), the five frozen model/helper source
files, reference LP results, pinned solver dependency, and a checksum manifest.
It runs independently of the working repository. Unzip it, enter `table3_lp/`,
and follow its README. Building the package from the repository is reproducible:

```bash
python3 BuildTable3LPPackage.py --output-dir reproducibility/table3_lp
```

The builder refuses to overwrite an existing directory or archive. Use a new
destination for a later package version.

## Replay from the working Code directory

Use Python 3.10 or later, `gurobipy==13.0.3`, and a Gurobi license that supports
these model sizes. The original LP batch used Gurobi 13.0.3, four independent
workers, and one solver thread per worker. The LP replay itself needs no NumPy,
OPL, DP data, or integer solver campaign.

```bash
python3 RunStaticLP.py \
  --input 'revision_R1/table3_bounds_2026-10-09/inputs/table2a_*.csv' \
  --source-dir revision_R1/table3_bounds_2026-10-09/source \
  --workers 4 --threads 1 \
  --check-against revision_R1/table3_bounds_2026-10-09/lp_results.csv \
  -f lp_reproduced.csv
```

To check eight cases first, add `--seeds 1 --escorts 3` and choose a separate
output file. Filters accept inclusive ranges and comma-separated values, such
as `--seeds 1,44` or `--escorts 3-8`. `--formulation loadflow` or
`--formulation escortflow` restricts the run to one method. For the future
four-target panel, supply the recorded four-target integer CSVs in the same way.

The default first solve has a 300-second solver limit. A `TIME_LIMIT` result is
retried using barrier without crossover, with a 600-second solver limit,
matching the fallback used in the Table 3 calculation. Override these with
`--time-limit` and `--retry-time-limit`, or disable the retry with `--no-retry`.
Only `OPTIMAL` LP results enter the CSV. Nonoptimal or failed cases are recorded
in diagnostics, the remaining independent cases continue, and the command exits
with a nonzero status if the batch is incomplete.

Repeat the identical command with `--resume` after an interruption or failure.
The runner checks the selected instances, source hashes, stored settings, status,
and objective components before reusing saved results. Changed existing instances
or model sources require a new output. Existing optimal records are reused.
`--resume --extend` additionally permits new instances, provided every previously
selected instance and all model sources remain unchanged. This is the mode
used automatically by `RunTable2bLP.py` as copied four-target input files grow.

Alongside the CSV, the runner saves a `.manifest.json` file and a `.logs/`
directory. Each result retains its coefficient, model and physical horizon,
target-load count, solver version, input fingerprint, source fingerprint, algorithm, fractional FT
and MV, LP objective, and elapsed wall time. Elapsed time includes construction,
extraction, and a retry when applicable; it is not the CPU time in Table 2.

`--check-against` verifies every selected LP objective against the saved
reference, with absolute tolerance `1e-6` and relative tolerance `1e-8`.
Fractional FT and MV can vary across alternative optimal LP solutions, so their
individual values are not required to match the reference decomposition. Each
new decomposition is checked against `FT + MV/R`.

The LP percentage gap in Table 2 uses the common best-known integer solution,
including either formulation's extension when applicable:

```text
Z_BK = FT_BK + MV_BK / R
gap (%) = 100 * (Z_BK - Z_LP) / Z_BK.
```

Percentages are calculated per instance and then averaged. A zero reference
objective contributes zero gap. This runner calculates the LP values; table
aggregation combines them with the best-known integer results.

## New two- and six-target experiments

`RunTable2Targets.py` runs all data needed for the additional Tables 2 and 3
panels in one command. Across the four layouts and seeds 1-100, it uses 8/12/16
escorts for two targets and 8/12/16/20 escorts for six targets. Each of the
5,600 integer searches uses the 300-second comparison cutoff and conditional
300-second extension, and is followed immediately by its matched LP solve.
The bound is saved in that same instance's CSV row. Both solves use the frozen
model source tree and the same coefficient/horizon.

After the existing four-target campaign has finished, in the Linux `Code` checkout:

```bash
git pull --ff-only
python -u RunTable2Targets.py --threads 16 --lp-threads 16
```

Use the same Python environment and Gurobi license as the current campaign, or
select it with `--python /path/to/python`. For a quick functional pilot, use a
fresh output directory and add `--layouts 13x7 --seeds 1`
`--weighted-time-limit 1 --extension-time-limit 1`.
These shortened budgets are for verification, not for reported experiments.

The new result CSV includes `loads` and `retrieval_mode` alongside layout,
escorts, seed, and method, so two- and six-target records cannot collide.
Legacy one-target references without these columns are interpreted as one
target in leave mode. Resume requires the schema and fingerprints from the
same runner version; an old completed CSV is still usable as a reference.

The integer CSV now includes `lp_relaxation_lower_bound` (the optimal value of
`FT+MV/R`), `lp_status`, `lp_elapsed_seconds`, fractional FT/MV, LP budgets and
threads, solver version, algorithm, and source/input fingerprints. Only an
`OPTIMAL` relaxation fills the bound column. Failed or unfinished LPs leave it
blank. Integer solve and CPU measurements exclude LP computation.

The integer result is saved before the LP starts. Integrated LPs run one at a
time, after the integer solver closes, and use no warm start or incumbent
cutoff. They do not strengthen the integer benchmark's search. `--no-lp`
disables them for an integer-only campaign. `ReproducePaper.py` uses the same
integrated sequence for both retrieval modes.

At the end, `RunStaticLP.py` collects the saved bounds into `lp_results.csv`,
checks their metadata, and solves only missing cases. `--fill-input-lp` fills
missing LP columns in new integrated integer CSVs from those successful
retries. Historical CSVs without the LP schema remain read-only. The launcher
also updates its per-configuration CSVs after recovery. `--lp-workers` applies
only to this recovery phase; `--lp-threads` applies to every LP solve.

The recovery phase saves its manifest, worker logs, and the complete restart
command `lp_commands.sh`. If only that phase fails or is interrupted,
use `bash /path/to/results_table2b_targets_.../lp_commands.sh`. This skips already
verified optimal LP values. It does not restart the integer campaign.
See [README.md](README.md) for settings, outputs, and progress monitoring.

For a standalone integer-plus-LP batch, use `--with-lp`:

```bash
python -u RunSafeWeightedStatic.py --formulation escortflow \
  -x 13 -y 7 -O 6 0 -l 2 -e 8,12,16 -r 1-100 -m leave \
  --threads 16 --with-lp --lp-threads 16 -f two_targets_escortflow.csv
```

`--lp` remains the LP-only switch; it cannot be combined with `--with-lp`.
The integrated integer-plus-LP workflow uses v5. Use the separate v4 replay
for the already running four-target campaign.

## Generate LPs directly from the usual experiment parameters

Add `--lp` to `EscortFlowStatic.py`, `LoadFlowStatic.py`, or
`RunSafeWeightedStatic.py`. Supply the usual grid, output cells, target count,
escort range, seed range, and retrieval mode. No integer CSV is required and
no integer optimization is performed. Both formulations use the same instance
generator as the integer campaign. The greedy reference supplies the data for
choosing each instance's coefficient and horizon.

Use `python` from the Conda environment that has NumPy, `gurobipy`, and the
Gurobi license. For the archived four-target campaign, use `--lp-protocol v4`:

```bash
python -u EscortFlowStatic.py -x 13 -y 7 -O 6 0 -l 4 -e 8,12,16 -r 1-100 \
  -m leave --lp --lp-protocol v4 --num_threads 16 -t 300 \
  -f lp_four_escortflow_13x7.csv

python -u LoadFlowStatic.py -x 13 -y 7 -O 6 0 -l 4 -e 8,12,16 -r 1-100 \
  -m leave --lp --lp-protocol v4 --num_threads 16 -t 300 \
  -f lp_four_loadflow_13x7.csv
```

These commands each generate all 300 selected instances, including seeds whose
integer results have not arrived. Change the grid and outputs to cover the
other configurations: `10x10` with `-O 0 0`, `16x10` with `-O 4 0 11 0`, and
`27x10` with `-O 4 0 13 0 22 0`. For a pilot, use `-r 1` and a separate filename.

The safe runner accepts the same switch:

```bash
python -u RunSafeWeightedStatic.py --formulation escortflow \
  -x 13 -y 7 -O 6 0 -l 4 -e 8,12,16 -r 1-100 -m leave \
  --lp --lp-protocol v4 --threads 16 -f lp_four_safe_escortflow_13x7.csv
```

`--lp-protocol v5` is the default for new campaigns, including the two- and
six-target campaigns. It uses the refined sufficient coefficient and common
physical horizon from the current safe integer runner. Change `-m leave` to
`-m continue` for the corresponding continue experiments. The LP time limit
is `-t` in the regular runners and `--weighted-time-limit` in the safe runner;
it defaults to 300 seconds. LPs run one at a time by default. Set `--lp-workers`
to change this, and `--lp-retry-time-limit` to change the 600-second barrier
retry. The integer extension budget is not used by `--lp`.

### Match the coefficient and horizon protocol

Let `D` be the sum of targets' nearest-output Manhattan distances, `d_max` the
largest such distance, `F_g` and `T_g` the greedy flow time and makespan, and
`K` the initial number of loads. Define `H_g=F_g-D+d_max` and `U_g=K*H_g`.

| Protocol | Coefficient R | EF physical horizon | LF physical horizon |
| --- | --- | --- | --- |
| v4, archived one-/four-target runs | `U_g+1` | `H_g+1` | `max(H_g,T_g+1)` |
| v5, current safe runner | `U_g-D+1` | `max(H_g,T_g+1)` | `max(H_g,T_g+1)` |

Both coefficients are sufficient for the same integer objective priorities.
Changing the coefficient or horizon can change the LP optimum and the progress
of a time-limited integer search. Use the protocol of the integer results when
calculating their LP gaps. The archived v4 formulation horizons can differ by
one period; retain that difference when reproducing the archived model.

An explicit `--flow-weight R` (alias `--flow_weight`) and `--horizon T` in the
regular runners override their generated values. `R` must satisfy the sufficient
coefficient bound. The physical horizon is `T+1` in EF and `T` in LF. For example,
the archived `13x7`, three-escort, one-target, seed-1 instance uses `R=881`,
EF `T=10`, and LF `T=11`. Use `RunStaticLP.py` with the archived inputs and
frozen model sources for exact replay of published results. Direct generation
uses the working repository's source at the start of the command.

### Output and restart

Choose a separate `-f` filename for each LP batch. An existing output is refused
unless `--resume` is supplied. The CSV contains only `OPTIMAL` LP results, with
the same schema as `RunStaticLP.py`. It is separate from the integer results.
Beside `name.csv`, the direct runner saves:

- `name.csv.instances.csv`: generated coordinates, `R`, horizons, and greedy reference.
- `name.csv.generation.json`: parameters and input/source checksums.
- `name.csv.source/`: frozen generation and LP model code.
- `name.csv.manifest.json`, `name.csv.logs/`, and, if needed, `name.csv.failures.json`: replay progress and diagnostics.

Repeat the same command with `--resume` to skip validated completed LPs and
retry missing cases. It uses the saved inputs and sources, including their
original `R`, even if the working code has changed. Changed instance parameters
or modified snapshots are refused; use a new filename for a different campaign.

For historical fixed-weight LP runs, the regular EF/LF runners accept
`--lp --legacy-lp --gamma 0.01`. This selects their older weights, heuristic
horizon, and CSV format. The historical shell scripts select that mode
explicitly. Integer commands without `--lp` retain their existing behavior.

## Verification

```bash
python -m unittest test_integrated_lp test_generated_lp test_static_lp test_table2b_lp
```

The tests cover malformed and duplicate input records, objective/status checks,
changed-source and changed-instance resume protection, retry behavior, and the
direct generation without integer CSVs, archived coefficient/horizon matching,
and completed-result resume. With Gurobi available, they also compare both
standard runners against the archived seed-1 LP objectives and solve small
direct LP batches. The separate
eight-case replay verifies both formulations on all four Table 3 layouts.
