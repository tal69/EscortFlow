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
Table 3 inputs match these physical horizons across formulations. Supplying the
same numerical `T` to both formulations would change the comparison.

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
and objective components before reusing saved results. Changed inputs or model
sources require a new output. Existing files are not overwritten.

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
escorts for two targets and 8/12/16/20 escorts for six targets. It first runs
5,600 integer searches with the 300-second comparison cutoff and conditional
300-second extension, checks paired instances, and then runs 5,600 LP solves
from the same frozen model source tree and recorded coefficients/horizons.

After the existing four-target campaign has finished, in the Linux `Code` checkout:

```bash
git pull --ff-only
python3 -u RunTable2Targets.py --threads 16 --lp-workers 1 --lp-threads 16
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

The LP stage saves `lp_results.csv`, its manifest, worker logs, and the complete
restart command `lp_commands.sh`. If only the LP stage fails or is interrupted,
use `bash /path/to/results_table2b_targets_.../lp_commands.sh`. This skips already
verified optimal LP values. It does not restart the integer campaign.
See [README.md](README.md) for settings, outputs, and progress monitoring.

## Switches in the standard static runners

Both `EscortFlowStatic.py` and `LoadFlowStatic.py` retain `--lp`. The new option
`--flow-weight R` (alias `--flow_weight`) sets `beta=1` and `gamma=1/R` explicitly.
It requires a positive integer `R` and a Gurobi LP run. It cannot be combined
with other objective weights or an objective cutoff. Existing integer runs and
the historical `--lp --gamma 0.01` commands retain their behavior.

For example, the frozen `13x7`, three-escort, one-target, seed-1 instance uses
`R=881` and physical horizon 11:

```bash
python3 EscortFlowStatic.py -x 13 -y 7 -O 6 0 -e 3 -l 1 -r 1 -m leave \
  --lp --flow-weight 881 --horizon 10 --num_threads 1 -t 300 \
  -f example_escort_lp.csv

python3 LoadFlowStatic.py -x 13 -y 7 -O 6 0 -e 3 -l 1 -r 1 -m leave \
  --lp --flow-weight 881 --horizon 11 --num_threads 1 -t 300 \
  -f example_load_lp.csv
```

Use the replay runner for exact table reproduction. Standard runners generate
instances from seeds and otherwise use their usual horizon selection; they do
not automatically load archived per-instance weights or horizons.

## Verification

```bash
python3 test_static_lp.py
```

The tests cover malformed and duplicate input records, objective/status checks,
changed-source and changed-instance resume protection, retry behavior, and the
new standard-runner options. With Gurobi available, they also compare both
standard runners against the archived seed-1 LP objectives. The separate
eight-case replay verifies both formulations on all four Table 3 layouts.
