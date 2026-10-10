# Escort Flow Optimization and Simulation

This repository accompanies two working papers by Tal Raviv and Yossi Bukchin:

- "Escort-Flow Formulation for Simultaneous Multi-Load Retrieval in Puzzle-Based Storage"
- a companion paper on rolling-horizon control for dynamic PBS retrieval

This repository contains optimization and simulation code for retrieval control in a puzzle-based storage (PBS) system with escorts. The main workflow is:

1. run a formulation-paper benchmark instance with `EscortFlowStatic.py` or `LoadFlowStatic.py`
2. run a dynamic rolling-horizon simulation with `EscortFlowSim_v8.py` and collect steady-state statistics directly into a CSV row
3. optionally save a raw simulation trace with `-a/--save_raw` and inspect it with `CI_Calculation.py` or `PBSAnimation.py`

The project currently uses `EscortFlowSim_v8.py` as its rolling-horizon dynamic simulator, with Gurobi accessed through the Python API.

## Main entry points

- `EscortFlowSim_v8.py`: primary simulator for dynamic request arrivals, rolling-horizon control, hybrid MILP/greedy policy, CSV reporting, and optional raw pickle export, using Gurobi directly from Python
- `EscortFlowStatic.py`: static escort-flow experiment runner for single-load and multi-load instances, defaulting to the Gurobi Python backend and also supporting greedy-only and naive-lower-bound modes
- `LoadFlowStatic.py`: static load-flow experiment runner; BM is the default, `--lm` switches to LM, and the Gurobi Python API is the default backend
- `RunStaticLP.py`: reproducible continuous LP replay from recorded coordinates, weights, and horizons, including archived Table 3 inputs
- `RunTable2bLP.py`: Mac Studio LP-bound runner for four-target Table 2(b), with automatic resume and relative-gap summaries
- `RunFourTargetLP.py`: all eight four-target LP batches directly from seeds, independently of unfinished integer results
- `ReproducePaper.py`: one command for all revised formulation-paper numerical tables, including leave/continue retrieval, both formulations, LP relaxations, and CSV/LaTeX summaries
- `EscortFlowStaticLex.py` and `LoadFlowStaticLex.py`: two-phase integer optimization, first flow time, then load movements at the best flow time found
- `CI_Calculation.py`: post-process one or more raw pickle files and compute steady-state means and confidence intervals using MSER-5 warmup deletion and batch selection
- `PBSAnimation.py`: animate PBS outputs from `EscortFlowStatic.py`, `LoadFlowStatic.py`, and `EscortFlowSim_v8.py`
- `OneStepHeuristic_v2.py`: greedy one-step escort heuristic used for static feasible solutions and horizon estimates, dynamic greedy control, and fallback when dynamic MILP results are rejected; see [One-step heuristic for static retrieval](#one-step-heuristic-for-static-retrieval)

`EscortFlowSim_v5.py` is retired, and the old `v5` file has been moved out of the repository into `Junk/` for local reference only.

## Requirements

Python:

- Python 3
- `numpy` via [requirements.txt](/Users/talraviv/Library/CloudStorage/Dropbox/research/PBS/EscrotsFlow/Code/requirements.txt)
- `tkinter` for animation
- standard-library modules such as `argparse`, `pickle`, `subprocess`, and `time` are used directly and do not need separate installation

Install the Python dependency with:

```bash
python3 -m pip install -r requirements.txt
```

Optimization:

- Gurobi with the Python API installed for all paper replication commands below
- IBM ILOG OPL / CPLEX with `oplrun` available on `PATH` only if you want to use the legacy `--opl` path

The legacy OPL paths expect the OPL model files in this repository, especially:

- `escort_flow_bm_rh_v3.mod`
- `escort_flow_bm_rh_static_v3.mod`

The static experiment scripts additionally depend on:

- `escort_flow_static_gurobi.py`, `escort_flow_static_lazy.py`, and `escort_flow_static_bnc.py` for the escort-flow Gurobi backends
- `load_flow_static_gurobi.py` for the load-flow Gurobi backend
- `pbs_escorts_bm_v3.mod` and related escort-flow OPL files only when the legacy `--opl` path is used
- `pbs_load_flow_multi.mod` and `pbs_load_flow_multi_lp.mod` only when the legacy `--opl` path is used
- `PBSCom.py` for common instance-generation and formatting helpers
- `PBS_DPHeuristic_bm.py` and `PBS_DPHeuristic_lm.py` when a DP-based upper bound is requested

For the DP-assisted single-load examples, large DP table files are also required:

- `BM_set10x10_k3_corner.p`
- `BM_set16x10_k3_2centers.p`

Those tables are used to derive an upper bound on the number of time steps for the single-load case. They are optional for the paper replication block below because the current formulation-paper commands use the greedy upper bound instead.

Download source for the large DP files:

- Dropbox archive: <https://www.dropbox.com/scl/fi/ug5eenojzkh4riv8ja1rc/ESCORTS.zip?rlkey=2r01my7aeetl5q0zzn6evzvtt&st=25y3rg7d&dl=0>

Place the extracted `.p` files in the project directory when using the DP-assisted single-load examples later in this README.

## Paper replication guide

Run all commands from the `Code/` directory:

```bash
cd /Users/talraviv/Library/CloudStorage/Dropbox/research/PBS/EscrotsFlow/Code
python3 -m pip install -r requirements.txt
```

The scripts append to CSV files. For a fresh replication, remove or rename existing CSV outputs before running the commands.

The paper experiments were developed with Python 3.11, Gurobi 13.0, and NumPy. Static benchmark runs in the formulation paper used a 300-second solver limit per instance. Dynamic runs use the request counts and seeds hard-coded in the shell scripts listed below.

### Formulation paper: retrieval benchmark tables

For the revised paper, use the complete runner below. The older shell scripts
in the following historical guide are retained for earlier experiments.

### Four-target LPs from seeds on the Mac Studio, while Linux is still running

For the currently running four-target campaign, activate your working Conda
environment and use **`python`**, which names its interpreter on your Mac Studio:

```bash
python -u RunFourTargetLP.py
```

Run this from `Code/`, including in a tmux window. No arguments are needed. It
runs **all eight parts** (13x7, 10x10, 16x10, 27x10, each with escort-flow and
load-flow), four target loads, 8/12/16 escorts, seeds 1-100, leave mode. That is
**2,400 LPs for 1,200 matched physical instances**. It uses the same seed and
greedy generator as the main runners' `--lp` option. It explicitly preserves
the current four-target campaign's **v4 coefficient and formulation-specific
horizons**, rather than the v5 settings of the next two-/six-target campaign.

Copied files in `Experiment Oct2026/table2b_*.csv` are optional verification
inputs. The script checks every available selected row's coordinates, R and
horizon against the generated instance before launching LPs. Missing files and
missing seeds do **not** reduce the run. It never writes to the integer CSVs.
Use `--input-dir /path/to/copied/results` to check another copied folder.

Results go into a separate **`results_four_target_lp/`** folder. It contains
`lp_results.csv` with all available optimal LP values, eight individual CSVs in
`parts/`, `coverage.json`, frozen generated inputs in `instances/`, and a frozen
`source/` tree with checksums in `campaign.json`. Each successful LP is saved
immediately. **Repeat the identical command to resume** completed parts and
retry missing LPs. Source updates do not change an existing run: it resumes
using its original frozen generator and models. Changing the seed/layout/escort
selection requires a separate `--output-dir`.

To check the plan, check the environment, or run a two-LP pilot:

```bash
python RunFourTargetLP.py --dry-run
python RunFourTargetLP.py --check-environment
python -u RunFourTargetLP.py --layouts 13x7 --escorts 16 --seeds 1 \
    --output-dir results_four_target_lp_pilot
```

For a seed range or explicit solver settings:

```bash
python -u RunFourTargetLP.py --seeds 1-100 --workers 1 --threads 16
```

The defaults run one LP at a time with the Mac's performance-core count for
threads, a 300-second solver limit, and a 600-second barrier retry after a time
limit. Only `OPTIMAL` results supply an LP lower bound. Python 3.10+, NumPy,
`gurobipy` and your full Gurobi license are required. The script prints the
actual Python path/version and uses that same interpreter for every child and
worker process. If NumPy is missing, install it with `python -m pip install numpy`
in this environment. LP values are generated now; percentage-gap table summaries
can be calculated when the complete integer results have been copied.

### Table 2(b): replay only the available copied integer results

Set up once on the Mac Studio, from `Code/`:

```bash
bash RunTable2bLP.sh --setup python3.13
```

This creates a dedicated environment at `~/.venvs/escortflow-table2b`, installs
`gurobipy==13.0.3`, checks your Gurobi license beyond the pip package's size limit,
and remembers the exact interpreter in `~/.config/escortflow/table2b-python`
(or `$XDG_CONFIG_HOME/escortflow/table2b-python`). It starts no paper experiments.
If `python3.13` is not on PATH, pass its absolute path instead. The LP runner
does not need NumPy or the older NumPy pin in `requirements.txt`.

Alternatively, keep an existing working environment without installing anything:

```bash
bash RunTable2bLP.sh --set-python /absolute/path/to/licensed/python
```

Thereafter, in every new terminal or tmux window, start or resume with:

```bash
bash RunTable2bLP.sh
```

No activation is needed. The launcher uses the saved interpreter, including for
worker processes and resumes that execute a frozen older runner. An installed
Python 3.13 does not guarantee that a command named `python3` uses it. The runner
now prints its actual interpreter path and version and distinguishes version,
module-import, and license failures. Check only the environment with
`bash RunTable2bLP.sh --check-environment`. A valid academic or commercial Gurobi
license must already be configured; installing `gurobipy` alone supplies only a
size-limited license. License failures must be resolved before setup is saved.

The underlying `python3 -u RunTable2bLP.py` command remains supported when
`python3` already names the correct environment. The script reads
`Experiment Oct2026/` beside the
script. It calculates both formulations' continuous LP relaxations for the
four-target leave-mode runs, seeds 1-100, 8/12/16 escorts, and all four paper
grids. A complete campaign has 2,400 LP solves. Missing load-flow files and
unfinished groups can be supplied later. To generate all requested LP instances
before those CSV rows arrive, use **`python -u RunFourTargetLP.py`** as described
above, or the individual parameter-based `--lp` commands in
[Direct LP runs from parameters and seeds](#direct-lp-runs-from-parameters-and-seeds).

Results go into the separate **`results_table2b_lp/`** folder. Each optimal LP
is saved immediately. **Repeat the same command to resume**, retry failed cases,
or process newly copied integer rows. Saved optimal values are reused. The
output folder preserves model sources and input snapshots for reproducibility.
Earlier instances must retain their recorded coordinates, coefficient, and
horizon. Use a separate output folder for a different campaign or seed/layout/
escort selection.

Each relaxation uses the exact recorded `R` and each formulation's recorded
horizon, including the archived four-target horizon difference. It minimizes
`FT + MV/R`, equivalent to `(R*FT + MV)/R`, with all integer variables continuous
and flow time free to optimize. Greedy solutions are read from the existing
files. The historical movement weight of 0.01 is not used.

The default runs **one LP at a time**, with the Mac's performance-core count as
the solver thread count, to limit memory use for the large load-flow LPs. Each
solve has a 300-second limit, followed by a barrier retry without crossover
for up to 600 seconds if it reaches the time limit. Only `OPTIMAL` results enter
the CSV. Settings can be changed on a resumed run:

```bash
bash RunTable2bLP.sh --threads 20 --time-limit 300 --retry-time-limit 600
```

Check the input coverage without solving or creating output files:

```bash
bash RunTable2bLP.sh --dry-run
```

For a two-LP pilot, one per formulation, use a separate output folder:

```bash
bash RunTable2bLP.sh --layouts 13x7 --escorts 16 --seeds 1 \
  --output-dir results_table2b_lp_pilot
```

Use `--input-dir "/path/to/copied/results"` for a different input folder and
`--output-dir "/path/to/LP/results"` for a different destination. The input
filenames are `table2b_escortflow_13x7.csv`, `table2b_loadflow_13x7.csv`, and their
counterparts for the other grids. Quote paths containing spaces.

The output folder contains:

- `lp_results.csv`: optimal LP objectives, fractional FT/MV, recorded `R` and
  horizons, instance fingerprints, and solver versions.
- `lp_gaps.csv`: per-instance percentages against the common lexicographically
  best feasible FT/MV pair from both formulations, including their extensions.
- `table2b_lp_summary.csv`: formulation-specific LP-gap means for Table 2(b).
  Incomplete groups stay blank.
- `table2b_lp_columns.tex`: a compact LaTeX fragment for the two LP-gap columns,
  with the selected seed count stated in a comment.
- `coverage.json`: available integer records, optimal LP records, and complete
  groups. `source/`, `inputs/`, manifests, and worker logs preserve the audit.

The gap is `100*(FT_BK + MV_BK/R - Z_LP)/(FT_BK + MV_BK/R)`. Percentages are
computed per instance before averaging over all selected seeds, including zero
gaps. A zero reference objective contributes zero. A summary cell requires all
100 seeds by default, or the explicit smaller count in a pilot. Best-known
references may improve as more integer results arrive. Run
`bash RunTable2bLP.sh --summarize-only` to refresh summaries using saved LPs
without further solves.

Run the command inside a tmux window on the Mac Studio to keep it running
through terminal disconnections. See [LP_REPRODUCIBILITY.md](LP_REPRODUCIBILITY.md)
for the replay conventions. Validation:

```bash
python3 -m unittest test_table2b_lp test_static_lp.LPReplayTests
```

### Revised paper: reproduce all numerical tables with one command

From the `Code/` directory, or an extracted reproducibility package containing
the files listed below, run:

```bash
python3 -u ReproducePaper.py
```

There are no required arguments. This runs **seeds 1-100**, all four paper
layouts, both load-flow and escort-flow, and all leave/continue experiments.
It calculates the integer results, conditional extensions, both continuous LP
relaxations, the method-comparison tables, and the bounds/solution/improvement
tables. The literature-review table is editorial content and needs no solver
experiment. No dynamic-paper experiments are included.

To reproduce the same tables over an inclusive subset of seeds:

```bash
python3 -u ReproducePaper.py 1-10
```

The equivalent form is `python3 -u ReproducePaper.py --seeds 1-10`. A single
seed, such as `1`, and a list, such as `1,7,20`, are also accepted. The script
uses the selected number of instances as the denominator in every percentage
and average, and records it in the output. It never assumes that a pilot has
100 instances.

**Requirements.** Use Python 3.10 or newer (3.11 recommended), NumPy, and
`gurobipy` with a working Gurobi license large enough for these models. Install
the packages in your solver environment with:

```bash
python3 -m pip install -r requirements.txt
python3 -m pip install gurobipy
```

The launcher checks the environment and license before solving. The models
use the Gurobi Python API; OPL, CPLEX, DP pickle files, and LaTeX are unnecessary
for running or generating the table fragments. For a separate solver environment,
use `--python /path/to/python`. The interpreter and solver versions are saved.

**Experimental coverage.** The four layouts and outputs are:

| Grid | Outputs | Escorts for approximately 70% occupancy |
| --- | --- | ---: |
| 13x7 | (6,0) | 27 |
| 10x10 | (0,0) | 30 |
| 16x10 | (4,0), (11,0) | 48 |
| 27x10 | (4,0), (13,0), (22,0) | 81 |

| Retrieval mode | Targets | Escort counts in each layout |
| --- | ---: | --- |
| leave | 1 | 3, 4, 5, 6, 7, 8 |
| leave | 2 | 8, 12, 16 |
| leave | 4 | 8, 12, 16, plus the approximately 70%-occupancy count |
| leave | 6 | 8, 12, 16, 20 |
| continue | 2, 4 | 8, 12, 16, plus the approximately 70%-occupancy count |
| continue | 6 | 8, 12, 16, 20, plus the approximately 70%-occupancy count |

This includes the existing leave benchmarks, `Run70Percent.py` and
`RunTable2Targets.py` configurations, and `RunContinue.py` configurations.
Six-target continue cases also include 20 escorts for comparison with leave
mode. There are **120 configuration rows**, **12,000 paired instances**,
**24,000 integer searches**, and **24,000 LP solves** at the default seed range.
Each integer instance is followed immediately by its LP relaxation. Both results
are saved in the same CSV row. All jobs run sequentially to avoid competing for memory.
Run the full campaign inside tmux; it is a substantial computation.

To inspect every command without creating files or requiring a solver license:

```bash
python3 ReproducePaper.py --dry-run
```

For a smaller pilot, use `python3 -u ReproducePaper.py 1 --layouts 13x7`.
The integer search defaults to up to 16 threads on Linux and the performance
core count on Apple Silicon. Set `--threads 16` explicitly for the Linux
benchmark. LPs default to one thread; `--lp-threads 16` changes this. Both
thread counts are recorded. Each integer search has a 300-second initial
solver cutoff and a conditional 300-second continuation of the same search
when FT remains unproved. Each LP has a 300-second limit and a 600-second
barrier retry after a time limit. Only `OPTIMAL` LP values enter the tables.
The integer CSV includes `lp_relaxation_lower_bound` in `FT+MV/R` units,
`lp_status`, and `lp_elapsed_seconds`, plus LP settings and source/input
fingerprints. Integer solver and CPU times exclude the separate LP solve;
`total_wall_time` records the integer portion, including greedy preparation.
The time-limit options shown by `--help` are intended for validation pilots;
changing them changes the experimental protocol.

**Results and restart.** Each launch creates a fresh
`results_paper_<timestamp>_<pid>/` directory beside the script. You can choose
its location with `--output-dir /path/to/new/results`. Existing directories
are rejected for a fresh launch. The directory contains:

- `source/`: the exact source files and documentation frozen before solving.
- `manifest.json` and `preflight.log`: seeds, configurations, settings, source
  hashes, environment, license check, run status, and restart sessions.
- `leave/` and `continue/`: `loadflow.csv`, `escortflow.csv`, `lp_results.csv`,
  per-configuration integer CSVs in `parts/`, solver logs, and LP manifests.
- `tables/instance_metrics.csv`: auditable metrics for every instance and method.
- `tables/method_comparison.csv`: the Table 2 metrics and new continue equivalents.
- `tables/bounds_and_solutions.csv`: the Table 3 metrics and new continue equivalents.
- `tables/table2_leave.tex` and `tables/table3_leave.tex`: single-target (a)
  and four-target (b) panels together, each pair with one shared caption.
- `tables/table{2,3}_{leave,continue}_{2,4,6}targets.tex` and single-target
  fragments: individual panels for inserting into the revised paper.
- `tables/audit.json`: input hashes, coverage checks, and exact metric definitions.

The LaTeX fragments preserve portrait headings, spell out Load-flow and
Escort-flow, use FT/MV and esc., omit utilization from the method comparison,
and bold the lower displayed time for each method pair. The single-target
Table 2 panel uses only 3-6 escorts; Table 3 retains 3-8 escorts. FT/MV gap
cells are empty at 100% corresponding optimality. LP gaps remain displayed.
Use `booktabs`, `makecell`, and `subfig` when inserting the fragments into LaTeX.
The script does not edit or compile the manuscript.

After interruption, resume with:

```bash
python3 -u ReproducePaper.py --resume /path/to/results_paper_...
```

Resume uses the saved seeds, settings, and source snapshot. It verifies their
hashes, keeps complete integer rows, runs only missing seeds, rebuilds merged
CSV inputs without duplicates, and skips previously saved optimal LPs. An
error row or truncated CSV is reported for inspection instead of silently
discarded. After a fully completed run, regenerate just the tables with:

```bash
python3 ReproducePaper.py --tables-only /path/to/results_paper_...
```

This last command uses only the recorded CSVs and Python's standard library;
it needs neither Gurobi nor a solver license.

**Objective and calculations.** New runs use the current sufficient integer
coefficient `R = K*(F_g-D+d_max)-D+1` in `Q = R*FT+MV`, with `K=N-e` initial
loads and the shared greedy warm start. Each LP reuses that run's actual
coordinates, coefficient, and physical horizon and minimizes the equivalent
scaled objective `FT+MV/R`.

Method-specific optimality rates, integer bounds, and times use only the
initial cutoff. Best-known FT and MV come from one lexicographically best
feasible pair across the greedy plan, both formulations, and both initial and
extension results. Component gaps are `100*(best-known - integer LB)/best-known`;
the MV bound is conditional on best-known FT. LP gaps are
`100*(FT_BK+MV_BK/R-LP)/(FT_BK+MV_BK/R)`. Shared FT/MV improvements are
`100*(greedy-best-known)/greedy`. These percentages are calculated per
instance and then averaged over **all** selected instances, including zeros.
The movement improvement can be negative when attaining a smaller FT requires
more movements. A zero reference contributes zero percent when the matching
component/bound is also zero.

To include previously recorded compatible integer runs in the best-known
reference, supply their CSVs on a fresh launch:

```bash
python3 -u ReproducePaper.py --reference-input old_loadflow.csv old_escortflow.csv
```

References must be safe-weighted v4/v5 integer CSVs with explicit feasible
FT/MV components. Their coordinates, target count, escorts, seed, retrieval
mode, and movement regime must match. They are copied into the results folder
and hashed. They affect shared best-known solutions and gap/improvement
references; they do not replace the new method-specific cutoff measurements.
Runs from another retrieval mode or target count are never pooled.

This is a fresh reproduction with the current refined coefficient. The
completed single-target tables used the earlier sufficient `R=K*H_g+1`;
their archived LP replay must retain that recorded coefficient, as described
in [LP_REPRODUCIBILITY.md](LP_REPRODUCIBILITY.md). Hardware, solver version,
thread count, and time-limited search can change runtimes and best-known
solutions, so a fresh run is not guaranteed to reproduce every printed digit
of an earlier run. Retain the recorded raw CSVs and their frozen sources in
the final reproducibility package for auditing the reported values.

**Files for the reproducibility package.** Include `ReproducePaper.py`,
`PaperTables.py`, `RunStaticCampaign.py`, `RunSafeWeightedStatic.py`,
`RunWeightedStatic.py`, `RunStaticLP.py`, `static_generated_lp.py`, `static_integrated_lp.py`,
`PBSCom.py`, `OneStepHeuristic_v2.py`,
`static_lexicographic.py`, `static_weighted_certification.py`,
`static_safe_weighted_search.py`, `escort_flow_static_gurobi.py`,
`load_flow_static_gurobi.py`, `requirements.txt`, `README.md`, and
`LP_REPRODUCIBILITY.md`. These exact files are collected automatically in
each results directory's `source/` folder. Keep any archived reference CSVs
with the pack when historical best-known solutions are used.

### Historical formulation-paper runners

The legacy weighted benchmark campaign in "Escort-Flow Formulation for Simultaneous Multi-Load Retrieval in Puzzle-Based Storage" is generated by `SingleLoadStatic.sh` and `FourLoadsStatic.sh`. These two scripts run both formulations and both ILP/LP-relaxation variants, so no separate LP-relaxation scripts are needed. These scripts retain an extra `9 x 5` layout that is not in the current manuscript's Table 2. Their scope is:

- Single target: five layouts, 3--8 escorts, 100 random instances per row.
- Four targets: five layouts, 8, 12, and 16 escorts, 100 random instances per row.

For the fixed `F+0.01*M` experiment and separate flow-time check matching
Table 2(a,b), use `RunTable2Weighted.sh`, described below. The additional
`RunTable2SafeWeighted.sh` version uses sufficient integer weights and continues
the same search only when the initial result cannot already prove minimum flow
time. `RunTable2Lex.sh` retains the two-phase lexicographic experiment.

Use this block for a clean formulation-paper replication:

```bash
cd /Users/talraviv/Library/CloudStorage/Dropbox/research/PBS/EscrotsFlow/Code
rm -f table1a_*.csv table1b_*.csv

bash SingleLoadStatic.sh
bash FourLoadsStatic.sh
```

Reference outputs from previous runs are stored under:

- `/Users/talraviv/Library/CloudStorage/Dropbox/research/PBS/EscrotsFlow/Experiment_Static_May2026/`
- `/Users/talraviv/Library/CloudStorage/Dropbox/research/PBS/EscrotsFlow/Experiment_Mar2026_take2/`

The historical LP commands in these shell scripts explicitly select
`--lp --legacy-lp` and use `F+0.01*M`. For the LP gaps in Table 2, use the
direct parameter-based `--lp` commands below, or `RunStaticLP.py` with the
recorded per-instance `R` and horizons. The separate
[LP reproducibility guide](LP_REPRODUCIBILITY.md) gives the full command and the
self-contained `reproducibility/table3_lp.zip` package.

### Two-phase lexicographic static retrieval

The new `EscortFlowStaticLex.py` and `LoadFlowStaticLex.py` entry points reuse the
existing formulations, instance seeds, and default horizon selection. The
original entry points retain weighted optimization by default; adding
`--lexicographic` selects the same two-phase procedure. The new entry points
default to separate `res_escort_flow_lex.csv` and `res_load_flow_lex.csv` files.

1. Minimize the unweighted, integer total flow time.
2. Add an equality fixing flow time to the best integer incumbent from phase 1,
   then minimize the unweighted number of load movements. An escort block move
   counts one movement per load shifted.

Both phases use `MIPGap=0` and `MIPGapAbs=0.999`. The small margin below one avoids
premature stopping at the numerical boundary: the original setting of exactly
one produced some Gurobi `OPTIMAL` results whose integer incumbent and reported
bound still differed by one. Certification separately requires the integer
incumbent-to-bound gap to be below one with an additional `1e-6` numerical margin.
This prevents a bound microscopically above an integer from falsely certifying a
one-unit gap, including on time-limited solves. No weighted objective cutoff is
carried into either phase. See the
[Gurobi gap parameter definitions](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#mipgapabs).

Use two time limits, in seconds:

- `--phase1_time_limit`: cap on the flow-time phase, also capped by the total limit.
  If omitted, phase 1 can use the entire total budget.
- `--time_limit` (alias `--total_time_limit`, or `-t`): total solver budget for both
  phases. Phase 2 receives this limit minus the actual Gurobi `Runtime` of phase 1,
  rather than minus its configured cap. For example, with a 100-second first-phase
  cap and a 300-second total limit, a 40-second first solve leaves 260 seconds.

The total defaults to 300 seconds, except when `--work_limit` is supplied without
a time limit. Work is also shared across both phases. Solver runtime excludes
instance generation, heuristic evaluation, model construction, and result export;
`Wall Clock Time` additionally includes construction and result extraction.
Gurobi time limits can have a small termination overhead.

If phase 1 stops at a limit with an incumbent, phase 2 still fixes that best
attained flow time. Such a run is only certified lexicographically optimal when
both objective gaps prove optimality. If phase 1 has no incumbent, or the shared
time/work budget is exhausted, phase 2 is skipped. The phase-one feasible plan is
retained if phase 2 returns no incumbent. A numerical failure or user interruption
in phase 1 also prevents starting phase 2.

Examples from `Code/`, using an interpreter with `gurobipy` and a working license:

```bash
python3 EscortFlowStaticLex.py -x 5 -y 5 -O 0 0 -e 4 -l 2 -r 1-10 \
  -m leave --phase1_time_limit 100 --time_limit 300 --num_threads 8 \
  -f escort_lex.csv

python3 LoadFlowStaticLex.py -x 5 -y 5 -O 0 0 -e 4 -l 2 -r 1-10 \
  -m leave --phase1_time_limit 100 --time_limit 300 --num_threads 8 \
  -f load_lex.csv
```

Escort-flow supports `stay`, `leave`, and `continue`, including `--lazy`, `--bnc`,
and the existing compatible `--warmstart` options. Load-flow supports `leave`,
with BM by default and `--lm` available. Lexicographic runs reject `--lp`, `--opl`,
`--cutoff`, and nondefault objective weights. Objective weights do not enter the
two optimizations; the load-flow CSV retains their default values as legacy
instance metadata.

The lexicographic CSV records each phase's status, bound, absolute gap, solver
runtime, work, and proof flag, plus both time caps and the computed phase-two
limit. `Lexicographic Optimal=1` requires proof in both phases. A movement optimum
at an unproven flow time has `Solver Status=FLOWTIME_NOT_PROVEN`. If Gurobi returns
`OPTIMAL` for phase 2 but its integer gap does not prove optimality, the overall
status is `MOVEMENTS_NOT_PROVEN`; the raw phase status remains available. Flow-time and
movement bounds are reported separately in their own units. Use separate CSVs
for weighted and lexicographic runs because their result columns differ.

Optimality refers to the selected finite planning horizon. The existing automatic
heuristic horizon is unchanged for comparison with earlier experiments; for
multiple targets it is not a proof that the unrestricted optimum is included.
Use `--horizon T` in either weighted or lexicographic runs to hold the horizon
fixed within a formulation. `T` is the last indexed decision period, with indices
`0,...,T`; escort-flow can arrive at time `T+1`, while load-flow records retrieval
through time `T`. Account for these conventions when comparing formulations.
The larger sufficient horizon discussed below remains relevant when unrestricted
optimality is required.

Regression checks cover small independently verified optima across all four
backends, first-step arrivals, zero targets, infeasibility, strict gap stopping,
budget allocation, callbacks, and retained incumbents:

```bash
python3 test_static_lexicographic.py
```

### Sufficient integer weights with flow-proof timing

Run the new version on Linux, including inside an existing tmux session:

```bash
bash RunTable2SafeWeighted.sh --threads 16
```

This runs the same Table 2(a,b) instances and both formulations sequentially.
Defaults are 16 threads, a 300-second first phase, and a conditional extension
of up to 300 more seconds. There is one optimization call and a complete common
greedy warm start. The solver records the first observed flow-time proof while
the first phase continues toward proving both objectives. If flow time remains
unproved at the cutoff, the extension attempts to certify that saved solution.

For each instance, let `F_g` be the greedy solution's flow time, `D=sum(d_i)`
the sum of the targets' nearest-output Manhattan distances, `N` the number of
cells, and `e` the initial number of escorts. The program sets

```text
K   = N-e
H_g = F_g-D+max(d_i)
U_g = K*H_g
R   = U_g-D+1
```

and minimizes the integer objective `R*F+M`. Every feasible plan has `M >= D`,
and a minimum-flow, minimum-movement plan has `M_star <= U_g`. Consequently
`R > U_g-D` guarantees that a weighted optimum is globally lexicographically
optimal. The feasible representative implies `U_g >= D`, so `R >= 1`.
The zero-flow case uses `D=H_g=U_g=0` and `R=1`. The weighted model covers the
sufficient physical horizon `H_g`, which can be longer than the old
greedy-makespan horizon. The same instance has the same coefficient in both
formulations. The coefficient and both movement bounds are saved per row:
`safe_movement_lower_bound=D` and `safe_movement_bound=U_g` (the latter retains
its meaning as an upper bound, not the range `U_g-D`).
Both formulations use physical horizon `H=max(H_g,C_g+1)`, retaining the
complete greedy trace. Their array indices differ: EF uses `T=H-1`, whereas
LF uses `T=H`. The CSV records both the index and the physical horizon.

The weighted solve uses `MIPGap=0` and `MIPGapAbs=0.999`:

1. An independently verified integer gap below one proves both objectives.
2. During the search, check flow optimality whenever the observed weighted
   lower bound improves or a new candidate has smaller flow. Unchanged bounds
   and movement-only improvements skip the proof comparison. Check the reporting
   cutoff and final result unconditionally.
   The incumbent-specific movement bound can prove flow time while movement
   optimality remains unresolved. A matching analytical distance bound also
   proves flow time.
3. Save the first observed proof time, flow value, node count, bound, and proof
   source. A flow proof alone does not stop the first phase. At the cutoff,
   freeze the incumbent observed by that time. Continue only if its flow is
   unproved, stopping on its proof, a lower-flow counterexample, or the total
   600-second cap. `MIPFocus=0` remains unchanged throughout.

The flow-only gap certificate retains the criterion
`gap < R+M_w-U_minus`, where `gap=R*F_w+M_w-L_Q`,
`H_minus=F_w-1-D+max(d_i)`, and `U_minus=K*H_minus`.
For v5 the threshold simplifies to `1+M_w-D+K*(F_g-F_w+1)`.
The model must cover `H_minus`; when `F_w-1 < D`, the distance bound already
proves minimum flow. The existing numerical margins and the independent
full-integer-gap certificate below one remain unchanged.

`--weighted-time-limit` controls the first-phase limit and reporting cutoff
(default 300 seconds). Proof checks are triggered by bound or flow changes,
with no time or node interval setting.
`--extension-time-limit` supplies up to 300 additional solver seconds by default,
used only when the saved candidate's flow remains unproved.
`--certification-time-limit` remains an alias for that option. The snapshot
uses incumbents and bounds observed by the reporting cutoff, never a later
solution retroactively. If no incumbent is available by that cutoff, the row
explicitly records that condition and stops without inventing a candidate to
certify. Final solutions remain separately reported.

Results are saved under a new `results_table2_safe_weighted_<timestamp>_<pid>/`
directory, with per-layout CSVs in `parts/`, solver logs in `logs/`, and merged
CSVs in the top level. The main solution columns describe the initial-cutoff
candidate. Separate final-solution columns show the result after any extension,
with separate bounds, proof flags, and elapsed times. A later improvement never
overwrites the main experimental result. The original candidate's flow proof
is reported as certified, disproved, or unresolved. `--dry-run` previews all
16 batch commands. `RunSafeWeightedStatic.py --help` describes individual batches.
See [the method notes](weighted_flow_certification_notes.md) for the guarantees.

The main `flowtime`, `movements`, `scaled_best_bound`, and `scaled_absolute_gap`
columns belong to the reporting cutoff. Their `final_` counterparts belong to
the end of the search. `flow_proven` and `lexicographic_proven` assess the initial
candidate using all available evidence; `phase1_flow_proven` records whether
flow was already proved by the cutoff. `final_flow_proven` and
`final_lexicographic_proven` assess the final candidate. `weighted_runtime` is
the initial-stage runtime; `final_runtime` is the entire solver runtime, not
the extension alone. `extension_runtime` measures elapsed solver time beyond
the cutoff, including termination overhead. `weighted_cpu_time` adds model
construction and is marked as an estimate when the snapshot was taken during
the ongoing search. These are elapsed times, not summed CPU time over threads.

`first_flow_proof_runtime` is the first observed proof time measured in solver
seconds. `first_flow_proof_cpu_time` adds model construction to that observation.
`first_flow_proof_flowtime` identifies the certified flow value, and
`first_flow_proof_node_count`, `first_flow_proof_scaled_bound`, and
`first_flow_proof_work` record its checkpoint. The source identifies a callback
or the final-result check; the method identifies the mathematical certificate.
These fields are blank when no reliable proof was observed. Gurobi controls when
callbacks run, so this time is an upper bound on the actual instant
when proof became possible. `final_runtime` still records total solver time.
CSV protocol `safe_integer_flow_timing_v5` and check mode `bound_or_flow_change`
identify this version. Archived v4 used the larger, still sufficient coefficient
`R=U_g+1`. V5 retains its timing and conditional-extension rules while reducing
the coefficient using the movement lower bound. Always start v5 in a fresh
results directory; preserve v4 outputs and do not append v5 rows to them.
Old weighted bounds, gaps, and proof times must not be reused under the new
coefficient or reported as v5 results.

Monitoring computes distance bounds once and caches the scalar proof criterion
by flow value. Movement improvements reuse that criterion. Bound-change checks
use the saved candidate and strongest observed bound, without additional
solution-vector reads. Once flow is proved, routine proof checks stop; cutoff
handling, contradictory-solution detection, and final validation remain active.

### Targeted 70%-occupancy experiment

`Run70Percent.py` runs four-target, leave-mode SBM cases on all four
Table 2 layouts and output locations. It compares EF and LF on the same 100
initial-state seeds per layout:

| Layout | Escorts | Occupancy | Outputs |
| --- | ---: | ---: | --- |
| 13x7 | 27 | 70.33% | (6,0) |
| 10x10 | 30 | 70% | (0,0) |
| 16x10 | 48 | 70% | (4,0), (11,0) |
| 27x10 | 81 | 70% | (4,0), (13,0), (22,0) |

There are 400 distinct instances and 800 solver runs, producing four additional
configuration rows for Table 2(b). Escort counts are the nearest integer to
30% of the cells, so 13x7 cannot have exactly 70% occupancy. Occupancy counts all
stored loads, including the four targets. Both formulations receive the same
complete greedy start and share the sufficient integer coefficient and
physical horizon. The main budget is 300 solver seconds, followed only when
flow remains unproved by a same-tree extension of up to 300 seconds. Main
solution quality, bound and proof rates use the initial-cutoff columns;
extension results remain separate. The event-driven flow-proof tracking is
the same as in the safe weighted Table 2 runner.

From this repository on Linux or macOS, with NumPy, gurobipy and a working
Gurobi license in the selected Python environment:

```bash
python3 Run70Percent.py --dry-run
python3 -u Run70Percent.py --threads 16
```

The default thread count uses the Mac's performance cores or up to 16 threads
on Linux. Use `--threads 16` on the numerical-experiment Linux box to retain
the current campaign setting. `--threads` overrides
it, and `--python /path/to/python3` selects another solver environment.
`NUM_THREADS` and `PYTHON` are also respected. Without `PYTHON` or `--python`,
the solver uses the interpreter running the launcher. `python` and `python3`
can select different installations; the launcher prints the selected executable
and includes its actual version in preflight diagnostics. To explicitly use
the active Conda `python` environment, run:

```bash
python -u Run70Percent.py --python python
```

For a short trial before the full campaign:

```bash
python3 -u Run70Percent.py --seeds 1-3
```

Each launch creates a fresh `results_70percent_<timestamp>_<pid>/` directory.
It saves partial per-layout CSVs in `parts/`, solver logs in `logs/`, and merged
`occupancy70_escortflow.csv` and `occupancy70_loadflow.csv` at the top level.
`environment.json` records hardware, solver version, budgets, threads, protocol,
Git state and source hashes. A frozen `source/` copy supplies every batch.
`commands.sh` records the exact solver commands, and `pairing.json` confirms
identical instances and common settings after all batches complete. Existing
result directories are never overwritten; interrupted runs retain written rows.

Run the expanded campaign on the same Linux box, with the same thread setting
and native Gurobi version as the Table 2 baseline. Preserve the earlier Mac
campaigns separately. Cross-machine runtime comparisons do not isolate the
effect of occupancy.

### Continue-mode and additional target-count campaigns

All three launchers share `RunStaticCampaign.py`, use both EF and LF, and retain
the same source snapshots, pairing checks, integer objective, greedy warm starts,
300-second main budget, and conditional 300-second extension. `RunContinue.py`
uses escort counts **8, 12, 16**, plus the layout-specific approximately
70%-occupancy counts **27, 30, 48, 81** above. `RunTable2Targets.py` uses
**8, 12, 16** escorts for **2 targets**, and **8, 12, 16, 20** for **6 targets**.
It calculates the matched LP relaxation immediately after each integer instance
and writes the bound into that instance's CSV row. This is enabled by default.

| Launcher | Retrieval mode | Targets | Configuration rows | Distinct instances | Integer runs | LP solves by default |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `Run70Percent.py` | leave | 4 | 4 | 400 | 800 | 0 |
| `RunContinue.py` | continue | 2, 4, 6 | 48 | 4,800 | 9,600 | 0 |
| `RunTable2Targets.py` | leave | 2, 6 | 28 | 2,800 | 5,600 | 5,600 |

Counts assume all four layouts and seeds 1-100. Commands for the Linux box:

```bash
python3 -u Run70Percent.py --threads 16
python3 -u RunContinue.py --threads 16
python -u RunTable2Targets.py --threads 16 --lp-threads 16
```

Run the campaigns separately so their processes do not compete for memory or
solver threads. Add `--dry-run` to inspect commands, `--seeds 1-3` for a pilot,
or `--layouts 13x7` to select one layout. All launchers accept `--python`,
`--weighted-time-limit`, `--extension-time-limit`, and `--output-dir`. In tmux,
no `nohup` is necessary. Fresh result directories are respectively
`results_70percent_*`, `results_continue_*`, and `results_table2b_targets_*`.
Merged files are `occupancy70_{escortflow,loadflow}.csv`,
`continue_{escortflow,loadflow}.csv`, and
`table2b_targets_{escortflow,loadflow}.csv`. Per-configuration files include
the target and escort counts whenever these vary.

After the current four-target run finishes, update the Linux checkout with
`git pull --ff-only`, then run `RunTable2Targets.py` in tmux using the solver's
Python environment. No new instance or model preparation is needed. To inspect
the complete plan without solving, use `python3 RunTable2Targets.py --dry-run`.
The result directory freezes the Python sources before any batch begins.

For Tables 2 and 3, the merged integer CSVs contain the common greedy FT/MV,
distance lower bound (`safe_movement_lower_bound`, also the naive FT bound),
the 300-second incumbent and certified bound, first FT proof time, solver and
CPU times, and final FT/MV/bound after any extension. Best-known FT/MV and
improvement over the greedy heuristic can therefore use both models' final
solutions, while the method comparison retains its initial 300-second cutoff.

For each instance, the sequence is greedy preparation, integer solve (including
any conditional extension), and then a separate continuous LP solve with the
same coordinates, `R`, and model horizon. The integer row is flushed before
starting the LP. Its `lp_relaxation_lower_bound` column is the optimal value of
`FT + MV/R`; `lp_status` and `lp_elapsed_seconds` record the LP outcome and time.
Additional columns retain the LP components, threads, budgets, solver version,
and input/source fingerprints. An unfinished or failed LP leaves the bound blank
and preserves the integer result. Integer solver/CPU times and `total_wall_time`
exclude LP computation. No LP bound is fed into the integer benchmark search.

At the end, the replay tool collects these saved optimal values into
`lp_results.csv` and only solves missing or failed cases. Completed relaxations
are reused after checking the instance, coefficient, horizon, and model-source
fingerprints. Successful retries fill missing LP columns in the new merged and
per-configuration integer CSVs. Target-load count remains part of each output key.
The model-specific LP gaps can then be calculated against the common best-known
integer solution. It never substitutes the obsolete fixed `0.01` coefficient.
See [LP_REPRODUCIBILITY.md](LP_REPRODUCIBILITY.md) for the gap formula.

Integrated LPs run sequentially and default to one solver thread to limit memory use.
`--lp-threads` changes this independently of integer `--threads`. `--lp-workers`
controls only parallel recovery of missing LP bounds at the end. The LP limit
is 300 seconds, followed by a 600-second barrier retry
after a time limit; change these using `--lp-time-limit` and
`--lp-retry-time-limit`. Add `--no-lp` for an integer-only pilot. The other two
launchers can also calculate LPs if explicitly given `--lp`.

LP failures leave valid integer rows in place and do not stop later integer
batches. The final recovery phase retries those LPs and reports an error if
coverage remains incomplete. To retry LPs after the integer campaign finishes:

```bash
RESULTS=/absolute/path/to/results_table2b_targets_...
bash "$RESULTS/lp_commands.sh"
```

This script uses the frozen sources and `--resume`, which verifies input and
source fingerprints before reusing results. The integer launcher requires a
fresh output directory and does not resume interrupted integer batches. Integrated
LP progress is in the per-configuration log alongside integer progress. Recovery
and collection are logged in `logs/table2b_targets_lp.log`, and
`lp_results.csv.manifest.json` records the completed and total counts.

For an individual integer-plus-LP batch, add `--with-lp` to the safe runner:

```bash
python -u RunSafeWeightedStatic.py --formulation loadflow \
  -x 13 -y 7 -O 6 0 -l 2 -e 8,12,16 -r 1-100 -m leave \
  --threads 16 --with-lp --lp-threads 16 -f two_targets_loadflow.csv
```

`--with-lp` means integer plus LP in one CSV. `--lp` still means an LP-only run.
Both integer-plus-LP workflows use the current v5 rule. The currently running
v4 campaign keeps its existing frozen code and can use the separate v4 replay.

In continue mode, a target is served when it reaches an output and immediately
becomes an ordinary blocking load. It does not create an escort. An initially
output-located target is already served at time zero. The greedy trace and both
Gurobi formulations use these same conventions. The LF target-to-blocker
conversion preserves occupancy and permits the retrieved blocker to move again;
there is no leave-mode output-service delay. The corrected greedy flow time sums
actual target arrival times without adding an extra iteration after retrieval.

EF and LF match initial states within every configuration. Equal seeds also
match leave and continue modes at a given target and escort count. Across target
counts, the existing generator gives nested target sets but different escort
locations because escorts follow the target prefix of the permutation. The
campaigns hold escort **counts** constant; they do not claim identical escort
locations across 2, 4, and 6 targets.

### Fixed-weight Table 2 experiment

`RunTable2Weighted.sh` runs both formulations on the current Table 2(a,b)
instances with **300 seconds for the weighted solve** and an independent
**300-second certification budget**. Both weighted formulations receive a complete
greedy warm start; each certification solve starts from its weighted incumbent.
The weighted objective is `100*F + M`, exactly equivalent to `F + 0.01*M`.
Both objectives are integer, so the solver uses `MIPGap=0`, `MIPGapAbs=0.999`,
followed by independent numerical certificate checks. Original weighted units
are used in the reported weighted objective, bound, and absolute gap.

From the repository directory on Linux:

```bash
bash RunTable2Weighted.sh --dry-run
nohup bash RunTable2Weighted.sh --python python3 --threads 12 > table2_weighted.log 2>&1 &
```

The script retains the Table 2 layouts, escort counts, seeds, leave retrieval,
BM movement, and sequential execution of `RunTable2Lex.sh`. It writes a new
timestamped result directory with `parts/`, `logs/`, merged CSVs, environment
metadata, and the executed commands. Each batch log reports per-instance
progress. Existing output directories are rejected to prevent mixing runs.

Certification runs for proven weighted optima, or for unproven weighted
candidates whose **absolute gap is below 0.1 in `F + 0.01*M` units**. The gate is
strict and configurable with `--certification-gap-threshold`. A larger gap or
missing incumbent skips certification. The gate is not itself an optimality
certificate. `--weighted-time-limit` and `--certification-time-limit` set the
two independent budgets. The runtime columns measure Gurobi optimization;
the corresponding elapsed-time columns also include model construction, which
does not consume either Gurobi time limit.

The check minimizes flow time in a separate model and stops when its lower bound
proves the saved weighted candidate's flow time, or when it finds a better flow
time. The weighted result is preserved in either case. Its full solution is saved
before model disposal and transferred into the certification model, with explicit
values for every variable. The check selects a sufficient horizon covering every potentially better-flow
schedule, shortening or extending the weighted horizon as needed.
The sufficient horizon is calculated from the weighted incumbent's flow time,
not the heuristic's. The certificate uses this sufficient horizon directly. The start is truncated
or padded with idle periods accordingly. Redundant post-retrieval movements
removed from the start are counted in the CSV; all retrieval times remain unchanged.
Weighted optimality remains scoped to its recorded horizon. A separate global
lexicographic flag also requires that the weighted horizon cover all candidates
with the saved flow time. See [the method notes](weighted_flow_certification_notes.md)
for the small-weight proof, horizon argument, and numerical safeguards.

`RunWeightedStatic.py --help` documents the standalone batch interface. Existing
weighted and two-phase CLI behavior remains available. `LoadFlowStatic.py` also
now accepts `--warmstart` for integer Gurobi BM leave runs; it uses the same
physical greedy plan as escort-flow, including blocking-load movements.

### Earlier lexicographic Linux batch run for Table 2(a) and (b)

`RunTable2Lex.sh` runs both lexicographic formulations with a **270-second phase-one
cap and a 300-second total solver limit per instance**. Phase 2 receives 300 seconds
minus the actual phase-one solver runtime. Keep this script with the current
`Code/` files on the Linux machine and use Python with NumPy, `gurobipy`, and an
active Gurobi license. Python 3.11 is recommended.

From `Code/`:

```bash
# Inspect the complete command plan without running Gurobi.
bash RunTable2Lex.sh --dry-run

# Run in the background with the active environment's Python.
nohup bash RunTable2Lex.sh --python python3 --threads 12 > table2_lex.log 2>&1 &
```

The script follows the current manuscript, including its exact four layouts:

| Grid | Output cells |
|---|---|
| 13 x 7 | (6, 0) |
| 10 x 10 | (0, 0) |
| 16 x 10 | (4, 0), (11, 0) |
| 27 x 10 | (4, 0), (13, 0), (22, 0) |

Table 2(a) uses one target and 3--8 escorts; Table 2(b) uses four targets and
8, 12, or 16 escorts. Each configuration uses seeds 1--100, `leave` mode, and
simultaneous block movement. This is 3,600 distinct instances, each solved by both
formulations, or **7,200 lexicographic instance solves**. Each solve may have two
optimization phases. Jobs run sequentially with 12 threads by default, balanced
MIP emphasis, and the existing automatic horizon selection. This batch generates
the new integer results for Table 2's configurations; it does not rerun the
historical weighted LP relaxation columns.

Each launch creates a fresh timestamped results directory beside the script,
containing:

- `table2a_escortflow_lex.csv` and `table2a_loadflow_lex.csv`: 2,400 rows each.
- `table2b_escortflow_lex.csv` and `table2b_loadflow_lex.csv`: 1,200 rows each.
- `parts/`: the original per-layout CSV batches, including any incomplete batch.
- `logs/`: per-layout solver logs.
- `environment.txt` and `commands.sh`: runtime versions, settings, and exact commands.

Completed batches are checked for missing/duplicate instance rows, solver errors,
and the requested limits before being merged into CSVs with one header each.
Ordinary time-limit results remain valid experiment observations. A failed batch
stops the script while preserving prior results and diagnostic files. An existing
output directory is rejected, preventing accidental duplicate appends.

Use `--part a` or `--part b` to run a single table part, `--seeds 1-10` for a smaller
subset, `--output-dir DIR` to choose a new destination, or `--python PATH` to select
a different interpreter. `PYTHON` and `NUM_THREADS` can also set their respective
defaults. Run `bash RunTable2Lex.sh --help` for all options.

### Dynamic rolling-horizon paper: simulation experiments

These commands reproduce the raw CSV files used for the dynamic rolling-horizon paper. The shell scripts contain the exact command lines for each experiment family.

```bash
cd /Users/talraviv/Library/CloudStorage/Dropbox/research/PBS/EscrotsFlow/Code
DYNAMIC_DIR=../replication/dynamic_paper
mkdir -p "$DYNAMIC_DIR"
rm -f FullFactorial9x5-o0.2.csv FullFactorial13x7.csv \
  Modular_vs_integrated.csv AttentionTest.csv HybridRatioTest.csv

bash FullFactor9x5.sh
bash FullFactor13x7.sh
bash TestIntegrated.sh
bash TestAtten.sh
bash TestHybridRatio.sh

mv FullFactorial9x5-o0.2.csv "$DYNAMIC_DIR/"
mv FullFactorial13x7.csv "$DYNAMIC_DIR/"
mv Modular_vs_integrated.csv "$DYNAMIC_DIR/"
mv AttentionTest.csv "$DYNAMIC_DIR/"
mv HybridRatioTest.csv "$DYNAMIC_DIR/"
```

The generated CSV files correspond to the paper experiments as follows:

- `FullFactorial9x5-o0.2.csv`: full-factorial RTRH evaluation on the `9 x 5` grid with escorts `4` and `6`, arrival rates `0.2` and `0.4`, and the full/surrogate/hybrid/PLPR/attention factors.
- `FullFactorial13x7.csv`: full-factorial RTRH evaluation on the `13 x 7` grid with escorts `8` and `12`, arrival rates `0.2` and `0.4`, and the same algorithmic factors.
- `Modular_vs_integrated.csv`: modular-versus-integrated layout comparison.
- `AttentionTest.csv`: sensitivity to the attention limit.
- `HybridRatioTest.csv`: sensitivity to the hybrid fallback ratio.

`TestWarmStart.sh` is a separate warm-start sandbox and is not part of the main dynamic-paper replication set unless the warm-start appendix or supplementary analysis is being regenerated.

## One-step heuristic for static retrieval

The static paper describes this method in the online supplement, Section "Greedy heuristic", Algorithm `OneStep` (source: `static_escort_flow_IISE_supplement.tex`). All retrieval requests are known initially. The heuristic constructs a feasible retrieval plan, supplies an objective upper bound and information for choosing the ILP horizon, and serves as a computational benchmark.

`OneStep()` plans **one time step, or takt, of simultaneous block movements**. `SolveGreedy()` repeatedly calls it until all requested loads reach output cells. A single call can select several escort movements. An escort is an empty cell: moving it along a row or column shifts every load on that segment by one cell in the opposite direction. For example, moving an escort from `(1, 2)` to a target at `(3, 2)` shifts the intervening load to `(1, 2)` and the target to `(2, 2)`, all in one takt. The target advances one cell even though the escort travels two cells.

### The rule described in the paper

Give each target a unique, fixed priority ID and retain that ID when the target moves. The supplement reports assigning priorities by increasing **initial** Manhattan distance to the closest output. For each target's current position, select its closest output, breaking equal-distance ties lexicographically by output coordinates. Partition the other grid cells relative to that target and output:

| Zone | Definition | Purpose of an escort there |
| --- | --- | --- |
| A | Same row or column as the target, in a direction toward its designated output. These are rays extending to the grid boundary. | Can advance the target directly if its path is available. |
| B | Outside A, sharing a row or column with an A cell, without the target lying between them. | Can be repositioned into A. |
| C | Outside A and B, sharing a row or column with a B cell. | Can be repositioned into B. |
| D | All remaining cells; present only when target and output share a row or column. | Can be repositioned into C. |

At the beginning of each takt, remove completed requests and clear the set of reserved cells. Scan the remaining targets in increasing ID order. Skip a target if an earlier block movement has already moved it in this takt. For each other target, attempt these four operations in order:

1. **Advance the target.** Move an available Zone A escort to the target's cell. This shifts the target one cell toward its designated output. Where possible, extend that same escort movement beyond the target to advance additional lower-priority targets too. Recompute the zones around the target's new position if it moved.
2. **Prepare Zone A.** Move at most one available escort from B to A.
3. **Prepare Zone B.** Move at most one available escort from C to B.
4. **Prepare Zone C.** Move at most one available escort from D to C.

All four operations are attempted, including the preparation operations after a successful target advance. The paper calls for the shortest eligible escort movement at each operation. Preparation moves may also be extended to help lower-priority targets, while keeping the escort endpoint in the intended destination zone.

Every selected escort path must be horizontal or vertical, contain no other escort, and share no cell, including endpoints, with a previously selected path. Reserve the current target's cell even if it did not move, so subsequent decisions cannot push it backward. Lower-priority targets may be shifted incidentally. An escort already assigned a movement is unavailable for the rest of this takt. Thus, the operations construct one simultaneous batch; they cannot move the same escort through D, C, B, and A within a single takt.

With fixed priorities, the supplement proves that the highest-priority remaining target never moves farther from its output and advances at least once in every four takts. In the worst case an escort first progresses through `D -> C -> B -> A` over three takts, then advances the target in the fourth. Consequently, the paper's algorithm retrieves all targets within `4 * n * d_max` takts, where `n` is the initial number of targets and `d_max` is the largest distance from any grid cell to its closest output. The termination and acyclicity arguments rely on retaining the priority order across takts.

### Current implementation and paper differences

The implementation is in `OneStepHeuristic_v2.py`. Its public inputs are the grid dimensions, output locations, target locations, and escort locations. For `OneStep`, targets are a dictionary from locations to IDs; `SolveGreedy` accepts target locations and assigns the IDs itself.

- **Priority:** `acyclic=True` scans targets by fixed ID. `SolveGreedy` assigns these IDs once by initial distance, with target-coordinate ties. The default, `acyclic=False`, reorders targets at every takt by their current distance to the closest output. Both static runners currently use this default and expose no `--acyclic` option. Their default runs therefore do not implement the paper's fixed-priority rule, and the paper's proof should not be attributed to that default.
- **Candidate selection:** the code orders escorts by distance to the current target, with coordinate ties, for every operation. For the direct advance this is also the escort movement length; for preparation moves it need not select the shortest movement specified in the paper.
- **Output ties:** `build_dist_map` retains the first equally near output encountered in its input iteration order. The static runners pass a set, so movement directions do not explicitly implement the paper's lexicographic output tie rule.
- **Protection and extensions:** the code prefers moves that do not increase any lower-priority target's distance to its closest output. Only the highest-priority target may fall back to a harmful base move to ensure progress. If an extension fails this guard, the unextended move is tried. Extensions of preparation moves must end in the intended destination zone.

These distinctions matter for reproducing exact trajectories. Selecting `acyclic=True` enables the fixed-priority variant but does not change the other implementation details above.

### Running and interpreting the heuristic

Run the current static greedy benchmark from `Code/`:

```bash
python3 EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 12 -l 4 -m leave -r 11 \
  --greedy -f greedy_static.csv
```

This generates one instance with seed `11`, runs the heuristic, and appends its results to `greedy_static.csv`. It requires NumPy but does not invoke Gurobi or CPLEX. Use `-r 1-100` for 100 instances and add `-a` to export animation traces. `--greedy` supports `leave` and `continue`, not `stay`.

To select fixed priorities directly through the Python API:

```python
from OneStepHeuristic_v2 import SolveGreedy

makespan, total_flowtime, movements = SolveGreedy(
    5, 5,
    {(0, 0)},                 # output cells
    {(2, 2), (4, 4)},         # target-load cells
    {(0, 1), (3, 1)},         # escort cells
    acyclic=True,
    retrieval_mode="leave",
    max_steps=1000,
)
print(makespan, total_flowtime, movements)
```

In static `leave` mode, `SolveGreedy` removes a request on arrival, blocks that output cell for the following full takt, and then makes it available as an escort. Arrival at the end of takt `t` therefore allows that escort to move from takt `t + 2`. In `continue` mode, an arrived load stops being a target and remains a blocking load; no new escort is created. The `leave` timing is enforced by the `SolveGreedy` wrapper; calling `OneStep(..., retrieval_mode="leave")` directly instead converts a target already at an output into an escort at the start of that call.

In the static `leave` results, flow time is the sum of target arrival times, makespan is the last arrival time, and movements count all individual one-cell load shifts, including blocking loads. The `Greedy UB` column in `EscortFlowStatic.py` is `beta * total_flowtime + gamma * movements`. It is an objective upper bound, not a mean flow time. The weights affect the reported objective, not the heuristic's move choices. `OneStep`'s returned `moves` list also includes blocking-load shifts; use `return_escort_moves=True` for escort paths or `return_target_moves=True` for the target-ID movement map.

Without `--greedy`, the normal static BM workflow still runs this heuristic in `leave` and `continue` modes. Unless a DP horizon is supplied, the runners choose their horizon from its makespan, with backend-specific time indexing. In the standard `EscortFlowStatic.py --gurobi` backend, `--warmstart` supplies its trace before optimization, and `--cutoff` enables its objective cutoff. `LoadFlowStatic.py --warmstart` now supplies the same physical plan for Gurobi integer BM leave runs. The standalone greedy benchmark performs no MILP solve.

The main paper's horizon theorem gives the sufficient horizon bound `sum(f_i) - sum(d_i) + max(d_i)`, where `f_i` are feasible arrival times and `d_i` are initial distances to the closest outputs, under lexicographic minimization of flow time and movements. The current runners use the greedy makespan instead of that expression. A feasible greedy makespan supplies a horizon containing a feasible plan; by itself it does not establish that the horizon contains an unrestricted flow-time optimum.

## Static optimization scripts

The repository contains two static experiment runners:

- `EscortFlowStatic.py`: escort-flow formulation
- `LoadFlowStatic.py`: load-flow formulation

These scripts are the main entry points for the static benchmark experiments in this repository.

### `EscortFlowStatic.py`

Purpose:

- runs the static escort-flow ILP
- supports single-load and multi-load experiments
- can use a DP-based upper bound in the single-load case
- can export animation traces

Key dependencies:

- `PBSCom.py`
- `PBS_DPHeuristic_lm.py`
- `PBS_DPHeuristic_bm.py`
- `escort_flow_static_gurobi.py`, `escort_flow_static_lazy.py`, and `escort_flow_static_bnc.py`
- escort-flow OPL model files and `oplrun` only for the legacy `--opl` path

Common arguments:

- `-x`, `-y`: PBS dimensions, required
- `-O`: output cells as coordinate pairs, required
- `-e`: escort-count range, default `5`
- `-r`: replication/seed range, default `1`
- `-l`: number of target loads, default `1`
- `-m`: retrieval mode, one of `stay`, `leave`, `continue`, default `leave`
- `-f`: CSV result file, default `res_escort_flow.csv`
- `--beta`: flowtime weight, default `1.0`
- `--gamma`: movement weight, default `0.01`
- `-T`: legacy horizon scaling factor, default `1.6`; retained in the CLI but retired from the normal static workflow
- `-t`: solver time limit, default `300`
- `--num_threads`: solver thread count, default `8` on macOS and `12` on Linux
- `--work_limit`: Gurobi work limit in work units, default none
- `--mip_emphasis`: Gurobi MIP emphasis, one of `balanced`, `feasibility`, `optimality`, or `bound`, default `balanced`
- `--dp_file`: DP table file for the single-load case, default empty
- `-k`: `k'` parameter used with the DP heuristic, default `0`
- `--lp`: LP relaxation option, default off
- `--greedy`: solve the static instance with the greedy heuristic instead of a MILP backend, default off; supported for `continue` and `leave` modes. In static `leave` mode, an arrived target occupies its output for one further takt before the cell becomes an available escort. See [One-step heuristic for static retrieval](#one-step-heuristic-for-static-retrieval) for the algorithm, priority settings, and examples
- `--gurobi`: explicitly select the Gurobi Python backend; accepted for clarity but now redundant because Gurobi is the default
- `--opl`: switch to the legacy `oplrun` / CPLEX path instead of the default Gurobi backend
- `--warmstart`: enable heuristic MIP start with the static Gurobi backend, default off
- `--naive`: skip optimization entirely, compute only the naive lower bound for each generated instance, and write a reduced CSV row
- `--lazy` or `--lazy N`: use the lazy-constraint Gurobi backend, default off; a bare `--lazy` means `0`, and `N` is the number of initial time steps kept in the master problem
- `--bnc` or `--bnc N`: use the branch-and-cut Gurobi backend, default off; a bare `--bnc` means a per-separated-node user-cut cap of `2*T`, and `N` sets that cap explicitly
- `-a`: export animation trace, default off

Range syntax:

- ranges accept `n`, `n1,n2,...`, `start-end`, `start:end`, `start-end-step`, and `start:step:end`
- example: `-e 8-16-4` means `8, 12, 16`

Notes:

- the default solver path is the Gurobi Python backend; `--opl` switches to the legacy `oplrun` / CPLEX path
- for `stay` mode, the number of output cells must be at least the number of target loads
- the DP-based upper bound currently applies only to the single-load case
- without `--dp_file`, the static horizon upper bound is taken from the greedy heuristic
- on the standard `--gurobi` path, the full static escort-flow model is built explicitly in the master problem
- with the standard Gurobi backend, `--warmstart` supplies a complete greedy start before optimization; legacy weighted lazy/BnC backends retain their fallback-restart behavior
- `--lazy` selects a separate Gurobi backend that keeps the flow/supply structure in the master and enforces the target-movement coupling constraints lazily; if you pass `--lazy N`, the first `N` time steps of that coupling family stay in the master
- `--bnc` selects a separate branch-and-cut backend; the cheap strong constraints and constraint family `(9)` stay in the master, and constraint family `(8)` also stays explicit for the first `T // 8` time steps, while the remaining later `(8)` constraints are separated by enumeration in callbacks
- in the current BnC implementation, user cuts are generated at the root and at a small number of early branch-and-bound nodes near the root; the root node is uncapped, while later separated nodes use the configurable `--bnc N` cap, which defaults to `2*T` when omitted, and under that cap the strongest violations are added first; incumbent violations are still rejected with lazy constraints using a cap of `4*T`
- `--lazy`, `--bnc`, and `--work_limit` all imply or apply only to the Gurobi backend
- `--lazy` and `--bnc` are integer-only MILP options; they cannot be combined with `--lp`, `--greedy`, or with each other
- `--naive` is a reporting-only mode: it cannot be combined with solver-selection flags, greedy mode, DP files, warmstarts, animation export, or `-k`

Static CSV output:

- each row now begins with `Machine Name`, `Time Stamp`, and `version`, mirroring the dynamic simulator
- regular static rows also include `Greedy UB` and `Naive LB` before the solver-result columns
- regular static rows also include `MIP Emphasis`; it is the selected Gurobi emphasis for MILP solves and `-` when not applicable
- `Wall Clock Time` is the elapsed solve time measured by Python around the backend call
- `Work` is the Gurobi work value when a Gurobi backend is used; it is blank for the OPL/CPLEX path
- `User Cut Time` is nonzero only for the BnC backend and measures time spent inside the callback separation logic
- for the BnC backend, `Model` is `ILP-Gurobi-BnC`, and the resolved cut cap is written separately in `Max User Cut Per Node`
- with `--greedy`, the CSV row stops after `Naive LB`; the final 9 solver-only columns are omitted
- with `--naive`, the script writes a reduced CSV containing only the instance description and the naive lower bound

Single-load examples:

```bash
python3 EscortFlowStatic.py -x 10 -y 10 -O 0 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set10x10_k3_corner.p
```

```bash
python3 EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set16x10_k3_2centers.p
```

Multi-load example:

```bash
python3 EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 8-16-4 -l 4
```

For multi-load experiments, drop `-k` and `--dp_file`. The horizon upper bound is then taken from the greedy heuristic.

Gurobi backend examples:

```bash
python3 EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 12 -l 4 -m leave -r 11 --gurobi
```

```bash
python3 EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 12 -l 4 -m leave -r 11 --lazy 2
```

```bash
python3 EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 12 -l 4 -m leave -r 11 --bnc 10 --work_limit 500
```

With bare `--bnc`, the per-node user-cut cap defaults to `2*T`, where `T` is the horizon selected for that instance.

Naive lower-bound example:

```bash
python3 EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -r 1-100 -l 1 --naive
```

### `LoadFlowStatic.py`

Purpose:

- runs the static load-flow formulation
- supports the same benchmark family as `EscortFlowStatic.py`
- can solve the ILP or LP relaxation
- defaults to BM and switches to LM only with `--lm` or `--LM`
- uses the Gurobi Python API by default; `--opl` switches to the legacy `oplrun` / CPLEX path
- can also export animation traces

Key dependencies:

- `PBSCom.py`
- `PBS_DPHeuristic_lm.py`
- `PBS_DPHeuristic_bm.py`
- `load_flow_static_gurobi.py`
- `pbs_load_flow_multi.mod`, `pbs_load_flow_multi_lp.mod`, and `oplrun` only for the legacy `--opl` path

Common arguments:

- `-x`, `-y`: PBS dimensions, required
- `-O`: output cells as coordinate pairs, required
- `-e`: escort-count range, default `5`
- `-r`: replication/seed range, default `1`
- `-l`: number of target loads, default `1`
- `-m`: retrieval mode, default `leave`
- `-f`: CSV result file, default `res_load_flow.csv`
- `--alpha`: makespan weight, default `0.0`
- `--warmstart`: supply a complete greedy MIP start in Gurobi integer BM leave mode, default off; unsupported mode combinations are rejected
- `--beta`: flowtime weight, default `1.0`
- `--gamma`: movement weight, default `0.01`
- `-T`: legacy horizon scaling factor, default `2.0`; retained in the CLI but retired from the documented BM workflow
- `-t`: solver time limit, default `300`
- `--num_threads`: solver thread count, default `8` on macOS and `12` on Linux
- `--mip_emphasis`: Gurobi MIP emphasis, one of `balanced`, `feasibility`, `optimality`, or `bound`, default `balanced`
- `--lm` or `--LM`: run LM instead of the default BM mode, default off
- `--dp_file`: DP table file for single-load upper bounds, default empty
- `-k`: `k'` parameter for the DP heuristic, default `0`
- `--lp`: solve the LP relaxation instead of the ILP, default off
- `--gurobi`: explicitly select the Gurobi Python backend; accepted for clarity but now redundant because Gurobi is the default
- `-a`: export animation trace, default off

Notes:

- the default solver path is the Gurobi Python backend; `--opl` switches to the legacy `oplrun` / CPLEX path
- the script header says only `leave` is supported at present; that is the safe mode to use
- as in `EscortFlowStatic.py`, the DP table route is for the single-load case
- without `--dp_file`, BM runs in `leave` and `continue` mode use `OneStepHeuristic_v2` to get an upper bound on `T`

Single-load examples:

```bash
python3 LoadFlowStatic.py -x 10 -y 10 -O 0 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set10x10_k3_corner.p
```

```bash
python3 LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set16x10_k3_2centers.p
```

Multi-load example:

```bash
python3 LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 8-16-4 -l 4
```

### Direct LP runs from parameters and seeds

Add `--lp` to either regular static runner, or to `RunSafeWeightedStatic.py`.
The same dimensions, outputs, targets, escorts, seed, and retrieval mode
generate the same initial instance as the integer run. Integer CSVs are not
required. The runner obtains the greedy reference and calculates each
instance's `R` and horizon before solving the continuous `F+M/R` relaxation.

Use `--lp-protocol v4` for the archived single-/four-target experiment. For
example, these commands cover all 300 instances of one four-target grid,
including instances not yet present in the copied integer results:

```bash
python -u EscortFlowStatic.py -x 13 -y 7 -O 6 0 -l 4 -e 8,12,16 -r 1-100 \
  -m leave --lp --lp-protocol v4 --num_threads 16 -t 300 \
  -f lp_four_escortflow_13x7.csv

python -u LoadFlowStatic.py -x 13 -y 7 -O 6 0 -l 4 -e 8,12,16 -r 1-100 \
  -m leave --lp --lp-protocol v4 --num_threads 16 -t 300 \
  -f lp_four_loadflow_13x7.csv
```

Run `python` from the Conda environment with NumPy, `gurobipy`, and the Gurobi
license. To cover the other grids, use `10x10` with `-O 0 0`, `16x10` with
`-O 4 0 11 0`, and `27x10` with `-O 4 0 13 0 22 0`. Use a different `-f`
filename for each model and grid. To resume, repeat the same command with
`--resume`. Completed optimal LPs are reused after checking the saved inputs
and frozen sources. Instance parameters cannot be changed during a resume.

The default `--lp-protocol v5` matches the current safe integer runner and the
new two-/six-target campaigns. `-m continue` selects continue mode. The safe
runner uses `--formulation escortflow` or `--formulation loadflow`, `--threads`
for solver threads, and `--weighted-time-limit` for the LP budget. In all three
runners, `--lp-workers` defaults to one and `--lp-retry-time-limit` defaults to
600 seconds. LP outputs use the `RunStaticLP.py` schema and only contain
`OPTIMAL` values. Coordinate/coefficient snapshots and frozen sources are saved
beside the output CSV.

The earlier v4 coefficient remains sufficient, but changing `R` or the horizon
changes the LP being solved and can affect integer search performance. Keep
the coefficient/horizon protocol matched to the integer campaign. For exact
published-result replay, use the archived input CSVs and frozen model sources
with `RunStaticLP.py`. The regular runners also accept explicit `--flow-weight R`
and `--horizon T`; EF's physical horizon is `T+1` and LF's is `T`. Historical
fixed-weight commands require `--lp --legacy-lp`. See
[LP_REPRODUCIBILITY.md](LP_REPRODUCIBILITY.md) for formulas, files, and restart checks.

### DP helper modules

`PBS_DPHeuristic_lm.py` and `PBS_DPHeuristic_bm.py` are helper modules used by the static scripts when `--dp_file` is supplied. They:

- load a precomputed DP table from a pickle file
- generate an upper-bound policy for the single-load case
- provide a feasible horizon estimate before the full OPL model is solved

These helpers depend on:

- `PBSCom.py`
- a compatible DP table pickle

They are not the main experiment entry points; normally you use them through `EscortFlowStatic.py` or `LoadFlowStatic.py`.

## Typical workflow

### 1. Run a simulation

#### Current version: `EscortFlowSim_v8.py`

`v8` is the maintained dynamic simulator. It uses an in-process persistent Gurobi model, writes the CSV/raw outputs used by the rest of this repository, and is the sandbox for the current rolling-horizon experiments.

Example greedy-only run:

```bash
python3 EscortFlowSim_v8.py -x 9 -y 5 -O 4 0 -e 4 -R 0.4 \
  -S 1600 --greedy -f results.csv -H
```

Example surrogate rolling-horizon run:

```bash
python3 EscortFlowSim_v8.py -x 9 -y 5 -O 4 0 -e 4 -S 2000 \
  -T 5 -I 1 -E 1 -R 0.2 -M 6 -q spt -L -f results.csv -H
```

Example full-model BnC hybrid run:

```bash
python3 EscortFlowSim_v8.py -x 9 -y 5 -O 4 0 -e 8 -R 0.4 -S 100 \
  --full --bnc --hybrid -I -f sim_escort_flow.csv -H
```

Key notes for `v8`:

- `v8` builds the rolling-horizon Gurobi model once per active horizon shape and reuses it across epochs by updating only state-dependent data.
- `-t` is the per-solve Gurobi time limit, and `--num_threads` controls Gurobi threads.
- `v8` uses balanced Gurobi search by default (`MIPFocus=0`) and switches to optimality emphasis (`MIPFocus=2`) when `--warmstart` is enabled.
- `--warmstart` accepts `greedy`, `ilp`, or `ilp greedy`; by default unused arc variables are set to `0`, and `nozero` disables that.
- the warmstart label is appended to the CSV `Algorithm Name` whenever a warmstart is active.
- if a solve is interrupted with `CTRL-C`, the script exits instead of continuing to later instances in the loop.
- the interactive progress line shows both arrivals and departures, and after the run finishes it prints the raw average lead time after cooldown trimming and warmup deletion.

Model variants:

- surrogate model: controlled directly by `-T`, `-I`, and `-E`
- full model: enabled by `--full`; the planning horizon is recomputed from the greedy completion horizon at each solve
- with `--full` and no `-I`: `v8` uses a partial LP relaxation (`PLPR`), where only the first `-E` periods are integer and the remaining greedy-based horizon is fractional
- with `--full` and any use of `-I` (including bare `-I`): the entire greedy-based full horizon is integer
- in the CSV `Algorithm Name`, `PLPR` is added when the active model contains a fractional tail

Branch-and-cut in `v8`:

- `--bnc` is available for both the surrogate and full MILP paths
- under `--bnc`, the target-movement coupling constraints are removed from the master and separated in callbacks
- user cuts are generated only at the root node
- in PLPR mode, user-cut separation scans the integer periods only
- lazy constraints still enforce the full horizon, including fractional periods when they exist
- the per-callback cut budget is `2*T`, where `T` is the active fractional horizon for that solve

Warmstart behavior in `v8`:

- warmstarts are allowed only on all-integer surrogate models, so for surrogate `--warmstart ilp` and `--warmstart ilp greedy`, `--integer_horizon` must equal `--fractional_horizon`
- `--warmstart` is recorded in both the dedicated warmstart CSV columns and in the `Algorithm Name`

Important options:

- `-x`, `-y`: PBS dimensions, required
- `-O`: output cells in pairs of coordinates; defaults to `0 0` if omitted, though in practice you normally set them explicitly
- `-e`: number of escorts, default `8`
- `-S`: number of requests in the simulation, default `1000`
- `-R`: Poisson request arrival rate per time step, default `0.1`
- `-E`, `--epoch`: execution epoch length, default `1`
- `-T`: surrogate fractional horizon; ignored by `--full`
- `-I`: integer horizon; bare `-I` means all periods are integer
- `-t`: per-solve Gurobi time limit in seconds; default `--epoch`
- `-M`, `--max_attention`: attention limit, i.e. the maximum number of target loads considered concurrently; if omitted, all cells are eligible
- `-q`: queue management, either `fifo` or `spt`, default `spt`
- `-m`, `--offline`: offline rolling horizon for pure ILP mode only, default off; in offline mode the MILP sees requests visible at the current decision time, and the flag cannot be combined with `--greedy` or `--hybrid`
- `--full`: use the full MILP instead of the surrogate model, default off
- `--bnc`: use branch-and-cut separation on the movement-coupling family, default off
- `--greedy`: greedy heuristic only, default off; currently requires `--epoch 1`
- `--hybrid`: switch to greedy epochs when the number of new open requests is large enough, default off
- `--hybrid_ratio`: ratio used by `--hybrid`, default `1.0`
- `--acyclic`: use seniority-based greedy priority instead of distance-based order, default off
- `--gamma`: movement weight, default `0.01`
- `--distance_penalty`: MILP distance penalty, default `1`
- `--time_penalty`: MILP time penalty, default `1`
- `--num_threads`: Gurobi thread count; default `8` on macOS and `12` on Linux
- `-L`: write a detailed log file, default off
- `-a`, `--save_raw`: save a raw pickle trace for post-processing and animation, default off
- `-f`: output CSV file, default `sim_escort_flow.csv`
- `-H`: write CSV header if needed, default off
- `-o`, `--max_opt_gap`: reject MILP solutions above this gap and fall back to greedy, default `0.4`
- `--seed`: random seed, default `0`
- `--warmstart`: warmstart mode list; default is no warmstart

`v8` timing and control notes:

- The simulator plans once per epoch and then executes the resulting move list for that epoch.
- MILP visibility is delayed by one epoch: requests must be visible by the beginning of the previous epoch to enter the next MILP solve.
- With `--offline`, the MILP instead uses the requests visible at the current decision time.
- In hybrid mode, `old requests` are the currently open requests that were already visible at the beginning of the previous epoch, and `new requests` are the currently open requests that were not visible then.
- If an ILP solution is unusable, `v8` falls back to greedy on all target loads visible at the current decision time and fills the whole epoch greedily.
- fallback greedy runs are now split in the CSV into `Fallback No Feasible Runs`, `Fallback Gap Too High Runs`, and `Fallback Same State Runs`.

## Simulator outputs

Each run of `EscortFlowSim_v8.py` produces the following outputs:

- one appended row in the requested CSV file
- optionally, when `-a/--save_raw` is set, one raw pickle trace such as `sim_escort_flow_rawYYYY-MM-DD_HHMMSS.p`
- optionally one log file such as `sim_escort_flow_logYYYY-MM-DD_HHMMSS.txt`

The CSV includes:

- instance and algorithm settings
- `--number_of_requests` and the actual simulation end time
- lead time, waiting time, flow time, and excess time estimates
- confidence-interval half widths
- batching diagnostics after MSER-5 warmup deletion
- `cpu_time`
- separate `Solver Time`, `Greedy Heuristic Time`, `Model Construction Time`, and `Warmstart Selection Time` columns
- `Attention Limit` and `Actual Attention`
- `Algorithm Name` enriched with model choice, `PLPR`, `BnC`, queue cap suffix, and warmstart method when applicable
- total and reason-specific fallback-greedy counters

The reported means and confidence intervals exclude the final cooldown tail using an arrival-based cutoff: find the first request whose departure occurs after the last arrival, then remove all requests that arrived after that request.

Before batching and confidence-interval calculation, the remaining request-level observations also go through warmup deletion using MSER-5. This removes an initial prefix of transient requests so the reported steady-state estimates are based on the post-warmup portion of the run.

`cpu_time` means accumulated algorithm compute time:

- for MILP runs: summed solver time reported by Gurobi
- for greedy runs: summed time spent inside the greedy heuristic
- for mixed runs: both combined

For `v8`, `Model Construction Time` means the persistent Gurobi model build work plus the per-iteration model reset / RHS-update / warmstart-application work needed before each solve. For `--full`, the model is rebuilt only when the greedy-based active horizon changes.

When running interactively, the terminal also prints a one-line progress bar with both arrivals and departures. After the progress bar completes, `v8` prints the raw lead-time average computed on the same request sample used for steady-state reporting after cooldown trimming and MSER-5 warmup deletion, but before batching.

If a simulation is too short for MSER-5 or batch-size selection, the run still completes. The affected steady-state fields are left blank in the CSV, while unrelated fields are still reported.

## Post-process raw traces

You can summarize existing raw pickle traces with:

```bash
python3 CI_Calculation.py -p "sim_escort_flow_raw*.p" -f block_res.csv -H
```

This script:

- reads one or more raw pickle files
- computes lead, waiting, and flow-time steady-state summaries
- applies MSER-5 warmup deletion
- chooses a feasible batch size subject to a lag-1 autocorrelation check

CLI defaults for `CI_Calculation.py`:

- `-p`, `--pickle-file`: one or more pickle files or glob patterns, required
- `-f`, `--csv`: summary CSV output file, default `block_res.csv`
- `-H`, `--header`: write CSV header row, default off

## Animate PBS outputs

`PBSAnimation.py` is the unified animation viewer for this repository. It is
compatible with:

- raw simulation pickles produced by `EscortFlowSim_v8.py`
- exported animation script pickles produced by `EscortFlowStatic.py` with `-a`
- exported animation script pickles produced by `LoadFlowStatic.py` with `-a`

To inspect a raw `EscortFlowSim_v8.py` trace visually:

```bash
python3 PBSAnimation.py sim_escort_flow_raw2026-03-15_143541.p
```

To inspect an exported static script visually:

```bash
python3 PBSAnimation.py script_BM_leave_16_10_8_4_1.p
```

`PBSAnimation.py` uses `tkinter` and opens a desktop GUI.

## Repository layout

- `EscortFlowSim_v8.py`: main simulator
- `escort_flow_gurobi_v8.py`: Gurobi backend used by the `v8` simulator
- `EscortFlowStatic.py`: static escort-flow experiment runner
- `LoadFlowStatic.py`: static load-flow experiment runner
- `CI_Calculation.py`: steady-state analysis from raw traces
- `mser5.py`: MSER-5 warmup deletion helper
- `OneStepHeuristic_v2.py`: greedy step heuristic
- `PBSAnimation.py`: unified animation tool for static and dynamic PBS traces
- `PBSCom.py`: common PBS parsing/helpers
- `PBS_DPHeuristic_lm.py`, `PBS_DPHeuristic_bm.py`: DP-based upper-bound helpers for single-load static runs
- `escort_flow_*.mod`, `pbs_*.mod`: OPL model files
- `run_*.txt`: example command lines used for experiments
- `Junk/`, `kit/`: older or auxiliary copies of scripts

The previous `load_flow_multi.py` script has been moved to `Junk/` for local legacy reference and is no longer part of the tracked repository.

## Notes

- This is a research codebase, not a packaged Python library.
- CSV files are appended to by default.
- Some script headers still refer to earlier versions or older filenames; prefer the current behavior in the code.
- The main simulator assumes the current working directory contains the files it needs.
- Python dependencies are listed in `requirements.txt`, but Gurobi and the large DP pickle files must still be installed or downloaded separately.
- `requirements.txt` pins NumPy to `1.24.4`, which is one of the versions known to have been used successfully for these experiments.

## Changelog

### 2026-07-08 - `OneStepHeuristic_v2.py`: move-extension fixes (termination-proof conformance)

Three related corrections to the "extend a move to promote lower-priority targets" refinement, aligning the implementation with the greedy-heuristic termination analysis in the appendix of the dynamic paper (and the supplement of the static paper):

1. **Direction sign.** An escort move shifts the loads on its path one cell *against* the escort's travel direction `dir`. The extension condition in `extend_move_for_lower_priority_targets` compared a load's preferred direction with `dir` instead of `-dir`, so extensions almost never fired, and when they did, they extended toward loads that the move would push *away* from their outputs (such moves were then typically rejected by the lower-priority-harm guard). The condition now uses `-dir`.
2. **Destination-zone guard.** In the zone-descent operations (B->A, C->B, D->C), an extension could carry the escort's final cell beyond the intended destination zone, which would invalidate the zone-descent step of the progress lemma. Extensions of descent moves are now applied only if the escort's final cell remains in the destination zone (`extension_endpoint_ok` predicate returned by `find_zone_escorts`). Extensions of the promotion move (operation 1, Zone A) remain unrestricted, matching the papers.
3. **Guard fallback.** If an extended move would hurt lower-priority targets, the unextended base move is now tried before the candidate is discarded (previously the entire candidate was skipped, and for the highest-priority load the *extended* harmful move was stored as the forced-progress fallback; the fallback is now the base move).

The public API of `OneStep`/`SolveGreedy` is unchanged; `EscortFlowSim_v8.py` and the static runners need no modification.

Verification (`test_onestep_heuristic.py`): on 100 random instances in each of the two priority modes, the fixed heuristic solves every instance, satisfies the theoretical makespan bound `4*n*d_max` in acyclic mode, and produces legal, conflict-free moves in every time step. Mean flow time, makespan, and movement counts all improved slightly (acyclic mode: mean flow time 32.41 -> 31.50, mean moves 77.25 -> 74.89 on the test battery); about half of the instances follow different trajectories than before, confirming that extensions now activate.

#### Rerunning the experiments affected by this fix

`RunAllDynamic.sh` runs the full battery below in one go (preflight checks, per-family logs, and CSV collection under a timestamped `results_*/` directory); start it inside `tmux` or with `nohup`, as it takes many hours.

The fix changes (slightly improves) the greedy moves, so any simulation in which greedy decisions were executed should be regenerated before quoting its numbers. Because all dynamic runs use the optimality-gap fallback (`-o 0.2`), greedy moves occur in essentially every dynamic experiment. Seeds are hard-coded in the scripts, so before/after differences are attributable to this fix. The CSV outputs append; rename or remove previous outputs first (see the replication guide above). In decreasing order of expected sensitivity:

1. `TestHybridRatio.sh` and `TestAtten.sh` - the hybrid rule applies the greedy heuristic directly, so these meta-parameter studies are the most affected.
2. `FullFactor9x5.sh` and `FullFactor13x7.sh` - the hybrid-factor rows use greedy directly; all other rows use it through the fallback.
3. `TestIntegrated.sh` - real-time modular-vs-integrated comparison; greedy enters through the fallback.
4. The greedy-vs-optimum gap quoted in the papers (mean flow times 30-70% above the optimum in static benchmarks) should be re-estimated with `EscortFlowStatic.py --greedy`.

The static formulation tables (`SingleLoadStatic.sh`, `FourLoadsStatic.sh`) can also be affected: the greedy solution determines the reported greedy bound and the selected planning horizon, and changing that horizon can change solver performance and potentially the best objective within it. Regenerate these results when exact replication with the corrected heuristic is required; see the horizon distinction in [One-step heuristic for static retrieval](#one-step-heuristic-for-static-retrieval).

## Known limitations

- `requirements.txt` is intentionally minimal and only covers Python packages imported by the checked scripts.
- The Python dependency list is pinned only for NumPy; solver and other system-level dependencies are still not captured by a full environment definition.
- Automated tests cover the greedy one-step heuristic and the static lexicographic, weighted certification, warm-start, and experiment-runner workflows. The dynamic simulation pipeline does not yet have comparable automated coverage.
- Gurobi/model failures currently stop the simulation rather than degrading gracefully.
