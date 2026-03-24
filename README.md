# Escort Flow Optimization and Simulation

The repository accompany the working paper by Tal Raviv and Yossi bukchin "Escort-Flow Approach for the Multi-Load Retrieval in Puzzle-Based Storage: Static and Dynamic Approach."  

This repository contains simulation and optimization code for retrieval control in a puzzle-based storage (PBS) system with escorts in static and dynamic environment. The main workflow is:

1. run a dynamic simulation with `EscortFlowSim_v8.py` and collect steady-state statistics directly into a CSV row
2. or run a static optimization instance with `EscortFlowStatic.py` or `LoadFlowStatic.py` to solve a fixed retrieval problem
3. optionally save a raw simulation trace with `-a/--save_raw` and inspect it with `CI_Calculation.py` or `PBSAnimation.py`

The project currently uses `EscortFlowSim_v8.py` as its rolling-horizon dynamic simulator, with Gurobi accessed through the Python API.

## Main entry points

- `EscortFlowSim_v8.py`: primary simulator for dynamic request arrivals, rolling-horizon control, hybrid MILP/greedy policy, CSV reporting, and optional raw pickle export, using Gurobi directly from Python
- `EscortFlowStatic.py`: static escort-flow experiment runner for single-load and multi-load instances, defaulting to the Gurobi Python backend and also supporting greedy-only and naive-lower-bound modes
- `LoadFlowStatic.py`: static load-flow experiment runner; BM is the default, `--lm` switches to LM, and `--gurobi` uses the Gurobi Python API
- `CI_Calculation.py`: post-process one or more raw pickle files and compute steady-state means and confidence intervals using MSER-5 warmup deletion and batch selection
- `PBSAnimation.py`: animate PBS outputs from `EscortFlowStatic.py`, `LoadFlowStatic.py`, and `EscortFlowSim_v8.py`
- `OneStepHeuristic_v2.py`: greedy one-step escort heuristic used in greedy mode and as fallback when MILP results are rejected

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

- Gurobi with the Python API installed for `EscortFlowSim_v8.py`
- IBM ILOG OPL / CPLEX with `oplrun` available on `PATH`

The simulator expects the OPL model files in this repository, especially:

- `escort_flow_bm_rh_v3.mod`
- `escort_flow_bm_rh_static_v3.mod`

The static experiment scripts additionally depend on:

- `pbs_escorts_bm_v3.mod` and related escort-flow OPL files used by `EscortFlowStatic.py`
- `pbs_load_flow_multi.mod` and `pbs_load_flow_multi_lp.mod` used by `LoadFlowStatic.py`
- `PBSCom.py` for common instance-generation and formatting helpers
- `PBS_DPHeuristic_bm.py` and `PBS_DPHeuristic_lm.py` when a DP-based upper bound is requested

For the single-load static experiments, large DP table files are also required:

- `BM_set10x10_k3_corner.p`
- `BM_set16x10_k3_2centers.p`

Those tables are used to derive an upper bound on the number of time steps for the single-load case.

Download source for the large DP files:

- Dropbox archive: <https://www.dropbox.com/scl/fi/ug5eenojzkh4riv8ja1rc/ESCORTS.zip?rlkey=2r01my7aeetl5q0zzn6evzvtt&st=25y3rg7d&dl=0>

Place the extracted `.p` files in the project directory when using the single-load static commands shown below.

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
- escort-flow OPL model files in this repository
- `oplrun`

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
- `--work_limit`: Gurobi work limit in work units, default none
- `--mip_emphasis`: Gurobi MIP emphasis, one of `balanced`, `feasibility`, `optimality`, or `bound`, default `balanced`
- `--dp_file`: DP table file for the single-load case, default empty
- `-k`: `k'` parameter used with the DP heuristic, default `0`
- `--lp`: LP relaxation option, default off
- `--greedy`: solve the static instance with the greedy heuristic instead of a MILP backend, default off; currently supported for `continue` and `leave` modes. In `leave` mode, a target load that reached an output cell becomes an escort at the beginning of the next step
- `--gurobi`: explicitly select the Gurobi Python backend; accepted for clarity but now redundant because Gurobi is the default
- `--opl`: switch to the legacy `oplrun` / CPLEX path instead of the default Gurobi backend
- `--warmstart`: enable heuristic MIP start with the static Gurobi backend, default off
- `--naive`: skip optimization entirely, compute only the naive lower bound for each generated instance, and write a reduced CSV row
- `--lazy` or `--lazy N`: use the lazy-constraint Gurobi backend, default off; a bare `--lazy` means `0`, and `N` is the number of initial time steps kept in the master problem
- `--bnc` or `--bnc N`: use the branch-and-cut Gurobi backend, default off; a bare `--bnc` means a per-separated-node user-cut cap of `2*T`, and `N` sets that cap explicitly
- `-a`: export animation trace, default off

Range syntax:

- ranges use the format `min-max-step`
- example: `-e 8-20-4` means `8, 12, 16, 20`

Notes:

- the default solver path is the Gurobi Python backend; `--opl` switches to the legacy `oplrun` / CPLEX path
- for `stay` mode, the number of output cells must be at least the number of target loads
- the DP-based upper bound currently applies only to the single-load case
- without `--dp_file`, the static horizon upper bound is taken from the greedy heuristic
- on the standard `--gurobi` path, the full static escort-flow model is built explicitly in the master problem
- with any Gurobi backend, `--warmstart` is treated as a fallback: the solver first tries to find an incumbent on its own, and only if that fails does it restart once with the greedy start
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
python3 EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 8:4:20 -l 4
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
- can use the Gurobi Python API with `--gurobi`
- can also export animation traces

Key dependencies:

- `PBSCom.py`
- `PBS_DPHeuristic_lm.py`
- `PBS_DPHeuristic_bm.py`
- `pbs_load_flow_multi.mod`
- `pbs_load_flow_multi_lp.mod`
- `oplrun`

Common arguments:

- `-x`, `-y`: PBS dimensions, required
- `-O`: output cells as coordinate pairs, required
- `-e`: escort-count range, default `5`
- `-r`: replication/seed range, default `1`
- `-l`: number of target loads, default `1`
- `-m`: retrieval mode, default `leave`
- `-f`: CSV result file, default `res_load_flow.csv`
- `--alpha`: makespan weight, default `0.0`
- `--beta`: flowtime weight, default `1.0`
- `--gamma`: movement weight, default `0.01`
- `-T`: legacy horizon scaling factor, default `2.0`; retained in the CLI but retired from the documented BM workflow
- `-t`: solver time limit, default `300`
- `--mip_emphasis`: Gurobi MIP emphasis, one of `balanced`, `feasibility`, `optimality`, or `bound`, default `balanced`
- `--lm` or `--LM`: run LM instead of the default BM mode, default off
- `--dp_file`: DP table file for single-load upper bounds, default empty
- `-k`: `k'` parameter for the DP heuristic, default `0`
- `--lp`: solve the LP relaxation instead of the ILP, default off
- `--gurobi`: solve with the Gurobi Python API instead of `oplrun`, default off
- `-a`: export animation trace, default off

Notes:

- if `--gurobi` is omitted, the default solver path is CPLEX through `oplrun`
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
python3 LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 8-20-4 -l 4
```

To solve the LP relaxation for either static model, add `--lp`. Using a larger horizon can make the LP lower bound slightly weaker.

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

## Known limitations

- `requirements.txt` is intentionally minimal and only covers Python packages imported by the checked scripts.
- The Python dependency list is pinned only for NumPy; solver and other system-level dependencies are still not captured by a full environment definition.
- There is no automated test suite in the repository.
- Gurobi/model failures currently stop the simulation rather than degrading gracefully.
