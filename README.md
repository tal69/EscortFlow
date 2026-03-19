# Escort Flow Simulation

This repository contains simulation and optimization code for retrieval control in a puzzle-based storage (PBS) system with escorts. The main workflow is:

1. run a dynamic simulation with `EscortFlowSim_v7.py`
2. collect steady-state statistics directly into a CSV row
3. optionally save a raw trace with `-a/--save_raw` and inspect it with `CI_Calculation.py` or `PBSAnimation.py`

The project currently uses `EscortFlowSim_v7.py` as its rolling-horizon dynamic simulator, with Gurobi accessed through the Python API.

## Main entry points

- `EscortFlowSim_v7.py`: primary simulator for dynamic request arrivals, rolling-horizon control, hybrid MILP/greedy policy, CSV reporting, and optional raw pickle export, using Gurobi directly from Python
- `EscortFlowStatic.py`: static escort-flow ILP experiment runner for single-load and multi-load instances
- `load_flow_multi.py`: static load-flow ILP experiment runner for the corresponding benchmark instances
- `CI_Calculation.py`: post-process one or more raw pickle files and compute steady-state means and confidence intervals using MSER-5 warmup deletion and batch selection
- `PBSAnimation.py`: animate PBS outputs from `EscortFlowStatic.py`, `load_flow_multi.py`, and `EscortFlowSim_v7.py`
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

- Gurobi with the Python API installed for `EscortFlowSim_v7.py`
- IBM ILOG OPL / CPLEX with `oplrun` available on `PATH`

The simulator expects the OPL model files in this repository, especially:

- `escort_flow_bm_rh_v3.mod`
- `escort_flow_bm_rh_static_v3.mod`

The static experiment scripts additionally depend on:

- `pbs_escorts_bm_v3.mod` and related escort-flow OPL files used by `EscortFlowStatic.py`
- `pbs_load_flow_multi.mod` and `pbs_load_flow_multi_lp.mod` used by `load_flow_multi.py`
- `PBSCom.py` for common instance-generation and formatting helpers
- `PBS_DPHeuristic_bm.py` and `PBS_DPHeuristic_lm.py` when a DP-based upper bound is requested

For the single-load experiments described in the paper, large DP table files are also required. The appendix in the linked paper source refers to:

- `BM_set10x10_k3_corner.p`
- `BM_set16x10_k3_2centers.p`

Those tables are used to derive an upper bound on the number of time steps for the single-load case.

Download source for the large DP files:

- Dropbox archive: <https://www.dropbox.com/scl/fi/ug5eenojzkh4riv8ja1rc/ESCORTS.zip?rlkey=2r01my7aeetl5q0zzn6evzvtt&st=25y3rg7d&dl=0>

Place the extracted `.p` files in the project directory when using the single-load static commands shown below.

## Static optimization scripts

The repository contains two older but still useful experiment runners for the static setting:

- `EscortFlowStatic.py`: escort-flow formulation
- `load_flow_multi.py`: load-flow formulation

These scripts are the ones used to reproduce the static benchmark experiments discussed in Appendix D of the paper source.

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

- `-x`, `-y`: PBS dimensions
- `-O`: output cells as coordinate pairs
- `-e`: escort-count range
- `-r`: replication/seed range
- `-l`: number of target loads
- `-m`: retrieval mode, one of `stay`, `leave`, `continue`
- `-f`: CSV result file
- `-T`: horizon scaling factor when no DP upper bound is provided
- `-t`: solver time limit
- `--dp_file`: DP table file for the single-load case
- `-k`: `k'` parameter used with the DP heuristic
- `--lp`: LP relaxation option
- `--greedy`: solve the static instance with the greedy heuristic instead of OPL; currently supported for `continue` and `leave` modes. In `leave` mode, a target load that reached an output cell becomes an escort at the beginning of the next step
- `-a`: export animation trace

Range syntax:

- ranges use the format `min-max-step`
- example: `-e 8-20-4` means `8, 12, 16, 20`

Notes:

- for `stay` mode, the number of output cells must be at least the number of target loads
- the DP-based upper bound currently applies only to the single-load case
- for multi-load runs, the script guesses the horizon length from the instance dimensions and a scaling factor

Single-load examples from the appendix:

```bash
python3 EscortFlowStatic.py -x 10 -y 10 -O 0 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set10x10_k3_corner.p -b
```

```bash
python3 EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set16x10_k3_2centers.p -b
```

Multi-load example:

```bash
python3 EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 3-8 -l 4
```

For multi-load experiments, drop `-k` and `--dp_file`. If the guessed time horizon is too small and the model becomes infeasible, increase `-T`.

### `load_flow_multi.py`

Purpose:

- runs the static load-flow formulation
- supports the same benchmark family as `EscortFlowStatic.py`
- can solve the ILP or LP relaxation
- can also export animation traces

Key dependencies:

- `PBSCom.py`
- `PBS_DPHeuristic_lm.py`
- `PBS_DPHeuristic_bm.py`
- `pbs_load_flow_multi.mod`
- `pbs_load_flow_multi_lp.mod`
- `oplrun`

Common arguments:

- `-x`, `-y`: PBS dimensions
- `-O`: output cells as coordinate pairs
- `-e`: escort-count range
- `-r`: replication/seed range
- `-l`: number of target loads
- `-m`: retrieval mode
- `-f`: CSV result file
- `-T`: horizon scaling factor
- `-t`: solver time limit
- `-b`: block-movement regime
- `--dp_file`: DP table file for single-load upper bounds
- `-k`: `k'` parameter for the DP heuristic
- `--lp`: solve the LP relaxation instead of the ILP
- `-a`: export animation trace

Notes:

- the script header says only `leave` is supported at present; that is the safe mode to use
- as in `EscortFlowStatic.py`, the DP table route is for the single-load case
- without a DP table, the script guesses a time horizon from the problem dimensions and `-T`

Single-load examples from the appendix:

```bash
python3 load_flow_multi.py -x 10 -y 10 -O 0 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set10x10_k3_corner.p -b
```

```bash
python3 load_flow_multi.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set16x10_k3_2centers.p -b
```

Multi-load example:

```bash
python3 load_flow_multi.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 3-8 -l 4
```

To solve the LP relaxation for either static model, add `--lp`. As noted in the appendix, using a larger horizon can make the LP lower bound slightly weaker.

### DP helper modules

`PBS_DPHeuristic_lm.py` and `PBS_DPHeuristic_bm.py` are helper modules used by the static scripts when `--dp_file` is supplied. They:

- load a precomputed DP table from a pickle file
- generate an upper-bound policy for the single-load case
- provide a feasible horizon estimate before the full OPL model is solved

These helpers depend on:

- `PBSCom.py`
- a compatible DP table pickle

They are not the main experiment entry points; normally you use them through `EscortFlowStatic.py` or `load_flow_multi.py`.

## Typical workflow

### 1. Run a simulation

#### Current version: `EscortFlowSim_v7.py`

`v7` is the maintained dynamic simulator. It uses an in-process Gurobi model and writes the CSV/raw outputs used by the rest of this repository.

Example greedy-only run:

```bash
python3 EscortFlowSim_v7.py -x 9 -y 5 -O 4 0 -e 4 -R 0.4 \
  -S 1600 --greedy -f results.csv -H
```

Example rolling-horizon Gurobi run:

```bash
python3 EscortFlowSim_v7.py -x 9 -y 5 -O 4 0 -e 4 -S 2000 \
  -T 5 -I 1 -E 1 -R 0.2 -M 6 -q spt -L -f results.csv -H
```

Example hybrid run with multi-step epochs:

```bash
python3 EscortFlowSim_v7.py -x 9 -y 5 -O 4 0 -e 4 -S 2000 \
  -T 6 -I 2 -E 2 -R 0.2 -M 6 --hybrid --hybrid_ratio 1.0 -L -f results.csv -H
```

Key notes for `v7`:

- `-t` is the per-solve Gurobi time limit
- `--num_threads` controls Gurobi threads
- `v7` now builds the rolling-horizon Gurobi model once and reuses it across epochs by updating only the dynamic state-dependent data; this removes the previous full model rebuild at every solve
- `v7` uses balanced Gurobi search by default (`MIPFocus=0`), and switches to optimality emphasis (`MIPFocus=2`) when `--warmstart` is enabled; the solver MIP gap target remains `0.01%`
- `--warmstart` now accepts explicit modes: `greedy`, `ilp`, or `ilp greedy`; by default all unused arc variables are set to `0`, and `nozero` disables that
- when a warm start is used, `v7` clears any previous Gurobi solution / MIP-start state before applying the requested start vector, so warmstart behavior remains explicit instead of leaking across solves
- `--warmstart` is allowed only when `--integer_horizon == --fractional_horizon`; if `--integer_horizon` is omitted, it defaults to `--fractional_horizon`
- no OPL temp files are created
- the interactive progress line now shows both arrivals and departures, and after the run finishes it prints the raw average lead time after cooldown trimming and warmup deletion

Important options:

- `-x`, `-y`: PBS dimensions
- `-O`: output cells in pairs of coordinates; ranges are supported through `PBSCom` formatting
- `-e`: number of escorts
- `-S`: number of requests in the simulation
- `-R`: Poisson request arrival rate per time step
- `-E`, `--epoch`: execution epoch length
- `-T`: fractional horizon for the MILP; defaults to `max(--epoch + 4, --integer_horizon)`
- `-I`: integer horizon for the MILP
- `-t`: per-solve Gurobi time limit in seconds
- `-M`: maximum number of target loads considered concurrently; if omitted it defaults internally to `Lx*Ly`
- `-q`: queue management, either `fifo` or `spt`
- `-m`, `--offline`: offline rolling horizon for pure ILP mode only; in offline mode the MILP sees requests visible at the current decision time, and the flag cannot be combined with `--greedy` or `--hybrid`
- `--full`: use the full MILP instead of the surrogate model
- `--greedy`: greedy heuristic only; currently requires `--epoch 1`
- `--hybrid`: at each decision epoch, split the currently open visible requests into `old` requests already visible at the beginning of the previous epoch and `new` requests that were not visible then; use greedy when `number_of_new_open_requests >= number_of_old_open_requests * --hybrid_ratio`
- `--hybrid_ratio`: non-negative real ratio used by `--hybrid` (default `1.0`)
- `--num_threads`: Gurobi thread count; the default is `8` on macOS and `12` on Linux
- `-L`: write a detailed log file
- `-a`, `--save_raw`: save a raw pickle trace for post-processing and animation
- `-f`: output CSV file
- `-H`: write CSV header if needed

`v7` timing notes:

- The simulator plans once per epoch and then executes the resulting move list for that entire epoch.
- MILP visibility is delayed by one epoch: requests must be visible by the beginning of the previous epoch to enter the next MILP solve.
- With `--offline`, the MILP instead uses the requests visible at the current decision time.
- In hybrid mode, `old requests` are the currently open requests that were already visible at the beginning of the previous epoch, and `new requests` are the currently open requests that were not visible then.
- If an ILP solution is unusable, `v7` falls back to greedy on all target loads visible at the current decision time and fills the whole epoch greedily.

## Simulator outputs

Each run of `EscortFlowSim_v7.py` produces the following outputs:

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
- separate `Solver Time`, `Greedy Heuristic Time`, `Model Construction Time`, and `Warmstart Selection Time` columns for `v7`

The reported means and confidence intervals exclude the final cooldown tail using an arrival-based cutoff: find the first request whose departure occurs after the last arrival, then remove all requests that arrived after that request.

Before batching and confidence-interval calculation, the remaining request-level observations also go through warmup deletion using MSER-5. This removes an initial prefix of transient requests so the reported steady-state estimates are based on the post-warmup portion of the run.

`cpu_time` means accumulated algorithm compute time:

- for MILP runs: summed solver time reported by Gurobi
- for greedy runs: summed time spent inside the greedy heuristic
- for mixed runs: both combined

For `v7`, `Model Construction Time` now means the one-time persistent Gurobi model build plus the per-iteration model reset / RHS-update / warmstart-application work needed before each solve; it no longer reflects a full from-scratch rebuild at every epoch.

When running interactively, the terminal also prints a one-line progress bar with both arrivals and departures. After the progress bar completes, `v7` prints the raw lead-time average computed on the same request sample used for steady-state reporting after cooldown trimming and MSER-5 warmup deletion, but before batching.

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

## Animate PBS outputs

`PBSAnimation.py` is the unified animation viewer for this repository. It is
compatible with:

- raw simulation pickles produced by `EscortFlowSim_v7.py`
- exported animation script pickles produced by `EscortFlowStatic.py` with `-a`
- exported animation script pickles produced by `load_flow_multi.py` with `-a`

To inspect a raw `EscortFlowSim_v7.py` trace visually:

```bash
python3 PBSAnimation.py sim_escort_flow_raw2026-03-15_143541.p
```

To inspect an exported static script visually:

```bash
python3 PBSAnimation.py script_BM_leave_16_10_8_4_1.p
```

`PBSAnimation.py` uses `tkinter` and opens a desktop GUI.

## Repository layout

- `EscortFlowSim_v7.py`: main simulator
- `EscortFlowStatic.py`: static escort-flow experiment runner
- `load_flow_multi.py`: static load-flow experiment runner
- `CI_Calculation.py`: steady-state analysis from raw traces
- `mser5.py`: MSER-5 warmup deletion helper
- `OneStepHeuristic_v2.py`: greedy step heuristic
- `PBSAnimation.py`: unified animation tool for static and dynamic PBS traces
- `PBSCom.py`: common PBS parsing/helpers
- `PBS_DPHeuristic_lm.py`, `PBS_DPHeuristic_bm.py`: DP-based upper-bound helpers for single-load static runs
- `escort_flow_*.mod`, `pbs_*.mod`: OPL model files
- `run_*.txt`: example command lines used for experiments
- `Junk/`, `kit/`: older or auxiliary copies of scripts

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
