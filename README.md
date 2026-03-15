# Escort Flow Simulation

This repository contains simulation and optimization code for retrieval control in a puzzle-based storage (PBS) system with escorts. The main workflow is:

1. run a dynamic simulation with `EscortFlowSim_v5.py`
2. collect steady-state statistics directly into a CSV row
3. optionally inspect the raw trace with `CI_Calculation.py` or animate it with `FlowAnimationSim.py`

The project mixes Python simulation code with IBM ILOG OPL models (`.mod` / `.dat`) used by the rolling-horizon MILP controller.

## Main entry points

- `EscortFlowSim_v5.py`: primary simulator for dynamic request arrivals, rolling-horizon control, greedy fallback, CSV reporting, and raw pickle export
- `EscortFlow_v3.py`: static escort-flow ILP experiment runner for single-load and multi-load instances
- `load_flow_multi.py`: static load-flow ILP experiment runner for the corresponding benchmark instances
- `CI_Calculation.py`: post-process one or more raw pickle files and compute steady-state means and confidence intervals using MSER-5 warmup deletion and batch selection
- `FlowAnimationSim.py`: visualize a raw simulation trace
- `OneStepHeuristic_v2.py`: greedy one-step escort heuristic used in greedy mode and as fallback when MILP results are rejected

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

- IBM ILOG OPL / CPLEX with `oplrun` available on `PATH`

The simulator also expects the OPL model files in this repository, especially:

- `escort_flow_bm_rh_v3.mod`
- `escort_flow_bm_rh_static_v3.mod`
- `escort_flow_sim.dat`

The static experiment scripts additionally depend on:

- `pbs_escorts_bm_v3.mod` and related escort-flow OPL files used by `EscortFlow_v3.py`
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

- `EscortFlow_v3.py`: escort-flow formulation
- `load_flow_multi.py`: load-flow formulation

These scripts are the ones used to reproduce the static benchmark experiments discussed in Appendix D of the paper source.

### `EscortFlow_v3.py`

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
python3 EscortFlow_v3.py -x 10 -y 10 -O 0 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set10x10_k3_corner.p -b
```

```bash
python3 EscortFlow_v3.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 3-8 -k 3 \
  --dp_file BM_set16x10_k3_2centers.p -b
```

Multi-load example:

```bash
python3 EscortFlow_v3.py -x 16 -y 10 -O 4 0 11 0 -r 1-100 -m leave -e 3-8 -l 4
```

For multi-load experiments, drop `-k` and `--dp_file`. If the guessed time horizon is too small and the model becomes infeasible, increase `-T`.

### `load_flow_multi.py`

Purpose:

- runs the static load-flow formulation
- supports the same benchmark family as `EscortFlow_v3.py`
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
- as in `EscortFlow_v3.py`, the DP table route is for the single-load case
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

They are not the main experiment entry points; normally you use them through `EscortFlow_v3.py` or `load_flow_multi.py`.

## Typical workflow

### 1. Run a simulation

Example greedy-only run:

```bash
python3 EscortFlowSim_v5.py \
  -x 9 -y 5 \
  -O 4 0 \
  -e 4 \
  -R 0.4 \
  -S 1600 \
  --greedy \
  -f results.csv \
  -H
```
This will run a simulation of the RTRH drove by the greedy heuristic.
The simulation is of 9x5 with 4 escorts PBS unit with an output cell at (0,4) which is the middle of the bottom wall. 
New random request will arrive at rate of 0.4 per tie step. The simulation runs until the 1600^th request is ejected at
an output cell. The summary statistics of the simulation results are written to "results.csv." 

Example rolling-horizon MILP run:

```bash
python3 EscortFlowSim_v5.py \
  -x 9 -y 5 \
  -O 4 0 \
  -e 4 \
  -S 2000 \
  -T 5 \
  -I 1 \
  -E 1 \
  -R 0.2 \
  -M 6 \
  -q spt \
  -L \
  -f results.csv \
  -H
```

Important options:

- `-x`, `-y`: PBS dimensions
- `-O`: output cells in pairs of coordinates; ranges are supported through `PBSCom` formatting
- `-e`: number of escorts
- `-S`: number of requests in the simulation
- `-R`: Poisson request arrival rate per time step
- `-E`: execution horizon
- `-T`: fractional horizon for the MILP
- `-I`: integer horizon for the MILP
- `-t`: per-solve OPL/CPLEX time limit in seconds
- `-M`: maximum number of target loads considered concurrently
- `-q`: queue management, either `fifo` or `spt`
- `-m`: offline mode
- `--static`: use the full static MILP; requires `-m`
- `--greedy`: greedy heuristic only; forces offline behavior and skips OPL
- `-L`: write a detailed log file
- `-f`: output CSV file
- `-H`: write CSV header if needed

## Simulator outputs

Each run of `EscortFlowSim_v5.py` produces:

- one appended row in the requested CSV file
- one raw pickle trace such as `sim_escort_flow_rawYYYY-MM-DD_HHMMSS.p`
- optionally one log file such as `sim_escort_flow_logYYYY-MM-DD_HHMMSS.txt`

The CSV includes:

- instance and algorithm settings
- lead time, waiting time, flow time, and excess time estimates
- confidence-interval half widths
- batching diagnostics after MSER-5 warmup deletion
- `cpu_time`

`cpu_time` currently means accumulated algorithm compute time:

- for MILP runs: summed solver time reported by OPL/CPLEX
- for greedy runs: summed time spent inside the greedy heuristic
- for mixed runs: both combined

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

## Animate a simulation

To inspect a raw trace visually:

```bash
python3 FlowAnimationSim.py sim_escort_flow_raw2026-03-15_143541.p
```

`FlowAnimationSim.py` uses `tkinter` and opens a desktop GUI.

## Repository layout

- `EscortFlowSim_v5.py`: main simulator
- `EscortFlow_v3.py`: static escort-flow experiment runner
- `load_flow_multi.py`: static load-flow experiment runner
- `CI_Calculation.py`: steady-state analysis from raw traces
- `mser5.py`: MSER-5 warmup deletion helper
- `OneStepHeuristic_v2.py`: greedy step heuristic
- `FlowAnimationSim.py`: animation for simulation traces
- `PBSCom.py`: common PBS parsing/helpers
- `PBS_DPHeuristic_lm.py`, `PBS_DPHeuristic_bm.py`: DP-based upper-bound helpers for single-load static runs
- `escort_flow_*.mod`, `pbs_*.mod`: OPL model files
- `run_*.txt`: example command lines used for experiments
- `Junk/`, `kit/`: older or auxiliary copies of scripts

## Notes

- This is a research codebase, not a packaged Python library.
- CSV files are appended to by default.
- Some script headers still refer to earlier versions or older filenames; prefer the current behavior in the code.
- The main simulator assumes the current working directory contains the OPL model and data files it needs.
- Python dependencies are listed in `requirements.txt`, but OPL/CPLEX and the large DP pickle files must still be installed or downloaded separately.
- `requirements.txt` pins NumPy to `1.24.4`, which is one of the versions known to have been used successfully for these experiments.

## Known limitations

- `requirements.txt` is intentionally minimal and only covers Python packages imported by the checked scripts.
- The Python dependency list is pinned only for NumPy; OPL/CPLEX and other system-level dependencies are still not captured by a full environment definition.
- There is no automated test suite in the repository.
- OPL failures currently stop the simulation rather than degrading gracefully.
