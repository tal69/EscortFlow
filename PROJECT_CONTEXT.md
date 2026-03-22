# Project Context

This file is a working context snapshot for future sessions on this repo.

## Active Tracks

- Dynamic simulator:
  - Active track is `EscortFlowSim_v8.py` with backend `escort_flow_gurobi_v8.py`.
  - `v7` stays in the repo but is not the documented active simulator.
- Static runners:
  - `EscortFlowStatic.py`
  - `LoadFlowStatic.py`

## Dynamic `v8` Status

- `v8` is Gurobi-only. There is no OPL/CPLEX path in `EscortFlowSim_v8.py`.
- `--time_limit` in `v8` is `float`.
- `--num_threads` defaults to `8` on macOS and `12` on Linux.
- User-facing terminology uses `attention` instead of `balls in the air`.

### Dynamic model behavior

- Surrogate model:
  - Controlled by `-T`, `-I`, `-E`.
  - `--bnc` is supported.
- Full model:
  - Full horizon is based on the greedy heuristic horizon.
  - With no `-I`, the model is `PLPR`: only the first `-E` periods are integer and the rest are fractional.
  - With any `-I`, all periods in the full model are integer.
  - `-T` is ignored in `--full`.
  - `--bnc` is supported.

### Dynamic BnC behavior

- Implemented in `escort_flow_gurobi_v8.py`.
- Root user cuts only.
- Lazy constraints enforce incumbent feasibility.
- In `PLPR`, user cuts are generated only for integer periods; lazy constraints still check the full horizon.
- Default per-callback cut budget is `2*T`.

### Dynamic warmstart behavior

- Warmstart is reported in CSV `Algorithm Name` as:
  - `ws_greedy`
  - `ws_ilp`
  - `ws_both`
- Warmstart is intended as fallback support rather than strong search guidance.

### Dynamic CSV/reporting notes

- `Algorithm Name` includes:
  - `BnC` when used
  - `PLPR` for full model without `-I`
  - warmstart label if active
- Fallback greedy runs are split into:
  - `Fallback No Feasible Runs`
  - `Fallback Gap Too High Runs`
  - `Fallback Same State Runs`

## Static Solver Status

- Static scripts now default to Gurobi.
- `--opl` switches to the legacy OPL/CPLEX path.
- `--cplex` is accepted as a hidden alias for `--opl`.
- Both static scripts set solver threads to:
  - `12` on Linux
  - `8` on macOS
- This is applied to:
  - Gurobi Python backends
  - OPL/CPLEX models via a `threads` data parameter

## Static BnC Status

Implemented in `escort_flow_static_bnc.py`.

### What stays in the master

- Flow conservation constraints for the selected retrieval mode
- `(3)` output stay
- `(4)` and `(5)` supply
- `(7)` escort conflict constraints
- `(9)` stay/move implication for all times
- `(8)` explicit only for the first `T // 8` time steps
- `(10)` retrieval completion
- `q` aggregation constraints

### What is separated

- Only the later-time part of family `(8)` is separated.

### When cuts are added

- `MIPNODE` user cuts:
  - root node: uncapped
  - early non-root nodes: allowed while explored node count `< 15`
  - non-root nodes: capped at `2*T` by default, using the most violated cuts first
- `MIPSOL` lazy constraints:
  - capped at `4*T`

### Static CSV note

- In `EscortFlowStatic.py`, `Model` is now `ILP-Gurobi-BnC`.
- The numeric per-node cap is recorded in a separate CSV column:
  - `Max User Cut Per Node`

## Constraint (6) Status

- Constraint `(6)` has been removed from all escort-flow Gurobi models, static and dynamic:
  - `escort_flow_static_gurobi.py`
  - `escort_flow_static_lazy.py`
  - `escort_flow_static_bnc.py`
  - `escort_flow_gurobi.py`
  - `escort_flow_gurobi_v8.py`
- It was removed because current testing indicated it is redundant for objective/bound purposes in the tested cases and may reduce solve time.
- OPL files were not changed for this removal.

## Documentation Status

- `README.md` tracks `v8` as the active dynamic simulator.
- `v7` references were removed from `README.md`.
- The README documents:
  - `v8` dynamic behavior
  - surrogate vs full
  - `PLPR`
  - `--bnc`
  - warmstart naming
  - `attention` terminology
  - fallback-reason CSV columns

## Useful Reminder For Next Session

If restarting work, first read:

- `README.md`
- `PROJECT_CONTEXT.md`
- `EscortFlowSim_v8.py`
- `escort_flow_gurobi_v8.py`
- `EscortFlowStatic.py`
- `escort_flow_static_bnc.py`
