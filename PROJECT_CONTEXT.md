# Project Context

This file is a working context snapshot for future sessions on this repo.

## Active Tracks

- Dynamic simulator:
  - Active track is `EscortFlowSim_v8.py` with backend `escort_flow_gurobi_v8.py`.
  - `v7` has been moved out of the tracked repo into local junk storage.
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
- Both static scripts expose `--num_threads` and default it to:
  - `12` on Linux
  - `8` on macOS
- This is applied to:
  - Gurobi Python backends
  - OPL/CPLEX models via a `threads` data parameter

## Two-Phase Static Version (October 7, 2026)

- New entry points: `EscortFlowStaticLex.py` and `LoadFlowStaticLex.py`. Existing
  static entry points also accept `--lexicographic`; weighted mode remains the default.
- Shared implementation: `static_lexicographic.py`, two explicit Gurobi solves.
  First minimize integer flow time, then fix its best incumbent by equality and
  minimize integer load movements. Both phases set `MIPGap=0`, `MIPGapAbs=1`.
- `--phase1_time_limit` caps phase 1. `--time_limit` / `--total_time_limit` caps
  total solver runtime. Phase 2 receives the total minus actual phase-one Gurobi
  `Runtime`. Work limits are shared in the same manner.
- Phase 2 may follow a time-limited phase 1 with an incumbent; CSV distinguishes
  this from proven lexicographic optimality, which requires proof in both phases.
  No incumbent or no remaining shared budget skips phase 2. The best phase-one
  plan is preserved if phase 2 obtains no solution.
- Standard, lazy, and BnC escort backends are supported; load-flow supports its
  existing leave-mode BM and LM models. LP/OPL and weighted cutoffs are excluded.
- Existing automatic horizons are preserved. `--horizon` overrides them in either
  weighted or lexicographic mode; optimality is for the chosen finite horizon.
- Fixed first-step arrival omission in escort-flow reporting and the BnC zero-cut
  cap at `T=0`. Load-flow nonstay arc detection is now independent of gamma.
- `test_static_lexicographic.py` covers independently checked small optima, all
  four backends, strict gaps, timing/work allocation, and incumbent retention.
  Local licensed interpreter: `/Users/talraviv/miniconda3/bin/python3.11`.
- This implements the solution-method part of planned R2.3. The full comparative
  experiment and corresponding manuscript/response-letter revision remain pending.
- `RunTable2Lex.sh` now prepares the Linux experiment for current Table 2(a,b):
  13x7, 10x10, 16x10, and 27x10 only; one target with 3--8 escorts and four targets
  with 8/12/16 escorts; seeds 1--100; both formulations, leave/BM. It sets phase 1
  to 270 seconds and the shared total to 300 seconds, defaults to 12 threads, and
  runs sequentially. It writes four merged CSVs plus per-layout CSVs/logs and an
  environment/command record in a new results directory. Full experiment not
  launched during script creation; command scope and orchestration were checked.
  The older static shell scripts retain an additional 9x5 layout outside Table 2.

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
- `v7` has been removed from the tracked repo and from `README.md`.
- The README documents:
  - `v8` dynamic behavior
  - surrogate vs full
  - `PLPR`
  - `--bnc`
  - warmstart naming
  - `attention` terminology
  - fallback-reason CSV columns

## Paper Source Location

- The paper LaTeX/Overleaf sync directory for this project is:
  `/Users/talraviv/Library/CloudStorage/Dropbox/.metadata/Apps/Overleaf/Escort flow formulation for the PBS retrieval problem`

## Useful Reminder For Next Session

If restarting work, first read:

- `README.md`
- `PROJECT_CONTEXT.md`
- `EscortFlowSim_v8.py`
- `escort_flow_gurobi_v8.py`
- `EscortFlowStatic.py`
- `escort_flow_static_bnc.py`

## IISE R1 Revision Plan (October 7, 2026)

Tal asked to remember the following plan for future work. Do not implement it until he instructs us to proceed. This records planned work, not completed revisions.

- **R1.1:** Justify the importance of puzzle-based storage (PBS) using the literature.
- **R1.2:** Explain that the claim that we studied only a single target load is incorrect. Consider larger instances.
- **R1.3:** Formally prove that the two formulations are equivalent and that a formulation solution prescribes a physically feasible solution.
- **R2.1:** Conduct and report an experiment with continue mode.
- **R2.2:** Attempt to formally prove LP relaxation dominance. If a proof cannot be established, soften the claim.
- **R2.3:** Soften the claim that alpha = 0.01 guarantees lexicographic optimization. Implement a two-phase solution method that enforces the lexicographic order and compare the solutions.
- **R2.4:** Add experiments with two and three target loads.
- **R2.5:** Tal's interpretation is that the referee misread the tables: with more escorts, both formulations perform well, but the load-flow formulation obtains optimal solutions in less time. Explain this comparison in the response, citing the relevant table entries when drafting.
- **R2.7:** Check whether the multi-load SBM load-flow formulation and the relevant strengthening were already included in the conference publication by Bukchin and Raviv (2023). Address their provenance explicitly in the main text either way, citing the earlier publication for material already present and identifying any extensions or strengthening introduced in the current paper.

### R2.2 Proof Attempt (October 7, 2026)

Tal subsequently authorized work on the LP dominance proof, allowing minimum grid-size or density assumptions if needed. An explicit affine projection proof for the intended SBM models is saved in `revision_R1/lp_projection_proof_analytical.txt`. It requires no size or density threshold and preserves the target arrival objective and total load-movement objective, including the load-flow equality A.13. It proves weak LP dominance, not a strict advantage for every instance or a guarantee of faster integer optimization.

The proof assumes that targets initially at outputs have been handled separately, idle escort arcs are included in the involved-cell inequalities, and target movement arcs originating at outputs are zero. The last restriction is implemented in the escort-flow code but is omitted from the printed main-text formulation. The printed load-flow continue-mode adaptation also incorrectly retains a leave-only output-capacity constraint; this must be removed or replaced consistently when targets become blocking loads. These output conventions must be reconciled before inserting a general theorem into the manuscript. The proof is for SBM and does not establish dominance for the separately formulated SLM relaxation. This proof attempt does not implement the other revision-plan items.

Exact arithmetic verification is saved in `revision_R1/lp_dominance_verification.py` and `.json`. The literal printed leave-mode EF admits a binary artificial-target counterexample on a 4x4 grid at 93.75% utilization, with objective 3.02 versus a load-flow lower bound of 6.06. All printed EF constraints pass; the omitted output-origin target restriction eliminates it. A separate strict-dominance certificate for the intended formulations at the same size and utilization has LF objective at most 4.58 and EF objective at least 5.04, with both models feasible over the same 12-step horizon. These values are feasible certificates and analytical bounds, not computed optimal LP values. A genuinely fractional leave-mode projection also passes exact constraint checks and preserves objective 4.045.

### R2.2 Manuscript Integration (October 7, 2026)

At Tal's request, the SBM LP-dominance result is now Theorem 2 in Section 3.3 of the canonical Overleaf source `static_escort_flow_IISE.tex`, and its full projection proof is in new Section E of `static_escort_flow_IISE_supplement.tex`. Both files are in the paper source directory listed above. The sources now explicitly include idle escort arcs in involved-cell sets, prohibit nonidle target arcs from outputs, and state the initial-target/output convention. The supplement explicitly absorbs target flow at outputs, clarifies leave-mode output service, removes the leave-only output inequalities when adapting to continue mode, and treats A.13 as a conservation identity under these conventions. The theorem states weak LP dominance for SBM in both retrieval modes, without size or density thresholds.

The built-in standalone compiler could not resolve this project's external figures and auxiliary files. Both actual project sources passed local `pdflatex -draftmode` and bibliography/reference checks with all new references resolved, without producing another PDF. Their existing auxiliary reference files were updated to resolve the new theorem and Section E. Existing citation duplication and layout warnings outside the added proof remain. The open response letter was not edited during this integration.

Tal then requested the corresponding response-letter update. The existing `revision_R1/response_letter_R1.tex` now reports Theorem 2 and Supplement Section E as completed in AE.2, R1.2, R1.3, R2.2, R2.7, and R3.1, with consistent introductory progress notes. It explains the output-boundary clarifications and A.13, distinguishes weak SBM LP dominance from computational performance and SLM, and keeps physical recovery and reverse integer equivalence pending. R2.7 retains the requested check of Bukchin and Raviv (2023). All 26 reviewer-comment quotations were preserved exactly. The built-in editor compiler confirmed successful compilation of the updated response letter.
