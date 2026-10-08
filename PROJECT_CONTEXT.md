# Project Context

This file is a working context snapshot for future sessions on this repo.

The remaining review tasks and current experiment/reporting decisions are
maintained in `tasks.md`. Read that checklist before resuming revision work;
it includes decisions that supersede older proposed plans in the response letter.

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

## Static Objective Protocols (October 7-8, 2026)

### Earlier fixed-weight strategy (October 8, 2026)

The user first superseded the two-phase experiment with a fixed-weight strategy.
That version uses the weighted objective `F + 0.01*M`; minimum flow time is
an instance-specific empirical certificate, not a claimed universal property.
The user authorized experiment code and method notes first, without further
manuscript or reviewer-response edits at this stage.

`RunTable2Weighted.sh` and `RunWeightedStatic.py` implement that separate protocol:
300 seconds for integer-scaled `100*F+M`, with a separate configurable 300-second
pure-flow certificate budget. Complete greedy starts are enabled for both
weighted formulations. Certification receives the saved weighted solution, with
truncation or idle padding to fit its horizon. Load-flow removes redundant terminal
blocker moves at the boundary before padding or truncation because terminal
destination constraints are absent from the shorter model. This preserves all
retrieval times; the original weighted result is retained. The CSV records the
start source and removed-movement count.
The certificate uses T=H directly, even if smaller than the weighted horizon.
H uses weighted flow time, not greedy flow time. The user requested skipping
unproven weighted candidates unless
their absolute weighted gap is significantly below one; the configurable default
gate is strictly below 0.1 in original weighted units. Solver gaps use 0.999 in
integer objective units, with independent guarded proof flags. Certification
sets a sufficient horizon covering all no-worse-flow candidates. Weighted results and
certificate results remain separate, including counterexamples and proof scope.
The protocol and sufficient-small-weight argument are documented in
`weighted_flow_certification_notes.md`. The full experiment is not yet run.

`LoadFlowStatic.py --warmstart` now supports the same physical BM leave greedy
trace as escort-flow, translated into every target/blocker arc, retrieval, and
makespan variable. Unsupported modes are rejected. Existing weighted defaults
and the earlier two-phase runners remain available.

Validation: 62 unit/integration tests passed across the lexicographic, LF
warm-start, weighted certificate, new runner, and weighted solution-transfer
suites. Starts were checked with every variable fixed, including shorter and
longer horizons. Six tiny two-target instances (three per formulation) accepted
both starts and passed all certificates; all three load-flow certification
horizons were shorter than their weighted horizons. A
one-second 10x10, three-escort, seed-23 smoke run retained the greedy incumbent
and correctly skipped certification because its weighted gap was too large.
The shell's 16-batch dry run and syntax checks passed. These checks are not the
new numerical campaign and must not be reported as its results.

### Additional safe-weight experiment (October 8, 2026)

Tal subsequently requested a separate version with sufficient integer objective
coefficients and conditional certification. `RunTable2SafeWeighted.sh` and
`RunSafeWeightedStatic.py` implement this version. The fixed `alpha=0.01`
experiment remains available, and this request does not authorize manuscript or
response-letter changes.

The current protocol is `safe_integer_flow_timing_v5`. For each instance, let
`D=sum(d_i)` be the sum of the targets' nearest-output Manhattan distances,
`K=N-e`, `H_g=F_g-D+max(d_i)`, and `U_g=K*H_g`. It uses
`R=U_g-D+1` in `R*F+M`, with the same coefficient for both formulations.
The weighted model covers this sufficient physical horizon, not merely the
greedy makespan. A minimum-flow, minimum-movement plan can stop all movements
after its final retrieval and has `M_star <= U_g`. Every feasible plan has
`M >= D`. A worse-flow plan therefore costs at least
`R*(F_star+1)+D > R*F_star+U_g`, which proves the coefficient is sufficient.
The feasible representative also implies `U_g >= D`, hence `R >= 1`; the
zero-flow case has `D=H_g=U_g=0` and `R=1`. Backend `weight_scale` defaults to
100, preserving the existing fixed-weight behavior.

Archived v4 used the larger coefficient `R=U_g+1`, which remains valid. V5
changes the sufficient coefficient, not the certificate or timing semantics.
New v5 runs must use a fresh results directory. Preserve archived v4 outputs;
do not reuse their objective-specific bounds, gaps, or proof times under the
new coefficient, and do not mix protocol versions in a campaign.

Earlier, Tal requested preserving the same branch-and-bound tree for the extra
time, rather than rebuilding a pure-flow model. The original defaults were 300 seconds for
the initial reporting cutoff, up to 300 seconds of additional search, and 16
threads. There is one continuous optimize call with a 600-second total cap,
unchanged integer objective and MIPFocus 0, and callback-controlled stopping.
Before the initial cutoff it attempts both objectives (`MIPGap=0`,
`MIPGapAbs=0.999`). A checked gap below one proves global lexicographic
optimality. Otherwise it tries a distance bound and an incumbent-specific proof:
`L_Q > R*(F_w-1)+(N-e)*H_minus`, with a 0.001 boundary margin, where
`H_minus=F_w-1-D+max(d_i)`. The weighted horizon must cover `H_minus`.
For `gap=R*F_w+M_w-L_Q`, the unchanged sufficient criterion is
`gap < R+M_w-U_minus`, where `U_minus=K*H_minus`; in v5 its right-hand side is
`1+M_w-D+K*(F_g-F_w+1)`. The full lexicographic integer-gap criterion remains
strictly below one with the existing numerical safeguards.
If this certifies flow time at the initial cutoff, it stops and movements can
remain unproved. Otherwise the same search continues until the frozen flow
time is certified, a lower-flow counterexample is found, or the total budget
expires. The previous gap gate does not apply to this version.

Tal explicitly requested reporting BOTH the initial-cutoff solution and final
solution separately, with the initial solution as the main experimental
result. Later improvements must never overwrite the main candidate or its
initial bound/gap. The extension certifies the original candidate's flow time.
The callback freezes the best incumbent observed at or before the cutoff,
before processing a later solution. Missing incumbents/bounds are recorded
explicitly, and callback observation times distinguish checkpoint delays from
the nominal cutoff. No candidate-specific extension is attempted when the
initial cutoff has no incumbent.

The CSV records coefficients, movement/horizon bounds, integer and normalized
gaps, proof source, extension timing, separate initial/final solutions, and
distinct flow/lexicographic proof flags. V5 adds `safe_movement_lower_bound=D`;
`safe_movement_bound` remains the upper bound `U_g`, not `U_g-D`. Results go in a
fresh `results_table2_safe_weighted_*` directory. Method notes contain the proofs.

Validation of v5: 154 unit/integration tests passed across the safe-weight, callback timing, runner, fixed-weight certification, warm-start, and lexicographic suites. Six tiny paired EF/LF solves (3x2, two targets, two escorts, seeds 1, 3, and 4) matched in coefficient, flow time, and movement count; all certified both objectives with one optimization call. Source/destination CSV validation, protocol separation, shell syntax, and whitespace checks passed. These are smoke checks, not the full numerical campaign. Evidence: `revision_R1/v5_tight_weight_smoke_rgn8ddbj/summary.json`.

Historical validation of the initial safe-weight version: 115 unit/integration
tests passed, including the
existing lexicographic/fixed-weight/warm-start suites, the new certificate
math, frozen-cutoff callback cases, and runner reporting. A live Gurobi test
asserted one optimize call, unchanged objective/focus, and the combined budget.
Six tiny EF/LF solves matched in coefficient and solution and proved both
objectives without extension. Short-budget 10x10, three-escort, seed-23 runs in
both formulations preserved the cutoff candidate `F=25,M=113` and reported the
later `F=25,M=112` separately; neither falsely certified flow time. These smoke
runs use shortened limits and are not numerical-campaign results. Shell syntax,
the 16-batch dry run, CSV merging, and diff whitespace checks passed.

On October 8, Tal revised this protocol to measure the first flow-time proof.
He briefly selected a 300-second overall limit, then corrected it: the first
phase has 300 seconds, followed by up to 300 extra seconds only if the saved
candidate's flow remains unproved. `RunTable2SafeWeighted.sh` and
`RunSafeWeightedStatic.py` then used protocol `safe_integer_flow_timing_v4`, with
a 300-second extension allowance and `stop_on_flow_proof=True` in
both backends. Tal then requested replacing time/node sampling with checks
only when the lower bound improves or a new candidate has smaller flow.
Unchanged/weaker bounds and movement-only improvements skip proof comparisons.
Mandatory cutoff and final checks always remain active. The
flow-specific threshold is cached and recalculated only when flow changes;
distance bounds and load count are computed once. Monitoring uses existing
candidate data and adds no solution-vector reads when checking bound changes.

New `first_flow_proof_*` fields record observed solver runtime, elapsed time
including model construction, certified flow, node count, bound, work,
observation source, and mathematical method. Unproved or invalidated proof
times are blank. `final_runtime` remains the total solve time. Callback sampling
can delay observation, so the time is an observed upper bound rather than an
exact latent proof instant. There is no additional time/node sampling delay.
The strongest valid observed bound is reused for a new candidate. Once a proof
is recorded, routine proof checks stop while contradiction detection and final
validation remain active. CSV check mode is `bound_or_flow_change`.
The cutoff candidate and final candidate remain
separate. Flow proof alone does not stop the first phase early, but the extension
stops when that candidate is certified, disproved by a lower-flow counterexample,
or the combined time limit expires. Existing fixed-weight entry points retain their separate certification
behavior, and the backends keep the old conditional stop as a compatible default
for callers that do not request the new tracking mode.

Historical validation of the v4 timing revision: 149 unit/integration tests
passed, including
event filtering, cached-threshold and no-extra-solution-read assertions,
bound validation against the smallest known feasible weighted objective, reliable proof timing,
continued first-phase search after flow proof, exact conditional-stop observation
times, and blank times for unproved instances.
Six live tiny EF/LF solves matched flow, movements, and coefficient, produced
valid first-proof timestamps no later than total runtime, and proved both
objectives. CSV merging, shell syntax, all 16 dry-run commands, and whitespace
checks passed. Two short-budget live 10x10, three-escort, seed-23 runs in EF/LF
verified conditional extension and time-limit reporting while preserving the
initial candidate `F=25,M=113`; LF's later `F=25,M=112` remained separate. No
flow proof was falsely reported. These are validation runs, not new paper
experiment results.

### Earlier two-phase implementation

- New entry points: `EscortFlowStaticLex.py` and `LoadFlowStaticLex.py`. Existing
  static entry points also accept `--lexicographic`; weighted mode remains the default.
- Shared implementation: `static_lexicographic.py`, two explicit Gurobi solves.
  First minimize integer flow time, then fix its best incumbent by equality and
  minimize integer load movements. Both phases set `MIPGap=0`, `MIPGapAbs=0.999`,
  with a separate strict integer-gap certificate below 1.
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
- Numerical stopping correction: the first Linux 13x7 batch contained three rows
  with phase-two `OPTIMAL` but integer gap exactly 1 (escorts/seed: 3/94, 5/57,
  6/14). Changed the solver tolerance from 1 to 0.999, retaining the strict <1
  certificate, and report `MOVEMENTS_NOT_PROVEN` if raw `OPTIMAL` lacks that
  certificate. Local reruns proved the same flow time 11 and respective movement
  counts 61, 57, 53. Original Linux runtimes are not replaced by local rerun times.
- A broader workbook audit also found false-positive proof flags when a raw bound
  lay microscopically above an integer. Certification now requires gap < 1 - 1e-6;
  raw bounds/gaps remain available, and the solver tolerance stays 0.999. In the
  original files, 13x7 has 4 phase-one and 5 phase-two such flags; the partial 10x10
  file has 20 phase-one and 3 phase-two flags (some phases share an instance).
  Six phase-one rows also have worse flow time than a known weighted solution.
  Corrected local 13x7 / 3 escorts / seed 42 rerun gives flow 10, movements 62,
  whereas its original Linux row reported flow 11 as proven. Original result
  workbooks are preserved and require reruns before use as certified optima.

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

### R3.M4 Follow-Up After Expanded Experiments (October 8, 2026)

Tal explicitly requested returning to Referee 3 minor comment M4 after the new
experiments are complete, including the different numbers of escorts and the
continue retrieval mode. Keep this item open until that expanded campaign has
been completed and checked; code implementation or smoke tests are not enough.

Then recompute the greedy-relative improvements on common matched instance
sets, examine the findings by escort count and retrieval mode, and update both
the manuscript discussion and the R3.M4 response. Distinguish incumbent quality
from feasibility, LP bounds, optimality certification, and runtime. State the
population and aggregation method explicitly, and do not assume the current
leave-mode findings transfer to continue mode.

The October 8 audit corrected the historical comparison to 21.9% for both
formulations on the same 32 rows (3,200 instances). This is an interim result,
not the answer for the expanded campaign. Evidence and the reproducible audit
are in `revision_R1/r3_revision_2026-10-08/audit_m3_m4.py` and `.json`.
Recheck the recorded historical objective-component discrepancies and revised
initial-state, objective, and horizon conventions when replacing those results.
The response letter is `revision_R1/response_letter_R1.tex`; the canonical main
paper remains in the linked Overleaf source directory.

### R3.3 Computational Movement Scope (October 8, 2026)

Tal chose to reinforce the paper's existing focus on simultaneous block
movement rather than add SLM experiments. The opening of Section 4 in the
canonical main manuscript now states that all reported computational results
and comparative performance conclusions concern SBM. The SLM formulation in
the supplement is presented as an extension of the modeling framework.
Response R3.3 reports this clarification, and its pending list no longer asks
whether to add SLM evidence. The continue-mode experiments and remaining
computational audit and reporting work are still open.

### Targeted 70%-Occupancy Launcher (October 8, 2026)

Tal authorized preparing the limited high-escort experiment for his 96 GB Mac
Studio. `Run70Percent.py` now launches leave-mode SBM runs with four targets,
seeds 1-100, and both formulations on 16x10 with 48 escorts and outputs (4,0),
(11,0), and 27x10 with 81 escorts and outputs (4,0), (13,0), (22,0). This is
exactly 70% occupancy, 200 distinct instances and 400 solves. Default budgets
are 300 solver seconds for the main result and a conditional additional 300
seconds in the same search only when its flow time remains unproved.

The launcher detects Mac performance cores by default, supports a thread and
Python override, checks the Gurobi license, and records the actual machine,
memory, versions, source hashes and protocol. It freezes the Python sources
for the campaign, preserves partial per-layout CSVs and logs, merges results
by formulation, and validates paired initial states and common settings.
Use `python3 -u Run70Percent.py` from Code on the Mac Studio. The README has
preview and short-trial commands and documents all result files. Interpret
cross-machine runtime comparisons separately from the paired comparison on
the Mac Studio.

The multi-target smoke check exposed an index mismatch in the safe runner's
horizon policy. Both formulations now share physical horizon
`H=max(H_g,C_g+1)`; EF uses `T=H-1` and LF uses `T=H`. This preserves the
complete common greedy start and the sufficient-horizon guarantee. Existing
archived v4 results are unchanged. The new campaign uses the current v5
coefficient `R=U_g-D+1`, with the actual protocol and horizon saved per row.

Validation passed 101 safe-runner/certificate/search tests, 16 warm-start
tests, CLI preview/override/overwrite checks, and four licensed Gurobi smoke
solves (one seed per layout and formulation, two threads, five-second budget).
All four accepted the starts and proved both objectives; paired instances,
coefficients, physical horizons and objective values matched. Smoke evidence
is in `revision_R1/occupancy70_launcher_smoke_matched_2026-10-08/`.
The initial one-second smoke directory is retained and records the mismatch
that was then corrected. The full 400-solve campaign has not been launched,
and task 4 in `tasks.md`, R2.5 reporting and the expanded R3.M4 analysis remain
open until those results are collected and audited.

### Expanded static campaigns and continue mode (October 9, 2026)

Tal requested rerunning approximately 70% occupancy on the same Linux box as
the main experiment and adding four rows to Table 2(b). `Run70Percent.py`
now covers 13x7 with 27 escorts, 10x10 with 30, 16x10 with 48, and 27x10 with
81, preserving the paper's output locations. The first occupancy is 64/91
(70.33%); the remaining three are exactly 70%. With four targets and seeds
1-100, this is 400 distinct instances and 800 solves. The earlier completed
two-layout Mac Gurobi 13.0.1 campaign is audited separately in
`revision_R1/occupancy70_analysis_2026-10-09/`; it does not replace the new Linux
campaign. Earlier launch-status statements above are historical.

New entry points share `RunStaticCampaign.py`:

- `RunContinue.py`: continue mode, 2/4/6 targets, four layouts, and 8/12/16
  escorts plus the approximately 70% category. Default: 4,800 instances and
  9,600 solves.
- `RunTable2Targets.py`: leave mode, 2/6 targets, the same layouts and escort
  categories. Default: 3,200 instances and 6,400 solves. This supersedes the
  earlier proposed 2/3-target plan in `tasks.md`.

All use the safe v5 objective, paired greedy starts, sufficient common physical
horizons, 300-second main results, and conditional 300-second extensions.
Linux defaults to at most 16 threads; explicitly pass `--threads 16` on the
experiment machine. Mac defaults to performance cores. `--layouts` and
`--seeds` allow a pilot; `--dry-run` creates no files or solver processes.
Every campaign freezes its Python sources, logs actual Gurobi/native versions
and machine settings, preserves partial files and validates paired settings.
The README lists merged CSV names and launch commands.

Continue-mode readiness required substantive corrections. `SolveGreedy` now
serves initial output targets at time zero and removes targets from the active
set upon arrival, so it does not add one extra flow-time unit per target or an
extra terminal iteration. Escort count is constant. It raises a catchable
error on a failed step cap and accepts success at the exact cap. The safe
runner uses a conservative step cap independent of high escort counts.
EF treats initially served output targets as ordinary blockers. LF now has an
explicit retrieval-mode config: continue retrieval transfers `q` units from
target flow to blocking flow, including time-zero supply, with no leave-mode
service delay. Both complete-start encoders implement these same semantics.
The legacy `LoadFlowStatic.py` CLI also passes its mode correctly and accepts
continue BM warm starts. Unsupported Gurobi stay requests fail explicitly.
CSV merging rejects mixing leave and continue modes. The coefficient and
certificate proofs remain valid because `K=N-e` is an upper bound on the
number of loads moving per period in either mode.

Validation evidence is in `revision_R1/continue_validation_2026-10-09/`.
All 4,800 planned continue greedy traces (48 configurations, seeds 1-100)
were physically replayed and checked for disjoint legal block shifts, constant
escort count, complete retrieval, exact flow time and movement count. Tiny
continue cases with 2/4/6 targets accepted complete fixed starts in both
models; optima agreed between EF/LF and with independent physical-state
Dijkstra searches where enumerated. Initial output targets, all targets at
outputs, repeated output use and blocker conversion are covered. These are
validation checks, not completed paper experiments.

The final static suite passed 173 unit/integration tests with native Gurobi
13.0.3. Three launcher smoke campaigns completed: eight solves across all four
70% layouts, 24 continue solves covering every target/escort category on 13x7,
and 16 additional-target leave solves on 13x7. All 48 accepted complete starts
and retained feasible incumbents; pairing and proof timestamps passed checks.
Short budgets (15 seconds for the 70% smoke, two seconds for the other two,
no extensions, one thread) deliberately allow incomplete objective proofs.
`validation_summary.json` records these distinctions and current source hashes.
