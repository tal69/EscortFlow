# Weighted retrieval and flow-time certification

Method notes for the revision, updated October 9, 2026. These notes describe the mathematical claims and implemented experiment protocol. They do not assert that the new experiments have already established those claims.

## Current version: sufficient integer weights and flow-proof timing (v5)

`RunSafeWeightedStatic.py` and `RunTable2SafeWeighted.sh` implement a separate
protocol from the fixed `F+0.01*M` experiment described below. They use an
instance-specific integer objective `Q=R*F+M` that guarantees lexicographic
optimality when solved to an integer absolute gap below one. Both formulations
receive the same greedy physical start and use the same coefficient per instance.

Let `F_g` be the greedy feasible flow time, `N` the number of cells, `e` the
initial number of escorts, and `d_i` target i's nearest-output Manhattan distance.
Set

```text
D   = sum_i(d_i)
K   = N-e
H_g = F_g - D + max_i(d_i)
U_g = K*H_g
R   = U_g - D + 1.
```

Every feasible retrieval plan has `M >= D`: each target requires at least its
nearest-output Manhattan distance in individual load movements, and blocking-load
movements can only increase the total.

The weighted horizon covers `H_g`, as well as the complete greedy trace. Both
formulations use the same physical horizon `H=max(H_g,C_g+1)`, where `C_g`
is the greedy makespan. Escort-flow uses movement index `T=H-1`, with arrivals
at `t+1`; load-flow uses retrieval index `T=H`. The extra greedy tail preserves
output service in the complete warm start. A horizon supplied only by the greedy
makespan would not be sufficient for the following global guarantee.

To prove the guarantee, choose a minimum-flow plan with the fewest movements.
Its flow is at most `F_g`, so its final retrieval occurs by `H_g`. Remove all
movements after its final retrieval. At most `N-e` loads can move in a period,
so this plan has at most `U_g` movements and is represented by the weighted
model. Write its objectives as `(F_star,M_star)`, with `D <= M_star <= U_g`.
Any plan with flow time at least one greater has weighted objective at least
`R*(F_star+1)+D`, whereas the selected plan's objective is at most
`R*F_star+U_g`. Their separation is at least

```text
R + D - U_g = 1.
```

Thus `R > U_g-D` suffices, and the integer choice `R=U_g-D+1` guarantees
lexicographic optimality. Among minimum-flow plans, minimizing the weighted
objective minimizes movements. This argument needs an upper bound on an optimal
representative and a lower bound on every feasible plan, not an upper bound on
every plan with redundant post-retrieval movements. The same feasible
representative establishes `U_g >= D`, so `R >= 1`. Invalid inputs that violate
this consistency must not be treated as a valid coefficient certificate.
The zero-flow case gives `D=H_g=U_g=0` and `R=1`.

### Retrieval modes and expanded campaigns

The same sufficient coefficient, horizon and lower-bound certificate apply to
leave and continue retrieval. In leave mode the number of stored loads decreases;
in continue mode it stays at `K=N-e`. In both modes, at most `K` loads move in
one period and all movement can stop after the final target retrieval. The
nearest-output distance bound and the last-retrieval bound are unchanged.

`RunSafeWeightedStatic.py --retrieval-mode continue` uses the common continue
greedy trace. A target is served immediately upon arrival and becomes an
ordinary blocking load, preserving escort count. Initially output-located
targets are already served at time zero. The greedy flow time sums arrival
times without an extra post-arrival iteration. EF represents served targets
implicitly as blockers; LF explicitly transfers their flow from commodity 1
to commodity 2 through `q`. Unlike leave mode, continue mode has no output
service takt blocking movement of the retrieved load. Complete starts preserve
these conventions, including initial output targets and repeated use of an output.

`Run70Percent.py` now prepares four additional four-target Table 2(b) rows on
all four layouts, using escort counts 27, 30, 48 and 81. Occupancy is 64/91 on
13x7 and exactly 70% on the other layouts. `RunContinue.py` combines 2, 4 and
6 targets with escort counts 8, 12, 16 and the approximately 70% category on
each layout. `RunTable2Targets.py` prepares the same categories in leave mode
with 2 and 6 targets. All use the shared `RunStaticCampaign.py`, frozen source
copies, common coefficients and physical horizons, and the existing v5
main/extension reporting rules. These campaigns are intended for the same Linux
machine, native solver version and thread count as the main numerical study.

Mode labels and source hashes distinguish runs; existing leave-mode files and
the earlier Mac 70% campaigns must remain separate. No stored result is rewritten.
The launchers reject mixed retrieval modes when merging CSVs. Their native
Gurobi preflight records the actual linked solver version, rather than relying
only on the Python package version.

The weighted search uses one uninterrupted optimization call with a 300-second
first-phase limit and a conditional extension of up to 300 seconds. It records
the first observed flow-time proof while the first phase continues toward
proving both objectives. Tal corrected the briefly selected 300-second overall
limit to restore this conditional extension.
The objective remains `R*F+M`, `MIPFocus=0` remains unchanged, and the same model,
incumbent, cuts, and branch-and-bound tree remain active. `TimeLimit` is the sum
of the allowances, and the callback implements conditional stopping. There
is no model rebuild, second optimization call, or change to a pure-flow
objective. The initial complete greedy warm start is applied once.

The solver retains the same weighted objective, `MIPGap=0`, and
`MIPGapAbs=0.999`. By default, flow optimality does not stop the first phase
early, but ends the extension when it certifies the saved first-phase candidate.
The optional `--stop-at-flow-proof` setting requests termination at the first
observed flow proof in either phase. It uses the same cached, event-driven
criterion and does not change the objective, coefficient, horizon, or proof
requirements. The completed solve still reports its movement count, weighted
gap, all proof flags, and all existing timing KPIs. Movement optimality is
claimed only if its separate certificate also succeeds. An early result fills
both main and final result columns when the initial cutoff has not been reached.
Callback-requested termination is labeled `FLOW_PROVEN_EARLY`; final validation
and Gurobi termination overhead can make total runtime slightly exceed the
recorded first-proof time.

The option defaults to false in both solver configurations and all campaign
launchers. Direct Python callers set `stop_at_flow_proof=True` in either
configuration with `objective_mode="weighted_integer"` and a positive finite
`time_limit`; omitting an extension limit gives a single-budget search.
Each integer CSV records `stop_at_flow_proof` as 0 or 1; campaign
metadata records the Boolean setting. Merge checks reject mixed settings within
or across batches, and pairing requires the same setting in both formulations.
Existing results without this schema field remain separate. Direct continuous
LP runs do not accept this integer flow-proof stopping option; integrated LP
bounds can still be computed after an integer solve stops early.
The reporting and stopping rules are:

1. A reliable, consistent weighted lower bound with reconstructed integer gap
   below `1-1e-6` proves global lexicographic optimality under the coefficient
   and horizon guarantee above. The solver can finish before the cutoff in
   this case.
2. At available MIP/MIPSOL callbacks, compare the proof criterion only when
   the weighted lower bound improves or a new candidate has smaller flow,
   including candidates that do not improve the weighted incumbent. Reuse the
   strongest valid observed bound when a later candidate arrives with a weaker
   bound. Record the first certified flow
   value and its observed runtime, node count, bound, work, observation source,
   and certificate method. The final result receives an unconditional check.
3. At the cutoff, freeze the solution and bound observed by that time. If its
   flow time is already proved, stop at the next supported checkpoint. Otherwise
   continue the same search until that flow is certified, a lower-flow feasible
   counterexample disproves it, or the combined 600-second cap is reached.
   Report the initial and final solutions separately. Missing cutoff candidates
   are recorded without filling them using later solutions; no candidate-specific
   extension is run in that case.

The callback observes each new feasible solution and keeps only improving
weighted incumbents for the initial snapshot. A first callback after the
cutoff must freeze the earlier cached candidate before processing any new
post-cutoff candidate. Later movement or flow improvements never overwrite the
initial experiment result. Any late feasible solution with smaller flow is a
counterexample, even if its weighted objective does not improve the current
incumbent. The final weighted incumbent is also reported in separate columns.
Raw callback termination status is distinguished from the proof outcome.

For the gap check, let `(F_w,M_w)` be the frozen initial incumbent and `L_Q` a valid
weighted lower bound. Put `B=F_w-1`. If `B < D`, the distance bound
already proves minimum flow. Otherwise a hypothetical minimum-flow plan with
flow at most `B` can be normalized to have

```text
H_minus = B - D + max_i(d_i)
U_minus = K*H_minus.
```

Provided the weighted model covers that physical horizon, such a plan would
have objective at most `R*B+U_minus`. Therefore

```text
L_Q > R*B + U_minus
```

proves minimum flow. Equivalently, with `gap=R*F_w+M_w-L_Q`, the sufficient
threshold remains `gap < R+M_w-U_minus`. Substituting the v5 coefficient gives

```text
R+M_w-U_minus = 1 + M_w - D + K*(F_g-F_w+1).
```

Only this substitution changes; the lower-bound proof and the independent
full-integer-gap criterion below one are unchanged. The implementation leaves a
`0.001` margin from the flow-certificate boundary, plus the common numerical
comparison margin.
Missing, nonfinite, inconsistent, or unreliable solver bounds cannot establish
this gap certificate. The threshold is derived for the frozen candidate and
continues to refer to that same candidate as the bound improves. Any separate
certificate for the final incumbent is recomputed using that final incumbent.
The initial a posteriori check requires no extra search when it succeeds. It
proves flow time alone; the CSV must not mark movement or full lexicographic
optimality unless separately established.

Target-to-output distance bounds and the initial load count are computed once.
The scalar bound threshold for a candidate flow value is cached and recalculated
only when that value changes. Movement improvements at the same flow value do
not change the threshold. Bound-change monitoring reuses the saved candidate and
compares an improved weighted bound to this cached criterion, without reading
another solution vector. Unchanged or weaker bounds and movement-only
improvements require no proof comparison. Mandatory cutoff and final checks
remain active, even when bounds do not change. Once a certificate is recorded,
routine checking ends while contradiction detection and final validation remain.

`first_flow_proof_runtime` is the observed solver time of the first proof.
`first_flow_proof_cpu_time` adds model construction; both are elapsed time,
not CPU time accumulated across threads. `first_flow_proof_flowtime` identifies
the flow value proved, which must match a reported candidate before attaching
that proof to the candidate. `first_flow_proof_source` distinguishes MIP,
MIPSOL, and final-result observations; `first_flow_proof_method` distinguishes
distance, weighted-gap, and full-weighted-optimum certificates. The timestamp
is an upper bound on the instant proof became possible. There is no artificial
time/node sampling delay, but a long operation without a supported callback can
delay observation. If proof is first established by the final-result check,
the timestamp is the final runtime. Unproved cases have blank timing fields.
Unreliable final solver statuses or a later contradiction invalidate recorded
proof timing. Protocol `safe_integer_flow_timing_v5`, with recorded check mode
`bound_or_flow_change`, separates these results from the earlier untimed and
interval-based campaigns. Gurobi still invokes callbacks at its own checkpoints;
the event filtering takes place inside the callback.

CSV rows identify the coefficient, horizon, movement bounds, integer and
normalized objectives/bounds/gaps, gap-certificate threshold and outcome,
proof source, time before/after the reporting cutoff, and any counterexample.
The new `safe_movement_lower_bound` field records `D`. The existing
`safe_movement_bound` field continues to record the upper bound `U_g`, not the
range `U_g-D`; therefore `R=safe_movement_bound-safe_movement_lower_bound+1`.
The main solution, bound, and gap columns belong to the initial reporting
cutoff; final columns belong to the end of the continuous search. A certificate
using a later bound must not be represented as a bound obtained within the
initial 300 seconds. Snapshot observation times make callback timing explicit.
Callbacks execute at solver checkpoints, so the phase transition and stopping
request can lag the nominal cutoff. Gurobi's termination overhead can also
exceed the time cap slightly. Missing initial incumbents or bounds are recorded
explicitly rather than filled with later observations. Existing fixed-weight
entry points and their separate pure-flow certification remain unchanged.
The implementation uses Gurobi's documented [MIP and MIPSOL callback data](https://docs.gurobi.com/projects/optimizer/en/current/reference/numericcodes/callbacks.html)
and [callback termination mechanism](https://docs.gurobi.com/projects/optimizer/en/current/reference/python/model.html#Model.terminate).

### Protocol history and result separation

Archived protocol `safe_integer_flow_timing_v4` used `R=U_g+1`. That larger
coefficient remains sufficient, and changing to v5 does not invalidate earlier
certificates established under their recorded coefficient and horizon. V5 uses
the movement lower bound to reduce the coefficient and records the new CSV field.
Timing, event filtering, frozen-cutoff reporting, and the conditional extension
retain their existing rules.

Start any v5 campaign in a fresh results directory. Preserve v4 outputs and their
protocol identifiers; do not append v5 rows to a v4 campaign or reinterpret old
weighted bounds, gaps, or proof times under the new coefficient. A solver bound
belongs to its objective, so an old bound cannot be reused as a bound on the new
`R*F+M`. Recomputing an incumbent's score from its `F` and `M` does not convert its
old solver bound or timing into a v5 result.

Validation of v5: 154 unit/integration tests passed across the safe-weight, callback timing, runner, fixed-weight certification, warm-start, and lexicographic suites. Six tiny paired EF/LF solves (3x2, two targets, two escorts, seeds 1, 3, and 4) matched in coefficient, flow time, and movement count; all certified both objectives with one optimization call. Source/destination CSV validation, protocol separation, shell syntax, and whitespace checks passed. These are smoke checks, not the full numerical campaign. Evidence: `revision_R1/v5_tight_weight_smoke_rgn8ddbj/summary.json`.

## Earlier fixed-weight objective and claims

The following sections document the separate fixed-weight workflow and its
certification rationale. They do not redefine the current v5 safe-weight protocol.

Let `F` denote total integer flow time and `M` the integer number of individual one-cell load movements, including movements of blocking loads. A block movement contributes its length to `M`. The earlier fixed-weight protocol uses

```text
minimize Z = F + alpha * M, with alpha = 0.01.
```

In the current command-line programs, the movement coefficient is named `gamma`; the paper's proposed `alpha` therefore corresponds to `--gamma 0.01`, with `--beta 1`. Load-flow's existing `--alpha` parameter is a makespan coefficient and must remain zero for this experiment.

The formulation defines a weighted optimization problem. An empirical finding that its solutions also minimize flow time is a separate, instance-specific result. A weighted optimum with minimum flow time also minimizes movements among all minimum-flow solutions: a solution with the same flow time and fewer movements would have a smaller weighted objective. Consequently it is lexicographically optimal on that instance. A merely feasible or loosely terminated weighted solution requires an additional movement certificate before making the full lexicographic claim.

## Why a sufficiently small weight exists

Suppose the feasible set under consideration satisfies `0 <= M <= U`, where `U` is finite. For two feasible plans with `F_2 >= F_1 + 1`,

```text
(F_2 + alpha*M_2) - (F_1 + alpha*M_1) >= 1 - alpha*U.
```

Thus every `0 < alpha < 1/U` makes all strict flow-time improvements preferable, regardless of movements. Among equal-flow plans, any positive weight prefers fewer movements. For example, `alpha = 1/(U+1)` suffices. When `U = 0`, every positive weight suffices. More generally, a valid bound on the movement range, `M_max - M_min`, can replace `U`.

For a rectangular PBS with `N = Lx*Ly` cells and `H` modeled movement periods, a safe bound for the present BM and LM formulations is

```text
U = N*H.
```

Each individual load moves at most one cell per period. In escort-flow, the nonoverlapping escort paths give the same bound because each selected path of length `d` occupies `d+1` cells and contributes `d` movements. The current implementations create movement variables for `t = 0,...,T`, so the conservative implementation-level bound is `U = N*(T+1)`. It is safe even if some final-period movement variables are redundant. No claim is made that this bound is tight or that `0.01` satisfies it for the tested instances.

If the problem has no prescribed finite horizon, the argument needs a finite horizon shown to contain a lexicographic optimum, or a valid bound on the movements of at least one such optimum. A bound derived only from an arbitrary restricted model does not establish unrestricted equivalence.

## Numerical implications

With a fixed flow time, a one-movement improvement changes the original weighted objective by only `alpha`. A guaranteed movement-optimality certificate therefore needs an absolute objective gap strictly below `alpha`, subject to numerical safeguards. A relative MIP-gap setting alone may be too loose, depending on the objective magnitude.

Very small weights can create a large objective coefficient range and make secondary improvements difficult to distinguish under finite precision and stopping tolerances. These are potential numerical and computational disadvantages, not a theorem that the particular models become unstable. Gurobi's [numerical guidance on hierarchical objectives](https://docs.gurobi.com/projects/optimizer/en/current/concepts/numericguide/tolerances_scaling.html#improving-ranges-for-variables-and-constraints) discusses the risks of aggregating objectives with widely separated weights.

For the experimental choice `alpha = 0.01`, use the exactly equivalent objective

```text
minimize Q = 100*F + M.
```

`Q` is integer valued. This scaling preserves the weighted ordering exactly and permits integer absolute-gap stopping. It does not turn `alpha = 0.01` into a universally lexicographic weight. The coefficients 100 and 1 are moderate; this transformation avoids needing a `0.00999` stopping threshold in the original units. Scaling can still change the solver's numerical search path, so the revision should state the implementation used.

Use `MIPGap = 0` and `MIPGapAbs = 0.999` for this integer objective. Independently reconstruct the incumbent as `Q_w = 100*F_w + M_w` from validated integer values and check that its gap to a valid global lower bound is less than `1 - epsilon`, using the project's conservative certificate margin. Do not derive a proof solely from Gurobi's `OPTIMAL` status or from rounding the bound upward. A bound materially above the reconstructed incumbent should trigger a consistency failure. Report both `Q` and the original `Z = F + 0.01*M`, including appropriately scaled bounds and gaps.

If `alpha` were much smaller, multiplying by `1/alpha` could avoid a tiny absolute-gap number but would introduce a large flow-time coefficient. Scaling changes the numerical representation, not the underlying separation of priorities. The distinction between MIP gap and LP dual-feasibility tolerance should remain explicit; they are different parameters.

## A certificate that the weighted candidate minimizes flow time

Let `(F_w, M_w)` be the saved feasible weighted candidate, and let `L_F` be a valid global lower bound for a separate minimization of `F` over the relevant feasible set. Since `F` is integer,

```text
L_F > F_w - 1  implies  F_w is minimum.
```

The proposed guarded target, `L_F > F_w - 0.999`, is sufficient. It leaves 0.001 between the stopping target and the mathematically critical boundary. This is a bound certificate; finding a new flow-optimal incumbent is unnecessary because the weighted candidate already attains `F_w` and is feasible. The implementation uses a further `1e-6` comparison margin and a bound-stop target just above that threshold. Each weighted solve receives a complete greedy MIP start. The certification solve receives the saved weighted incumbent itself, including all variable values and explicit zeros, truncated or padded with idle periods to fit the certification horizon without changing retrieval times. The saved weighted candidate also supplies the external feasible flow-time upper bound.

For a minimization model, Gurobi's [`BestBdStop`](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#parameterBestBdStop) can terminate once the global bound reaches `F_w - 0.999`. [`BestObjStop`](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#parameterBestObjStop) can alternatively stop upon finding a strictly better integer flow time, using a target near `F_w - 1` with an explicit numerical margin. Both may return `USER_OBJ_LIMIT`; the implementation must inspect the bound and incumbent to determine what was established. The final certificate must be recomputed from saved candidate values and the returned bound, independent of the termination status.

A pure-flow LP relaxation or a valid analytical lower bound may certify some candidates without a MIP search. Bound-focused MIP settings can be evaluated, but faster verification is an empirical possibility, not a guarantee.

### Preserve the right feasible set

Do not fix `F` to `F_w` for the flow-time certification solve. Do not constrain movements to `M <= M_w`. Either restriction could hide a better-flow plan. Remove any explicit weighted-objective cutoff constraint and reset objective-specific solver parameters before changing the objective. A weighted sublevel restriction may exclude a lower-flow plan with many more movements.

Certification on the same finite horizon would prove minimum flow time only within that horizon. The automatic greedy-makespan horizon supplies a feasible plan but, by itself, is not a proof that an unrestricted optimum fits.

The separate certification model uses a horizon covering every candidate with `F <= F_w`. If `d_i` is a valid individual lower bound on target `i`'s arrival time, then such a candidate satisfies

```text
C_max <= F_w - sum_i(d_i) + max_i(d_i).
```

This follows from `f_i <= F_w - sum_{j != i}(d_j)`. The implementation uses nearest-output Manhattan distances for `d_i`. For proving only the absence of a better flow time, `F_w - 1` can replace `F_w`. Write the displayed bound as `H`. Current escort-flow arrivals can occur at `t+1` for movement indices `0,...,T`, so `T >= max(0,H-1)` covers arrival times through `H`. Load-flow retrieval variables are indexed `0,...,T`, so it requires `T >= H`. The new runner sets certification `T = H` for both, even when this is shorter than `weighted_T`, while covering all no-worse-flow schedules. Every arrival of the weighted incumbent fits this horizon, so only its post-retrieval suffix can be discarded. The value of `H` is computed from the weighted incumbent's flow time, never from the greedy flow time. It records both horizons and their proof scope.

If movement optimality is also claimed for a larger feasible set, the weighted lower bound must be valid for that same set, or a separate fixed-flow movement certificate is needed. A weighted proof from the original restricted horizon cannot automatically be combined with a larger-horizon flow proof to establish unrestricted lexicographic optimality. In particular, a sufficient horizon check for a global movement claim at `F_w` is `weighted_T + 1 >= H` for escort-flow, or `weighted_T >= H` for load-flow. Otherwise distinguish global flow certification from lexicographic certification within the original weighted horizon.

When changing the horizon of a load-flow incumbent, the boundary period `min(weighted_T, H)` occurs after every target has been retrieved. Any redundant blocking-load movements in that period are replaced with stationary arcs before padding or truncation. This is necessary because the shorter model does not constrain their final destinations or opposite-direction conflicts. It preserves all retrieval times and can only reduce movements in the supplied start; the original weighted result remains unchanged. The CSV records `certification_warmstart_source=weighted_solution` and the number of removed post-retrieval movements. Escort-flow preserves the incumbent arcs through the new horizon. When extending it, the start adds the required output-service stays before converting newly retrieved targets into stationary escorts.

## Outcome rules and reporting

1. Solve the weighted model with a complete greedy MIP start and a 300-second solver-time cap. Record the actual time, work, objective, bound, status, candidate `(F_w,M_w)`, and warm-start outcome. Retain a validated greedy fallback separately if no solver incumbent is returned.
2. Preserve the weighted candidate and its results before deciding whether to verify it. Every proven weighted optimum is eligible. For a suboptimal weighted candidate, verify it only when its absolute weighted gap is below a configurable eligibility threshold that is significantly less than one, initially `0.1` in the original `F + 0.01*M` units. Divide the scaled objective gap by 100 before applying this gate. Missing, nonfinite, or inconsistent bounds do not pass the gate. Mark a candidate that fails it `FLOW_SKIPPED_WEIGHTED_GAP`, not unresolved or counterexample. Verification has a separate configurable budget, initially 300 seconds. Report verification time separately from the weighted benchmark time; also report their sum when describing the complete certification procedure.
3. If the pure-flow lower bound passes the guarded target, mark the original weighted candidate `FLOW_CERTIFIED`.
4. If verification finds an integer flow time smaller than `F_w`, mark the original candidate `FLOW_COUNTEREXAMPLE` and preserve both plans. The weighted candidate is then not flow-optimal. Do not silently replace it and retain the original claim about what the weighted run produced.
5. If neither event occurs before the verification budget, mark it `FLOW_UNRESOLVED`. A timeout is not evidence of a counterexample or evidence of optimality.
6. Mark `LEX_CERTIFIED` only when flow optimality and movement optimality at that flow are both justified on the same feasible set. A weighted global optimum plus the flow certificate suffices. More generally, a valid weighted gap below `alpha` suffices for movements at fixed `F_w`: any one-movement improvement at the same flow would reduce the weighted objective by at least `alpha`. With `alpha = 0.01`, a guarded gap below one for `Q = 100F+M` supplies this certificate and also certifies weighted optimality on the integer objective lattice.

Apply the same physical greedy plan to escort-flow and load-flow, with formulation-specific encodings of all model variables. Verify that both encodings have matching flow time and movement count and are feasible under the same retrieval convention. Report warm starts as a revision improvement, while making clear that the previous static experiment scripts left them disabled. For fair runtime comparisons, use the same machine, solver version, thread count, horizon policy, instance data, and benchmark budgets for both formulations.

Implement the revised workflow in a separate experiment runner. Keep the existing weighted command-line behavior available for reproducing the previous experiment; the new workflow explicitly selects the integer-scaled weighted objective and the subsequent pure-flow certificate model.

The paper can report the number of weighted candidates with certified minimum flow time, certified lexicographic optimality, counterexamples, and unresolved verification. It should not extrapolate an all-tested-instances finding to a general guarantee for `alpha = 0.01`.
