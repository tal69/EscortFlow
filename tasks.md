# Remaining R1 Revision Tasks

Updated October 9, 2026. This checklist records the remaining review work and
the decisions made during the revision. An implemented runner or a successful
smoke test does not complete an experimental task; the results must also be
validated and incorporated into the manuscript and response letter.

## Decisions to preserve

- Compare escort flow (EF) and load flow (LF) using the best feasible solution
  and valid lower bound available within the first **300 solver seconds**.
  Use the saved first-phase fields and certificates. Improvements or
  certificates obtained during an extension belong in separate reporting.
- Match initial states, output cells, movement and retrieval modes, physical
  horizons, coefficients, solver settings, and time budgets between methods.
  Use the common complete greedy warm start in both formulations.
- Preserve the coefficient and protocol actually used in every run. The
  archived v4 coefficient remains sufficient. The smaller v5 coefficient does
  not invalidate v4 results or by itself require restarting the campaign.
  Keep the protocols separate and do not reinterpret stored bounds under a
  different coefficient. A performance pilot comparing the coefficients is
  optional future work.
- Computational claims concern **simultaneous block movement (SBM)**. The
  supplementary SLM formulation is an extension of the modeling framework;
  adding SLM experiments is no longer a pending revision task.
- For the utilization issue, explain the crossover already visible in the
  existing results and add a limited approximately 70%-occupancy check on all four layouts.
  A full 50%, 60%, and 70% sweep is not part of the current plan.

## Open major tasks

### 1. Complete and incorporate the replacement numerical campaign

Review references: AE.3, R1.2, R2.3.

- [ ] Finish matched EF/LF runs for the selected single-target and multi-target
  configurations, using sufficient physical horizons and integer coefficients.
- [ ] Validate the complete run records and distinguish flow-only certification
  from certification of both lexicographic objectives at the 300-second cutoff.
- [ ] Replace the historical tables and corresponding discussion. State the
  coefficient rule, protocol, warm starts, solver version, hardware, and limits
  actually used. Keep any extension results separate.
- [ ] Update R1.2 with the final counts of configurations, initial states,
  distinct instances, and solver runs across both formulations.

### 2. Add multi-target continue-mode evidence

Review references: R2.1, R3.3.

- [x] Validate both implementations, including the treatment of retrieved loads,
  common warm starts, and physical time indexing in continue mode.
- [ ] Run and analyze matched multi-target continue-mode experiments.
  `RunContinue.py` prepares 2, 4, and 6 targets, all four layouts, and escort
  counts 8, 12, 16 plus approximately 70% occupancy. There are 4,800 distinct
  instances and 9,600 solver runs with seeds 1-100. Run on the numerical-study
  Linux box with 16 threads and preserve the initial 300-second comparisons.
- [ ] Incorporate the results and scope the abstract and contributions according
  to the operating modes actually tested.

### 3. Add intermediate target counts

Review reference: R2.4.

- [ ] Add two- and six-target cases under matched conditions (Tal's October 9
  instruction supersedes the earlier two-/three-target plan).
  `RunTable2Targets.py` prepares leave-mode Table 2(b) cases with the same four
  escort categories as the continue campaign: 3,200 instances and 6,400 solves.
- [ ] Report how solution quality, bounds, certification rates, and computation
  time change with target count, while separating effects of layout, utilization,
  retrieval mode, and horizon.

### 4. Explain the existing crossover and add targeted 70% cases

Review reference: R2.5.

- [ ] Explain how LF becomes relatively more competitive as escort count
  increases. Distinguish the crossover in flow-certification time from the
  crossover in time to prove both objectives, and avoid claiming a universal
  threshold across layouts.
- [ ] Add approximately 70%-occupancy cases on **13x7 with 27 escorts**,
  **10x10 with 30**, **16x10 with 48**, and **27x10 with 81**, preserving the
  existing outputs and matched protocol. The 13x7 occupancy is 64/91=70.33%;
  the others are exactly 70%. Four targets and seeds 1-100 give 400 distinct
  instances and 800 solves, adding four rows to Table 2(b).
  `Run70Percent.py` prepares the expanded campaign on the same Linux box. It saves the
  actual machine, thread count, solver version, coefficient protocol, source
  snapshot, and separate cutoff/final results. The complete experiment and its
  incorporation remain pending; launcher validation does not close this task.
  The earlier two-layout Mac Gurobi 13.0.1 run is complete and audited, with LF
  faster and all 400 solves optimal. Keep its hardware/version separate from
  the forthcoming Linux results.
- [ ] Update the manuscript and R2.5 reply around this narrower scope. The
  current response-letter draft still proposes a broader 50%/60%/70% sweep;
  replace that wording before submission.

### 5. Audit benchmark provenance and implementation consistency

Review references: R2.7, R3.1, R3.2, R3.3.

- [ ] Check Bukchin and Raviv (2023) and identify which parts of the multi-target
  SBM load-flow model are previously published and which are new or strengthened.
  State the provenance explicitly in the main text.
- [ ] Reconcile the supplement's statement that its load-flow model differs from
  the model used in the experiments. Trace results to the implemented constraints,
  including equality (A.13) and the output conventions.
- [ ] Audit both computational implementations against the clarified physical
  model: initial targets at outputs, first-step conflict indexing, output service,
  terminal layers, and physical horizon conventions.
- [ ] If relevant benchmark strengthening was omitted, assess its effect before
  retaining a claim about comparison with the strongest available formulation.

### 6. Rebuild the computational reporting

Review references: R2.6, R3.3, R3.M4.

- [ ] Organize Section 4 by escort count/utilization, target count, layout, and
  the EF/LF comparison. Explain the practical implications of the observations
  without inferring that tighter LP bounds alone cause shorter runtimes.
- [ ] Report medians, quantiles or interquartile ranges, timeout rates, and a
  performance profile or suitable distribution plot alongside averages.
  Distinguish time-limited runs from completed optimal solves.
- [ ] Define LP bounds, root-node bounds, MIP gaps, feasible-solution rates,
  first flow-proof times, and full-objective proof times precisely. Explain
  whether elapsed times include model construction; the CSV `cpu_time` fields
  represent elapsed wall-clock time rather than summed processor time.
- [ ] Audit per-instance selection and aggregation: common populations,
  denominators, missing incumbents, greedy fallback, LP-gap formulas, and the
  order in which best values and averages are calculated. A certificate from
  EF must not be presented as a certificate obtained independently by LF.
- [ ] Collect variables, constraints, and nonzero counts for both formulations
  under matched conditions.
- [ ] Resolve archived objective/component discrepancies and refresh the
  greedy-relative improvements using the expanded matched campaign, including
  continue mode. Update the discussion and R3.M4 reply; the corrected historical
  percentages are interim results.

### 7. Clarify the operational horizon and practical claims

Review references: R3.2, R2.8.

- [ ] Explain in the main text how the static mathematical horizon differs from
  a retrieval batch, service window, prediction horizon, and execution interval
  in rolling-horizon control.
- [ ] Align the abstract and conclusion with a static exact optimization
  contribution and a basis for future dynamic methods. Do not imply that the
  reported experiments demonstrate real-time readiness or rolling-horizon
  performance.

## Final editorial and submission pass

- [ ] Standardize cell-set notation in the main paper and supplement. Verify
  sources for the layouts and applications; identify synthetic stress tests
  accurately. (R2.9)
- [ ] Check consistency of every table, figure, theorem, supplement reference,
  metric definition, and numerical claim after inserting the final results.
- [ ] Replace proposed-work language in the response letter with the work
  actually completed, remove drafting notes and resolved pending items, and
  finalize the Associate Editor's summary.
- [ ] Check the final page count and retain the overlength explanation only if
  needed. Compile and inspect the final manuscript, supplement, and response.

## Completed substantive items

- [x] Constructive physical-feasibility and integer-equivalence proof.
- [x] General weak LP-dominance proof for SBM under the stated conventions.
- [x] Sufficient physical-horizon bound, integer objective coefficient, and
  flow-time certificate, with proofs in the supplement.
- [x] PBS motivation paragraph distinguishing retrieval optimization from
  system-level technology selection.
- [x] Explicit SBM scope at the beginning of Section 4 and in response R3.3;
  the pending decision about SLM experiments has been removed.
- [x] R3.M1 citation formatting, R3.M2 figure/caption clarification, and R3.M3
  Manhattan-distance terminology and historical bound verification.
- [x] Audit and correction of the historical greedy-improvement calculation.
  Refreshing that comparison after the expanded campaign remains open above.

## Reference files

- Review comments, current replies, and pending details:
  `revision_R1/response_letter_R1.tex`.
- Method and protocol details: `weighted_flow_certification_notes.md`.
- Project status and manuscript source location: `PROJECT_CONTEXT.md`.

This checklist incorporates later decisions from the conversation. If an older
proposed plan in the response letter conflicts with it, reconcile the reply
before submission rather than treating the older proposal as a new requirement.
