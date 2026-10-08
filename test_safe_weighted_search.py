"""Regression checks for one uninterrupted safe-weighted search.

The first-budget candidate is immutable while the same search tree is used
to establish or refute its flow-time optimality during the extension.
"""

import contextlib
import io
import math
import time
import unittest
from unittest import mock

import gurobipy as gp
from gurobipy import GRB

from static_safe_weighted_search import SafeWeightedContinuation
from static_weighted_certification import solve_weighted_or_certificate
from escort_flow_static_gurobi import StaticEscortFlowGurobiSolver


class CallbackExpression:
    def __init__(self, terms, constant=0):
        self.terms = list(terms)
        self.constant = constant

    def size(self):
        return len(self.terms)

    def getVar(self, index):
        return self.terms[index][0]

    def getCoeff(self, index):
        return self.terms[index][1]

    def getConstant(self):
        return self.constant


class CallbackModel:
    """Only expose legal callback data for the scripted callback location."""

    def __init__(self):
        self.values = {}
        self.queries = []
        self.solution = {}
        self.termination_count = 0
        self.solution_error = None

    def set_event(self, where, runtime, *, flow=5, movements=7,
                  bound=180, best=212, objective=None, nodes=0, work=None):
        self.where = where
        self.solution = {"flow": flow, "movement": movements}
        self.values = {
            GRB.Callback.RUNTIME: runtime,
            GRB.Callback.WORK: runtime / 10 if work is None else work,
        }
        if where == GRB.Callback.MIP:
            self.values.update({
                GRB.Callback.MIP_OBJBST: best,
                GRB.Callback.MIP_OBJBND: bound,
                GRB.Callback.MIP_NODCNT: nodes,
            })
        elif where == GRB.Callback.MIPSOL:
            self.values.update({
                GRB.Callback.MIPSOL_OBJ: 41 * flow + movements if objective is None else objective,
                GRB.Callback.MIPSOL_OBJBST: best,
                GRB.Callback.MIPSOL_OBJBND: bound,
                GRB.Callback.MIPSOL_NODCNT: nodes,
            })
        return self

    def cbGet(self, what):
        self.queries.append((self.where, what))
        if what not in self.values:
            raise AssertionError(f"Illegal query {what} at callback {self.where}")
        return self.values[what]

    def cbGetSolution(self, variables):
        if self.where != GRB.Callback.MIPSOL:
            raise AssertionError("Solution vector queried outside MIPSOL")
        if self.solution_error is not None:
            raise self.solution_error
        if isinstance(variables, (list, tuple)):
            return [self.solution[variable] for variable in variables]
        return self.solution[variables]

    def terminate(self):
        self.termination_count += 1


class ContinuousSearchTests(unittest.TestCase):
    CONTEXT = dict(
        targets=((1, 0), (2, 0)), outputs=((0, 0),), cell_count=12,
        escort_count=2, physical_horizon=4,
        weighted_time_limit=300, extension_time_limit=300,
    )

    def setUp(self):
        self.callback = SafeWeightedContinuation(
            CallbackExpression([("flow", 1)]),
            CallbackExpression([("movement", 1)]), 41, dict(self.CONTEXT))
        self.model = CallbackModel()

    def event(self, where, runtime, **kwargs):
        self.model.set_event(where, runtime, **kwargs)
        with contextlib.redirect_stdout(io.StringIO()):
            self.callback(self.model, where)

    @staticmethod
    def result(*, flow=5, movements=7, bound=195, runtime=400,
               status="INTERRUPTED", proven=False, **overrides):
        objective = 41 * flow + movements
        result = dict(
            has_solution=True, flowtime=flow, movements=movements,
            scaled_objective=objective, scaled_best_bound=bound,
            scaled_absolute_gap=abs(objective - bound),
            objective=objective / 41, best_bound=bound / 41,
            absolute_gap=abs(objective - bound) / 41,
            weight_scale=41, bound_consistent=True,
            weighted_proven=proven, status_name=status, runtime=runtime,
            work=runtime / 10, node_count=10,
        )
        result.update(overrides)
        return result

    def seed_candidate(self, *, bound=180, runtime=20, flow=5, movements=7):
        objective = 41 * flow + movements
        self.event(GRB.Callback.MIPSOL, runtime, flow=flow, movements=movements,
                   bound=bound, best=objective)
        self.event(GRB.Callback.MIP, runtime + 1, bound=bound, best=objective)

    def test_flow_proof_does_not_end_the_weighted_budget_early(self):
        self.seed_candidate(bound=195)
        self.event(GRB.Callback.MIP, 299, bound=195)
        self.assertEqual(self.model.termination_count, 0)
        self.event(GRB.Callback.MIP, 300.1, bound=195)
        self.assertGreater(self.model.termination_count, 0)
        extra = self.callback.finalize(self.result(runtime=300.1))
        self.assertTrue(extra["flow_proven"])
        self.assertEqual(extra["phase1_snapshot"]["flowtime"], 5)

    def test_extension_continues_and_stops_when_fixed_candidate_is_certified(self):
        self.seed_candidate()
        self.event(GRB.Callback.MIP, 299, bound=180, nodes=10)
        self.event(GRB.Callback.MIP, 300.1, bound=180, nodes=11)
        self.assertEqual(self.model.termination_count, 0)
        self.event(GRB.Callback.MIP, 400, bound=195, nodes=30)
        self.assertGreater(self.model.termination_count, 0)
        extra = self.callback.finalize(self.result())
        self.assertTrue(extra["flow_proven"])
        self.assertTrue(extra["extension_used"])
        self.assertEqual(extra["optimization_calls"], 1)
        snapshot = extra["phase1_snapshot"]
        self.assertEqual((snapshot["flowtime"], snapshot["movements"]), (5, 7))
        self.assertEqual(snapshot["scaled_best_bound"], 180)

    def test_nonimproving_predeadline_mipsol_cannot_replace_best_candidate(self):
        self.seed_candidate()
        self.event(GRB.Callback.MIPSOL, 100, flow=5, movements=30, best=212)
        self.event(GRB.Callback.MIP, 300.1, bound=195)
        extra = self.callback.finalize(self.result(runtime=300.1))
        snapshot = extra["phase1_snapshot"]
        self.assertEqual((snapshot["flowtime"], snapshot["movements"]), (5, 7))

    def test_late_equal_flow_movement_improvement_does_not_replace_snapshot(self):
        self.seed_candidate(bound=180)
        # The first callback after the cutoff is itself an improved solution.
        self.event(GRB.Callback.MIPSOL, 301, flow=5, movements=2, best=212, bound=180)
        self.assertEqual(self.model.termination_count, 0)
        self.event(GRB.Callback.MIP, 310, best=207, bound=195)
        extra = self.callback.finalize(self.result(movements=2, runtime=310))
        snapshot = extra["phase1_snapshot"]
        self.assertEqual((snapshot["flowtime"], snapshot["movements"]), (5, 7))
        self.assertEqual(snapshot["scaled_objective"], 212)
        self.assertTrue(extra["flow_proven"])

    def test_late_lower_flow_is_counterexample_even_with_worse_weighted_objective(self):
        self.seed_candidate(bound=130)
        self.event(GRB.Callback.MIPSOL, 301, flow=4, movements=100,
                   best=212, bound=130)
        self.assertGreater(self.model.termination_count, 0)
        extra = self.callback.finalize(self.result(bound=130, runtime=301))
        snapshot = extra["phase1_snapshot"]
        self.assertEqual((snapshot["flowtime"], snapshot["movements"]), (5, 7))
        self.assertFalse(extra["flow_proven"])
        self.assertTrue(extra["counterexample"])

    def test_bound_after_deadline_is_not_backfilled_into_phase1_snapshot(self):
        self.seed_candidate(bound=180)
        self.event(GRB.Callback.MIP, 299, bound=181)
        self.event(GRB.Callback.MIP, 301, bound=195)
        extra = self.callback.finalize(self.result(runtime=301))
        self.assertEqual(extra["phase1_snapshot"]["scaled_best_bound"], 181)
        self.assertTrue(extra["flow_proven"])

    def test_strict_bound_threshold_and_invalid_bound_do_not_stop_search(self):
        for bound in (194, 194.001, 194.0010005, math.nan, math.inf, 213):
            with self.subTest(bound=bound):
                self.setUp()
                self.seed_candidate()
                self.event(GRB.Callback.MIP, 301, bound=bound)
                self.assertEqual(self.model.termination_count, 0)

    def test_callback_exceptions_terminate_and_are_raised_to_the_caller(self):
        self.model.solution_error = RuntimeError("scripted solution retrieval failure")
        self.event(GRB.Callback.MIPSOL, 20)
        self.assertGreater(self.model.termination_count, 0)
        with self.assertRaisesRegex(RuntimeError, "callback failed") as raised:
            self.callback.raise_if_failed()
        self.assertEqual(str(raised.exception.__cause__), "scripted solution retrieval failure")

    def test_polling_callback_does_not_query_unavailable_information(self):
        self.event(GRB.Callback.POLLING, 400)
        self.assertFalse(self.model.queries)
        self.assertEqual(self.model.termination_count, 0)

    def test_solution_observed_exactly_at_cutoff_belongs_to_first_phase(self):
        self.seed_candidate(bound=130)
        self.event(GRB.Callback.MIPSOL, 300, flow=4, movements=12, best=212, bound=130)
        self.event(GRB.Callback.MIP, 301, best=176, bound=145)
        extra = self.callback.finalize(self.result(flow=4, movements=12, bound=145, runtime=301))
        snapshot = extra["phase1_snapshot"]
        self.assertEqual((snapshot["flowtime"], snapshot["movements"]), (4, 12))
        self.assertEqual(snapshot["incumbent_runtime"], 300)
        self.assertTrue(extra["flow_proven"])
        self.assertFalse(extra["counterexample"])

    def test_natural_early_optimum_uses_actual_final_result_as_snapshot(self):
        final = self.result(bound=212, runtime=12, status="OPTIMAL", proven=True)
        extra = self.callback.finalize(final)
        self.assertEqual(extra["phase1_snapshot"]["runtime"], 12)
        self.assertEqual(extra["phase1_snapshot"]["snapshot_source"], "SOLVE_FINISHED")
        self.assertFalse(extra["phase1_snapshot_missing"])
        self.assertFalse(extra["extension_used"])
        self.assertTrue(extra["flow_proven"])

    def test_missing_checkpoint_does_not_invent_predeadline_incumbent(self):
        final = self.result(runtime=601, status="TIME_LIMIT")
        extra = self.callback.finalize(final)
        self.assertTrue(extra["phase1_snapshot_missing"])
        self.assertFalse(extra["phase1_snapshot"]["has_solution"])
        self.assertIsNone(extra["phase1_snapshot"]["flowtime"])
        self.assertIsNone(extra["phase1_snapshot"]["scaled_best_bound"])
        self.assertTrue(final["has_solution"])
        self.assertEqual(final["flowtime"], 5)
        self.assertFalse(extra["flow_proven"])

    def test_late_first_solution_is_not_used_as_the_primary_candidate(self):
        self.event(GRB.Callback.MIP, 299, bound=130, best=GRB.INFINITY)
        self.event(GRB.Callback.MIPSOL, 301, bound=130)
        self.assertGreater(self.model.termination_count, 0)
        final = self.result(bound=130, runtime=301)
        extra = self.callback.finalize(final)
        self.assertFalse(extra["phase1_snapshot"]["has_solution"])
        self.assertEqual(extra["stop_reason"], "NO_PHASE1_SOLUTION")
        self.assertTrue(final["has_solution"])

    def test_unreliable_final_status_invalidates_saved_bound_certificates(self):
        for status in ("NUMERIC", "SUBOPTIMAL", "UNKNOWN"):
            for bound in (195, 211.5):
                with self.subTest(status=status, bound=bound):
                    self.setUp()
                    self.seed_candidate(bound=bound)
                    self.event(GRB.Callback.MIP, 301, bound=bound)
                    extra = self.callback.finalize(self.result(bound=bound, runtime=301, status=status))
                    self.assertFalse(extra["flow_proven"])
                    self.assertFalse(extra["final_flow_proven"])
                    self.assertFalse(extra["phase1_snapshot"]["flow_proven"])
                    self.assertFalse(extra["phase1_snapshot"]["weighted_proven"])

    def test_final_improved_solution_is_separate_and_finalize_does_not_mutate_it(self):
        self.seed_candidate(bound=180)
        self.event(GRB.Callback.MIP, 301, bound=180)
        self.event(GRB.Callback.MIPSOL, 320, flow=5, movements=2, best=212, bound=195)
        final = self.result(movements=2, bound=195, runtime=320)
        saved = dict(final)
        extra = self.callback.finalize(final)
        self.assertEqual(final, saved)
        self.assertEqual(extra["phase1_snapshot"]["movements"], 7)
        self.assertEqual(final["movements"], 2)
        self.assertEqual(extra["phase1_snapshot"]["scaled_best_bound"], 180)
        self.assertEqual(final["scaled_best_bound"], 195)

    def test_nonfinite_later_checkpoint_does_not_erase_valid_phase1_bound(self):
        self.seed_candidate(bound=180)
        self.event(GRB.Callback.MIP, 290, bound=181)
        self.event(GRB.Callback.MIP, 299, bound=-GRB.INFINITY)
        self.event(GRB.Callback.MIP, 301, bound=195)
        extra = self.callback.finalize(self.result(runtime=301))
        self.assertEqual(extra["phase1_snapshot"]["scaled_best_bound"], 181)
        self.assertEqual(extra["phase1_snapshot"]["bound_checkpoint_runtime"], 290)

    def test_expression_coefficients_and_constants_are_used_for_both_objectives(self):
        self.callback = SafeWeightedContinuation(
            CallbackExpression([("flow", 3)], constant=2),
            CallbackExpression([("movement", 2)], constant=3), 41, dict(self.CONTEXT))
        self.event(GRB.Callback.MIPSOL, 20, flow=1, movements=2, objective=212, bound=180)
        self.event(GRB.Callback.MIP, 301, bound=195)
        extra = self.callback.finalize(self.result(runtime=301))
        snapshot = extra["phase1_snapshot"]
        self.assertEqual((snapshot["flowtime"], snapshot["movements"]), (5, 7))
        self.assertEqual(snapshot["scaled_objective"], 212)

    def test_components_are_rounded_before_large_weight_amplifies_residuals(self):
        self.event(GRB.Callback.MIPSOL, 20, flow=5.00009, movements=7.00009,
                   objective=41 * 5.00009 + 7.00009, bound=180)
        self.event(GRB.Callback.MIP, 301, bound=195)
        self.callback.raise_if_failed()
        extra = self.callback.finalize(self.result(runtime=301))
        snapshot = extra["phase1_snapshot"]
        self.assertEqual((snapshot["flowtime"], snapshot["movements"]), (5, 7))
        self.assertEqual(snapshot["scaled_objective"], 212)

    def test_early_flow_proof_is_timed_but_initial_weighted_budget_continues(self):
        self.seed_candidate(bound=180)
        self.event(GRB.Callback.MIP, 40, bound=195, nodes=10)
        self.event(GRB.Callback.MIP, 299, bound=200, nodes=20)
        self.assertEqual(self.model.termination_count, 0)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 40)
        self.event(GRB.Callback.MIP, 300.1, bound=200, nodes=21)
        self.assertEqual(self.model.termination_count, 1)
        extra = self.callback.finalize(self.result(bound=200, runtime=300.2))
        self.assertEqual(extra["first_flow_proof_runtime"], 40)
        self.assertEqual(extra["stop_reason"], "FLOW_PROVEN_AT_PHASE_LIMIT")
        self.assertFalse(extra["extension_used"])
        self.assertFalse(extra["phase1_snapshot"]["weighted_proven"])

    def test_extension_stopping_proof_records_subsecond_bound_improvement(self):
        self.seed_candidate(bound=180)
        self.event(GRB.Callback.MIP, 299, bound=180, nodes=10)
        self.event(GRB.Callback.MIP, 300.1, bound=180, nodes=11)
        self.assertIsNone(self.callback.first_flow_proof)
        # A changed bound is checked immediately, even just 0.1s after the
        # preceding callback and with only one new node.
        self.event(GRB.Callback.MIP, 300.2, bound=195, nodes=12, work=31)
        self.assertEqual(self.model.termination_count, 1)
        extra = self.callback.finalize(self.result(bound=195, runtime=300.25))
        self.assertEqual(extra["first_flow_proof_runtime"], 300.2)
        self.assertEqual(extra["first_flow_proof_source"], "CALLBACK_MIP")
        self.assertEqual(extra["first_flow_proof_work"], 31)
        self.assertEqual(extra["stop_reason"], "FLOW_PROVEN_DURING_EXTENSION")
        self.assertTrue(extra["extension_used"])

    def test_phase1_bound_proof_keeps_its_first_time_when_later_bound_is_weaker(self):
        self.seed_candidate(bound=180)
        self.event(GRB.Callback.MIP, 298, bound=180, nodes=10)
        self.event(GRB.Callback.MIP, 299, bound=195, nodes=11)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 299)
        # A MIPSOL location may expose a weaker bound than the last MIP
        # observation. The phase snapshot retains its valid earlier bound.
        self.event(GRB.Callback.MIP, 300.1, bound=180, nodes=12)
        self.assertEqual(self.model.termination_count, 1)
        extra = self.callback.finalize(self.result(bound=195, runtime=300.2))
        self.assertEqual(extra["first_flow_proof_runtime"], 299)
        self.assertEqual(extra["first_flow_proof_scaled_bound"], 195)
        self.assertEqual(extra["first_flow_proof_source"], "CALLBACK_MIP")
        self.assertTrue(extra["flow_proven"])


class FlowProofTimingTests(unittest.TestCase):
    """Proof timing describes solver observations without ending optimization."""

    CONTEXT = ContinuousSearchTests.CONTEXT
    event = ContinuousSearchTests.event
    seed_candidate = ContinuousSearchTests.seed_candidate
    result = staticmethod(ContinuousSearchTests.result)

    def setUp(self):
        self.model = CallbackModel()
        self.callback = SafeWeightedContinuation(
            CallbackExpression([("flow", 1)]),
            CallbackExpression([("movement", 1)]), 41,
            dict(self.CONTEXT, stop_on_flow_proof=False))

    def test_proof_time_is_first_observed_bound_and_total_search_continues(self):
        self.seed_candidate(bound=180)
        self.event(GRB.Callback.MIP, 40, bound=195, nodes=123, work=4.5)
        self.event(GRB.Callback.MIP, 301, bound=197, nodes=234)
        self.event(GRB.Callback.MIP, 450, bound=209, nodes=345)
        self.assertEqual(self.model.termination_count, 0)
        extra = self.callback.finalize(self.result(bound=209, runtime=600, status="TIME_LIMIT"))
        self.assertEqual(extra["first_flow_proof_runtime"], 40)
        self.assertEqual(extra["first_flow_proof_node_count"], 123)
        self.assertEqual(extra["first_flow_proof_work"], 4.5)
        self.assertEqual(extra["first_flow_proof_flowtime"], 5)
        self.assertEqual(extra["first_flow_proof_scaled_bound"], 195)
        self.assertEqual(extra["first_flow_proof_source"], "CALLBACK_MIP")
        self.assertEqual(extra["first_flow_proof_method"], "weighted_gap")
        self.assertEqual(extra["stop_reason"], "TIME_LIMIT")
        self.assertTrue(extra["flow_proven"])
        self.assertFalse(extra["phase1_snapshot"]["weighted_proven"])

    def test_bound_improvement_is_checked_immediately_without_time_sampling(self):
        self.seed_candidate(bound=180, runtime=20)
        count = self.callback.flow_proof_check_count
        self.event(GRB.Callback.MIP, 20.01, bound=195, nodes=0)
        self.assertEqual(self.callback.flow_proof_check_count, count + 1)
        extra = self.callback.finalize(self.result(runtime=30, status="TIME_LIMIT"))
        self.assertEqual(extra["first_flow_proof_runtime"], 20.01)

    def test_unchanged_bounds_skip_checks_even_after_large_time_and_node_jumps(self):
        self.seed_candidate(bound=180, runtime=20)
        count = self.callback.flow_proof_check_count
        for runtime, nodes in ((21, 100), (100, 1000), (299, 100000), (450, 1000000)):
            self.event(GRB.Callback.MIP, runtime, bound=180, nodes=nodes)
        self.assertEqual(self.callback.flow_proof_check_count, count)
        self.assertIsNone(self.callback.first_flow_proof)

    def test_new_flow_rechecks_unchanged_bound_at_candidate_observation(self):
        self.seed_candidate(bound=145, runtime=20)
        self.event(GRB.Callback.MIPSOL, 20.1, flow=4, movements=12,
                   bound=145, best=212, nodes=1)
        extra = self.callback.finalize(self.result(flow=4, movements=12, bound=145,
                                                  runtime=30, status="TIME_LIMIT"))
        self.assertEqual(extra["first_flow_proof_runtime"], 20.1)
        self.assertEqual(extra["first_flow_proof_flowtime"], 4)
        self.assertEqual(extra["first_flow_proof_source"], "CALLBACK_MIPSOL")

    def test_nonweighted_incumbent_can_establish_analytical_flow_proof(self):
        self.seed_candidate(bound=130)
        self.event(GRB.Callback.MIPSOL, 100, flow=3, movements=100,
                   bound=130, best=212)
        final = self.result(bound=130, runtime=300, status="TIME_LIMIT")
        extra = self.callback.finalize(final)
        self.assertEqual(extra["phase1_snapshot"]["flowtime"], 5)
        self.assertEqual(extra["first_flow_proof_flowtime"], 3)
        self.assertEqual(extra["first_flow_proof_method"], "distance_bound")
        self.assertEqual(extra["first_flow_proof_runtime"], 100)
        self.assertTrue(extra["counterexample"])
        self.assertFalse(extra["flow_proven"])
        self.assertEqual(self.model.termination_count, 0)

    def test_missing_cutoff_incumbent_does_not_end_search_or_invent_snapshot(self):
        self.event(GRB.Callback.MIP, 299, bound=130, best=GRB.INFINITY)
        self.event(GRB.Callback.MIPSOL, 301, flow=3, movements=12, bound=130)
        self.event(GRB.Callback.MIP, 450, bound=134)
        self.assertEqual(self.model.termination_count, 0)
        extra = self.callback.finalize(self.result(flow=3, movements=12, bound=134,
                                                  runtime=600, status="TIME_LIMIT"))
        self.assertFalse(extra["phase1_snapshot"]["has_solution"])
        self.assertEqual(extra["first_flow_proof_runtime"], 301)
        self.assertTrue(extra["final_flow_proven"])
        self.assertEqual(extra["stop_reason"], "TIME_LIMIT")

    def test_lower_flow_after_cutoff_is_reported_without_ending_search(self):
        self.seed_candidate(bound=130)
        self.event(GRB.Callback.MIPSOL, 301, flow=4, movements=100,
                   best=212, bound=130)
        self.event(GRB.Callback.MIP, 450, bound=145)
        self.assertEqual(self.model.termination_count, 0)
        extra = self.callback.finalize(self.result(bound=145, runtime=600, status="TIME_LIMIT"))
        self.assertTrue(extra["counterexample"])
        self.assertFalse(extra["flow_proven"])
        self.assertEqual(extra["first_flow_proof_flowtime"], 4)
        self.assertEqual(extra["first_flow_proof_runtime"], 450)
        self.assertEqual(extra["stop_reason"], "TIME_LIMIT")

    def test_final_fallback_reports_final_runtime_without_backdating(self):
        self.seed_candidate(bound=180)
        self.event(GRB.Callback.MIP, 30, bound=190, nodes=10)
        extra = self.callback.finalize(self.result(bound=212, runtime=30.25,
                                                  status="OPTIMAL", proven=True))
        self.assertEqual(extra["first_flow_proof_runtime"], 30.25)
        self.assertEqual(extra["first_flow_proof_source"], "SOLVE_FINISHED")
        self.assertEqual(extra["first_flow_proof_node_count"], 10)
        self.assertEqual(extra["stop_reason"], "WEIGHTED_OPTIMUM")

    def test_time_limit_with_unproved_flow_leaves_timing_empty(self):
        self.seed_candidate(bound=180)
        extra = self.callback.finalize(self.result(bound=180, runtime=300, status="TIME_LIMIT"))
        for field in ("runtime", "node_count", "work", "flowtime", "scaled_bound", "source", "method"):
            self.assertIsNone(extra["first_flow_proof_" + field])
        self.assertFalse(extra["flow_proven"])

    def test_unreliable_final_status_clears_observed_proof_timing(self):
        for status in ("NUMERIC", "SUBOPTIMAL", "UNKNOWN"):
            with self.subTest(status=status):
                self.setUp()
                self.seed_candidate(bound=195)
                self.assertIsNotNone(self.callback.first_flow_proof)
                extra = self.callback.finalize(self.result(runtime=300, status=status))
                self.assertIsNone(extra["first_flow_proof_runtime"])
                self.assertTrue(extra["flow_proof_invalidated"])
                self.assertEqual(extra["flow_proof_invalidation_reason"], "UNRELIABLE_FINAL_STATUS")
                self.assertFalse(extra["flow_proven"])
                self.assertFalse(extra["final_flow_proven"])

    def test_better_flow_after_recorded_proof_invalidates_inconsistent_evidence(self):
        self.seed_candidate(bound=195)
        self.event(GRB.Callback.MIPSOL, 100, flow=4, movements=100, bound=130, best=212)
        extra = self.callback.finalize(self.result(flow=4, movements=100, bound=145,
                                                  runtime=300, status="TIME_LIMIT"))
        self.assertIsNone(extra["first_flow_proof_runtime"])
        self.assertTrue(extra["flow_proof_invalidated"])
        self.assertEqual(extra["flow_proof_invalidation_reason"], "BETTER_FLOW_AFTER_RECORDED_PROOF")
        self.assertFalse(extra["final_flow_proven"])

    def test_zero_extension_keeps_full_budget_behavior_after_cutoff(self):
        self.callback.extension_limit = 0
        self.seed_candidate(bound=195)
        self.event(GRB.Callback.MIP, 300.1, bound=195)
        self.assertEqual(self.model.termination_count, 0)
        extra = self.callback.finalize(self.result(runtime=300.1, status="TIME_LIMIT"))
        self.assertFalse(extra["extension_used"])
        self.assertEqual(extra["first_flow_proof_runtime"], 20)

    def test_bound_checks_reuse_precalculated_criterion_without_solution_reads(self):
        with mock.patch.object(self.callback, "_build_flow_criterion", wraps=self.callback._build_flow_criterion) as build:
            self.seed_candidate(bound=130)
            for runtime in (22, 24, 26, 28, 30):
                with mock.patch.object(self.model, "cbGetSolution", side_effect=AssertionError("Periodic proof read a solution vector")):
                    self.event(GRB.Callback.MIP, runtime, bound=runtime + 130)
            self.event(GRB.Callback.MIPSOL, 31, flow=5, movements=6, bound=130, best=212)
            self.event(GRB.Callback.MIPSOL, 32, flow=5, movements=5, bound=130, best=211)
            self.assertEqual(build.call_count, 1)
            self.assertEqual(build.call_args.args, (5,))
            self.assertEqual(set(self.callback.flow_criteria), {5})
            self.assertEqual(self.callback.flow_proof_check_count, 6)
            self.event(GRB.Callback.MIPSOL, 33, flow=4, movements=12, bound=130, best=210)
            self.assertEqual(build.call_count, 2)
            self.assertEqual(build.call_args.args, (4,))
            self.assertEqual(set(self.callback.flow_criteria), {4, 5})

    def test_same_flow_movement_improvements_do_not_recheck_criterion(self):
        self.seed_candidate(bound=180)
        checks = self.callback.flow_proof_check_count
        with mock.patch.object(self.callback, "_cached_flow_proof_method", wraps=self.callback._cached_flow_proof_method) as criterion:
            self.event(GRB.Callback.MIPSOL, 30, flow=5, movements=6, bound=180, best=212)
            self.event(GRB.Callback.MIPSOL, 40, flow=5, movements=5, bound=180, best=211)
            self.assertEqual(criterion.call_count, 0)
        self.assertEqual(self.callback.flow_proof_check_count, checks)
        self.assertEqual(self.callback.best_observed_weighted_objective, 210)

    def test_nonfinite_weaker_and_inconsistent_bounds_do_not_poison_cached_bound(self):
        self.seed_candidate(bound=180)
        checks = self.callback.flow_proof_check_count
        for index, bound in enumerate((math.nan, math.inf, -GRB.INFINITY, 179, 180, 213)):
            self.event(GRB.Callback.MIP, 30 + index, bound=bound)
        self.assertEqual(self.callback.flow_proof_check_count, checks)
        self.assertEqual(self.callback.best_flow_proof_bound, 180)
        self.event(GRB.Callback.MIP, 40, bound=195)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 40)
        self.assertEqual(self.callback.first_flow_proof["scaled_bound"], 195)

    def test_new_candidate_uses_stronger_bound_observed_before_any_candidate(self):
        self.event(GRB.Callback.MIP, 10, bound=195, best=GRB.INFINITY)
        self.assertEqual(self.callback.best_flow_proof_bound, 195)
        self.assertEqual(self.callback.flow_proof_check_count, 0)
        self.event(GRB.Callback.MIPSOL, 20, bound=180)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 20)
        self.assertEqual(self.callback.first_flow_proof["scaled_bound"], 195)
        self.assertEqual(self.callback.first_flow_proof["source"], "CALLBACK_MIPSOL")

    def test_inconsistent_provisional_bound_is_discarded_on_first_candidate(self):
        self.event(GRB.Callback.MIP, 10, bound=250, best=GRB.INFINITY)
        self.event(GRB.Callback.MIPSOL, 20, bound=180)
        self.assertIsNone(self.callback.first_flow_proof)
        self.assertEqual(self.callback.best_flow_proof_bound, 180)
        extra = self.callback.finalize(self.result(bound=195, runtime=30, status="TIME_LIMIT"))
        self.assertEqual(extra["first_flow_proof_runtime"], 30)
        self.assertEqual(extra["first_flow_proof_scaled_bound"], 195)
        self.assertEqual(extra["first_flow_proof_source"], "SOLVE_FINISHED")

    def test_bounds_are_validated_against_best_weighted_incumbent_not_old_flow_candidate(self):
        self.event(GRB.Callback.MIPSOL, 20, flow=5, movements=100, bound=180)
        self.event(GRB.Callback.MIPSOL, 21, flow=5, movements=7, bound=180, best=305)
        self.assertEqual(self.callback.best_flow_candidate["scaled_objective"], 305)
        self.assertEqual(self.callback.best_observed_weighted_objective, 212)
        self.event(GRB.Callback.MIP, 22, bound=250)
        self.assertIsNone(self.callback.first_flow_proof)
        self.assertEqual(self.callback.best_flow_proof_bound, 180)
        extra = self.callback.finalize(self.result(bound=195, runtime=30, status="TIME_LIMIT"))
        self.assertEqual(extra["first_flow_proof_runtime"], 30)
        self.assertEqual(extra["first_flow_proof_scaled_bound"], 195)

    def test_same_flow_better_movements_invalidate_inconsistent_recorded_proof_bound(self):
        self.event(GRB.Callback.MIPSOL, 20, flow=5, movements=100, bound=305)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 20)
        self.event(GRB.Callback.MIPSOL, 21, flow=5, movements=7, bound=180, best=305)
        self.assertIsNone(self.callback.first_flow_proof)
        extra = self.callback.finalize(self.result(bound=195, runtime=30, status="TIME_LIMIT"))
        self.assertIsNone(extra["first_flow_proof_runtime"])
        self.assertTrue(extra["flow_proof_invalidated"])
        self.assertEqual(extra["flow_proof_invalidation_reason"], "INCONSISTENT_RECORDED_FLOW_BOUND")
        self.assertFalse(extra["flow_proven"])
        self.assertFalse(extra["final_flow_proven"])

    def test_final_candidate_revalidates_recorded_bound_against_its_weighted_objective(self):
        self.event(GRB.Callback.MIPSOL, 20, flow=5, movements=100, bound=305)
        extra = self.callback.finalize(self.result(bound=195, runtime=30, status="TIME_LIMIT"))
        self.assertIsNone(extra["first_flow_proof_runtime"])
        self.assertEqual(extra["flow_proof_invalidation_reason"], "INCONSISTENT_RECORDED_FLOW_BOUND")
        self.assertFalse(extra["final_flow_proven"])

    def test_final_weighted_upper_bound_is_checked_even_when_best_flow_witness_has_smaller_flow(self):
        self.event(GRB.Callback.MIPSOL, 20, flow=4, movements=100, bound=264)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 20)
        final = self.result(flow=5, movements=7, bound=130, runtime=30, status="TIME_LIMIT")
        extra = self.callback.finalize(final)
        self.assertEqual(self.callback.best_flow_candidate["flowtime"], 4)
        self.assertEqual(self.callback.best_observed_weighted_objective, 212)
        self.assertIsNone(extra["first_flow_proof_runtime"])
        self.assertEqual(extra["flow_proof_invalidation_reason"], "INCONSISTENT_RECORDED_FLOW_BOUND")
        self.assertFalse(extra["final_flow_proven"])

    def test_every_strict_bound_improvement_is_seen_including_tiny_threshold_crossing(self):
        self.seed_candidate(bound=194.0010005)
        self.assertIsNone(self.callback.first_flow_proof)
        self.event(GRB.Callback.MIP, 21.01, bound=194.0010011)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 21.01)

    def test_unchanged_bound_after_cutoff_does_not_repeat_frozen_candidate_check(self):
        self.callback.stop_on_flow_proof = True
        self.seed_candidate(bound=180)
        with mock.patch.object(self.callback, "_cached_flow_proof_method", wraps=self.callback._cached_flow_proof_method) as criterion:
            self.event(GRB.Callback.MIP, 300.1, bound=180)
            self.assertEqual(criterion.call_count, 1)  # Required cutoff assessment.
            self.assertIsNotNone(self.callback.phase1_snapshot)
            self.event(GRB.Callback.MIP, 350, bound=180)
            self.event(GRB.Callback.MIP, 500, bound=180)
            self.assertEqual(criterion.call_count, 1)
        self.assertEqual(self.model.termination_count, 0)
        self.event(GRB.Callback.MIP, 500.01, bound=195)
        self.assertEqual(self.model.termination_count, 1)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 500.01)

    def test_instance_geometry_is_read_once_and_public_helpers_are_not_in_callback(self):
        from static_weighted_certification import _safe_weight_instance
        with mock.patch("static_safe_weighted_search._safe_weight_instance", wraps=_safe_weight_instance) as geometry:
            self.setUp()
            with mock.patch("static_safe_weighted_search.gap_flow_certificate", side_effect=AssertionError("Expensive certificate helper in callback")), \
                    mock.patch("static_safe_weighted_search.safe_weight_parameters", side_effect=AssertionError("Instance parameter helper in callback")):
                self.seed_candidate(bound=130)
                self.event(GRB.Callback.MIP, 50, bound=180)
                self.event(GRB.Callback.MIP, 299, bound=190)
                self.event(GRB.Callback.MIP, 301, bound=195)
                self.event(GRB.Callback.MIP, 450, bound=200)
            self.assertEqual(geometry.call_count, 1)
        self.assertEqual(self.callback.first_flow_proof["runtime"], 301)
        extra = self.callback.finalize(self.result(bound=200, runtime=600, status="TIME_LIMIT"))
        self.assertTrue(extra["flow_proven"])

    def test_cached_scalar_criterion_agrees_with_authoritative_certificate(self):
        for flow, movements in ((3, 2), (4, 12), (5, 7), (5, 100)):
            objective = 41 * flow + movements
            incumbent = dict(has_solution=True, flowtime=flow, movements=movements,
                             scaled_objective=objective)
            for bound in (None, math.nan, math.inf, 130, 144.001, 144.0010005, 145,
                          194.001, 194.0010005, 195, objective - .999999, objective, objective + 1):
                with self.subTest(flow=flow, movements=movements, bound=bound):
                    method = self.callback._cached_flow_proof_method(incumbent, bound)
                    authoritative = self.callback._proof(incumbent, bound)
                    self.assertEqual(bool(method), authoritative["flow_proven"])
                    self.assertEqual(method or "", authoritative["proof_source"])


class CountingModel:
    def __init__(self, model):
        self.model = model
        self.optimize_calls = []
        self.objectives = []
        self.disposed = False

    def __getattr__(self, key):
        return getattr(self.model, key)

    def setObjective(self, expression, sense):
        self.objectives.append((expression, sense))
        self.model.setObjective(expression, sense)

    def optimize(self, callback):
        self.optimize_calls.append({
            "time_limit": self.model.Params.TimeLimit,
            "focus": self.model.Params.MIPFocus,
            "relative_gap": self.model.Params.MIPGap,
            "absolute_gap": self.model.Params.MIPGapAbs,
            "callback": callback,
        })
        self.model.optimize(callback)

    def dispose(self):
        self.disposed = True
        self.model.dispose()


class LiveContinuousSearchTests(unittest.TestCase):
    def test_real_gurobi_uses_one_optimize_with_full_budget_and_fixed_objective(self):
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            model = CountingModel(gp.Model(env=env))
            model.Params.Threads = 1
            model.Params.MIPFocus = 0
            flow = model.addVar(lb=3, vtype=GRB.INTEGER)
            movement = model.addVar(lb=2, vtype=GRB.INTEGER)
            flow_expression, movement_expression = gp.LinExpr(flow), gp.LinExpr(movement)
            context = dict(ContinuousSearchTests.CONTEXT,
                           weighted_time_limit=15, extension_time_limit=2)
            result = solve_weighted_or_certificate(
                model, flow_expression, movement_expression,
                lambda: dict(makespan=2, animation_moves=[]),
                status_name=StaticEscortFlowGurobiSolver._status_name,
                solve_start=time.perf_counter(), mode="weighted_integer",
                weight_scale=41, flow_proof_context=context)
        self.assertTrue(model.disposed)
        self.assertEqual(len(model.optimize_calls), 1)
        call = model.optimize_calls[0]
        self.assertEqual(call["time_limit"], 17)
        self.assertEqual(call["focus"], 0)
        self.assertEqual(call["relative_gap"], 0)
        self.assertEqual(call["absolute_gap"], .999)
        self.assertIsInstance(call["callback"], SafeWeightedContinuation)
        self.assertEqual(len(model.objectives), 1)
        self.assertEqual([model.objectives[0][0].getCoeff(i) for i in range(2)], [41, 1])
        self.assertTrue(result["weighted_proven"])
        self.assertTrue(result["flow_proven"])
        self.assertTrue(result["final_flow_proven"])
        self.assertFalse(result["extension_used"])
        self.assertEqual(result["optimization_calls"], 1)
        self.assertEqual((result["flowtime"], result["movements"]), (3, 2))
        self.assertEqual(result["phase1_snapshot"]["flowtime"], 3)
        self.assertEqual(result["phase1_snapshot"]["snapshot_source"], "SOLVE_FINISHED")


if __name__ == "__main__":
    unittest.main()
