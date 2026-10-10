"""Regression checks for weighted solves and independent flow certificates.

Run with a licensed Gurobi Python interpreter:
    python3 -m unittest test_static_weighted_certification.py
"""

import contextlib
import io
import math
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import gurobipy as gp
from gurobipy import GRB

from escort_flow_static_gurobi import StaticEscortFlowGurobiSolver, StaticGurobiConfig
from load_flow_static_gurobi import LoadFlowStaticGurobiSolver, LoadFlowStaticGurobiConfig
from static_weighted_certification import (
    certification_eligibility, solve_weighted_or_certificate, sufficient_flow_horizon,
)


class Expression:
    def __init__(self, values, coefficients):
        self.values, self.coefficients = values, coefficients

    def __rmul__(self, coefficient):
        return Expression(self.values, {key: coefficient * value
                                       for key, value in self.coefficients.items()})

    def __add__(self, other):
        coefficients = dict(self.coefficients)
        for key, value in other.coefficients.items():
            coefficients[key] = coefficients.get(key, 0) + value
        return Expression(self.values, coefficients)

    def getValue(self):
        return sum(value * self.values[key] for key, value in self.coefficients.items())


class ScriptedModel:
    def __init__(self, *, bound=1012, status=GRB.OPTIMAL, solutions=1,
                 time_limit=300, work_limit=20, error=None):
        self.Params = SimpleNamespace(TimeLimit=time_limit, WorkLimit=work_limit,
                                      Cutoff=3, BestBdStop=2, BestObjStop=1)
        self.Status, self.SolCount, self.ObjBound = status, solutions, bound
        self.Runtime, self.Work = 7.0, 2.0
        self.error = error
        self.calls = []
        self.disposed = False

    def setObjective(self, expression, sense):
        self.expression, self.sense = expression, sense

    def optimize(self):
        self.calls.append(vars(self.Params).copy())
        if self.error:
            raise self.error

    def dispose(self):
        self.disposed = True


class ScriptedCertificateTests(unittest.TestCase):
    def run_model(self, mode="weighted_integer", *, flow=10, movements=12,
                  target=None, extract_error=None, **kwargs):
        model = ScriptedModel(**kwargs)
        values = dict(flow=flow, movements=movements)

        def extract():
            if extract_error:
                raise extract_error
            return dict(makespan=5, animation_moves=[[]])

        try:
            result = solve_weighted_or_certificate(
                model, Expression(values, {"flow": 1}),
                Expression(values, {"movements": 1}), extract,
                status_name=StaticEscortFlowGurobiSolver._status_name,
                solve_start=time.perf_counter(), mode=mode, target=target,
            )
        finally:
            self.assertTrue(model.disposed)
        return result, model

    def test_weighted_scaling_and_reported_units(self):
        result, model = self.run_model(bound=1011.5, status=GRB.TIME_LIMIT)
        self.assertEqual(model.expression.coefficients, {"flow": 100, "movements": 1})
        self.assertEqual(model.sense, GRB.MINIMIZE)
        self.assertAlmostEqual(result["objective"], 10.12)
        self.assertAlmostEqual(result["best_bound"], 10.115)
        self.assertAlmostEqual(result["absolute_gap"], .005)
        self.assertTrue(result["weighted_proven"])
        self.assertEqual(model.Params.MIPGap, 0)
        self.assertEqual(model.Params.MIPGapAbs, .999)
        self.assertEqual(model.Params.Cutoff, GRB.INFINITY)
        self.assertEqual(model.Params.BestBdStop, GRB.INFINITY)
        self.assertEqual(model.Params.BestObjStop, -GRB.INFINITY)

    def test_raw_optimal_does_not_override_scaled_integer_gap(self):
        for bound, expected in ((1011, False), (1011 + 1e-12, False),
                                (1011 + 5e-7, False), (1011.01, True)):
            with self.subTest(bound=bound):
                result, _ = self.run_model(bound=bound)
                self.assertEqual(result["status_name"], "OPTIMAL")
                self.assertEqual(result["weighted_proven"], expected)

    def test_weighted_bound_above_incumbent_cannot_certify_or_pass_gate(self):
        for bound in (1012 + 2e-6, 1012.5):
            with self.subTest(bound=bound):
                result, _ = self.run_model(bound=bound)
                self.assertFalse(result["bound_consistent"])
                self.assertFalse(result["weighted_proven"])
                self.assertEqual(certification_eligibility(result),
                                 (False, "INCONSISTENT_WEIGHTED_BOUND"))

    def test_scaled_gap_is_converted_before_eligibility_check(self):
        for bound, expected_gap, expected_eligible in ((1002, .1, False),
                                                       (1002.1, .099, True)):
            with self.subTest(bound=bound):
                result, _ = self.run_model(bound=bound, status=GRB.TIME_LIMIT)
                self.assertFalse(result["weighted_proven"])
                self.assertAlmostEqual(result["absolute_gap"], expected_gap)
                self.assertEqual(certification_eligibility(result)[0], expected_eligible)

    def test_bound_certifies_flow_without_any_incumbent(self):
        result, model = self.run_model("flow_certificate", target=10, bound=9.5,
                                       status=GRB.USER_OBJ_LIMIT, solutions=0)
        self.assertTrue(result["flow_proven"])
        self.assertFalse(result["has_solution"])
        self.assertFalse(result["counterexample"])
        self.assertIsNone(result["objective"])
        self.assertEqual(model.expression.coefficients, {"flow": 1})
        self.assertGreater(model.Params.BestBdStop, 10 - .999)
        self.assertGreater(model.Params.BestObjStop, 9)
        self.assertLess(model.Params.BestObjStop, 10)
        self.assertEqual(model.Params.Cutoff, GRB.INFINITY)

    def test_flow_proof_gap_boundary_includes_numerical_guard(self):
        for bound, expected in ((9, False), (9.001, False),
                                (9.001 + 5e-7, False), (9.001 + 2e-6, True),
                                (10, True), (10.01, False),
                                (math.inf, False), (math.nan, False), (None, False)):
            with self.subTest(bound=bound):
                result, _ = self.run_model("flow_certificate", target=10, bound=bound,
                                           status=GRB.TIME_LIMIT, solutions=0)
                self.assertEqual(result["flow_proven"], expected)

    def test_better_flow_is_counterexample_even_if_reported_bound_is_high(self):
        result, _ = self.run_model("flow_certificate", target=10, flow=9, bound=10,
                                   status=GRB.USER_OBJ_LIMIT)
        self.assertTrue(result["counterexample"])
        self.assertFalse(result["flow_proven"])
        self.assertEqual(result["flowtime"], 9)

    def test_unreliable_status_cannot_certify_or_claim_counterexample(self):
        for status in (GRB.NUMERIC, GRB.SUBOPTIMAL, GRB.INFEASIBLE,
                       GRB.INF_OR_UNBD, GRB.UNBOUNDED):
            with self.subTest(status=status):
                weighted, _ = self.run_model(status=status)
                self.assertFalse(weighted["weighted_proven"])
                flow, _ = self.run_model("flow_certificate", target=10,
                                         bound=10, status=status, solutions=0)
                self.assertFalse(flow["flow_proven"])
                counterexample, _ = self.run_model("flow_certificate", target=10,
                                                   flow=9, bound=9, status=status)
                self.assertFalse(counterexample["counterexample"])

    def test_no_weighted_incumbent_is_not_proven_even_with_optimal_status(self):
        result, _ = self.run_model(solutions=0)
        self.assertFalse(result["weighted_proven"])
        self.assertFalse(result["has_solution"])

    def test_independent_time_and_work_budgets_are_preserved(self):
        _, weighted = self.run_model(time_limit=300, work_limit=20)
        _, certificate = self.run_model("flow_certificate", target=10, bound=10,
                                        time_limit=270, work_limit=30)
        self.assertEqual(weighted.calls[0]["TimeLimit"], 300)
        self.assertEqual(certificate.calls[0]["TimeLimit"], 270)
        self.assertEqual(weighted.calls[0]["WorkLimit"], 20)
        self.assertEqual(certificate.calls[0]["WorkLimit"], 30)
        self.assertEqual(len(weighted.calls), 1)
        self.assertEqual(len(certificate.calls), 1)

    def test_model_is_disposed_on_optimization_and_extraction_failures(self):
        for kwargs in (dict(error=RuntimeError("optimize failed")),
                       dict(extract_error=RuntimeError("extract failed")),
                       dict(flow=10.25)):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises((RuntimeError, ValueError)):
                    self.run_model(**kwargs)


class EligibilityTests(unittest.TestCase):
    def test_gate_uses_strict_absolute_gap_in_original_units(self):
        for gap, expected in ((0, True), (.00999, True), (.099999, True),
                              (.1, False), (.11, False), (1, False),
                              (None, False), (math.inf, False), (math.nan, False)):
            with self.subTest(gap=gap):
                eligible, reason = certification_eligibility(dict(
                    has_solution=True, status_name="TIME_LIMIT", absolute_gap=gap,
                ))
                self.assertEqual(eligible, expected)
                self.assertEqual(reason, "" if expected else "WEIGHTED_GAP_TOO_LARGE")

    def test_proven_weighted_solution_is_eligible_and_missing_solution_is_not(self):
        self.assertEqual(certification_eligibility(dict(
            has_solution=True, status_name="OPTIMAL", weighted_proven=True,
        )), (True, ""))
        self.assertEqual(certification_eligibility(dict(
            has_solution=False, status_name="TIME_LIMIT", absolute_gap=0,
        )), (False, "NO_WEIGHTED_SOLUTION"))
        self.assertEqual(certification_eligibility(dict(
            has_solution=True, status_name="NUMERIC", weighted_proven=True,
        )), (False, "UNRELIABLE_WEIGHTED_STATUS"))

    def test_invalid_thresholds_are_rejected(self):
        for threshold in (0, -1, 1, math.inf, math.nan):
            with self.subTest(threshold=threshold), self.assertRaises(ValueError):
                certification_eligibility({}, threshold)

    def test_common_horizon_covers_each_targets_possible_arrival(self):
        targets, outputs = {(0, 0), (2, 0), (2, 2)}, {(0, 0)}
        # With F <= 12 and minimum distances 0, 2, 4, no arrival can
        # exceed 12 - 6 + 4 = 10. LF q[10] represents this arrival;
        # EF permits a target movement at t=9, so T=10 also covers it.
        self.assertEqual(sufficient_flow_horizon(targets, outputs, 12), 10)
        self.assertEqual(sufficient_flow_horizon(set(), outputs, 0), 0)


class WeightedIntegrationTests(unittest.TestCase):
    BACKENDS = (("escort", StaticEscortFlowGurobiSolver, StaticGurobiConfig),
                ("load", LoadFlowStaticGurobiSolver, LoadFlowStaticGurobiConfig))

    def configuration(self, backend, mode="weighted_integer", **overrides):
        name, _, config_class = backend
        config = dict(Lx=3, Ly=2, output_cells=((0, 0),), beta=1, gamma=.01,
                      time_limit=15, threads=1, objective_mode=mode)
        config.update(dict(alpha=0, move_method="BM") if name == "load"
                      else dict(retrieval_mode="leave"))
        config.update(overrides)
        return config_class(**config)

    def solve_case(self, backend, mode="weighted_integer", **overrides):
        with contextlib.redirect_stdout(io.StringIO()):
            solver = backend[1](self.configuration(backend, mode, **overrides))
            try:
                return solver.solve({(1, 1), (2, 0)}, {(0, 0), (2, 1)}, 8)
            finally:
                solver.close()

    def test_scaled_weighted_preserves_legacy_optimum_for_both_formulations(self):
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]):
                legacy = self.solve_case(backend, "legacy")
                weighted = self.solve_case(backend)
                self.assertTrue(weighted["weighted_proven"])
                self.assertEqual((weighted["flowtime"], weighted["movements"]), (8, 9))
                self.assertAlmostEqual(weighted["objective"], legacy["objective"])
                self.assertAlmostEqual(weighted["objective"], 8.09)
                self.assertLess(weighted["absolute_gap"], .01)

    def test_pure_flow_can_certify_optimal_candidate_or_find_counterexample(self):
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]):
                certified = self.solve_case(backend, "flow_certificate", certification_target=8)
                self.assertTrue(certified["flow_proven"])
                self.assertFalse(certified["counterexample"])
                counterexample = self.solve_case(backend, "flow_certificate", certification_target=9)
                self.assertTrue(counterexample["counterexample"])
                self.assertEqual(counterexample["flowtime"], 8)
                self.assertFalse(counterexample["flow_proven"])

    def test_backend_preserves_independent_configured_budgets(self):
        original_optimize = gp.Model.optimize
        calls = []

        def capture(model, *args, **kwargs):
            calls.append((model.Params.TimeLimit, model.Params.WorkLimit,
                          model.Params.MIPGap, model.Params.MIPGapAbs))
            return original_optimize(model, *args, **kwargs)

        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]), patch.object(gp.Model, "optimize", capture):
                self.solve_case(backend, time_limit=15, work_limit=20)
                self.solve_case(backend, "flow_certificate", certification_target=8,
                                time_limit=12, work_limit=30)
                self.assertEqual(calls[-2:], [(15, 20, 0, .999), (12, 30, 0, .999)])

    def test_flow_stop_retains_nonminimal_movement_start_in_both_modes(self):
        # Both independent row shifts are feasible; the second is unnecessary.
        # F=1 reaches the distance bound, but M=2 exceeds the optimum M=1.
        targets, escorts = {(1, 0)}, {(0, 0), (0, 1)}
        target_moves = [{0: ((1, 0), (0, 0))}]
        escort_moves = [[(0, 0, 1, 0), (0, 1, 1, 1)]]
        for backend in self.BACKENDS:
            for retrieval_mode in ("leave", "continue"):
                results = []
                for enabled in (True, False):
                    with self.subTest(backend=backend[0], mode=retrieval_mode, enabled=enabled):
                        with contextlib.redirect_stdout(io.StringIO()):
                            solver = backend[1](self.configuration(
                                backend, retrieval_mode=retrieval_mode, weight_scale=4,
                                gamma=.25, stop_at_flow_proof=enabled))
                            try:
                                horizon = 1 if backend[0] == "escort" else 2
                                start = solver.build_warmstart_from_trace(
                                    targets, escorts, horizon, target_moves, escort_moves)
                                result = solver.solve(targets, escorts, horizon, warmstart=start)
                            finally:
                                solver.close()
                        results.append(result)
                stopped, complete = results
                self.assertEqual((stopped["flowtime"], stopped["movements"], stopped["makespan"]), (1, 2, 1))
                self.assertTrue(stopped["has_solution"])
                self.assertTrue(stopped["flow_proven"])
                self.assertTrue(stopped["final_flow_proven"])
                self.assertFalse(stopped["weighted_proven"])
                self.assertEqual(stopped["stop_reason"], "FLOW_PROVEN_EARLY")
                self.assertEqual(stopped["status_name"], "INTERRUPTED")
                self.assertFalse(stopped["extension_used"])
                self.assertEqual(stopped["optimization_calls"], 1)
                self.assertEqual(stopped["phase1_snapshot"]["movements"], 2)
                self.assertEqual(stopped["phase1_snapshot"]["runtime"], stopped["runtime"])
                self.assertIsNotNone(stopped["first_flow_proof_runtime"])
                self.assertIsNotNone(stopped["first_flow_proof_cpu_time"])
                self.assertIsNotNone(stopped["animation_moves"])
                self.assertIn("solution_warmstart", stopped)
                self.assertEqual((complete["flowtime"], complete["movements"]), (1, 1))
                self.assertTrue(complete["weighted_proven"])

    def test_invalid_early_stop_options_fail_before_solver_creation(self):
        for backend in self.BACKENDS:
            for overrides in (dict(stop_at_flow_proof=1),
                              dict(stop_at_flow_proof=True, objective_mode="legacy"),
                              dict(stop_at_flow_proof=True, objective_mode="flow_certificate", certification_target=8),
                              dict(stop_at_flow_proof=True, time_limit=None),
                              dict(stop_at_flow_proof=True, time_limit=0),
                              dict(stop_at_flow_proof=True, time_limit=math.inf)):
                with self.subTest(backend=backend[0], overrides=overrides), self.assertRaises(ValueError):
                    backend[1](self.configuration(backend, **overrides))

    def test_unsupported_modes_weights_and_targets_fail_before_solver_creation(self):
        invalid = [dict(lp=True), dict(lexicographic=True), dict(beta=2),
                   dict(gamma=.1), dict(certification_target=10),
                   dict(objective_mode="other"), dict(objective_mode="flow_certificate"),
                   dict(objective_mode="flow_certificate", certification_target=-1),
                   dict(objective_mode="flow_certificate", certification_target=8.5)]
        for backend in self.BACKENDS:
            special = [dict(alpha=1), dict(move_method="LM")] if backend[0] == "load" else [dict(retrieval_mode="stay")]
            for overrides in invalid + special:
                with self.subTest(backend=backend[0], overrides=overrides), self.assertRaises(ValueError):
                    backend[1](self.configuration(backend, **overrides))

    def test_cutoffs_are_rejected_for_weighted_and_certificate_modes(self):
        for backend in self.BACKENDS:
            for mode in ("weighted_integer", "flow_certificate"):
                with self.subTest(backend=backend[0], mode=mode), contextlib.redirect_stdout(io.StringIO()):
                    config = self.configuration(backend, mode, certification_target=8 if mode == "flow_certificate" else None)
                    solver = backend[1](config)
                    try:
                        with self.assertRaisesRegex(ValueError, "cutoff"):
                            solver.solve({(1, 1)}, {(0, 0)}, 8, objective_cutoff=10)
                    finally:
                        solver.close()


if __name__ == "__main__":
    unittest.main()
