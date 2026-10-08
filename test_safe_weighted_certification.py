"""Mathematical and backend checks for safely weighted flow certificates."""

import contextlib
import io
import math
import time
import unittest

from gurobipy import GRB

from escort_flow_static_gurobi import StaticEscortFlowGurobiSolver, StaticGurobiConfig
from load_flow_static_gurobi import LoadFlowStaticGurobiSolver, LoadFlowStaticGurobiConfig
from static_weighted_certification import (
    gap_flow_certificate, safe_weight_parameters, solve_weighted_or_certificate,
)
from test_static_weighted_certification import Expression, ScriptedModel


class SafeWeightMathematicsTests(unittest.TestCase):
    TARGETS = ((1, 0), (2, 0))
    OUTPUTS = ((0, 0),)

    def certificate(self, *, flow=5, movements=7, bound=180, horizon=4,
                    weight=41, **overrides):
        weighted = dict(has_solution=True, status_name="TIME_LIMIT", flowtime=flow,
                        movements=movements, weight_scale=weight,
                        scaled_best_bound=bound, bound_consistent=True)
        weighted.update(overrides)
        return gap_flow_certificate(weighted, self.TARGETS, self.OUTPUTS, 12, 2,
                                    horizon, weight)

    def test_safe_weight_uses_flow_upper_bound_and_initial_load_count(self):
        # d=(1,2); H=5-3+2=4, at most ten loads move in each period.
        self.assertEqual(safe_weight_parameters(self.TARGETS, self.OUTPUTS, 12, 2, 5),
                         dict(flow_weight=41, movement_bound=40, flow_horizon=4))
        self.assertEqual(safe_weight_parameters(((0, 0),), self.OUTPUTS, 12, 2, 0),
                         dict(flow_weight=1, movement_bound=0, flow_horizon=0))

    def test_safe_weight_separates_every_worse_flow_from_bounded_lex_optimum(self):
        # This exhausts a small abstract objective space, including plans with
        # redundant movement counts exceeding U. Only a lex-optimal plan needs
        # to satisfy the movement upper bound used in the proof.
        for upper_movements in range(12):
            weight = upper_movements + 1
            for optimal_flow in range(5):
                for optimal_movements in range(upper_movements + 1):
                    optimal_value = weight * optimal_flow + optimal_movements
                    for worse_flow in range(optimal_flow + 1, optimal_flow + 4):
                        for worse_movements in range(2 * upper_movements + 3):
                            self.assertLess(optimal_value, weight * worse_flow + worse_movements)

    def test_safe_parameters_reject_impossible_or_noninteger_inputs(self):
        for cell_count, escorts, flow in ((0, 0, 5), (12, 13, 5), (12, 11, 5),
                                          (12, 2, 2), (12, 2, 5.5), (12, 2, math.inf),
                                          (12, 2, True), (12.5, 2, 5)):
            with self.subTest(cell_count=cell_count, escorts=escorts, flow=flow), self.assertRaises(ValueError):
                safe_weight_parameters(self.TARGETS, self.OUTPUTS, cell_count, escorts, flow)
        with self.assertRaises(ValueError):
            safe_weight_parameters(self.TARGETS, (), 12, 2, 5)

    def test_gap_above_one_can_certify_flow_without_proving_weighted_optimum(self):
        result = self.certificate(bound=195)
        self.assertTrue(result["flow_proven"])
        self.assertEqual(result["proof_source"], "weighted_gap")
        self.assertEqual(result["flow_lower_bound"], 5)
        self.assertEqual(result["better_flow_horizon"], 3)
        self.assertEqual(result["better_flow_movement_bound"], 30)
        self.assertAlmostEqual(result["weighted_bound_threshold"], 194.001)
        self.assertAlmostEqual(result["weighted_gap_threshold"], 17.999)
        self.assertEqual(result["scaled_absolute_gap"], 17)

    def test_strict_gap_test_keeps_numerical_margin(self):
        for bound, proven in ((194, False), (194 + 1e-6, False),
                               (194.001, False), (194.0010005, False),
                               (194.001002, True), (195, True)):
            with self.subTest(bound=bound):
                self.assertEqual(self.certificate(bound=bound)["flow_proven"], proven)

    def test_gap_below_flow_weight_alone_does_not_certify(self):
        result = self.certificate(bound=180)
        self.assertLess(result["scaled_absolute_gap"], 41)
        self.assertFalse(result["flow_proven"])
        self.assertEqual(result["reason"], "WEIGHTED_GAP_TOO_LARGE")

    def test_proof_requires_horizon_containing_hypothetical_better_flow(self):
        # H_minus=3; the incumbent itself can fit in an even shorter horizon,
        # but that does not license a global lower-bound certificate.
        for horizon, proven in ((2, False), (3, True), (4, True)):
            with self.subTest(horizon=horizon):
                result = self.certificate(bound=195, horizon=horizon)
                self.assertEqual(result["flow_proven"], proven)
                self.assertEqual(result["horizon_sufficient"], proven)

    def test_distance_bound_proves_without_a_solver_bound(self):
        result = self.certificate(flow=3, bound=None, horizon=0)
        self.assertTrue(result["flow_proven"])
        self.assertEqual(result["proof_source"], "distance_bound")
        self.assertEqual(result["flow_lower_bound"], 3)

    def test_invalid_results_cannot_create_a_certificate(self):
        cases = [dict(has_solution=False), dict(flow=5.5), dict(movements=7.1),
                 dict(flow=-1), dict(flow=2), dict(bound=math.inf), dict(bound=math.nan),
                 dict(bound=None), dict(bound=212.01), dict(bound_consistent=False),
                 dict(scaled_objective=213), dict(weight_scale=100)]
        for status in ("NUMERIC", "SUBOPTIMAL", "INFEASIBLE", "UNBOUNDED", "INF_OR_UNBD",
                       "LOADED", "INPROGRESS", "UNKNOWN", None):
            cases.append(dict(status_name=status))
        for kwargs in cases:
            with self.subTest(kwargs=kwargs):
                self.assertFalse(self.certificate(**kwargs)["flow_proven"])

    def test_saved_gap_fields_are_not_trusted_without_recomputing(self):
        result = self.certificate(bound=180, scaled_absolute_gap=0, absolute_gap=0,
                                  weighted_proven=True)
        self.assertFalse(result["flow_proven"])
        self.assertEqual(result["scaled_absolute_gap"], 32)

    def test_normalized_bound_fallback_preserves_legacy_data_compatibility(self):
        weighted = dict(has_solution=True, status_name="TIME_LIMIT", flowtime=5,
                        movements=7, best_bound=195 / 41)
        result = gap_flow_certificate(weighted, self.TARGETS, self.OUTPUTS, 12, 2, 4, 41)
        self.assertTrue(result["flow_proven"])
        self.assertEqual(result["scaled_absolute_gap"], 17)


class DynamicWeightBackendTests(unittest.TestCase):
    def test_unfinished_or_unknown_solver_status_never_proves_a_tight_gap(self):
        values = dict(flow=10, movements=12)
        for status in (GRB.LOADED, GRB.INPROGRESS, 12345):
            for mode in ("weighted_integer", "flow_certificate"):
                with self.subTest(status=status, mode=mode):
                    model = ScriptedModel(status=status, bound=1012 if mode == "weighted_integer" else 10)
                    result = solve_weighted_or_certificate(
                        model, Expression(values, {"flow": 1}), Expression(values, {"movements": 1}),
                        lambda: dict(makespan=5, animation_moves=[]),
                        status_name=StaticEscortFlowGurobiSolver._status_name,
                        solve_start=time.perf_counter(), mode=mode, target=10 if mode == "flow_certificate" else None)
                    self.assertFalse(result["weighted_proven"])
                    self.assertFalse(result["flow_proven"])

    def test_helper_optimizes_integer_scale_and_preserves_normalized_units(self):
        model = ScriptedModel(bound=110061.5, status=GRB.TIME_LIMIT)
        values = dict(flow=10, movements=12)
        result = solve_weighted_or_certificate(
            model, Expression(values, {"flow": 1}), Expression(values, {"movements": 1}),
            lambda: dict(makespan=5, animation_moves=[]),
            status_name=StaticEscortFlowGurobiSolver._status_name,
            solve_start=time.perf_counter(), mode="weighted_integer", weight_scale=11005,
        )
        self.assertEqual(model.expression.coefficients, {"flow": 11005, "movements": 1})
        self.assertEqual(result["scaled_objective"], 110062)
        self.assertEqual(result["scaled_best_bound"], 110061.5)
        self.assertEqual(result["scaled_absolute_gap"], .5)
        self.assertEqual(result["objective"], 110062 / 11005)
        self.assertEqual(result["best_bound"], 110061.5 / 11005)
        self.assertEqual(result["absolute_gap"], .5 / 11005)
        self.assertTrue(result["weighted_proven"])
        self.assertTrue(model.disposed)

    def test_scale_must_be_a_positive_integer(self):
        for scale in (0, -1, 1.1, 100.0, math.nan, math.inf, True, "100"):
            with self.subTest(scale=scale), self.assertRaises(ValueError):
                model = ScriptedModel()
                solve_weighted_or_certificate(
                    model, None, None, None, status_name=lambda _: "OPTIMAL",
                    solve_start=time.perf_counter(), mode="weighted_integer", weight_scale=scale)
            self.assertTrue(model.disposed)

    def test_both_backends_use_dynamic_coefficients_and_find_the_lex_optimum(self):
        for name, solver_class, config_class in (
            ("escort", StaticEscortFlowGurobiSolver, StaticGurobiConfig),
            ("load", LoadFlowStaticGurobiSolver, LoadFlowStaticGurobiConfig),
        ):
            with self.subTest(backend=name), contextlib.redirect_stdout(io.StringIO()):
                config = dict(Lx=3, Ly=2, output_cells=((0, 0),), beta=1,
                              gamma=1 / 11005, weight_scale=11005, time_limit=15,
                              threads=1, objective_mode="weighted_integer")
                config.update(dict(alpha=0, move_method="BM") if name == "load"
                              else dict(retrieval_mode="leave"))
                solver = solver_class(config_class(**config))
                try:
                    result = solver.solve({(1, 1), (2, 0)}, {(0, 0), (2, 1)}, 8)
                finally:
                    solver.close()
                self.assertTrue(result["weighted_proven"])
                self.assertEqual((result["flowtime"], result["movements"]), (8, 9))
                self.assertEqual(result["scaled_objective"], 11005 * 8 + 9)
                self.assertLess(result["scaled_absolute_gap"], 1)


if __name__ == "__main__":
    unittest.main()
