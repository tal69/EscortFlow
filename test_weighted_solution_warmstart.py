"""Verify actual weighted incumbents remain complete, feasible MIP starts.

Fixing every variable to its supplied Start is intentional: merely optimizing
with a MIP start would let Gurobi repair or ignore an incorrect transfer.
"""

from copy import deepcopy
import contextlib
import io
import math
import unittest
from unittest.mock import patch

import gurobipy as gp
from gurobipy import GRB

from escort_flow_static_gurobi import StaticEscortFlowGurobiSolver, StaticGurobiConfig
from load_flow_static_gurobi import LoadFlowStaticGurobiSolver, LoadFlowStaticGurobiConfig
from static_weighted_certification import sufficient_flow_horizon


class WeightedSolutionWarmstartTests(unittest.TestCase):
    BACKENDS = (
        ("escort", StaticEscortFlowGurobiSolver, StaticGurobiConfig),
        ("load", LoadFlowStaticGurobiSolver, LoadFlowStaticGurobiConfig),
    )

    def make_solver(self, backend, *, mode="weighted_integer", target=None,
                    Lx=3, Ly=2):
        name, solver_type, config_type = backend
        config = dict(Lx=Lx, Ly=Ly, output_cells=((0, 0),), beta=1,
                      gamma=.01, time_limit=20, threads=1,
                      objective_mode=mode, certification_target=target)
        config.update(dict(retrieval_mode="leave") if name == "escort"
                      else dict(alpha=0, move_method="BM"))
        solver = solver_type(config_type(**config))
        self.addCleanup(solver.close)
        return solver

    def solve_fixed_start(self, solver, targets, escorts, horizon, start):
        original_optimize = gp.Model.optimize
        checked = []

        def optimize_fixed(model, *args, **kwargs):
            model.update()
            variables = model.getVars()
            self.assertTrue(variables)
            for variable in variables:
                value = variable.Start
                self.assertTrue(math.isfinite(value))
                self.assertNotEqual(value, GRB.UNDEFINED,
                                    "Every variable must receive an explicit MIP start")
                variable.LB = value
                variable.UB = value
            # Require an actual feasible fixed solution, rather than stopping
            # solely because a bound happens to meet the certification target.
            model.Params.BestBdStop = GRB.INFINITY
            model.Params.BestObjStop = -GRB.INFINITY
            checked.append(len(variables))
            return original_optimize(model, *args, **kwargs)

        with patch.object(gp.Model, "optimize", optimize_fixed):
            result = solver.solve(targets, escorts, horizon, warmstart=start)
        self.assertEqual(len(checked), 1)
        self.assertTrue(result["has_solution"], result)
        self.assertEqual(result["status_name"], "OPTIMAL")
        return result

    def test_real_weighted_solution_is_feasible_at_shorter_same_and_larger_horizons(self):
        targets, escorts = {(1, 1), (2, 0)}, {(0, 0), (2, 1)}
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]), contextlib.redirect_stdout(io.StringIO()):
                weighted_solver = self.make_solver(backend)
                weighted = weighted_solver.solve(targets, escorts, 8)
                self.assertTrue(weighted["weighted_proven"])
                self.assertEqual(weighted["solution_horizon"], 8)
                original = deepcopy(weighted)
                certificate_solver = self.make_solver(
                    backend, mode="flow_certificate", target=weighted["flowtime"])
                safe_horizon = sufficient_flow_horizon(targets, {(0, 0)}, weighted["flowtime"])
                self.assertLess(safe_horizon, weighted["solution_horizon"])
                for horizon in (safe_horizon, 8, 11):
                    start = certificate_solver.build_warmstart_from_solution(weighted, horizon)
                    fixed = self.solve_fixed_start(certificate_solver, targets, escorts,
                                                   horizon, start)
                    self.assertEqual(fixed["flowtime"], weighted["flowtime"])
                    self.assertEqual(fixed["movements"], weighted["movements"])
                    self.assertEqual(start["removed_post_retrieval_movements"], 0)
                    self.assertTrue(fixed["flow_proven"])
                self.assertEqual(weighted, original, "Transfer must not mutate the weighted result")

    def test_initial_output_retrieval_with_zero_horizon_extends_for_both_models(self):
        targets, escorts = {(0, 0)}, {(2, 0)}
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]), contextlib.redirect_stdout(io.StringIO()):
                solver = self.make_solver(backend, Lx=3, Ly=1)
                weighted = solver.solve(targets, escorts, 0)
                self.assertEqual((weighted["flowtime"], weighted["movements"]), (0, 0))
                start = solver.build_warmstart_from_solution(weighted, 3)
                fixed = self.solve_fixed_start(solver, targets, escorts, 3, start)
                self.assertEqual((fixed["flowtime"], fixed["movements"]), (0, 0))
                if backend[0] == "escort":
                    self.assertEqual(start["x_a"][((0, 0, 0, 0), 0)], 1)
                    for t in range(1, 4):
                        self.assertEqual(start["x_e"][((0, 0, 0, 0), t)], 1)

    def test_last_period_arrival_is_completed_before_escort_padding(self):
        targets, escorts = {(1, 0)}, {(0, 0)}
        for backend in self.BACKENDS:
            horizon = 0 if backend[0] == "escort" else 1
            with self.subTest(backend=backend[0]), contextlib.redirect_stdout(io.StringIO()):
                solver = self.make_solver(backend, Lx=3, Ly=1)
                weighted = solver.solve(targets, escorts, horizon)
                self.assertEqual((weighted["flowtime"], weighted["movements"]), (1, 1))
                start = solver.build_warmstart_from_solution(weighted, horizon + 3)
                fixed = self.solve_fixed_start(solver, targets, escorts, horizon + 3, start)
                self.assertEqual((fixed["flowtime"], fixed["movements"]), (1, 1))
                if backend[0] == "escort":
                    self.assertEqual(start["x_a"][((0, 0, 0, 0), 1)], 1)
                    self.assertEqual(start["x_e"].get(((0, 0, 0, 0), 1), 0), 0)
                    for t in (2, 3):
                        self.assertEqual(start["x_e"][((0, 0, 0, 0), t)], 1)

    def test_load_terminal_swap_is_removed_only_when_horizon_grows(self):
        # The original LF model imposes no opposite-direction constraint on
        # its last period. This feasible but suboptimal incumbent swaps two
        # blocking loads after the target has already been retrieved at t=0.
        # Copying those arcs into a longer model would make the start invalid.
        with contextlib.redirect_stdout(io.StringIO()):
            solver = self.make_solver(self.BACKENDS[1], Lx=3, Ly=1)
            targets, escorts = {(0, 0)}, set()
            terminal_swap = {"x": {((1, 0, 2, 0), 0, 2): 1.0,
                                     ((2, 0, 1, 0), 0, 2): 1.0},
                             "q": {((0, 0), 0): 1.0}, "z": 0.0}
            weighted = self.solve_fixed_start(solver, targets, escorts, 0, terminal_swap)
            self.assertEqual((weighted["flowtime"], weighted["movements"]), (0, 2))
            original = deepcopy(weighted)
            same_horizon = solver.build_warmstart_from_solution(weighted, 0)
            unchanged = self.solve_fixed_start(solver, targets, escorts, 0, same_horizon)
            self.assertEqual(unchanged["movements"], 2)
            self.assertEqual(same_horizon["removed_post_retrieval_movements"], 0)
            start = solver.build_warmstart_from_solution(weighted, 3)
            fixed = self.solve_fixed_start(solver, targets, escorts, 3, start)
            self.assertEqual((fixed["flowtime"], fixed["movements"]), (0, 0))
            self.assertEqual(start["removed_post_retrieval_movements"], 2)
            self.assertEqual(weighted, original)

    def test_horizon_before_last_retrieval_cannot_truncate_weighted_schedule(self):
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]), contextlib.redirect_stdout(io.StringIO()):
                solver = self.make_solver(backend)
                weighted = solver.solve({(1, 1), (2, 0)}, {(0, 0), (2, 1)}, 8)
                minimum_horizon = weighted["makespan"] - (backend[0] == "escort")
                self.assertGreater(minimum_horizon, 0)
                start = solver.build_warmstart_from_solution(weighted, minimum_horizon)
                fixed = self.solve_fixed_start(solver, {(1, 1), (2, 0)}, {(0, 0), (2, 1)},
                                               minimum_horizon, start)
                self.assertEqual(fixed["flowtime"], weighted["flowtime"])
                with self.assertRaises(ValueError):
                    solver.build_warmstart_from_solution(weighted, minimum_horizon - 1)

    def test_load_terminal_destination_collision_is_removed_on_extension(self):
        # At t=T, one horizontal and one vertical blocker can enter the same
        # terminal cell because its destination-capacity row is outside the
        # original horizon. The next period must begin from their distinct
        # origins instead, after eliminating these redundant final moves.
        with contextlib.redirect_stdout(io.StringIO()):
            solver = self.make_solver(self.BACKENDS[1], Lx=2, Ly=2)
            targets, escorts = {(0, 0)}, {(1, 1)}
            terminal_collision = {
                "x": {((1, 0, 1, 1), 0, 2): 1.0,
                      ((0, 1, 1, 1), 0, 2): 1.0},
                "q": {((0, 0), 0): 1.0}, "z": 0.0,
            }
            weighted = self.solve_fixed_start(solver, targets, escorts, 0,
                                               terminal_collision)
            self.assertEqual((weighted["flowtime"], weighted["movements"]), (0, 2))
            start = solver.build_warmstart_from_solution(weighted, 3)
            fixed = self.solve_fixed_start(solver, targets, escorts, 3, start)
            self.assertEqual((fixed["flowtime"], fixed["movements"]), (0, 0))
            self.assertEqual(start["removed_post_retrieval_movements"], 2)

    def test_shortening_discards_only_moves_after_retrieval(self):
        targets, escorts = {(0, 0)}, {(2, 0)}
        source_starts = {
            "escort": {
                "x_a": {((0, 0, 0, 0), 0): 1.0},
                "x_e": {((2, 0, 2, 0), 0): 1.0,
                        ((0, 0, 0, 0), 1): 1.0, ((2, 0, 1, 0), 1): 1.0,
                        ((0, 0, 0, 0), 2): 1.0, ((1, 0, 2, 0), 2): 1.0,
                        ((0, 0, 0, 0), 3): 1.0, ((2, 0, 2, 0), 3): 1.0},
                "q": {(0, 0): 0.0},
            },
            "load": {
                "x": {((1, 0, 1, 0), 0, 2): 1.0,
                      ((1, 0, 2, 0), 1, 2): 1.0,
                      ((2, 0, 1, 0), 2, 2): 1.0,
                      ((1, 0, 1, 0), 3, 2): 1.0},
                "q": {((0, 0), 0): 1.0}, "z": 0.0,
            },
        }
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]), contextlib.redirect_stdout(io.StringIO()):
                solver = self.make_solver(backend, Lx=3, Ly=1)
                weighted = self.solve_fixed_start(solver, targets, escorts, 3,
                                                   source_starts[backend[0]])
                self.assertEqual((weighted["flowtime"], weighted["movements"]), (0, 2))
                original = deepcopy(weighted)
                safe_horizon = sufficient_flow_horizon(targets, {(0, 0)}, weighted["flowtime"])
                self.assertEqual(safe_horizon, 0)
                start = solver.build_warmstart_from_solution(weighted, safe_horizon)
                fixed = self.solve_fixed_start(solver, targets, escorts, safe_horizon, start)
                self.assertEqual((fixed["flowtime"], fixed["movements"]), (0, 0))
                self.assertEqual(start["removed_post_retrieval_movements"], 2)
                self.assertEqual(weighted, original)


if __name__ == "__main__":
    unittest.main()
