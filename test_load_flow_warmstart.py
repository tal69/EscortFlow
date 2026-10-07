"""Validate complete load-flow starts against the actual Gurobi constraints.

Run with: python3 -m unittest test_load_flow_warmstart.py
"""

import contextlib
import io
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import gurobipy as gp
import numpy as np

from OneStepHeuristic_v2 import SolveGreedy
from escort_flow_static_gurobi import StaticEscortFlowGurobiSolver, StaticGurobiConfig
from load_flow_static_gurobi import LoadFlowStaticGurobiConfig, LoadFlowStaticGurobiSolver


class LoadFlowWarmStartTests(unittest.TestCase):
    def new_solver(self, width=3, height=2, outputs=((0, 0),), **overrides):
        config = dict(Lx=width, Ly=height, output_cells=outputs, move_method="BM",
                      alpha=0.25, beta=1.0, gamma=0.01, time_limit=15, threads=1)
        config.update(overrides)
        solver = LoadFlowStaticGurobiSolver(LoadFlowStaticGurobiConfig(**config))
        self.addCleanup(solver.close)
        return solver

    def assert_trace_feasible(self, targets, escorts, *, width=3, height=2,
                              outputs=((0, 0),)):
        trace = SolveGreedy(width, height, set(outputs), targets, escorts,
                            retrieval_mode="leave", return_trace=True)
        makespan, flowtime, movements, _, escort_moves, target_moves = trace
        solver = self.new_solver(width, height, outputs)
        # Use the minimum LF horizon: the final q event occurs at makespan.
        start = solver.build_warmstart_from_trace(targets, escorts, makespan,
                                                  target_moves, escort_moves)
        original_optimize = gp.Model.optimize
        observed_start = []

        def fix_complete_start(model, *args, **kwargs):
            model.update()
            for var in model.getVars():
                value = var.Start
                self.assertNotEqual(value, gp.GRB.UNDEFINED)
                self.assertIn(value, (0.0, 1.0, float(makespan)))
                var.LB = value
                var.UB = value
                observed_start.append(value)
            return original_optimize(model, *args, **kwargs)

        with patch.object(gp.Model, "optimize", fix_complete_start):
            result = solver.solve(targets, escorts, makespan, warmstart=start)
        self.assertTrue(observed_start)
        self.assertIn(0.0, observed_start)
        self.assertTrue(result["has_solution"])
        self.assertEqual(result["status_name"], "OPTIMAL")
        self.assertAlmostEqual(result["makespan"], makespan)
        self.assertAlmostEqual(result["flowtime"], flowtime)
        self.assertAlmostEqual(result["movements"], movements)
        self.assertAlmostEqual(result["objective"], .25 * makespan + flowtime + .01 * movements)

        # Both formulations receive this exact trace and therefore the same
        # target arrivals and the same number of physical load shifts.
        escort_solver = StaticEscortFlowGurobiSolver(StaticGurobiConfig(
            Lx=width, Ly=height, output_cells=outputs, retrieval_mode="leave",
            beta=1.0, gamma=.01, time_limit=15, threads=1,
        ))
        self.addCleanup(escort_solver.close)
        escort_start = escort_solver.build_warmstart_from_trace(
            targets, escorts, max(0, makespan - 1), target_moves, escort_moves,
        )
        summary = escort_solver.summarize_warmstart(escort_start)
        self.assertAlmostEqual(summary["flowtime"], flowtime)
        self.assertAlmostEqual(summary["movements"], movements)

    def test_adjacent_targets_share_output_with_service_delay(self):
        self.assert_trace_feasible({(1, 0), (2, 0)}, {(0, 0), (0, 1)})

    def test_initial_target_at_output_and_blocking_load_movements(self):
        self.assert_trace_feasible({(0, 0), (1, 0)}, {(2, 1)})

    def test_multiple_outputs_and_simultaneous_retrievals(self):
        self.assert_trace_feasible({(1, 0), (2, 0), (2, 1)},
                                   {(0, 0), (3, 0), (0, 1)}, width=4,
                                   outputs=((0, 0), (3, 0)))

    def test_all_targets_initially_at_outputs(self):
        self.assert_trace_feasible({(0, 0)}, {(2, 1)})

    def test_no_targets(self):
        self.assert_trace_feasible(set(), {(2, 1)})

    def test_numpy_coordinates_from_generated_instances(self):
        targets = {(np.int64(1), np.int64(1)), (np.int64(2), np.int64(0))}
        escorts = {(np.int64(0), np.int64(0)), (np.int64(2), np.int64(1))}
        self.assert_trace_feasible(targets, escorts)

    def test_gurobi_accepts_start_in_weighted_and_lexicographic_solves(self):
        targets, escorts = {(1, 0), (2, 0)}, {(0, 0), (0, 1)}
        trace = SolveGreedy(3, 2, {(0, 0)}, targets, escorts,
                            retrieval_mode="leave", return_trace=True)
        for lexicographic in (False, True):
            with self.subTest(lexicographic=lexicographic):
                solver = self.new_solver(alpha=0.0, lexicographic=lexicographic)
                start = solver.build_warmstart_from_trace(targets, escorts, 4, trace[5], trace[4])
                log = io.StringIO()
                with contextlib.redirect_stdout(log):
                    result = solver.solve(targets, escorts, 4, warmstart=start)
                self.assertIn("Loaded user MIP start", log.getvalue())
                self.assertTrue(result["has_solution"])

    def test_short_horizon_and_inconsistent_trace_fail_explicitly(self):
        solver = self.new_solver()
        targets, escorts = {(1, 0)}, {(0, 0)}
        trace = SolveGreedy(3, 2, {(0, 0)}, targets, escorts,
                            retrieval_mode="leave", return_trace=True)
        with self.assertRaisesRegex(ValueError, "horizon at least"):
            solver.build_warmstart_from_trace(targets, escorts, 0, trace[5], trace[4])
        with self.assertRaisesRegex(ValueError, "traces disagree"):
            solver.build_warmstart_from_trace(targets, escorts, 1, [{}], trace[4])

    def test_unsupported_cli_combinations_fail_before_solving(self):
        script = Path(__file__).with_name("LoadFlowStatic.py")
        for extra in (["--lm"], ["--lp"], ["--opl"],
                      ["--retrieval_mode", "continue"], ["--dp_file", "missing.p"]):
            with self.subTest(extra=extra):
                result = subprocess.run(
                    [sys.executable, str(script), "-x", "3", "-y", "2", "-O", "0", "0",
                     "--warmstart", *extra], capture_output=True, text=True,
                )
                self.assertEqual(result.returncode, 2)
                self.assertIn("--warmstart", result.stderr)


if __name__ == "__main__":
    unittest.main()
