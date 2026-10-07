"""Regression checks for the two-phase static retrieval solvers.

Run with a Python interpreter that has gurobipy and a working Gurobi license::

    python3 test_static_lexicographic.py

The small leave-mode optima below were independently obtained by exhaustive
search over vertex-disjoint straight escort movements, including the output's
one-takt service delay. Scripted model runs make time/work-limit checks
deterministic without waiting for a real solver to reach a limit.
"""

import contextlib
import io
from types import SimpleNamespace
import time
import unittest

from gurobipy import GRB

from escort_flow_static_gurobi import StaticEscortFlowGurobiSolver, StaticGurobiConfig
from escort_flow_static_lazy import (
    LazyStaticEscortFlowGurobiSolver,
    StaticGurobiConfig as LazyConfig,
)
from escort_flow_static_bnc import (
    BnCStaticEscortFlowGurobiSolver,
    StaticGurobiConfig as BnCConfig,
)
from load_flow_static_gurobi import LoadFlowStaticGurobiConfig, LoadFlowStaticGurobiSolver
from static_lexicographic import solve_lexicographic


class RetrievalIntegrationTests(unittest.TestCase):
    BACKENDS = (
        ("escort", StaticEscortFlowGurobiSolver, StaticGurobiConfig),
        ("lazy", LazyStaticEscortFlowGurobiSolver, LazyConfig),
        ("bnc", BnCStaticEscortFlowGurobiSolver, BnCConfig),
        ("load", LoadFlowStaticGurobiSolver, LoadFlowStaticGurobiConfig),
    )

    def solve_case(self, backend, targets, escorts, *, width=3, height=2, horizon=8):
        name, solver_class, config_class = backend
        kwargs = dict(
            Lx=width, Ly=height, output_cells=((0, 0),),
            beta=1.0, gamma=0.01, time_limit=20, threads=1,
            lexicographic=True, phase1_time_limit=10,
        )
        if name == "load":
            kwargs.update(move_method="BM", alpha=0.0)
        else:
            kwargs.update(retrieval_mode="leave")
        with contextlib.redirect_stdout(io.StringIO()):
            solver = solver_class(config_class(**kwargs))
            try:
                return solver.solve(set(targets), set(escorts), horizon)
            finally:
                solver.close()

    def assert_optimum(self, result, flowtime, movements):
        self.assertTrue(result["has_solution"])
        self.assertTrue(result["lexicographic_optimal"])
        self.assertAlmostEqual(result["flowtime"], flowtime)
        self.assertAlmostEqual(result["movements"], movements)
        for key, objective in (("phase1", flowtime), ("phase2", movements)):
            phase = result[key]
            self.assertTrue(phase["proven_optimal"])
            self.assertAlmostEqual(phase["objective"], objective)
            self.assertLess(phase["absolute_gap"], 1.0)

    def test_flow_priority_and_movement_tie_break(self):
        # Minimizing movements first gives (flowtime, movements) = (10, 7).
        # The required flow-first optimum instead has objective pair (8, 9).
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]):
                result = self.solve_case(backend, {(1, 1), (2, 0)}, {(0, 0), (2, 1)})
                self.assert_optimum(result, 8, 9)

    def test_output_service_delay_and_first_transition(self):
        # Arrivals at times 1 and 3 must both contribute to total flow time.
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]):
                result = self.solve_case(backend, {(1, 0), (2, 0)}, {(0, 0), (0, 1)})
                self.assert_optimum(result, 4, 3)

    def test_single_transition_uses_fewest_load_shifts(self):
        # Moving the escort two cells also retrieves the target at time 1,
        # but unnecessarily shifts a blocking load. One cell is optimal.
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]):
                horizon = 1 if backend[0] == "load" else 0
                result = self.solve_case(
                    backend, {(1, 0)}, {(0, 0)}, height=1, horizon=horizon,
                )
                self.assert_optimum(result, 1, 1)

    def test_empty_and_already_retrieved_targets(self):
        for backend in self.BACKENDS:
            for targets in (set(), {(0, 0)}):
                with self.subTest(backend=backend[0], targets=targets):
                    result = self.solve_case(backend, targets, {(2, 1)}, horizon=2)
                    self.assert_optimum(result, 0, 0)

    def test_infeasible_phase_one_does_not_run_phase_two(self):
        for backend in self.BACKENDS:
            with self.subTest(backend=backend[0]):
                result = self.solve_case(backend, {(1, 0)}, set(), horizon=2)
                self.assertFalse(result["has_solution"])
                self.assertFalse(result["lexicographic_optimal"])
                self.assertEqual(result["phase1"]["status_name"], "INFEASIBLE")
                self.assertEqual(result["phase2"]["status_name"], "NOT_RUN")
                self.assertTrue(result["phase2_skip_reason"])


class ScriptedExpression:
    def __init__(self, model, name):
        self.model = model
        self.name = name

    def getValue(self):
        return self.model.current[self.name]

    def __eq__(self, value):
        return self.name, "==", value


class ScriptedModel:
    """Small Gurobi protocol double: each optimize call consumes one run."""

    def __init__(self, runs):
        self.runs = list(runs)
        self.current = {}
        self.Params = SimpleNamespace(
            TimeLimit=GRB.INFINITY, WorkLimit=GRB.INFINITY,
            MIPGap=1e-4, MIPGapAbs=1e-10, Cutoff=GRB.INFINITY,
        )
        self.calls = []
        self.constraints = []
        self.variables = [SimpleNamespace(X=1.0, Start=GRB.UNDEFINED)]
        self.disposed = False
        self.SolCount = 0
        self.Status = GRB.LOADED

    def setObjective(self, expression, sense):
        self.expression = expression

    def addConstr(self, expression, name=""):
        self.constraints.append(expression)

    def update(self):
        pass

    def getVars(self):
        return self.variables

    def getAttr(self, attribute, variables=None):
        if variables is None:
            return getattr(self, attribute)
        return [getattr(variable, attribute) for variable in variables]

    def setAttr(self, attribute, variables, values):
        for variable, value in zip(variables, values):
            setattr(variable, attribute, value)

    def optimize(self, callback=None):
        self.calls.append(dict(
            objective=self.expression.name, time_limit=self.Params.TimeLimit,
            work_limit=self.Params.WorkLimit, mip_gap=self.Params.MIPGap,
            absolute_gap=self.Params.MIPGapAbs, callback=callback,
        ))
        self.current = self.runs[len(self.calls) - 1]
        self.Status = self.current.get("status", GRB.OPTIMAL)
        self.SolCount = self.current.get("solutions", 1)
        self.Runtime = self.current.get("runtime", 0.5)
        self.Work = self.current.get("work", 0.1)
        self.ObjBound = self.current.get("bound", self.current[self.expression.name])
        self.ObjVal = self.current[self.expression.name]

    def dispose(self):
        self.disposed = True


def scripted_run(*, flowtime=10, movements=12, status=GRB.OPTIMAL,
                 bound=None, runtime=0.5, work=0.1, solutions=1):
    run = dict(flowtime=flowtime, movements=movements, status=status,
               runtime=runtime, work=work, solutions=solutions)
    if bound is not None:
        run["bound"] = bound
    return run


class PhaseBudgetTests(unittest.TestCase):
    def solve_runs(self, runs, **limits):
        model = ScriptedModel(runs)

        def extract_solution():
            return dict(
                makespan=3, flowtime=model.current["flowtime"],
                movements=model.current["movements"], animation_moves=[[]],
            )

        result = solve_lexicographic(
            model, ScriptedExpression(model, "flowtime"),
            ScriptedExpression(model, "movements"), extract_solution,
            status_name=StaticEscortFlowGurobiSolver._status_name,
            solve_start=time.perf_counter(), **limits,
        )
        self.assertTrue(model.disposed)
        for call in model.calls:
            self.assertEqual(call["mip_gap"], 0)
            self.assertEqual(call["absolute_gap"], 1)
        return result, model

    def test_phase_one_limit_continues_at_best_attained_flow(self):
        result, model = self.solve_runs([
            scripted_run(status=GRB.TIME_LIMIT, bound=8, runtime=3, work=2),
            scripted_run(movements=4, runtime=5, work=3),
        ], time_limit=10, phase1_time_limit=3, work_limit=10)
        self.assertEqual([c["time_limit"] for c in model.calls], [3, 7])
        self.assertEqual([c["work_limit"] for c in model.calls], [10, 8])
        self.assertIn(("flowtime", "==", 10), model.constraints)
        self.assertEqual(result["movements"], 4)
        self.assertFalse(result["phase1"]["proven_optimal"])
        self.assertTrue(result["phase2"]["proven_optimal"])
        self.assertFalse(result["lexicographic_optimal"])

    def test_unused_phase_one_time_is_available_to_phase_two(self):
        result, model = self.solve_runs([
            scripted_run(runtime=0.5), scripted_run(movements=4),
        ], time_limit=10, phase1_time_limit=3)
        self.assertEqual([c["time_limit"] for c in model.calls], [3, 9.5])
        self.assertTrue(result["lexicographic_optimal"])

    def test_total_budget_caps_phase_one_and_prevents_phase_two(self):
        result, model = self.solve_runs([
            scripted_run(status=GRB.TIME_LIMIT, runtime=2),
        ], time_limit=2, phase1_time_limit=3)
        self.assertEqual(len(model.calls), 1)
        self.assertEqual(model.calls[0]["time_limit"], 2)
        self.assertEqual(result["phase2"]["status_name"], "NOT_RUN")
        self.assertEqual(result["movements"], 12)
        self.assertTrue(result["phase2_skip_reason"])

    def test_shared_work_exhaustion_preserves_phase_one_solution(self):
        result, model = self.solve_runs([
            scripted_run(status=GRB.WORK_LIMIT, work=2),
        ], time_limit=10, phase1_time_limit=3, work_limit=2)
        self.assertEqual(len(model.calls), 1)
        self.assertEqual(result["phase2"]["status_name"], "NOT_RUN")
        self.assertTrue(result["has_solution"])
        self.assertEqual(result["movements"], 12)

    def test_no_incumbent_skips_phase_two(self):
        result, model = self.solve_runs([
            scripted_run(status=GRB.TIME_LIMIT, solutions=0, bound=8),
        ], time_limit=10, phase1_time_limit=3)
        self.assertEqual(len(model.calls), 1)
        self.assertFalse(result["has_solution"])
        self.assertEqual(result["phase2"]["status_name"], "NOT_RUN")
        self.assertFalse(result["lexicographic_optimal"])

    def test_phase_two_without_incumbent_preserves_phase_one(self):
        result, model = self.solve_runs([
            scripted_run(),
            scripted_run(movements=0, status=GRB.TIME_LIMIT, solutions=0, bound=2),
        ], time_limit=10, phase1_time_limit=3)
        self.assertEqual(len(model.calls), 2)
        self.assertTrue(result["has_solution"])
        self.assertEqual(result["flowtime"], 10)
        self.assertEqual(result["movements"], 12)
        self.assertFalse(result["lexicographic_optimal"])

    def test_absolute_gap_boundary_is_strict(self):
        for bound, proven in ((9.0, False), (9.0001, True)):
            with self.subTest(bound=bound):
                result, _ = self.solve_runs([
                    scripted_run(status=GRB.TIME_LIMIT, bound=bound),
                    scripted_run(movements=4),
                ], time_limit=10, phase1_time_limit=3)
                self.assertEqual(result["phase1"]["proven_optimal"], proven)
                self.assertEqual(result["lexicographic_optimal"], proven)

    def test_phase_one_cap_without_total_limit(self):
        _, model = self.solve_runs([
            scripted_run(), scripted_run(movements=4),
        ], phase1_time_limit=3)
        self.assertEqual(model.calls[0]["time_limit"], 3)
        self.assertGreaterEqual(model.calls[1]["time_limit"], GRB.INFINITY)

    def test_numeric_phase_two_cannot_certify_an_optimum(self):
        # A reported zero gap does not override a numeric/suboptimal status.
        for status in (GRB.NUMERIC, GRB.SUBOPTIMAL):
            with self.subTest(status=status):
                result, _ = self.solve_runs([
                    scripted_run(),
                    scripted_run(movements=4, status=status, bound=4),
                ], time_limit=10, phase1_time_limit=3)
                self.assertTrue(result["has_solution"])
                self.assertFalse(result["phase2"]["proven_optimal"])
                self.assertFalse(result["lexicographic_optimal"])
                self.assertEqual(
                    result["status_name"], StaticEscortFlowGurobiSolver._status_name(status),
                )

    def test_interrupted_phase_one_does_not_restart_optimization(self):
        result, model = self.solve_runs([
            scripted_run(status=GRB.INTERRUPTED, bound=8),
        ], time_limit=10, phase1_time_limit=3)
        self.assertEqual(len(model.calls), 1)
        self.assertTrue(result["has_solution"])
        self.assertEqual(result["phase2"]["status_name"], "NOT_RUN")
        self.assertEqual(result["status_name"], "INTERRUPTED")
        self.assertFalse(result["lexicographic_optimal"])

    def test_separation_callback_and_incumbent_reach_both_phases(self):
        callback = lambda model, where: None
        _, model = self.solve_runs([
            scripted_run(), scripted_run(movements=4),
        ], time_limit=10, phase1_time_limit=3, callback=callback)
        self.assertEqual([call["callback"] for call in model.calls], [callback, callback])
        self.assertEqual(model.variables[0].Start, model.variables[0].X)


if __name__ == "__main__":
    unittest.main(verbosity=2)
