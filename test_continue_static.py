"""Validate continue-mode starts and optima against a physical SBM search."""
from contextlib import redirect_stdout
import heapq
import io
import itertools
import unittest
from unittest.mock import patch

import gurobipy as gp

from OneStepHeuristic_v2 import SolveGreedy
from escort_flow_static_gurobi import StaticEscortFlowGurobiSolver, StaticGurobiConfig
from load_flow_static_gurobi import LoadFlowStaticGurobiSolver, LoadFlowStaticGurobiConfig
from static_weighted_certification import safe_weight_parameters


def physical_optimum(width, height, outputs, targets, escorts):
    """Dijkstra on actual disjoint block shifts, independent of either ILP."""
    outputs = frozenset(outputs)
    initial = (frozenset(targets) - outputs, frozenset(escorts))
    best = {initial: (0, 0)}
    queue, serial = [(0, 0, 0, initial)], itertools.count(1)
    while queue:
        flow, movements, _, state = heapq.heappop(queue)
        if best[state] != (flow, movements):
            continue
        targets, escorts = state
        if not targets:
            return flow, movements
        options = []
        for origin in sorted(escorts):
            moves = [(frozenset(), (), origin)]
            for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                current, path, shifts = origin, {origin}, []
                while True:
                    source = (current[0] + dx, current[1] + dy)
                    if not (0 <= source[0] < width and 0 <= source[1] < height) or source in escorts:
                        break
                    path.add(source); shifts.append((source, current)); current = source
                    moves.append((frozenset(path), tuple(shifts), source))
            options.append(moves)
        for combination in itertools.product(*options):
            used, shifts, next_escorts = set(), {}, set()
            for cells, load_shifts, destination in combination:
                if used & cells:
                    break
                used.update(cells); shifts.update(load_shifts); next_escorts.add(destination)
            else:
                if not shifts:
                    continue
                next_targets = frozenset(shifts.get(loc, loc) for loc in targets) - outputs
                next_state = (next_targets, frozenset(next_escorts))
                cost = (flow + len(targets), movements + len(shifts))
                if cost < best.get(next_state, (float("inf"), float("inf"))):
                    best[next_state] = cost
                    heapq.heappush(queue, (*cost, next(serial), next_state))
    raise AssertionError("No physical retrieval schedule")


class ContinueTests(unittest.TestCase):
    def solver(self, method, width, height, outputs, weight):
        common = dict(Lx=width, Ly=height, output_cells=tuple(outputs), retrieval_mode="continue",
                      beta=1, gamma=1 / weight, time_limit=20, threads=1,
                      weight_scale=weight, objective_mode="weighted_integer")
        if method == "escortflow":
            solver = StaticEscortFlowGurobiSolver(StaticGurobiConfig(**common))
        else:
            solver = LoadFlowStaticGurobiSolver(LoadFlowStaticGurobiConfig(move_method="BM", alpha=0, **common))
        self.addCleanup(solver.close)
        return solver

    def check_case(self, width, height, outputs, targets, escorts, oracle=True):
        trace = SolveGreedy(width, height, set(outputs), targets, escorts,
                            retrieval_mode="continue", return_trace=True)
        makespan, flow, movements, _, escort_moves, target_moves = trace
        active = set(targets) - set(outputs)
        expected_flow = 0
        for step in target_moves:
            expected_flow += len(active)
            by_source = dict(step.values())
            active = {by_source.get(loc, loc) for loc in active} - set(outputs)
        self.assertFalse(active)
        self.assertEqual(expected_flow, flow)
        parameters = safe_weight_parameters(targets, outputs, width * height, len(escorts), flow)
        weight = parameters["flow_weight"]
        physical_horizon = max(parameters["flow_horizon"], makespan + 1)
        solved = []
        for method in ["escortflow", "loadflow"]:
            with self.subTest(method=method):
                solver = self.solver(method, width, height, outputs, weight)
                horizon = physical_horizon - (method == "escortflow")
                start = solver.build_warmstart_from_trace(targets, escorts, horizon, target_moves, escort_moves)
                original = gp.Model.optimize
                def fix_start(model, *args, **kwargs):
                    model.update()
                    for variable in model.getVars():
                        self.assertNotEqual(variable.Start, gp.GRB.UNDEFINED)
                        variable.LB = variable.Start; variable.UB = variable.Start
                    return original(model, *args, **kwargs)
                with patch.object(gp.Model, "optimize", fix_start), redirect_stdout(io.StringIO()):
                    fixed = solver.solve(targets, escorts, horizon, warmstart=start)
                self.assertEqual(fixed["status_name"], "OPTIMAL")
                self.assertEqual((fixed["flowtime"], fixed["movements"]), (flow, movements))
                self.assertFalse(any(destination == (None, None)
                                     for step in fixed["animation_moves"] for _, destination in step))
                log = io.StringIO()
                with redirect_stdout(log):
                    result = solver.solve(targets, escorts, horizon, warmstart=start)
                self.assertIn("Loaded user MIP start", log.getvalue())
                self.assertTrue(result["weighted_proven"])
                solved.append((result["flowtime"], result["movements"]))
        self.assertEqual(solved[0], solved[1])
        if oracle:
            self.assertEqual(solved[0], physical_optimum(width, height, outputs, targets, escorts))

    def test_two_targets_single_output(self):
        self.check_case(3, 2, {(0, 0)}, {(1, 0), (2, 0)}, {(0, 0), (0, 1)})

    def test_initial_target_at_output_remains_a_movable_blocker(self):
        self.check_case(3, 2, {(0, 0)}, {(0, 0), (2, 0)}, {(1, 0), (0, 1)})

    def test_consecutive_arrivals_can_shift_retrieved_blocker(self):
        self.check_case(3, 2, {(1, 0)}, {(0, 0), (2, 0)}, {(1, 0), (1, 1)})

    def test_four_targets(self):
        self.check_case(4, 2, {(0, 0), (3, 0)}, {(1, 0), (2, 0), (1, 1), (2, 1)}, {(0, 0), (3, 0)})

    def test_six_targets(self):
        self.check_case(4, 3, {(0, 0), (3, 0)}, {(1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (3, 1)},
                        {(0, 0), (3, 0), (0, 2), (3, 2)}, oracle=False)

    def test_all_targets_already_at_outputs(self):
        self.check_case(3, 2, {(0, 0)}, {(0, 0)}, {(2, 1)})


if __name__ == "__main__":
    unittest.main()
