import csv
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import RunWeightedStatic as runner


class FakeSolver:
    def __init__(self, result):
        self.result = result
        self.starts = []
        self.solution_starts = []
        self.solves = []
        self.closed = False

    def build_warmstart_from_trace(self, targets, escorts, horizon, target_history, escort_history):
        start = (targets, escorts, horizon, target_history, escort_history)
        self.starts.append(start)
        return start

    def solve(self, targets, escorts, horizon, warmstart):
        self.solves.append((targets, escorts, horizon, warmstart))
        return self.result

    def build_warmstart_from_solution(self, result, horizon):
        start = dict(weighted_solution=result["solution_warmstart"], horizon=horizon,
                     removed_post_retrieval_movements=0)
        self.solution_starts.append((result, horizon, start))
        return start

    def close(self):
        self.closed = True


class WeightedRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.args = SimpleNamespace(
            formulation="escortflow", Lx=2, Ly=2, loads=1, outputs=[(0, 0)],
            threads=1, weighted_time_limit=300, certification_time_limit=71,
            certification_gap_threshold=0.1,
            output=Path(self.temporary.name) / "batch.csv",
        )
        self.weighted = dict(has_solution=True, status_name="TIME_LIMIT", cpu_time=301,
                             runtime=300, work=42, best_bound=4.13, objective=4.18,
                             absolute_gap=0.05, weighted_proven=False, flowtime=4,
                             movements=18, makespan=4, bound_consistent=True,
                             animation_moves=[[[1, 0, 0, 0]]],
                             solution_horizon=4, solution_warmstart={"saved": "weighted incumbent"})
        self.certificate = dict(has_solution=True, status_name="USER_OBJ_LIMIT", cpu_time=12,
                                runtime=11, work=1, best_bound=4, flowtime=4,
                                movements=23, makespan=4, flow_proven=True,
                                counterexample=False, animation_moves=[])

    def run_fake(self, weighted=None, certificate=None):
        solver = FakeSolver(dict(self.weighted if weighted is None else weighted))
        certificate_solver = FakeSolver(dict(self.certificate if certificate is None else certificate))
        calls = []

        def factory(args, **kwargs):
            calls.append(kwargs)
            return certificate_solver

        # The fake greedy trace makes the weighted physical horizon five.
        trace = (5, 5, 22, [], [{"escort": 1}], [{"target": 1}])
        with patch.object(runner.OneStepHeuristic_v2, "SolveGreedy", return_value=trace), redirect_stdout(io.StringIO()):
            row = runner.run_instance(self.args, 7, 1, solver, factory)
        return row, solver, certificate_solver, calls

    def test_independent_budget_and_original_candidate_preserved(self):
        row, solver, certificate_solver, calls = self.run_fake()
        self.assertEqual(calls, [dict(objective_mode="flow_certificate", time_limit=71, certification_target=4)])
        self.assertEqual(row["weighted_runtime"], 300)
        self.assertEqual(row["certification_runtime"], 11)
        self.assertEqual(row["movements"], 18)
        self.assertEqual(row["certification_best_movements"], 23)
        self.assertEqual(row["flow_proven"], 1)
        self.assertEqual(row["lexicographic_proven"], 0)
        self.assertEqual(row["scaled_objective"], 418)
        self.assertEqual(certificate_solver.starts, [])
        self.assertIs(certificate_solver.solution_starts[0][0], solver.result)
        self.assertEqual(certificate_solver.solves[0][3]["weighted_solution"],
                         self.weighted["solution_warmstart"])
        self.assertEqual(row["certification_warmstart_source"], "weighted_solution")
        self.assertTrue(certificate_solver.closed)

    def test_strict_gap_gate_and_no_incumbent_skip(self):
        for changes in (dict(absolute_gap=0.1), dict(absolute_gap=0.9),
                        dict(has_solution=False), dict(bound_consistent=False)):
            with self.subTest(changes=changes):
                result = dict(self.weighted, **changes)
                row, _, certificate_solver, calls = self.run_fake(weighted=result)
                self.assertEqual(row["certification_eligible"], 0)
                self.assertEqual(row["certification_status"], "NOT_RUN")
                self.assertEqual(calls, [])
                self.assertEqual(certificate_solver.starts, [])
                self.assertEqual(certificate_solver.solution_starts, [])

    def test_missing_saved_solution_reports_error_instead_of_using_greedy(self):
        weighted = dict(self.weighted)
        del weighted["solution_warmstart"]
        row, _, certificate_solver, _ = self.run_fake(weighted=weighted)
        self.assertIn("solution_warmstart", row["error"])
        self.assertEqual(certificate_solver.starts, [])
        self.assertEqual(certificate_solver.solves, [])
        self.assertTrue(certificate_solver.closed)

    def test_proven_solution_runs_and_horizon_scope_is_precise(self):
        weighted = dict(self.weighted, weighted_proven=True, absolute_gap=0)
        row, _, _, _ = self.run_fake(weighted=weighted)
        self.assertEqual(row["global_flow_horizon_requirement"], weighted["flowtime"])
        self.assertNotEqual(row["global_flow_horizon_requirement"], row["greedy_flowtime"])
        self.assertEqual(row["certification_horizon"], weighted["flowtime"])
        self.assertEqual(row["lexicographic_proven"], 1)
        self.assertEqual(row["lexicographic_proven_within_weighted_horizon"], 1)
        weighted.update(flowtime=10, objective=10.18)
        row, _, _, _ = self.run_fake(weighted=weighted)
        self.assertEqual(row["global_flow_horizon_requirement"], 10)
        self.assertEqual(row["lexicographic_proven"], 0)
        self.assertEqual(row["lexicographic_proven_within_weighted_horizon"], 1)
        self.assertGreater(row["certification_horizon"], row["weighted_horizon"])

    def test_certificate_horizon_shrinks_using_weighted_flow_not_greedy(self):
        weighted = dict(self.weighted, flowtime=2, makespan=2, objective=2.18)
        row, _, certificate_solver, _ = self.run_fake(weighted=weighted)
        self.assertEqual(row["weighted_horizon"], 4)
        self.assertEqual(row["certification_horizon"], 2)
        self.assertEqual(certificate_solver.solution_starts[0][1], 2)

    def test_counterexample_saves_both_plans_without_changing_weighted_result(self):
        certificate = dict(self.certificate, flowtime=3, movements=30,
                           flow_proven=False, counterexample=True,
                           animation_moves=[[[0, 1, 0, 0]]])
        row, _, _, _ = self.run_fake(certificate=certificate)
        self.assertEqual((row["flowtime"], row["movements"]), (4, 18))
        self.assertEqual(row["counterexample"], 1)
        self.assertEqual(row["flow_proven"], 0)
        witness = json.loads(Path(row["counterexample_file"]).read_text())
        self.assertEqual(witness["weighted_flowtime"], 4)
        self.assertEqual(witness["certificate_flowtime"], 3)
        self.assertEqual(witness["weighted_moves"], self.weighted["animation_moves"])
        self.assertEqual(witness["certificate_moves"], certificate["animation_moves"])

    def test_merge_rejects_duplicate_or_missing_instances(self):
        row, _, _, _ = self.run_fake()
        source = self.args.output
        with source.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=runner.FIELDNAMES)
            writer.writeheader()
            writer.writerow(row)
        destination = source.with_name("merged.csv")
        with redirect_stdout(io.StringIO()):
            runner.merge_batch(source, destination, "7", "1", 300, 71, 0.1)
        with destination.open(newline="") as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 1)
        with self.assertRaisesRegex(ValueError, "Missing or duplicate"):
            runner.merge_batch(source, destination, "7-8", "1", 300, 71, 0.1)
        with self.assertRaisesRegex(ValueError, "Unexpected certification_gap_threshold"):
            runner.merge_batch(source, destination, "7", "1", 300, 71, 0.2)

    def test_shell_dry_run_contains_sixteen_sequential_batches(self):
        script = Path(runner.__file__).with_name("RunTable2Weighted.sh")
        result = subprocess.run(["bash", str(script), "--dry-run", "--seeds", "1", "--output-dir", "/tmp/fake_results"],
                                capture_output=True, text=True, check=True)
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 16)
        self.assertEqual(sum("--formulation escortflow" in line for line in lines), 8)
        self.assertEqual(sum("--formulation loadflow" in line for line in lines), 8)
        self.assertTrue(all("--weighted-time-limit 300" in line for line in lines))
        self.assertTrue(all("--certification-time-limit 300" in line for line in lines))
        self.assertTrue(all("--certification-gap-threshold 0.1" in line for line in lines))
        self.assertEqual(sum("-l 1" in line for line in lines), 8)
        self.assertEqual(sum("-l 4" in line for line in lines), 8)


if __name__ == "__main__":
    unittest.main()
