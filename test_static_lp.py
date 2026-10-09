import csv
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import RunStaticLP as runner

HERE = Path(__file__).resolve().parent


def input_row(**changes):
    row = {"Lx x Ly": "2x2", "formulation": "escortflow", "retrieval_mode": "leave",
           "IOs": "[(0, 0)]", "Target Loads": "[(1, 0)]", "Escorts": "[(1, 1)]",
           "#Loads": "1", "# Escorts": "1", "seed": "7", "flow_weight": "5",
           "movement_weight": "0.2", "weighted_horizon": "1", "weighted_physical_horizon": "2"}
    row.update(changes)
    return row


def write_inputs(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def result_row(problem, source_hash):
    row = {key: problem[key] for key in ("layout", "escorts", "loads", "seed", "method", "flow_weight",
                                       "horizon", "physical_horizon", "retrieval_mode", "problem_sha256")}
    row.update(lp_objective=1.2, lp_flow=1, lp_movements=1, status="OPTIMAL", elapsed_seconds=0.1,
               source_sha256=source_hash, solver_version="13.0.3", algorithm="automatic")
    return row


class LPReplayTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def test_recorded_coefficients_and_horizon_are_preserved(self):
        escort = runner.parse_instance(input_row())
        load = runner.parse_instance(input_row(formulation="loadflow", weighted_horizon="2"))
        self.assertEqual(escort["flow_weight"], 5)
        self.assertEqual(escort["physical_horizon"], load["physical_horizon"])
        self.assertNotEqual(escort["horizon"], load["horizon"])
        self.assertNotEqual(escort["problem_sha256"], load["problem_sha256"])

    def test_malformed_or_inconsistent_model_inputs_are_rejected(self):
        for changes in ({"flow_weight": "0"}, {"weighted_physical_horizon": "1"},
                        {"movement_weight": "0.01"}, {"# Escorts": "2"},
                        {"Target Loads": "[(1, 1)]"}, {"Escorts": "[(2, 0)]"},
                        {"IOs": "[(0, 0), (0, 0)]"}, {"Moves": "LM"}, {"movement_mode": "LM"}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                runner.parse_instance(input_row(**changes))

    def test_duplicate_labels_cannot_hide_different_coefficients(self):
        path = self.directory / "inputs.csv"
        write_inputs(path, [input_row(), input_row(flow_weight="10", movement_weight="0.1")])
        with self.assertRaisesRegex(ValueError, "Duplicate instance"):
            runner.read_instances([path])

    def test_two_and_six_target_cases_have_distinct_result_and_reference_keys(self):
        path = self.directory / "inputs.csv"
        shared = {"Lx x Ly": "3x3", "Escorts": "[(2, 2)]"}
        two = input_row(**shared, **{"#Loads": "2", "Target Loads": "[(1, 0), (2, 0)]"})
        six = input_row(**shared, **{"#Loads": "6", "Target Loads":
                                    "[(1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2)]"})
        write_inputs(path, [two, six])
        problems = runner.read_instances([path])
        self.assertEqual([p["loads"] for p in problems], [2, 6])
        rows = [result_row(p, "source") for p in problems]
        self.assertNotEqual(runner.reference_key(rows[0]), runner.reference_key(rows[1]))
        reference_path = self.directory / "reference.csv"
        write_inputs(reference_path, rows)
        reference = runner.read_reference(reference_path)
        self.assertEqual(len(reference), 2)
        for row in rows:
            runner.check_reference(row, reference)
        changed_mode = dict(rows[0], retrieval_mode="continue")
        self.assertNotEqual(runner.reference_key(rows[0]), runner.reference_key(changed_mode))

    def test_legacy_single_target_reference_is_accepted(self):
        row = result_row(runner.parse_instance(input_row()), "source")
        old_reference = {key: value for key, value in row.items() if key not in {"loads", "retrieval_mode"}}
        path = self.directory / "reference.csv"
        write_inputs(path, [old_reference])
        runner.check_reference(row, runner.read_reference(path))

    def test_only_optimal_consistent_finite_values_are_accepted(self):
        valid = result_row(runner.parse_instance(input_row()), "source")
        runner.validate_value(valid)
        for changes in ({"status": "TIME_LIMIT"}, {"lp_objective": 2},
                        {"lp_flow": float("nan")}, {"lp_movements": -1}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                runner.validate_value(dict(valid, **changes))

    def test_reference_accepts_alternative_optimal_decompositions(self):
        row = result_row(runner.parse_instance(input_row()), "source")
        reference = {runner.reference_key(row): dict(row)}
        row.update(lp_flow=1.1, lp_movements=0.5)
        runner.validate_value(row)
        runner.check_reference(row, reference)
        row["lp_objective"] = 1.2001
        with self.assertRaisesRegex(ValueError, "differs from reference"):
            runner.check_reference(row, reference)

    def prepare_completed_run(self):
        path = self.directory / "inputs.csv"
        output = self.directory / "lp.csv"
        write_inputs(path, [input_row()])
        args = runner.parse_args(["--input", str(path), "--source-dir", str(HERE),
                                  "-f", str(output), "--resume"])
        problem = runner.parse_instance(input_row())
        sources = {name: runner.sha256(HERE / name) for name in runner.MODEL_FILES}
        manifest = dict(identity=dict(objective="F + M/R", model_sources=sources,
                                      problem_sha256=[problem["problem_sha256"]]), sessions=[])
        manifest_path = Path(str(output) + ".manifest.json")
        manifest_path.write_text(json.dumps(manifest))
        with output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=runner.FIELDS)
            writer.writeheader()
            writer.writerow(result_row(problem, runner.fingerprint(sources)))
        return args, manifest_path, path

    def test_verified_complete_resume_does_not_resolve(self):
        args, _, _ = self.prepare_completed_run()
        with patch.object(runner.concurrent.futures, "ProcessPoolExecutor") as pool, redirect_stdout(io.StringIO()):
            self.assertEqual(runner.run(args), 0)
        pool.assert_not_called()

    def test_resume_refuses_changed_sources(self):
        args, manifest_path, _ = self.prepare_completed_run()
        manifest = json.loads(manifest_path.read_text())
        manifest["identity"]["model_sources"][runner.MODEL_FILES[0]] = "changed"
        manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "Resume refused"):
            runner.run(args)

    def test_resume_refuses_changed_instance(self):
        args, _, path = self.prepare_completed_run()
        write_inputs(path, [input_row(weighted_horizon="2", weighted_physical_horizon="3")])
        with self.assertRaisesRegex(ValueError, "Resume refused"):
            runner.run(args)

    def test_resume_refuses_corrupted_result_components(self):
        args, _, _ = self.prepare_completed_run()
        with args.output.open(newline="") as handle:
            row = next(csv.DictReader(handle))
        row["lp_objective"] = "2"
        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=runner.FIELDS)
            writer.writeheader()
            writer.writerow(row)
        with self.assertRaisesRegex(ValueError, "objective disagrees"):
            runner.run(args)

    def test_resume_refuses_wrong_target_count_in_saved_metadata(self):
        args, _, _ = self.prepare_completed_run()
        with args.output.open(newline="") as handle:
            row = next(csv.DictReader(handle))
        row["loads"] = "6"
        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=runner.FIELDS)
            writer.writeheader()
            writer.writerow(row)
        with self.assertRaisesRegex(ValueError, "metadata disagrees"):
            runner.run(args)

    def test_retry_and_no_retry_never_accept_time_limited_lp_values(self):
        problem = runner.parse_instance(input_row())
        settings = dict(threads=1, time_limit=300, retry_time_limit=600, source_sha256="source")
        configs, parameters = [], []
        results = [dict(status_name="TIME_LIMIT", has_solution=True),
                   dict(status_name="OPTIMAL", has_solution=True,
                        objective=1.2, flowtime=1.1, movements=0.5)]

        class FakeSolver:
            def __init__(self, config):
                configs.append(config)
                self.env = SimpleNamespace(setParam=lambda *args: parameters.append(args))

            def solve(self, targets, escorts, horizon):
                return results.pop(0)

            def close(self):
                pass

        module = SimpleNamespace(StaticEscortFlowGurobiSolver=FakeSolver, StaticGurobiConfig=SimpleNamespace)
        gp = SimpleNamespace(gurobi=SimpleNamespace(version=lambda: (13, 0, 3)))
        with patch.object(runner.importlib, "import_module", side_effect=lambda name: gp if name == "gurobipy" else module):
            row = runner.solve_instance(problem, settings)
            self.assertEqual(row["algorithm"], "barrier_retry")
            self.assertEqual([c.time_limit for c in configs], [300, 600])
            self.assertTrue(all(c.lp and c.gamma == 0.2 for c in configs))
            self.assertEqual(parameters, [("Method", 2), ("Crossover", 0)])
            results[:] = [dict(status_name="TIME_LIMIT", has_solution=True)]
            with self.assertRaisesRegex(RuntimeError, "LP status TIME_LIMIT"):
                runner.solve_instance(problem, dict(settings, retry_time_limit=None))


class StandardLPOptionTests(unittest.TestCase):
    def test_invalid_flow_weight_combinations_fail_before_solving(self):
        base = ["-x", "2", "-y", "2", "-O", "0", "0"]
        for script in ("EscortFlowStatic.py", "LoadFlowStatic.py"):
            for flags in (["--flow-weight", "5"], ["--lp", "--flow-weight", "0"],
                          ["--lp", "--flow-weight", "5", "--gamma", "0.3"],
                          ["--lp", "--flow-weight", "5", "--cutoff"]):
                with self.subTest(script=script, flags=flags):
                    result = subprocess.run([sys.executable, str(HERE / script), *base, *flags],
                                            capture_output=True, text=True, cwd=HERE)
                    self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                    self.assertIn("--flow-weight", result.stderr)

    @unittest.skipUnless(importlib.util.find_spec("gurobipy"), "Gurobi Python API is unavailable")
    def test_both_standard_runners_match_archived_lp_objectives(self):
        archive = HERE / "revision_R1" / "table3_bounds_2026-10-09"
        if not (archive / "lp_results.csv").exists():
            self.skipTest("Archived Table 3 references are unavailable")
        references = runner.read_reference(archive / "lp_results.csv")
        with tempfile.TemporaryDirectory() as temporary:
            for method, script, horizon in (("escortflow", "EscortFlowStatic.py", 10),
                                             ("loadflow", "LoadFlowStatic.py", 11)):
                with self.subTest(method=method):
                    output = Path(temporary) / f"{method}.csv"
                    command = [sys.executable, str(HERE / script), "-x", "13", "-y", "7", "-O", "6", "0",
                               "-e", "3", "-l", "1", "-r", "1", "-m", "leave", "--lp",
                               "--flow-weight", "881", "--horizon", str(horizon), "--num_threads", "1",
                               "-t", "300", "-f", str(output)]
                    result = subprocess.run(command, capture_output=True, text=True, cwd=HERE, timeout=330)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    with output.open() as handle:
                        rows = list(csv.DictReader((line for line in handle if line.strip()), skipinitialspace=True))
                    self.assertEqual(len(rows), 1)
                    row = rows[0]
                    self.assertEqual(row["Model"].strip(), "LP-Gurobi")
                    self.assertAlmostEqual(float(row["gamma"]), 1 / 881, places=12)
                    self.assertAlmostEqual(float(row["obj"]),
                                           float(references[("13x7", 3, 1, 1, method, "leave")]["lp_objective"]), places=6)


if __name__ == "__main__":
    unittest.main()
