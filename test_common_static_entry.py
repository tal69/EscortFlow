"""Check the public reproduction commands and a small actual paired solve."""
import csv
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent


class CommonStaticEntryTests(unittest.TestCase):
    def dry_run(self, destination, **overrides):
        # Keep the test independent of campaign overrides in the user's shell.
        env = dict(os.environ)
        for key in ("PYTHON", "SEEDS", "THREADS", "LAYOUTS", "TARGET_COUNTS",
                    "ESCORT_COUNTS", "TIME_LIMIT", "EXTENSION_LIMIT", "WITH_LP",
                    "LP_THREADS", "LP_TIME_LIMIT", "LP_RETRY_LIMIT"):
            env.pop(key, None)
        env.update(overrides)
        result = subprocess.run(["bash", str(HERE / "RunContinue.sh"), "--dry-run",
                                 str(destination)], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return [shlex.split(line) for line in result.stdout.splitlines()]

    def test_default_plan_has_complete_paired_scope_and_correct_settings(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "results"
            commands = self.dry_run(destination)
            self.assertFalse(destination.exists())
        self.assertEqual(len(commands), 72)
        outputs = set()
        configurations = set()
        for lf, ef in zip(commands[::2], commands[1::2]):
            self.assertEqual(lf[2], ef[2])
            self.assertEqual(Path(lf[2]).name, "SolveStatic.py")
            self.assertEqual(lf[lf.index("--formulation") + 1], "loadflow")
            self.assertEqual(ef[ef.index("--formulation") + 1], "escortflow")
            for command in (lf, ef):
                for option, value in (("-r", "1-100"), ("-m", "continue"),
                                      ("--threads", "16"), ("--weighted-time-limit", "300"),
                                      ("--extension-time-limit", "300"), ("--lp-threads", "1"),
                                      ("--lp-time-limit", "300"), ("--lp-retry-time-limit", "600")):
                    self.assertEqual(command[command.index(option) + 1], value)
                self.assertIn("--with-lp", command)
                self.assertNotIn("--stop-at-flow-proof", command)
                outputs.add(command[command.index("-f") + 1])
            # Apart from formulation and output name the paired options are identical.
            def common(command):
                copy = command[3:]
                for option in ("--formulation", "-f"):
                    index = copy.index(option)
                    del copy[index:index + 2]
                return copy
            self.assertEqual(common(lf), common(ef))
            configurations.add(tuple(lf[lf.index(flag) + 1] for flag in ("-x", "-y", "-l", "-e")))
        self.assertEqual(len(outputs), 72)
        self.assertEqual(len(configurations), 36)
        self.assertEqual([commands[0][commands[0].index(f) + 1] for f in ("-x", "-y", "-l")],
                         ["13", "7", "4"])
        output_coordinates = {("13", "7"): ["6", "0"], ("10", "10"): ["0", "0"],
                              ("16", "10"): ["4", "0", "11", "0"],
                              ("27", "10"): ["4", "0", "13", "0", "22", "0"]}
        for command in commands:
            layout = tuple(command[command.index(f) + 1] for f in ("-x", "-y"))
            self.assertEqual(command[command.index("-O") + 1:command.index("-l")],
                             output_coordinates[layout])

    def test_subset_and_no_lp_are_explicit_and_write_nothing(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "results"
            commands = self.dry_run(destination, SEEDS="7-9", THREADS="1", LAYOUTS="16x10",
                                    TARGET_COUNTS="2", ESCORT_COUNTS="12", WITH_LP="0")
            self.assertFalse(destination.exists())
        self.assertEqual(len(commands), 2)
        for command in commands:
            self.assertEqual(command[command.index("-r") + 1], "7-9")
            self.assertEqual(command[command.index("--threads") + 1], "1")
            self.assertNotIn("--with-lp", command)

    def test_existing_campaign_directory_is_refused_without_modification(self):
        with tempfile.TemporaryDirectory() as temporary:
            marker = Path(temporary) / "keep.txt"
            marker.write_text("existing results")
            result = subprocess.run(["bash", str(HERE / "RunContinue.sh"), temporary],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn("already exists", result.stderr)
            self.assertEqual(marker.read_text(), "existing results")
            self.assertEqual(list(Path(temporary).iterdir()), [marker])

    def test_public_entry_uses_matching_corrected_defaults_and_accepted_starts(self):
        rows = []
        with tempfile.TemporaryDirectory() as temporary:
            for method in ("loadflow", "escortflow"):
                output = Path(temporary) / (method + ".csv")
                command = [sys.executable, str(HERE / "SolveStatic.py"), "--formulation", method,
                           "-x", "3", "-y", "2", "-O", "0", "0", "-l", "2", "-e", "2",
                           "-r", "1", "-m", "continue", "--threads", "1",
                           "--weighted-time-limit", "10", "--extension-time-limit", "0",
                           "--with-lp", "--lp-time-limit", "10", "--lp-retry-time-limit", "10",
                           "-f", str(output)]
                result = subprocess.run(command, capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("Loaded user MIP start", result.stdout)
                with output.open() as handle:
                    batch = list(csv.DictReader(handle))
                self.assertEqual(len(batch), 1)
                row = batch[0]
                self.assertEqual(row["error"], "")
                self.assertEqual(row["protocol"], "safe_integer_flow_timing_v6")
                self.assertEqual(row["warmstart"], "1")
                self.assertEqual(row["lexicographic_proven"], "1")
                self.assertEqual(row["lp_status"], "OPTIMAL")
                for field in ("model_num_variables", "model_num_constraints", "model_num_nonzeros"):
                    self.assertGreater(int(row[field]), 0)
                rows.append(row)
            lf, ef = rows
            for field in ("IOs", "Escorts", "Target Loads", "flow_weight", "weighted_physical_horizon",
                          "greedy_flowtime", "greedy_movements", "flowtime", "movements"):
                self.assertEqual(lf[field], ef[field], field)
            self.assertEqual(int(lf["weighted_horizon"]), int(ef["weighted_horizon"]) + 1)


if __name__ == "__main__":
    unittest.main()
