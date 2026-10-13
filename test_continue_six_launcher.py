"""Verify revised campaign scope and budgets without launching experiments."""

import json
import importlib.util
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

import RunContinueSixTargets as launcher


class SixTargetLauncherTests(unittest.TestCase):
    def test_exact_scope_and_method_specific_budgets(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "new-results"
            result = subprocess.run([sys.executable, launcher.__file__, "--dry-run", str(output)],
                                    capture_output=True, text=True, check=True)
            self.assertFalse(output.exists())
            commands = [shlex.split(line) for line in result.stdout.splitlines()]
        self.assertEqual(len(commands), 24)
        observed = set()
        for command in commands:
            value = lambda flag: command[command.index(flag)+1]
            method = value("--formulation")
            observed.add((value("-x"), value("-y"), value("-e"), method))
            self.assertEqual(value("-l"), "6")
            self.assertEqual(value("-r"), "1-100")
            self.assertEqual(value("-m"), "continue")
            self.assertEqual(value("--weighted-time-limit"), "300")
            self.assertEqual(value("--extension-time-limit"), "0" if method == "loadflow" else "300")
            self.assertEqual(value("--threads"), "16")
            self.assertEqual(value("--lp-time-limit"), "300")
            self.assertEqual(value("--lp-retry-time-limit"), "600")
            self.assertIn("--with-lp", command)
            self.assertNotIn("--stop-at-flow-proof", command)
        self.assertEqual(len(observed), 24)
        self.assertEqual({x[2] for x in observed}, {"12", "16", "20"})

    def test_refuses_existing_output_without_modifying_it(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            sentinel = output / "kept.txt"
            sentinel.write_text("keep")
            result = subprocess.run([sys.executable, launcher.__file__, str(output)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("already exists", result.stderr)
            self.assertEqual(sentinel.read_text(), "keep")
            self.assertEqual(list(output.iterdir()), [sentinel])

    def test_source_validation_detects_changes_before_freezing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            source.mkdir()
            names = ("SolveStatic.py", "RunSafeWeightedStatic.py",
                     "escort_flow_static_gurobi.py", "load_flow_static_gurobi.py")
            for name in names:
                (source / name).write_text("original")
            env = dict(protocol=launcher.PROTOCOL, retrieval_mode="continue",
                       source_sha256={name:launcher.sha(source/name) for name in names})
            (root/"environment.json").write_text(json.dumps(env))
            self.assertEqual(launcher.validate_source(source), env)
            (source/"SolveStatic.py").write_text("modified")
            with self.assertRaisesRegex(ValueError, "integrity"):
                launcher.validate_source(source)

    def test_seed_subsets_and_duplicates(self):
        self.assertEqual(launcher.seeds("1-3,8"), [1,2,3,8])
        for text in ("1,1", "0", "101", "4-2", "nan"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                launcher.seeds(text)

    def test_table_plan_uses_revised_six_target_source_and_order(self):
        path = Path(launcher.__file__).parent / "revision_R1/continue_tables_2026-10-11/build_tables.py"
        spec = importlib.util.spec_from_file_location("continue_builder", path)
        builder = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(builder)
        self.assertEqual(builder.ESCORTS_BY_TARGET[6], (12,16,20))
        self.assertEqual(builder.ESCORTS_BY_TARGET[4], (8,12,16))
        self.assertEqual(builder.campaign_root(6).name, "results_continue_v6_six_targets")
        self.assertEqual(builder.campaign_root(2).name, "results_continue_v6")
        self.assertEqual(builder.extension_budget(6, "loadflow"), 0)
        self.assertEqual(builder.extension_budget(4, "loadflow"), 300)
        template = "\n".join(builder.tables.tabular_header(2)) + "\n13$\\times$7 (1) & 8 \\\\\n"
        fragment = builder.render_fragment(2, template, {2,4,6}, [], [], [])
        rows = builder.table_rows(fragment)
        self.assertEqual(len(rows), 36)
        self.assertEqual([k[2] for k in rows if k[:2] == (6,"13x7")], [12,16,20])


if __name__ == "__main__":
    unittest.main()
