"""Check campaign coverage, pairing and the Linux launch commands without solving."""
import csv
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

import RunStaticCampaign as campaign


class CampaignTests(unittest.TestCase):
    def dry_run(self, script):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "not-created"
            result = subprocess.run([sys.executable, str(Path(__file__).with_name(script)),
                                     "--dry-run", "--threads", "16", "--output-dir", str(output)],
                                    capture_output=True, text=True, check=True)
            self.assertFalse(output.exists())
            return [shlex.split(line) for line in result.stdout.splitlines()]

    def test_all_campaigns_have_exact_configurations_and_budgets(self):
        for script, count, mode, loads in [("Run70Percent.py", 8, "leave", {4}),
                                           ("RunContinue.py", 96, "continue", {2, 4, 6}),
                                           ("RunTable2Targets.py", 64, "leave", {2, 6})]:
            with self.subTest(script=script):
                commands = self.dry_run(script)
                self.assertEqual(len(commands), count)
                observed, destinations = set(), set()
                for command in commands:
                    option = lambda name: command[command.index(name) + 1]
                    self.assertEqual(option("--retrieval-mode"), mode)
                    self.assertEqual(option("--threads"), "16")
                    self.assertEqual(option("--weighted-time-limit"), "300")
                    self.assertEqual(option("--extension-time-limit"), "300")
                    self.assertEqual(option("-r"), "1-100")
                    observed.add((option("-x"), option("-y"), option("-e"), int(option("-l")), option("--formulation")))
                    destinations.add(option("-f"))
                self.assertEqual(len(observed), count)
                self.assertEqual(len(destinations), count)
                self.assertEqual({x[3] for x in observed}, loads)
                self.assertEqual({(x[0], x[1]) for x in observed}, {("13", "7"), ("10", "10"), ("16", "10"), ("27", "10")})
                for lx, ly, occupancy_escorts, _ in campaign.LAYOUTS:
                    counts = {int(x[2]) for x in observed if (x[0], x[1]) == (str(lx), str(ly))}
                    expected = {occupancy_escorts} if script == "Run70Percent.py" else {8, 12, 16, occupancy_escorts}
                    self.assertEqual(counts, expected)

    def test_exact_occupancy_is_recorded(self):
        layouts = ["{}x{}".format(x, y) for x, y, _, _ in campaign.LAYOUTS]
        configs = list(campaign.configurations("occupancy70", layouts))
        self.assertEqual(len(configs), 4)
        self.assertAlmostEqual(configs[0]["occupancy"], 64 / 91)
        self.assertTrue(all(abs(c["occupancy"] - .7) < .004 for c in configs))

    def test_pairing_rejects_missing_or_different_instances(self):
        config = dict(lx=3, ly=2, escorts=1, loads=2)
        row = {field: "shared" for field in (
            "protocol", "threads", "warmstart", "retrieval_mode", "movement_mode", "flow_weight",
            "movement_integer_weight", "safe_movement_bound", "safe_movement_lower_bound",
            "safe_flow_horizon", "weighted_physical_horizon", "weighted_global_scope",
            "greedy_flowtime", "greedy_movements", "greedy_makespan", "weighted_time_limit",
            "extension_time_limit", "flow_proof_check_mode")}
        row.update({"Lx x Ly": "3x2", "# Escorts": "1", "#Loads": "2", "seed": "1",
                    "IOs": "[(0, 0)]", "Escorts": "[(0, 1)]", "Target Loads": "[(1, 0), (2, 0)]"})
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            def write(method, records):
                with (p / ("continue_" + method + ".csv")).open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerows(records)
            write("escortflow", [row]); write("loadflow", [row])
            with redirect_stdout(io.StringIO()):
                campaign.validate_pairs(p, [1], [config], "continue")
            self.assertEqual(json.loads((p / "pairing.json").read_text())["matched_instances"], 1)
            write("loadflow", [dict(row, weighted_physical_horizon="different")])
            with self.assertRaisesRegex(ValueError, "Different paired"):
                campaign.validate_pairs(p, [1], [config], "continue")
            write("loadflow", [])
            with self.assertRaisesRegex(ValueError, "Missing or unexpected"):
                campaign.validate_pairs(p, [1], [config], "continue")


if __name__ == "__main__":
    unittest.main()
