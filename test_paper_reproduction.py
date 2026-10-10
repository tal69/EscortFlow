"""Coverage, metric definitions, source integrity, and restart checks for paper reproduction."""

import csv
import io
import importlib.util
from contextlib import redirect_stdout
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import PaperTables as tables
import ReproducePaper as paper
from RunStaticLP import parse_instance


def example_row(method="loadflow"):
    # Two targets, D=4, d_max=2, K=10, R=37, common physical horizon 6.
    return {"Lx x Ly": "4x3", "# Escorts": "2", "#Loads": "2", "seed": "1",
        "formulation": method, "retrieval_mode": "continue", "movement_mode": "BM",
        "IOs": "[(0,0)]", "Escorts": "[(0,0),(3,2)]", "Target Loads": "[(2,0),(1,1)]",
        "flow_weight": "37", "movement_weight": str(1/37), "weighted_horizon": "5" if method == "escortflow" else "6",
        "weighted_physical_horizon": "6", "safe_flow_horizon": "4", "safe_movement_bound": "40",
        "safe_movement_lower_bound": "4", "greedy_flowtime": "6", "greedy_movements": "10",
        "has_solution": "1", "flowtime": "6", "movements": "8", "scaled_objective": "230",
        "scaled_best_bound": "179", "weighted_bound_consistent": "1", "phase1_flow_proven": "0",
        "weighted_proven": "0", "weighted_cpu_time": "12", "weighted_cpu_time_is_estimate": "1",
        "weighted_time_limit": "300", "extension_time_limit": "300", "threads": "16", "warmstart": "1",
        "final_has_solution": "1", "final_flowtime": "5", "final_movements": "7",
        "final_lexicographic_proven": "0", "protocol": "safe_integer_flow_timing_v5", "error": ""}


def example_lp(row):
    problem = parse_instance(row)
    return dict(problem, lp_objective=4 + 4/37, lp_flow=4, lp_movements=4, status="OPTIMAL")


class ReproductionTests(unittest.TestCase):
    def test_default_plan_covers_both_models_modes_and_all_requested_rows(self):
        configs = list(paper.paper_configurations([f"{x}x{y}" for x, y, _, _ in paper.LAYOUTS]))
        self.assertEqual(len(configs), 120)
        self.assertEqual(len({(c["mode"], c["lx"], c["ly"], c["loads"], c["escorts"]) for c in configs}), 120)
        for lx, ly, occupancy, _ in paper.LAYOUTS:
            counts = lambda mode, loads: {c["escorts"] for c in configs if
                (c["lx"], c["ly"], c["mode"], c["loads"]) == (lx, ly, mode, loads)}
            self.assertEqual(counts("leave", 1), set(range(3,9)))
            self.assertEqual(counts("leave", 2), {8,12,16})
            self.assertEqual(counts("leave", 4), {8,12,16,occupancy})
            self.assertEqual(counts("leave", 6), {8,12,16,20})
            self.assertEqual(counts("continue", 1), set())
            self.assertEqual(counts("continue", 2), {8,12,16,occupancy})
            self.assertEqual(counts("continue", 4), {8,12,16,occupancy})
            self.assertEqual(counts("continue", 6), {8,12,16,20,occupancy})

    def test_parameter_free_and_seed_subset_dry_runs_require_no_solver(self):
        for seed_args, expected in (([], "1-100"), (["5-7"], "5-7"), (["--seeds", "2-4"], "2-4")):
            with tempfile.TemporaryDirectory() as directory:
                destination = Path(directory) / "not-created"
                result = subprocess.run([sys.executable, paper.__file__, *seed_args,
                    "--dry-run", "--output-dir", str(destination)], text=True, capture_output=True, check=True)
                self.assertFalse(destination.exists())
                commands = [shlex.split(line) for line in result.stdout.splitlines() if "RunSafeWeightedStatic.py" in line]
                self.assertEqual(len(commands), 240)
                self.assertTrue(all(c[c.index("--seeds")+1] == expected for c in commands))
                self.assertEqual(sum("RunStaticLP.py" in line for line in result.stdout.splitlines()), 2)

    def test_invalid_seed_and_output_settings_rejected(self):
        for arguments in (["1,1"], ["-1"], ["--threads", "0"], ["--weighted-time-limit", "0"],
                          ["1", "--seeds", "2"], ["--layouts", "13x7", "13x7"]):
            with self.subTest(arguments=arguments), redirect_stdout(io.StringIO()), patch("sys.stderr", io.StringIO()):
                with self.assertRaises(SystemExit):
                    paper.parse_args(arguments)

    def test_multi_target_integer_bounds_and_shared_extension_reference(self):
        row = example_row()
        metrics = tables.component_metrics(row, (5,7), example_lp(row))
        # ceil((179 + 10*(4-2)) / (37+10)) = 5, whereas single-target formula gives 4.
        self.assertEqual(metrics["ft_lower_bound"], 5)
        self.assertEqual(metrics["mv_lower_bound_at_best_ft"], 4)
        self.assertEqual(metrics["ft_gap_pct"], 0)
        self.assertAlmostEqual(metrics["mv_gap_pct"], 100*3/7)
        self.assertAlmostEqual(metrics["ft_improvement_pct"], 100/6)
        self.assertEqual(metrics["mv_improvement_pct"], 30)
        self.assertEqual(metrics["ft_cpu_seconds"], 12)

    def test_cutoff_certificates_and_blank_integer_gaps(self):
        row = example_row()
        row.update(flowtime="5", movements="7", scaled_objective="192", scaled_best_bound="192",
                   phase1_flow_proven="1", weighted_proven="1", first_flow_proof_cpu_time="3",
                   first_flow_proof_runtime="2", first_flow_proof_flowtime="5")
        records = []
        for method in paper.METHODS:
            r = dict(row, formulation=method, weighted_horizon="5" if method == "escortflow" else "6")
            records.append(tables.component_metrics(r, (5,7), example_lp(r)))
        comparisons, _ = tables.aggregate(records)
        rendered = tables.render_fragment(2, comparisons)
        data = next(line for line in rendered.splitlines() if "4$\\times$3" in line).split(" & ")
        self.assertEqual(len(data), 16)
        for offset in (2,9):
            self.assertEqual(data[offset:offset+4], ["100","100","",""])
            self.assertNotEqual(data[offset+4], "")
        self.assertEqual(rendered.count(r"\textbf{"), 4)  # both timing ties are bold
        self.assertNotIn("landscape", rendered)

    def test_instance_percentages_are_averaged_not_ratios_of_means(self):
        first = tables.component_metrics(example_row(), (5,7), example_lp(example_row()))
        second = dict(first, seed=2, ft_gap_pct=50, best_ft=100)
        comparisons, _ = tables.aggregate([first, second])
        self.assertEqual(comparisons[0]["ft_gap_pct"], 25)
        self.assertEqual(comparisons[0]["instances"], 2)

    def test_lp_replay_cannot_mix_modes_horizons_or_nonoptimal_values(self):
        row = example_row()
        for change in ({"retrieval_mode":"leave"}, {"physical_horizon":7}, {"flow_weight":38}, {"status":"TIME_LIMIT"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                tables.component_metrics(row, (5,7), dict(example_lp(row), **change))

    def test_movement_improvement_can_be_negative_and_zero_reference_is_defined(self):
        self.assertEqual(tables.improvement(4,5), -25)
        self.assertEqual(tables.improvement(0,0), 0)
        self.assertEqual(tables.percentage_gap(0,0), 0)
        with self.assertRaises(ValueError):
            tables.percentage_gap(0,1)

    def test_resume_rejects_modified_frozen_sources_before_solving(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/"source").mkdir()
            path = root/"source"/"model.py"
            path.write_text("original")
            manifest = dict(source_sha256={"model.py":paper.sha256(path)}, reference_inputs=[])
            (root/"manifest.json").write_text(json.dumps(manifest))
            path.write_text("modified")
            with self.assertRaisesRegex(ValueError, "Frozen source changed"):
                paper.main(["--resume", str(root)])

    @unittest.skipUnless(importlib.util.find_spec("gurobipy") and importlib.util.find_spec("numpy"),
                         "integer CSV validation needs the solver environment, but performs no solve")
    def test_partial_restart_recovers_pending_rows_and_runs_only_missing_seeds(self):
        from RunSafeWeightedStatic import FIELDNAMES
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            configs = [dict(lx=4, ly=3, escorts=2, loads=2, outputs=[0,0], mode=mode)
                       for mode in ("leave", "continue")]
            manifest = dict(configurations=configs, seeds="1-3", seed_values=[1,2,3],
                            settings=dict(threads=16, weighted_time_limit=300, extension_time_limit=300))

            def row_for(method, mode, seed):
                row = dict.fromkeys(FIELDNAMES, "")
                row.update(example_row(method), seed=str(seed), retrieval_mode=mode, total_time_limit="600",
                           **{"Solver Status":"BUDGET_REACHED", "movement_integer_weight":"1",
                              "weighted_global_scope":"1", "optimization_calls":"1",
                              "flow_proof_check_mode":"bound_or_flow_change"})
                return row

            def write_rows(path, rows):
                with path.open("w", newline="") as handle:
                    writer=csv.DictWriter(handle,fieldnames=FIELDNAMES)
                    writer.writeheader();writer.writerows(rows)

            for config in configs:
                for name in ("parts", "logs"):
                    (root/config["mode"]/name).mkdir(parents=True,exist_ok=True)
                for method in paper.METHODS:
                    path=root/config["mode"]/"parts"/(paper.stem(config,method)+".csv")
                    seed_values=[1] if config["mode"]=="leave" and method=="loadflow" else [1,2,3]
                    write_rows(path,[row_for(method,config["mode"],s) for s in seed_values])
                    if len(seed_values)==1:
                        write_rows(path.with_suffix(".pending.csv"),[row_for(method,config["mode"],2)])
            launched=[]

            def fake_solve(command, **kwargs):
                destination=Path(command[command.index("-f")+1])
                self.assertFalse(destination.exists())
                self.assertEqual(command[command.index("--seeds")+1],"3")
                launched.append(command)
                write_rows(destination,[row_for("loadflow","leave",3)])
                return SimpleNamespace(returncode=0)

            with patch("ReproducePaper.subprocess.run",side_effect=fake_solve),redirect_stdout(io.StringIO()):
                paper.run_integer_batches(sys.executable,root,manifest)
            self.assertEqual(len(launched),1)
            rows=tables.read_csv(root/"leave"/"loadflow.csv")
            self.assertEqual([r["seed"] for r in rows],["1","2","3"])
            self.assertEqual(rows[0],row_for("loadflow","leave",1))
            self.assertEqual(rows[1],row_for("loadflow","leave",2))
            self.assertFalse(list(root.glob("*/parts/*.pending.csv")))


if __name__ == "__main__":
    unittest.main()
