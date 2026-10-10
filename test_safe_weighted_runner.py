import csv
from contextlib import redirect_stdout
import io
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import RunSafeWeightedStatic as runner


class FakeSolver:
    def __init__(self, result):
        self.result = result
        self.starts = []
        self.solves = []
        self.closed = False

    def build_warmstart_from_trace(self, targets, escorts, horizon, target_history, escort_history):
        start = (targets, escorts, horizon, target_history, escort_history)
        self.starts.append(start)
        return start

    def solve(self, targets, escorts, horizon, warmstart):
        self.solves.append((targets, escorts, horizon, warmstart))
        return self.result

    def close(self):
        self.closed = True


class SafeWeightedRunnerTests(unittest.TestCase):
    def test_continue_mode_reaches_both_backends(self):
        self.args.retrieval_mode = "continue"
        for formulation, constructor in [
            ("escortflow", "escort_flow_static_gurobi.StaticEscortFlowGurobiSolver"),
            ("loadflow", "load_flow_static_gurobi.LoadFlowStaticGurobiSolver")]:
            self.args.formulation = formulation
            with patch(constructor) as make:
                runner.make_solver(self.args, flow_weight=34)
                self.assertEqual(make.call_args.args[0].retrieval_mode, "continue")

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.args = SimpleNamespace(
            formulation="escortflow", Lx=3, Ly=3, loads=1, outputs=[(0, 0)],
            threads=16, weighted_time_limit=300, extension_time_limit=71,
            output=Path(self.temporary.name) / "batch.csv",
        )
        self.targets = [(2, 0)]
        self.escorts = [(2, 1), (2, 2)]
        # Seven initial loads, H=5, U=35, D=2, R=34. Incumbent Z=146.
        self.initial = self.solution(status_name="BUDGET_REACHED", runtime=300,
                                     bound_checkpoint_runtime=299.9, node_count=100)
        self.final = self.solution(status_name="TIME_LIMIT", runtime=371, cpu_time=372, work=65)
        self.trace = (5, 5, 22, [], [{"escort": 1}], [{"target": 1}])

    def solution(self, **changes):
        result = dict(has_solution=True, status_name="TIME_LIMIT", runtime=300, cpu_time=301,
                      work=42, weighted_proven=False, flowtime=4, movements=10, makespan=4,
                      bound_consistent=True, weight_scale=34, scaled_best_bound=100)
        result.update(changes)
        if result["has_solution"]:
            objective = result["weight_scale"] * result["flowtime"] + result["movements"]
            result.update(scaled_objective=objective, objective=objective / result["weight_scale"])
            bound = result["scaled_best_bound"]
            result.update(best_bound=None if bound is None else bound / result["weight_scale"],
                          scaled_absolute_gap=None if bound is None else abs(objective - bound),
                          absolute_gap=None if bound is None else abs(objective - bound) / result["weight_scale"])
        else:
            for key in ("flowtime", "movements", "makespan", "scaled_objective", "objective",
                        "scaled_absolute_gap", "absolute_gap"):
                result[key] = None
            bound = result["scaled_best_bound"]
            result["best_bound"] = None if bound is None else bound / result["weight_scale"]
        return result

    def run_fake(self, initial=None, final=None, **metadata):
        search = dict(self.final if final is None else final)
        search.update(phase1_snapshot=dict(self.initial if initial is None else initial),
                      extension_used=True, extension_runtime=71, phase_transition_runtime=300.01,
                      optimization_calls=1, stop_reason="TOTAL_TIME_LIMIT")
        search.update(metadata)
        solver = FakeSolver(search)
        calls = []

        def factory(args, **kwargs):
            calls.append(kwargs)
            return solver

        with patch.object(runner.OneStepHeuristic_v2, "SolveGreedy", return_value=self.trace), \
                patch.object(runner, "GeneretaeRandomInstance", return_value=(self.targets, self.escorts)), \
                redirect_stdout(io.StringIO()):
            row = runner.run_instance(self.args, 7, 2, factory)
        return row, solver, calls

    def test_exact_early_weighted_proof_proves_both_without_extension(self):
        exact = self.solution(weighted_proven=True, status_name="OPTIMAL", scaled_best_bound=146, runtime=2)
        row, solver, calls = self.run_fake(initial=exact, final=exact, extension_used=False,
                                           extension_runtime=0, stop_reason="WEIGHTED_OPTIMUM")
        self.assertEqual(row["error"], "")
        self.assertEqual((row["flow_weight"], row["safe_movement_bound"], row["safe_flow_horizon"]), (34, 35, 5))
        self.assertEqual(row["safe_movement_lower_bound"], 2)
        self.assertEqual(row["protocol"], "safe_integer_flow_timing_v6")
        self.assertEqual(row["weighted_horizon"], 5)
        self.assertEqual(row["weighted_global_scope"], 1)
        self.assertEqual(row["flow_proven"], 1)
        self.assertEqual(row["lexicographic_proven"], 1)
        self.assertEqual(row["flow_proof_source"], "safe_weighted_optimum")
        self.assertEqual(row["extension_used"], 0)
        self.assertEqual(row["weighted_runtime"], 2)
        self.assertEqual(row["final_runtime"], 2)
        self.assertEqual(len(calls), 1)
        self.assertEqual(len(solver.solves), 1)
        self.assertTrue(solver.closed)

    def test_gap_proof_at_first_budget_preserves_phase_one_metrics(self):
        initial = self.solution(status_name="BUDGET_REACHED", scaled_best_bound=124)
        final = self.solution(status_name="INTERRUPTED", scaled_best_bound=124, runtime=300.01)
        row, _, calls = self.run_fake(initial=initial, final=final, extension_used=False,
                                      extension_runtime=0, stop_reason="FLOW_PROVEN_AT_PHASE1_LIMIT")
        self.assertEqual(row["gap_flow_proven"], 1)
        self.assertEqual(row["flow_proof_source"], "weighted_gap")
        self.assertEqual(row["lexicographic_proven"], 0)
        self.assertEqual(row["better_flow_movement_bound"], 21)
        self.assertAlmostEqual(row["weighted_bound_threshold"], 123.001)
        self.assertAlmostEqual(row["weighted_gap_threshold"], 22.999)
        self.assertEqual(row["scaled_absolute_gap"], 22)
        self.assertEqual(row["extension_used"], 0)
        self.assertEqual(len(calls), 1)

    def test_early_flow_stop_reports_all_incumbent_kpis_for_both_formulations(self):
        self.args.stop_at_flow_proof = True
        stopped = self.solution(status_name="INTERRUPTED", scaled_best_bound=124,
                                runtime=12.5, cpu_time=13.2, flow_proven=True)
        for formulation in ("escortflow", "loadflow"):
            with self.subTest(formulation=formulation):
                self.args.formulation = formulation
                row, _, _ = self.run_fake(initial=stopped, final=stopped,
                    extension_used=False, extension_runtime=0, stop_reason="FLOW_PROVEN_EARLY",
                    first_flow_proof_runtime=12.5, first_flow_proof_cpu_time=13.2,
                    first_flow_proof_flowtime=4, first_flow_proof_method="weighted_gap")
                self.assertEqual(row["error"], "")
                self.assertEqual(row["stop_at_flow_proof"], 1)
                self.assertEqual(row["search_stop_reason"], "FLOW_PROVEN_EARLY")
                self.assertEqual((row["flowtime"], row["movements"], row["scaled_objective"]), (4, 10, 146))
                self.assertEqual((row["final_flowtime"], row["final_movements"], row["final_scaled_objective"]), (4, 10, 146))
                self.assertEqual(row["scaled_absolute_gap"], 22)
                self.assertEqual(row["flow_proven"], 1)
                self.assertEqual(row["lexicographic_proven"], 0)
                self.assertEqual(row["final_lexicographic_proven"], 0)
                self.assertEqual(row["weighted_runtime"], 12.5)
                self.assertEqual(row["final_runtime"], 12.5)
                self.assertEqual(row["first_flow_proof_runtime"], 12.5)

    def test_extension_bound_certifies_original_without_overwriting_original_gap(self):
        final = self.solution(status_name="INTERRUPTED", scaled_best_bound=124, runtime=320)
        row, solver, calls = self.run_fake(final=final, extension_runtime=20, stop_reason="FLOW_PROVEN")
        self.assertEqual(row["error"], "")
        self.assertEqual(calls, [dict(flow_weight=34)])
        self.assertEqual(len(solver.solves), 1)
        self.assertEqual(len(solver.starts), 1)
        self.assertEqual(row["weighted_runtime"], 300)
        self.assertEqual(row["final_runtime"], 320)
        self.assertEqual(row["extension_runtime"], 20)
        self.assertEqual(row["scaled_best_bound"], 100)
        self.assertEqual(row["scaled_absolute_gap"], 46)
        self.assertEqual(row["proof_scaled_best_bound"], 124)
        self.assertEqual(row["proof_scaled_absolute_gap"], 22)
        self.assertEqual(row["flow_proven"], 1)
        self.assertEqual(row["flow_proof_source"], "weighted_gap")
        self.assertEqual(row["lexicographic_proven"], 0)
        self.assertEqual(row["optimization_calls"], 1)

    def test_movement_improvement_is_reported_only_in_final_columns(self):
        final = self.solution(status_name="INTERRUPTED", movements=8, scaled_best_bound=124, runtime=320)
        row, _, _ = self.run_fake(final=final, extension_runtime=20)
        self.assertEqual((row["flowtime"], row["movements"], row["scaled_objective"]), (4, 10, 146))
        self.assertEqual((row["final_flowtime"], row["final_movements"], row["final_scaled_objective"]), (4, 8, 144))
        self.assertEqual(row["flow_proven"], 1)
        self.assertEqual(row["lexicographic_proven"], 0)
        self.assertEqual(row["counterexample"], 0)

    def test_flow_improvement_is_a_counterexample_not_a_replacement_candidate(self):
        final = self.solution(status_name="INTERRUPTED", flowtime=3, makespan=3, movements=12,
                              scaled_best_bound=110, runtime=320)
        row, _, _ = self.run_fake(final=final, extension_runtime=20, stop_reason="FLOW_COUNTEREXAMPLE")
        self.assertEqual((row["flowtime"], row["movements"], row["makespan"]), (4, 10, 4))
        self.assertEqual((row["final_flowtime"], row["final_movements"], row["final_makespan"]), (3, 12, 3))
        self.assertEqual(row["counterexample"], 1)
        self.assertEqual(row["flow_proven"], 0)
        self.assertEqual(row["lexicographic_proven"], 0)
        self.assertEqual(row["search_stop_reason"], "FLOW_COUNTEREXAMPLE")
        self.assertEqual(row["counterexample_flowtime"], 3)
        self.assertEqual(row["counterexample_movements"], 12)

    def test_nonincumbent_lower_flow_witness_is_still_reported(self):
        witness = dict(has_solution=True, flowtime=3, movements=60, scaled_objective=162)
        row, _, _ = self.run_fake(flow_counterexample=witness, stop_reason="BETTER_FLOW_FOUND")
        self.assertEqual(row["final_flowtime"], 4)
        self.assertEqual(row["flowtime"], 4)
        self.assertEqual(row["counterexample"], 1)
        self.assertEqual(row["counterexample_flowtime"], 3)
        self.assertEqual(row["counterexample_scaled_objective"], 162)
        self.assertEqual(row["flow_proven"], 0)

    def test_later_exact_bound_proves_original_without_replacing_first_budget_proof_flag(self):
        final = self.solution(status_name="OPTIMAL", weighted_proven=True, scaled_best_bound=146)
        row, _, _ = self.run_fake(final=final)
        self.assertEqual(row["error"], "")
        self.assertEqual(row["weighted_proven"], 0)
        self.assertEqual(row["scaled_absolute_gap"], 46)
        self.assertEqual(row["proof_scaled_absolute_gap"], 0)
        self.assertEqual(row["lexicographic_proven"], 1)
        self.assertEqual(row["final_lexicographic_proven"], 1)

    def test_snapshot_provenance_and_estimated_phase_timing_are_preserved(self):
        initial = dict(self.initial, flow_proven=False, flow_proof_source="", incumbent_runtime=285.4,
                       snapshot_source="CALLBACK_OBSERVATIONS", cpu_time=301.5, cpu_time_is_estimate=True,
                       statistics_checkpoint_runtime=299.99)
        row, _, _ = self.run_fake(initial=initial, phase_transition_delay=0.01, phase1_snapshot_missing=False)
        self.assertEqual(row["weighted_cpu_time"], 301.5)
        self.assertEqual(row["weighted_cpu_time_is_estimate"], 1)
        self.assertEqual(row["phase1_incumbent_runtime"], 285.4)
        self.assertEqual(row["phase1_bound_checkpoint_runtime"], 299.9)
        self.assertEqual(row["phase1_statistics_checkpoint_runtime"], 299.99)
        self.assertEqual(row["phase1_snapshot_source"], "CALLBACK_OBSERVATIONS")
        self.assertEqual(row["phase_transition_delay"], 0.01)

    def test_distance_bound_does_not_require_a_weighted_bound(self):
        initial = self.solution(status_name="BUDGET_REACHED", flowtime=2, makespan=2, scaled_best_bound=None)
        final = dict(initial, status_name="INTERRUPTED")
        row, _, _ = self.run_fake(initial=initial, final=final, extension_used=False, extension_runtime=0)
        self.assertEqual(row["flow_proof_source"], "distance_bound")
        self.assertEqual(row["flow_lower_bound"], 2)
        self.assertEqual(row["flow_proven"], 1)
        self.assertEqual(row["lexicographic_proven"], 0)

    def test_unproved_extension_reports_total_budget_without_false_proof(self):
        row, solver, calls = self.run_fake()
        self.assertEqual(row["error"], "")
        self.assertEqual(row["total_time_limit"], 371)
        self.assertEqual(row["extension_used"], 1)
        self.assertEqual(row["weighted_runtime"], 300)
        self.assertEqual(row["final_runtime"], 371)
        self.assertEqual(row["final_status"], "TIME_LIMIT")
        self.assertEqual(row["Solver Status"], "BUDGET_REACHED")
        self.assertEqual(row["flow_proven"], 0)
        self.assertEqual(row["flow_proof_source"], "")
        self.assertEqual(row["lexicographic_proven"], 0)
        self.assertEqual(len(calls), 1)
        self.assertEqual(len(solver.solves), 1)

    def test_no_phase_one_incumbent_cannot_be_replaced_by_a_late_one(self):
        initial = self.solution(has_solution=False, status_name="BUDGET_REACHED")
        row, _, _ = self.run_fake(initial=initial, stop_reason="NO_PHASE1_SOLUTION")
        self.assertEqual(row["has_solution"], 0)
        self.assertIsNone(row["flowtime"])
        self.assertEqual(row["final_has_solution"], 1)
        self.assertEqual(row["final_flowtime"], 4)
        self.assertEqual(row["gap_certificate_reason"], "NO_WEIGHTED_SOLUTION")
        self.assertEqual(row["flow_proven"], 0)
        self.assertEqual(row["counterexample"], 0)

    def test_multi_target_horizon_and_coefficient_are_global_for_both_formulations(self):
        self.targets = [(2, 0), (1, 0)]
        self.args.loads = 2
        self.trace = (6, 10, 20, [], [{"escort": 1}], [{"target": 1}])
        initial = self.solution(weight_scale=61, flowtime=8, movements=12, scaled_best_bound=0)
        for formulation in ("escortflow", "loadflow"):
            with self.subTest(formulation=formulation):
                self.args.formulation = formulation
                row, solver, calls = self.run_fake(initial=initial, final=initial)
                self.assertEqual(row["error"], "")
                self.assertEqual(row["flow_weight"], 61)
                self.assertEqual(row["safe_movement_lower_bound"], 3)
                self.assertEqual(row["safe_movement_bound"], 63)
                self.assertEqual(row["safe_flow_horizon"], 9)
                self.assertEqual(row["weighted_physical_horizon"], 9)
                self.assertEqual(row["weighted_global_scope"], 1)
                # Both models can retrieve through period 9. EF's arrival
                # index is offset by one relative to LF's retrieval index.
                expected_index = 8 if formulation == "escortflow" else 9
                self.assertEqual(row["weighted_horizon"], expected_index)
                self.assertEqual(solver.starts[0][2], expected_index)
                self.assertEqual(solver.solves[0][2], expected_index)
                self.assertEqual(calls, [dict(flow_weight=61)])

    def test_common_horizon_preserves_complete_greedy_start(self):
        for formulation in ("escortflow", "loadflow"):
            with self.subTest(formulation=formulation):
                self.args.formulation = formulation
                row, solver, _ = self.run_fake()
                self.assertEqual(row["weighted_physical_horizon"], 6)
                expected_index = 5 if formulation == "escortflow" else 6
                self.assertEqual(row["weighted_horizon"], expected_index)
                self.assertEqual(solver.starts[0][2], expected_index)

    def test_invalid_weighted_proof_flag_cannot_prove_global_optimality(self):
        changes = (dict(status_name="NUMERIC"), dict(scaled_best_bound=145),
                   dict(scaled_best_bound=147), dict(weight_scale=100),
                   dict(bound_consistent=False), dict(scaled_objective=147))
        for changed in changes:
            with self.subTest(changed=changed):
                exact = self.solution(weighted_proven=True, status_name="OPTIMAL", scaled_best_bound=146)
                initial = dict(exact, **changed)
                row, _, _ = self.run_fake(initial=initial, final=exact)
                self.assertIn("Invalid weighted optimality certificate", row["error"])
                self.assertEqual(row["flow_proven"], 0)
                self.assertEqual(row["lexicographic_proven"], 0)

    def test_proven_result_that_contradicts_greedy_is_rejected(self):
        exact = self.solution(flowtime=6, weighted_proven=True, status_name="OPTIMAL", scaled_best_bound=214)
        row, _, _ = self.run_fake(initial=exact, final=exact)
        self.assertIn("contradicts", row["error"])
        self.assertEqual(row["flow_proven"], 0)
        self.assertEqual(row["lexicographic_proven"], 0)

    def test_minimum_safe_coefficient_accepted_but_one_below_rejected(self):
        exact = self.solution(weighted_proven=True, status_name="OPTIMAL", scaled_best_bound=146)
        runner._check_weighted_proof(exact, 34, 5, 192, 35, 2, True)
        unsafe = self.solution(weight_scale=33, weighted_proven=True, status_name="OPTIMAL",
                               scaled_best_bound=142)
        with self.assertRaisesRegex(ValueError, "safe-weight scope"):
            runner._check_weighted_proof(unsafe, 33, 5, 187, 35, 2, True)

    def test_weighted_proof_rejects_flow_or_movements_below_distance_bound(self):
        for changes in ({"flowtime":1}, {"movements":1}, {"flowtime":None},
                        {"movements":True}, {"has_solution":False}):
            with self.subTest(changes=changes):
                exact = self.solution(weighted_proven=True, status_name="OPTIMAL", scaled_best_bound=146)
                exact.update(changes)
                with self.assertRaisesRegex(ValueError, "Invalid weighted optimality certificate"):
                    runner._check_weighted_proof(exact, 34, 5, 192, 35, 2, True)

    def test_weighted_proof_rejects_invalid_movement_bounds(self):
        exact = self.solution(weighted_proven=True, status_name="OPTIMAL", scaled_best_bound=146)
        for upper, lower in ((1,2), (35,-1), (35,2.5), (float("nan"),2)):
            with self.subTest(upper=upper, lower=lower), \
                    self.assertRaisesRegex(ValueError, "Invalid safe-weight proof parameters"):
                runner._check_weighted_proof(exact, 34, 5, 192, upper, lower, True)

    def test_impossible_greedy_movement_count_is_rejected_before_solving(self):
        self.trace = (5, 5, 1, [], [{"escort":1}], [{"target":1}])
        row, solver, calls = self.run_fake()
        self.assertIn("Greedy movements violate", row["error"])
        self.assertEqual(calls, [])
        self.assertEqual(solver.solves, [])

    def test_multiple_optimization_calls_are_rejected(self):
        row, _, _ = self.run_fake(optimization_calls=2)
        self.assertIn("exactly one optimization call", row["error"])
        self.assertEqual(row["flow_proven"], 0)

    def test_configs_preserve_weighted_objective_focus_and_request_extension(self):
        classes = (("escortflow", "escort_flow_static_gurobi.StaticEscortFlowGurobiSolver"),
                   ("loadflow", "load_flow_static_gurobi.LoadFlowStaticGurobiSolver"))
        for formulation, constructor in classes:
            self.args.formulation = formulation
            with self.subTest(formulation=formulation), patch(constructor) as make:
                runner.make_solver(self.args, flow_weight=34)
                config = make.call_args.args[0]
                self.assertEqual(config.weight_scale, 34)
                self.assertEqual(config.gamma, 1 / 34)
                self.assertEqual(config.beta, 1)
                self.assertEqual(config.mip_focus, 0)
                self.assertEqual(config.objective_mode, "weighted_integer")
                self.assertEqual(config.time_limit, 300)
                self.assertEqual(config.flow_proof_extension_time_limit, 71)
                self.assertTrue(config.stop_on_flow_proof)
                self.assertFalse(config.stop_at_flow_proof)
                self.assertFalse(hasattr(config, "flow_proof_check_seconds"))
                self.assertFalse(hasattr(config, "flow_proof_check_nodes"))
                self.assertEqual(config.threads, 16)

    def test_early_stop_option_is_forwarded_to_both_solver_configs(self):
        self.args.stop_at_flow_proof = True
        for formulation, constructor in (
                ("escortflow", "escort_flow_static_gurobi.StaticEscortFlowGurobiSolver"),
                ("loadflow", "load_flow_static_gurobi.LoadFlowStaticGurobiSolver")):
            with self.subTest(formulation=formulation), patch(constructor) as make:
                self.args.formulation = formulation
                runner.make_solver(self.args, flow_weight=34)
                config = make.call_args.args[0]
                self.assertTrue(config.stop_at_flow_proof)
                self.assertTrue(config.stop_on_flow_proof)
                self.assertEqual(config.flow_proof_extension_time_limit, 71)

    def test_merge_validates_coverage_coefficients_and_budgets(self):
        row, _, _ = self.run_fake()
        source = self.args.output

        def write(row):
            with source.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=runner.FIELDNAMES)
                writer.writeheader()
                writer.writerow(row)

        write(row)
        destination = source.with_name("merged.csv")
        with redirect_stdout(io.StringIO()):
            runner.merge_batch(source, destination, "7", "2", 300, 71)
        with destination.open(newline="") as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 1)
        with self.assertRaisesRegex(ValueError, "Missing or duplicate"):
            runner.merge_batch(source, destination, "7-8", "2", 300, 71)
        with self.assertRaisesRegex(ValueError, "Unexpected extension_time_limit"):
            runner.merge_batch(source, destination, "7", "2", 300, 300)
        write(dict(row, flow_weight=100))
        with self.assertRaisesRegex(ValueError, "safe objective"):
            runner.merge_batch(source, destination, "7", "2", 300, 71)
        write(dict(row, retrieval_mode="continue"))
        with self.assertRaisesRegex(ValueError, "mixed retrieval modes"):
            runner.merge_batch(source, destination, "7", "2", 300, 71)
        write(dict(row, flow_proof_check_mode="time_or_nodes"))
        with self.assertRaisesRegex(ValueError, "safe objective"):
            runner.merge_batch(source, destination, "7", "2", 300, 71)

    def test_merge_rejects_old_protocol_schema_and_coefficient_on_both_sides(self):
        row, _, _ = self.run_fake()
        source = self.args.output
        destination = source.with_name("merged.csv")

        def write(path, value, fields=runner.FIELDNAMES):
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow({key:value[key] for key in fields})

        old_fields = [key for key in runner.FIELDNAMES if key != "safe_movement_lower_bound"]
        incompatible = ((dict(row, protocol="safe_integer_flow_timing_v4"), runner.FIELDNAMES, "protocol"),
                        (dict(row, protocol="safe_integer_flow_timing_v5"), runner.FIELDNAMES, "protocol"),
                        (row, old_fields, "schema"),
                        (dict(row, flow_weight=36, movement_weight=1/36), runner.FIELDNAMES, "safe objective"))
        for side in ("source", "destination"):
            for invalid, fields, message in incompatible:
                with self.subTest(side=side, message=message):
                    write(source, row)
                    write(destination, dict(row, seed=8))
                    write(source if side == "source" else destination, invalid, fields)
                    before = destination.read_bytes()
                    with self.assertRaisesRegex(ValueError, message):
                        runner.merge_batch(source, destination, "7", "2", 300, 71)
                    self.assertEqual(destination.read_bytes(), before)

    def test_merge_rejects_invalid_distance_bounds_and_duplicate_appends(self):
        row, _, _ = self.run_fake()
        source = self.args.output
        destination = source.with_name("merged.csv")

        def write(value):
            with source.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=runner.FIELDNAMES)
                writer.writeheader()
                writer.writerow(value)

        for invalid in (dict(row, safe_movement_lower_bound=-1),
                        dict(row, safe_movement_lower_bound=36),
                        dict(row, flowtime=1), dict(row, movements=1),
                        dict(row, final_flowtime=1), dict(row, final_movements=1)):
            with self.subTest(invalid=invalid):
                write(invalid)
                with self.assertRaises(ValueError):
                    runner.merge_batch(source, destination, "7", "2", 300, 71)
                self.assertFalse(destination.exists())
        write(row)
        with redirect_stdout(io.StringIO()):
            runner.merge_batch(source, destination, "7", "2", 300, 71)
        before = destination.read_bytes()
        with self.assertRaisesRegex(ValueError, "duplicate instance"):
            runner.merge_batch(source, destination, "7", "2", 300, 71)
        self.assertEqual(destination.read_bytes(), before)

    def test_merge_rejects_mixed_or_unexpected_flow_stop_settings(self):
        row, _, _ = self.run_fake()
        source = self.args.output
        destination = source.with_name("merged.csv")

        def write(path, records):
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=runner.FIELDNAMES)
                writer.writeheader()
                writer.writerows(records)

        self.assertEqual(row["stop_at_flow_proof"], 0)
        write(source, [dict(row, stop_at_flow_proof=1)])
        with self.assertRaisesRegex(ValueError, "Unexpected stop_at_flow_proof"):
            runner.merge_batch(source, destination, "7", "2", 300, 71, False)
        self.assertFalse(destination.exists())
        with redirect_stdout(io.StringIO()):
            runner.merge_batch(source, destination, "7", "2", 300, 71, True)
        before = destination.read_bytes()
        write(source, [dict(row, seed=8)])
        with self.assertRaisesRegex(ValueError, "mixed flow-proof stopping"):
            runner.merge_batch(source, destination, "8", "2", 300, 71)
        self.assertEqual(destination.read_bytes(), before)
        write(source, [row, dict(row, seed=8, stop_at_flow_proof=1)])
        with self.assertRaisesRegex(ValueError, "mixed flow-proof stopping"):
            runner.merge_batch(source, destination, "7-8", "2", 300, 71)
        write(source, [dict(row, stop_at_flow_proof=2)])
        with self.assertRaisesRegex(ValueError, "Invalid stop_at_flow_proof"):
            runner.merge_batch(source, destination, "7", "2", 300, 71)

    def test_wrapper_dry_run_creates_no_results_and_existing_directory_is_protected(self):
        script = Path(runner.__file__).with_name("RunTable2SafeWeighted.sh")
        results = Path(self.temporary.name) / "new_results"
        dry = subprocess.run(["bash", str(script), "--dry-run", "--output-dir", str(results)],
                             capture_output=True, text=True, check=True)
        self.assertEqual(dry.stdout.count("RunSafeWeightedStatic.py"), 16)
        self.assertNotIn("--stop-at-flow-proof", dry.stdout)
        early = subprocess.run(["bash", str(script), "--dry-run", "--stop-at-flow-proof",
                                "--output-dir", str(results)], capture_output=True, text=True, check=True)
        self.assertEqual(early.stdout.count("--stop-at-flow-proof"), 16)
        self.assertFalse(results.exists())
        results.mkdir()
        sentinel = results / "existing.csv"
        sentinel.write_text("preserve\n")
        attempted = subprocess.run(["bash", str(script), "--python", sys.executable,
                                    "--output-dir", str(results)], capture_output=True, text=True)
        self.assertNotEqual(attempted.returncode, 0)
        self.assertIn("Output directory already exists", attempted.stderr)
        self.assertEqual(sentinel.read_text(), "preserve\n")
        self.assertEqual(list(results.iterdir()), [sentinel])

    def test_cli_defaults_alias_and_existing_result_protection(self):
        arguments = ["--formulation", "escortflow", "-x", "3", "-y", "3", "-O", "0", "0",
                     "-e", "2", "-f", str(self.args.output)]
        args = runner.parse_args(arguments)
        self.assertEqual((args.threads, args.weighted_time_limit, args.extension_time_limit), (16, 300, 300))
        self.assertFalse(args.stop_at_flow_proof)
        self.assertTrue(runner.parse_args(arguments + ["--stop-at-flow-proof"]).stop_at_flow_proof)
        with self.assertRaises(SystemExit), patch("sys.stderr", new=io.StringIO()):
            runner.parse_args(arguments + ["--stop-at-flow-proof", "--lp"])
        self.assertFalse(hasattr(args, "flow_proof_check_seconds"))
        self.assertFalse(hasattr(args, "flow_proof_check_nodes"))
        alias = runner.parse_args(arguments + ["--certification-time-limit", "71"])
        self.assertEqual(alias.extension_time_limit, 71)
        self.assertFalse(hasattr(args, "certification_gap_threshold"))
        self.args.output.write_text("preserve\n")
        with self.assertRaises(SystemExit), patch("sys.stderr", new=io.StringIO()):
            runner.parse_args(arguments)
        self.assertEqual(self.args.output.read_text(), "preserve\n")

    def test_first_flow_proof_timing_remains_distinct_from_total_runtime(self):
        row, _, _ = self.run_fake(
            first_flow_proof_runtime=12.5, first_flow_proof_cpu_time=13.2,
            first_flow_proof_flowtime=4, first_flow_proof_source="weighted_gap",
            first_flow_proof_method="weighted_gap",
            first_flow_proof_node_count=24, first_flow_proof_scaled_bound=124,
            first_flow_proof_work=2.5)
        self.assertEqual(row["first_flow_proof_runtime"], 12.5)
        self.assertEqual(row["first_flow_proof_cpu_time"], 13.2)
        self.assertEqual(row["first_flow_proof_flowtime"], 4)
        self.assertEqual(row["first_flow_proof_source"], "weighted_gap")
        self.assertEqual(row["first_flow_proof_method"], "weighted_gap")
        self.assertEqual(row["final_runtime"], 371)
        self.assertEqual(row["weighted_runtime"], 300)

    def test_unobserved_flow_proof_time_is_missing(self):
        row, _, _ = self.run_fake()
        self.assertIsNone(row["first_flow_proof_runtime"])
        self.assertIsNone(row["first_flow_proof_cpu_time"])

    def test_contradictory_recorded_proof_clears_claims_but_keeps_metrics(self):
        exact = self.solution(weighted_proven=True, status_name="OPTIMAL", scaled_best_bound=146, runtime=2)
        row, _, _ = self.run_fake(initial=exact, final=exact,
                                  first_flow_proof_runtime=1, first_flow_proof_flowtime=3)
        self.assertIn("contradicts", row["error"])
        self.assertEqual(row["flowtime"], 4)
        self.assertEqual(row["final_runtime"], 2)
        for key in ["weighted_proven", "final_weighted_proven", "flow_proven", "lexicographic_proven",
                    "final_flow_proven", "final_lexicographic_proven"]:
            self.assertEqual(row[key], 0)
        self.assertIsNone(row["first_flow_proof_runtime"])

    def test_invalidated_proof_is_a_recorded_error_with_no_certificates(self):
        row, _, _ = self.run_fake(flow_proof_invalidated=True,
                                  flow_proof_invalidation_reason="BETTER_FLOW_AFTER_RECORDED_PROOF")
        self.assertIn("Flow proof invalidated", row["error"])
        self.assertEqual(row["flow_proof_invalidated"], 1)
        self.assertEqual(row["final_runtime"], 371)
        self.assertEqual(row["flow_proven"], 0)

    def test_cli_accepts_zero_extension_and_rejects_obsolete_check_intervals(self):
        arguments = ["--formulation", "loadflow", "-x", "3", "-y", "3", "-O", "0", "0",
                     "-e", "2", "-f", str(self.args.output)]
        args = runner.parse_args(arguments + ["--extension-time-limit", "0"])
        self.assertEqual(args.extension_time_limit, 0)
        for option, value in [("--extension-time-limit", "-1"), ("--extension-time-limit", "nan"),
                              ("--flow-proof-check-seconds", "0"), ("--flow-proof-check-nodes", "0"),
                              ("--flow-proof-check-nodes", "1.5")]:
            with self.subTest(option=option, value=value), self.assertRaises(SystemExit), \
                    patch("sys.stderr", new=io.StringIO()):
                runner.parse_args(arguments + [option, value])

    def test_report_identifies_bound_or_flow_change_checks(self):
        row, _, _ = self.run_fake()
        self.assertEqual(row["flow_proof_check_mode"], "bound_or_flow_change")
        self.assertNotIn("flow_proof_check_seconds", row)
        self.assertNotIn("flow_proof_check_nodes", row)


if __name__ == "__main__":
    unittest.main()
