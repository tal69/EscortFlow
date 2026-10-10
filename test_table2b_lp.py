"""Checks for Table 2(b) input fidelity, common references, and incomplete means."""

import csv
from contextlib import redirect_stdout
import io
from pathlib import Path
import tempfile
import subprocess
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import RunStaticLP as lp
import RunTable2bLP as runner


def integer_row(method='escortflow', seed=1):
    row = {'Lx x Ly': '13x7', 'formulation': method, 'retrieval_mode': 'leave',
           'IOs': '[(6, 0)]', 'Target Loads': '[(1, 0), (2, 0), (3, 0), (4, 0)]',
           'Escorts': '[(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (1, 1)]',
           '#Loads': '4', '# Escorts': '8', 'seed': str(seed), 'flow_weight': '1246',
           'movement_weight': str(1/1246), 'weighted_horizon': '15',
           'weighted_physical_horizon': '16' if method == 'escortflow' else '15',
           'protocol': 'safe_integer_flow_timing_v4', 'movement_mode': 'BM',
           'warmstart': '1', 'has_solution': '1', 'final_has_solution': '1',
           'weighted_global_scope': '1', 'error': '', 'safe_movement_bound': '1245',
           'safe_flow_horizon': '15', 'greedy_flowtime': '24', 'greedy_movements': '60',
           'flowtime': '20' if method == 'escortflow' else '21',
           'movements': '50' if method == 'escortflow' else '40',
           'final_flowtime': '19', 'final_movements': '49' if method == 'escortflow' else '45'}
    return row


def write_inputs(directory, rows):
    for method in runner.METHODS:
        selected = [r for r in rows if r['formulation'] == method]
        if selected:
            path = directory / f'table2b_{method}_13x7.csv'
            with path.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(selected[0]))
                writer.writeheader()
                writer.writerows(selected)


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


class Table2bLPTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.inputs = self.root / 'inputs'
        self.inputs.mkdir()
        self.output = self.root / 'lp_output'

    def args(self, seeds='1-100'):
        return runner.parse_args(['--input-dir', str(self.inputs), '--output-dir', str(self.output),
                                  '--layouts', '13x7', '--escorts', '8', '--seeds', seeds, '--threads', '1'])

    def test_archived_weights_and_different_physical_horizons_are_preserved(self):
        write_inputs(self.inputs, [integer_row(m) for m in runner.METHODS])
        _, problems, _, _ = runner.read_inputs(self.args())
        self.assertEqual({p['flow_weight'] for p in problems.values()}, {1246})
        self.assertEqual({p['physical_horizon'] for p in problems.values()}, {15, 16})

    def test_truncated_rows_are_not_silently_dropped(self):
        write_inputs(self.inputs, [integer_row()])
        path = self.inputs / 'table2b_escortflow_13x7.csv'
        path.write_text(path.read_text() + '13x7,escortflow\n')
        with self.assertRaisesRegex(ValueError, 'truncated CSV row'):
            runner.read_inputs(self.args())

    def test_duplicate_and_wrong_target_count_inputs_are_rejected(self):
        write_inputs(self.inputs, [integer_row(), integer_row()])
        with self.assertRaisesRegex(ValueError, 'Duplicate instance'):
            runner.read_inputs(self.args())
        row = integer_row()
        row['#Loads'], row['Target Loads'] = '1', '[(1, 0)]'
        write_inputs(self.inputs, [row])
        with self.assertRaisesRegex(ValueError, 'four-target'):
            runner.read_inputs(self.args())

    def test_v5_coefficient_is_not_converted_back_to_v4(self):
        row = integer_row()
        row.update(protocol='safe_integer_flow_timing_v5', flow_weight='1232', movement_weight=str(1/1232))
        write_inputs(self.inputs, [row])
        _, problems, _, _ = runner.read_inputs(self.args())
        self.assertEqual(next(iter(problems.values()))['flow_weight'], 1232)

    def test_dry_run_writes_nothing_and_starts_no_solver(self):
        write_inputs(self.inputs, [integer_row()])
        with patch.object(runner, 'automatic_threads', return_value=1), \
                patch.object(runner.subprocess, 'run', side_effect=AssertionError('solver started')), \
                redirect_stdout(io.StringIO()):
            self.assertEqual(runner.main(['--input-dir', str(self.inputs), '--output-dir', str(self.output),
                                         '--threads', '1', '--dry-run']), 0)
        self.assertFalse(self.output.exists())

    def prepare_results(self, seeds):
        write_inputs(self.inputs, [integer_row(m) for m in runner.METHODS])
        args = self.args(seeds)
        records, problems, _, _ = runner.read_inputs(args)
        source = self.output / 'source'
        source.mkdir(parents=True)
        for name in lp.MODEL_FILES:
            (source / name).write_text('frozen model ' + name)
        source_hash = lp.fingerprint({n: lp.sha256(source / n) for n in lp.MODEL_FILES})
        results = []
        for problem in problems.values():
            row = {f: problem[f] for f in ('layout', 'escorts', 'loads', 'seed', 'method', 'flow_weight',
                                           'horizon', 'physical_horizon', 'retrieval_mode', 'problem_sha256')}
            row.update(lp_objective=16+30/problem['flow_weight'], lp_flow=16, lp_movements=30,
                       status='OPTIMAL', elapsed_seconds=1, source_sha256=source_hash,
                       solver_version='test', algorithm='test')
            results.append(row)
        runner.save_csv(self.output / 'lp_results.csv', results, lp.FIELDS)
        return args, records, problems

    def test_common_reference_uses_one_lexicographic_pair_including_extensions(self):
        args, records, problems = self.prepare_results('1')
        with redirect_stdout(io.StringIO()):
            runner.export_gaps(self.output, records, problems, args)
        gaps = read_csv(self.output / 'lp_gaps.csv')
        self.assertEqual({(r['best_ft'], r['best_mv']) for r in gaps}, {('19', '45')})
        expected = 100*(19+45/1246-(16+30/1246))/(19+45/1246)
        self.assertTrue(all(abs(float(r['lp_gap_pct'])-expected) < 1e-10 for r in gaps))
        summary = read_csv(self.output / 'table2b_lp_summary.csv')
        self.assertTrue(all(r['complete_group'] == '1' for r in summary))

    def test_partial_groups_never_enter_the_100_instance_table_means(self):
        args, records, problems = self.prepare_results('1-100')
        with redirect_stdout(io.StringIO()):
            runner.export_gaps(self.output, records, problems, args)
        summary = read_csv(self.output / 'table2b_lp_summary.csv')
        self.assertTrue(all(r['complete_group'] == '0' and r['lp_gap_pct'] == '' for r in summary))
        self.assertTrue(all(r['expected_instances'] == '100' and r['lp_instances_complete'] == '1' for r in summary))

    def test_mismatched_cached_problem_cannot_enter_a_gap(self):
        args, records, problems = self.prepare_results('1')
        key = next(iter(problems))
        problems[key]['problem_sha256'] = 'changed'
        with self.assertRaisesRegex(ValueError, 'metadata differs'):
            runner.export_gaps(self.output, records, problems, args)

    def test_snapshot_round_trip_is_stable_on_mac_and_linux(self):
        text = runner.csv_text([integer_row()], list(integer_row()))
        path = self.root / 'snapshot.csv'
        path.write_text(text)
        self.assertEqual(path.read_text(), text)

    def test_environment_failure_identifies_interpreter_and_preserves_cause(self):
        failure = SimpleNamespace(returncode=1, stdout='', stderr='Python 3.9 is too old')
        with patch.object(runner.subprocess, 'run', return_value=failure):
            with self.assertRaisesRegex(ValueError, 'Python 3.9 is too old') as error:
                runner.check_environment()
        self.assertIn(sys.executable, str(error.exception))
        self.assertIn('RunTable2bLP.sh --setup', str(error.exception))

    def test_license_failure_is_distinct_from_version_failure(self):
        failure = SimpleNamespace(returncode=1, stdout='', stderr='Gurobi\'s size-limited license is active.')
        with patch.object(runner.subprocess, 'run', return_value=failure):
            with self.assertRaisesRegex(ValueError, 'size-limited license'):
                runner.check_environment()

    def test_environment_check_needs_no_integer_inputs_or_output_folder(self):
        with patch.object(runner, 'check_environment', return_value={'gurobi': '13.0.3'}), \
                patch.object(runner, 'read_inputs', side_effect=AssertionError('inputs read')), \
                redirect_stdout(io.StringIO()):
            self.assertEqual(runner.main(['--input-dir', str(self.inputs/'missing'),
                                         '--output-dir', str(self.output), '--threads', '1',
                                         '--check-environment']), 0)
        self.assertFalse(self.output.exists())

    def test_launcher_pins_python_with_spaces_independently_of_path(self):
        fake = self.root/'licensed python'
        fake.write_text('#!/bin/bash\n'
                        'if [[ "$1" == "-c" ]]; then printf "%s\\n" "$0"; exit 0; fi\n'
                        'if [[ "$3" == "--check-environment" ]]; then exit 0; fi\n'
                        'printf "%s\\n" "$@"\n')
        fake.chmod(0o755)
        import os
        environment = dict(os.environ, XDG_CONFIG_HOME=str(self.root/'config'))
        launcher = str(Path(runner.__file__).with_suffix('.sh'))
        setup = subprocess.run(['bash', launcher, '--set-python', str(fake)],
                               env=environment, capture_output=True, text=True)
        self.assertEqual(setup.returncode, 0, setup.stderr)
        saved = self.root/'config'/'escortflow'/'table2b-python'
        self.assertEqual(saved.read_text().strip(), str(fake))
        resumed = subprocess.run(['bash', launcher, '--dry-run', '--seeds', '1-2'],
                                 env=environment, capture_output=True, text=True)
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertIn('--dry-run\n--seeds\n1-2', resumed.stdout)

    def test_launcher_rejects_bad_environment_without_replacing_saved_python(self):
        fake = self.root/'bad python'
        fake.write_text('#!/bin/bash\nprintf "Wrong Python environment\\n" >&2\nexit 1\n')
        fake.chmod(0o755)
        config = self.root/'config'/'escortflow'
        config.mkdir(parents=True)
        saved = config/'table2b-python'
        saved.write_text('/previous/licensed/python\n')
        import os
        result = subprocess.run(['bash', str(Path(runner.__file__).with_suffix('.sh')),
                                 '--set-python', str(fake)],
                                env=dict(os.environ, XDG_CONFIG_HOME=str(config.parent)),
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(saved.read_text(), '/previous/licensed/python\n')


if __name__ == '__main__':
    unittest.main()
