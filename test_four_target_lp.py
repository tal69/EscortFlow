"""The independent Mac LP campaign must cover absent MIPs and preserve v4 inputs."""

import csv
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import RunFourTargetLP as runner
import RunStaticLP as lp


class FourTargetLPTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.argv = ['--input-dir', str(self.root/'missing_integer_results'),
                     '--output-dir', str(self.root/'lp'), '--layouts', '13x7',
                     '--escorts', '16', '--seeds', '1-2', '--threads', '1']
        self.args = runner.parse_args(self.argv)

    def test_no_arguments_select_exactly_all_eight_parts_and_2400_lps(self):
        args = runner.parse_args(['--threads', '1'])
        parts = list(runner.batches(args))
        self.assertEqual(len(parts), 8)
        self.assertEqual(len(args.seed_values)*len(args.escort_values)*len(parts), 2400)
        self.assertEqual(runner.selection(args)['protocol'], 'v4')
        self.assertEqual({b['method'] for b in parts}, {'escortflow', 'loadflow'})

    def test_dry_run_needs_no_result_files_and_writes_nothing(self):
        with patch.object(runner, 'environment', side_effect=AssertionError('solver checked')), \
                redirect_stdout(io.StringIO()) as output:
            self.assertEqual(runner.main(self.argv+['--dry-run']), 0)
        self.assertIn('Total: 4 LPs', output.getvalue())
        self.assertFalse(self.args.output_dir.exists())

    def prepare(self, references=None):
        saved = runner.freeze(self.args)
        with redirect_stdout(io.StringIO()):
            prepared, checks = runner.prepare(self.args, saved, references or {})
        return saved, prepared, checks

    def test_generation_includes_all_seeds_without_any_integer_csv(self):
        _, prepared, checks = self.prepare()
        self.assertEqual(len(prepared), 2)
        self.assertEqual([len(rows) for _, _, rows in prepared], [2, 2])
        self.assertTrue(all(c['integer_rows_verified'] == 0 for c in checks))
        for _, _, rows in prepared:
            self.assertEqual({r['seed'] for r in rows}, {1, 2})
            self.assertEqual({r['protocol'] for r in rows}, {'generated_lp_v4'})
        first = [lp.parse_instance(rows[0]) for _, _, rows in prepared]
        self.assertEqual({p['flow_weight'] for p in first}, {1351})
        self.assertEqual({p['physical_horizon'] for p in first}, {18, 19})

    def test_resume_reuses_generated_inputs_without_recomputing_greedy(self):
        self.prepare()
        saved = runner.freeze(self.args)
        with patch.object(runner, 'generated_row', side_effect=AssertionError('regenerated')), \
                redirect_stdout(io.StringIO()):
            _, checks = runner.prepare(self.args, saved, {})
        self.assertEqual(sum(c['generated'] for c in checks), 4)

    def test_available_integer_row_is_verified_but_does_not_filter_other_seeds(self):
        row = runner.generated_row('escortflow', 13, 7, [(6, 0)], 16, 4, 1, 'leave', 'v4')
        row['protocol'] = 'safe_integer_flow_timing_v4'
        self.args.input_dir.mkdir()
        runner.save_csv(self.args.input_dir/'table2b_escortflow_13x7.csv', [row], list(row))
        references, files = runner.read_references(self.args)
        _, prepared, checks = self.prepare(references)
        self.assertEqual(len(files), 1)
        self.assertEqual(sum(c['integer_rows_verified'] for c in checks), 1)
        self.assertEqual(sum(len(rows) for _, _, rows in prepared), 4)
        row['flow_weight'] += 1
        row['movement_weight'] = 1/row['flow_weight']
        changed = lp.parse_instance(row)
        with self.assertRaisesRegex(ValueError, 'flow_weight'):
            runner.compare_references(prepared[0][2], {lp.reference_key(changed): changed})

    def test_mismatched_seed_selection_and_modified_snapshots_are_rejected(self):
        _, prepared, _ = self.prepare()
        changed = runner.parse_args(self.argv+['--seeds', '1-3'])
        with self.assertRaisesRegex(ValueError, 'selection'):
            runner.freeze(changed)
        prepared[0][1].write_text(prepared[0][1].read_text()+'\n')
        with self.assertRaisesRegex(ValueError, 'input snapshot changed'):
            runner.freeze(self.args)

    def test_frozen_source_contains_all_lp_engine_dependencies_and_is_checked(self):
        self.prepare()
        source = self.args.output_dir/'source'
        self.assertTrue((source/'static_integrated_lp.py').is_file())
        self.assertTrue((source/'OneStepHeuristic_v2.py').is_file())
        module = source/'static_generated_lp.py'
        module.write_text(module.read_text()+'\n# altered\n')
        with self.assertRaisesRegex(ValueError, 'Frozen source changed'):
            runner.freeze(self.args)

    def test_unrecognized_output_folder_is_not_overwritten(self):
        self.args.output_dir.mkdir()
        marker = self.args.output_dir/'integer_results.csv'
        marker.write_text('preserve this\n')
        with self.assertRaisesRegex(ValueError, 'manifest'):
            runner.freeze(self.args)
        self.assertEqual(marker.read_text(), 'preserve this\n')

    def test_v5_reference_is_rejected_instead_of_silently_mixed(self):
        row = runner.generated_row('escortflow', 13, 7, [(6, 0)], 16, 4, 1, 'leave', 'v5')
        row['protocol'] = 'safe_integer_flow_timing_v5'
        self.args.input_dir.mkdir()
        runner.save_csv(self.args.input_dir/'table2b_escortflow_13x7.csv', [row], list(row))
        with self.assertRaisesRegex(ValueError, 'another protocol'):
            runner.read_references(self.args)

    def test_both_native_lps_match_archived_v4_values_and_resume_without_resolving(self):
        command = [sys.executable, '-u', str(runner.HERE/'RunFourTargetLP.py'),
                   *self.argv, '--seeds', '1', '--time-limit', '30']
        completed = subprocess.run(command, capture_output=True, text=True, timeout=90)
        self.assertEqual(completed.returncode, 0, completed.stdout+'\n'+completed.stderr)
        result = lp.read_reference(self.args.output_dir/'lp_results.csv')
        self.assertEqual(len(result), 2)
        for row in result.values():
            expected = 36.02673201574465 if row['method'] == 'escortflow' else 36.022434880688024
            self.assertAlmostEqual(float(row['lp_objective']), expected, places=6)
            self.assertEqual(row['status'], 'OPTIMAL')
        coverage = json.loads((self.args.output_dir/'coverage.json').read_text())
        self.assertTrue(coverage['complete'])
        hashes = {p.name: lp.sha256(p) for p in (self.args.output_dir/'parts').glob('*.csv')}
        resumed = subprocess.run(command, capture_output=True, text=True, timeout=90)
        self.assertEqual(resumed.returncode, 0, resumed.stdout+'\n'+resumed.stderr)
        self.assertEqual(resumed.stdout.count('0 to solve'), 2)
        self.assertEqual(hashes, {p.name: lp.sha256(p) for p in (self.args.output_dir/'parts').glob('*.csv')})


if __name__ == '__main__':
    unittest.main()
