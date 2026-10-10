"""Integrated LP bounds preserve integer measurements and survive interruption."""

import csv
from contextlib import redirect_stdout
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import RunSafeWeightedStatic as integer
import RunStaticLP as replay
import RunStaticCampaign as campaign
import static_integrated_lp as integrated

HERE = Path(__file__).resolve().parent


def example():
    row = dict.fromkeys(integer.FIELDNAMES, '')
    row.update({'Lx x Ly':'2x2', 'formulation':'escortflow', 'retrieval_mode':'leave',
        'movement_mode':'BM', 'IOs':'[(0, 0)]', 'Target Loads':'[(1, 0)]', 'Escorts':'[(1, 1)]',
        '#Loads':1, '# Escorts':1, 'seed':7, 'flow_weight':5, 'movement_weight':0.2,
        'weighted_horizon':1, 'weighted_physical_horizon':2, 'greedy_flowtime':1,
        'greedy_movements':1, 'has_solution':1, 'flowtime':1, 'movements':1,
        'final_has_solution':1, 'final_flowtime':1, 'final_movements':1,
        'weighted_runtime':300, 'weighted_cpu_time':301, 'final_runtime':600,
        'first_flow_proof_runtime':12, 'total_wall_time':605, 'error':'', 'lp_requested':0, 'lp_status':'NOT_RUN'})
    return row


def result(row):
    problem = replay.parse_instance(row)
    value = {key: problem[key] for key in ('layout','escorts','loads','seed','method','flow_weight',
                'horizon','physical_horizon','retrieval_mode','problem_sha256')}
    value.update(lp_objective=1.2, lp_flow=1, lp_movements=1, status='OPTIMAL', elapsed_seconds=0.1,
                 source_sha256=integrated.source_fingerprint(HERE), solver_version='13.0.3', algorithm='automatic')
    return value


class IntegratedLPTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.args = SimpleNamespace(lp_threads=1, lp_time_limit=300, lp_retry_time_limit=600)

    def write(self, path, rows):
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)

    def test_attached_bound_uses_exact_problem_and_preserves_every_integer_field(self):
        row = example()
        before = {key:value for key,value in row.items() if key not in integrated.LP_FIELDS}
        with patch.object(replay, 'solve_instance', return_value=result(row)) as solve, redirect_stdout(io.StringIO()):
            integrated.attach_lp(row, self.args, HERE)
        self.assertEqual(solve.call_args.args[0], replay.parse_instance(row))
        self.assertEqual(row['lp_relaxation_lower_bound'], 1.2)
        self.assertEqual(integrated.embedded_value(row)['lp_objective'], 1.2)
        self.assertEqual(before, {key:value for key,value in row.items() if key not in integrated.LP_FIELDS})

    def test_nonoptimal_or_invalid_lp_cannot_fill_lower_bound_column(self):
        bad = result(example()); bad.update(lp_objective=2, lp_flow=2, lp_movements=0)
        for outcome in (replay.LPNotOptimalError('LP status TIME_LIMIT', 'TIME_LIMIT'), bad):
            row = example()
            options = {'side_effect':outcome} if isinstance(outcome, Exception) else {'return_value':outcome}
            with patch.object(replay, 'solve_instance', **options), redirect_stdout(io.StringIO()):
                integrated.attach_lp(row, self.args, HERE)
            self.assertIn(row['lp_status'], ('TIME_LIMIT','ERROR'))
            self.assertEqual(row['lp_relaxation_lower_bound'], '')
            self.assertTrue(row['lp_error'])
            self.assertEqual(row['flowtime'], 1)

    def test_reuse_requires_matching_coordinates_weight_horizon_and_sources(self):
        row = example(); row.update(integrated.result_fields(result(row)))
        for change in ({'flow_weight':6,'movement_weight':1/6}, {'weighted_horizon':2,'weighted_physical_horizon':3},
                       {'seed':8}, {'lp_source_sha256':'changed'}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                integrated.embedded_value(dict(row, **change), integrated.source_fingerprint(HERE))

    def test_replay_exports_embedded_bounds_without_any_solver_work(self):
        row = example(); row.update(integrated.result_fields(result(row)))
        path = self.root/'integer.csv'; self.write(path,[row])
        output = self.root/'lp.csv'
        args = replay.parse_args(['--input',str(path),'-f',str(output),'--source-dir',str(HERE)])
        with patch.object(replay.concurrent.futures, 'ProcessPoolExecutor', side_effect=AssertionError('unexpected solve')), \
             redirect_stdout(io.StringIO()) as log:
            self.assertEqual(replay.run(args),0)
        self.assertIn('0 to solve', log.getvalue())
        with output.open(newline='') as handle:
            saved = next(csv.DictReader(handle))
        self.assertEqual(saved['status'],'OPTIMAL')
        self.assertEqual(float(saved['lp_objective']),1.2)

    def test_retry_backfill_changes_only_lp_fields_and_keeps_historical_csv_read_only(self):
        row = example(); row.update(lp_requested=1,lp_status='TIME_LIMIT',lp_error='time limit')
        path = self.root/'integer.csv'; self.write(path,[row])
        old = self.root/'old.csv'
        self.write(old,[{key:value for key,value in row.items() if key not in integrated.LP_FIELDS}])
        old_bytes = old.read_bytes()
        with path.open(newline='') as handle:
            before = next(csv.DictReader(handle))
        value = result(row)
        self.assertEqual(integrated.fill_inputs([path,old],{value['problem_sha256']:value},value['source_sha256'],
                         dict(threads=1,time_limit=300,retry_time_limit=600)),1)
        with path.open(newline='') as handle:
            after = next(csv.DictReader(handle))
        self.assertEqual({k:v for k,v in before.items() if k not in integrated.LP_FIELDS},
                         {k:v for k,v in after.items() if k not in integrated.LP_FIELDS})
        self.assertEqual(after['lp_status'],'OPTIMAL')
        self.assertEqual(old.read_bytes(),old_bytes)

    def test_integer_checkpoint_is_saved_before_lp_interruption(self):
        path = self.root/'interrupted.csv'
        command = ['--formulation','escortflow','-x','2','-y','2','-O','0','0','-e','1',
                   '-r','7','--with-lp','-f',str(path)]
        with patch.object(integer,'run_instance',return_value=example()), \
             patch.object(integer,'attach_lp',side_effect=KeyboardInterrupt), self.assertRaises(KeyboardInterrupt):
            integer.main(command)
        with path.open(newline='') as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['flowtime'],'1')
        self.assertEqual(rows[0]['weighted_runtime'],'300')
        self.assertEqual(rows[0]['lp_status'],'PENDING')

    def test_lp_failure_keeps_all_integer_rows_and_reports_separate_exit_code(self):
        path = self.root/'partial.csv'
        command = ['--formulation','escortflow','-x','2','-y','2','-O','0','0','-e','1',
                   '-r','7-8','--with-lp','-f',str(path)]
        def fake_integer(args,seed,count):
            return dict(example(),seed=seed)
        def fail(row,args,source):
            row.update(lp_status='ERROR',lp_error='failure',lp_elapsed_seconds=1)
        with patch.object(integer,'run_instance',side_effect=fake_integer),patch.object(integer,'attach_lp',side_effect=fail):
            self.assertEqual(integer.main(command),3)
        with path.open(newline='') as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([r['seed'] for r in rows],['7','8'])
        self.assertTrue(all(r['lp_relaxation_lower_bound']=='' and r['flowtime']=='1' for r in rows))


@unittest.skipUnless(importlib.util.find_spec('gurobipy') and importlib.util.find_spec('numpy'), 'solver environment unavailable')
class IntegratedLPSolverTests(unittest.TestCase):
    def test_campaign_reuses_and_recovers_lp_bounds_from_frozen_sources(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/'campaign'
            settings = dict(campaign.CAMPAIGNS['target_counts'], loads=(2,), escorts_by_loads={2:(2,)})
            with patch.object(campaign,'LAYOUTS',((3,2,2,(0,0)),)), \
                 patch.dict(campaign.CAMPAIGNS,{'target_counts':settings}),redirect_stdout(io.StringIO()):
                self.assertEqual(campaign.main('target_counts',['--layouts','3x2','--seeds','1-2',
                    '--threads','1','--weighted-time-limit','10','--extension-time-limit','0',
                    '--lp-threads','1','--lp-time-limit','10','--lp-retry-time-limit','10',
                    '--output-dir',str(root)]),0)
            paths = list(root.glob('table2b_targets_*.csv')) + list((root/'parts').glob('*.csv'))
            before = {}
            for path in paths:
                with path.open(newline='') as handle:
                    rows = list(csv.DictReader(handle))
                self.assertTrue(all(row['lp_status']=='OPTIMAL' for row in rows))
                before[path] = [{k:v for k,v in row.items() if k not in integrated.LP_FIELDS} for row in rows]
            collected = root/'lp_results.csv'
            with collected.open(newline='') as handle:
                results = list(csv.DictReader(handle))
            self.assertEqual(len(results),4)
            missing = results.pop(0)['problem_sha256']
            with collected.open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=replay.FIELDS)
                writer.writeheader();writer.writerows(results)
            for path in paths:
                with path.open(newline='') as handle:
                    reader=csv.DictReader(handle);fields=reader.fieldnames;rows=list(reader)
                for row in rows:
                    if row['lp_problem_sha256']==missing:
                        row.update({key:'' for key in integrated.VALUE_FIELDS})
                        row.update(lp_status='TIME_LIMIT',lp_error='simulated interrupted LP')
                with path.open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=fields)
                    writer.writeheader();writer.writerows(rows)
            recovered=subprocess.run(['bash',str(root/'lp_commands.sh')],capture_output=True,text=True,timeout=60)
            self.assertEqual(recovered.returncode,0,recovered.stdout+recovered.stderr)
            self.assertIn('1 to solve',recovered.stdout)
            for path in paths:
                with path.open(newline='') as handle:
                    rows=list(csv.DictReader(handle))
                self.assertTrue(all(row['lp_status']=='OPTIMAL' for row in rows))
                self.assertEqual(before[path],[{k:v for k,v in row.items() if k not in integrated.LP_FIELDS} for row in rows])

    def test_both_models_and_modes_save_optimal_bounds_and_replay_without_resolving(self):
        with tempfile.TemporaryDirectory() as temporary:
            for method in ('escortflow','loadflow'):
                for mode in ('leave','continue'):
                    path = Path(temporary)/(method+'_'+mode+'.csv')
                    command = [sys.executable,str(HERE/'RunSafeWeightedStatic.py'),'--formulation',method,
                         '-x','3','-y','2','-O','0','0','-l','2','-e','2','-r','1','-m',mode,
                         '--threads','1','--weighted-time-limit','10','--extension-time-limit','0',
                         '--with-lp','--lp-threads','1','--lp-time-limit','10','--lp-retry-time-limit','10','-f',str(path)]
                    completed = subprocess.run(command,capture_output=True,text=True,timeout=60)
                    self.assertEqual(completed.returncode,0,completed.stdout+completed.stderr)
                    with path.open(newline='') as handle:
                        row = next(csv.DictReader(handle))
                    value = integrated.embedded_value(row,integrated.source_fingerprint(HERE))
                    self.assertEqual(value['status'],'OPTIMAL')
                    self.assertLessEqual(float(value['lp_objective']),float(row['objective'])+1e-6)
                    collected = subprocess.run([sys.executable,str(HERE/'RunStaticLP.py'),'--input',str(path),
                        '--source-dir',str(HERE),'-f',str(path)+'.lp.csv'],capture_output=True,text=True,timeout=60)
                    self.assertEqual(collected.returncode,0,collected.stdout+collected.stderr)
                    self.assertIn('0 to solve',collected.stdout)


if __name__ == '__main__':
    unittest.main()
