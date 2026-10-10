"""Direct LP generation must match integer inputs and preserve reproducible resumes."""
import ast
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

import RunStaticLP as lp
import static_generated_lp as direct


class GeneratedLPTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.args = SimpleNamespace(Lx=3, Ly=2, output_cells=[0,0], load_num=2,
                                    escorts_range='2', reps_range='1-2', retrieval_mode='leave',
                                    csv=str(self.root/'lp.csv'), num_threads=1, time_limit=20,
                                    lp_protocol='v5', resume=False)

    def test_both_models_generate_identical_initial_states_and_current_physical_horizons(self):
        rows = [direct.generated_row(method,3,2,[(0,0)],2,2,7,'leave','v5') for method in ('escortflow','loadflow')]
        problems = [lp.parse_instance(r) for r in rows]
        for key in ('outputs','targets','escort_cells','flow_weight','physical_horizon'):
            self.assertEqual(problems[0][key],problems[1][key])
        self.assertEqual(problems[0]['horizon']+1,problems[1]['horizon'])

    def test_archived_four_target_coefficient_and_horizons_match_saved_inputs(self):
        for method in ('escortflow','loadflow'):
            path=direct.HERE/'Experiment Oct2026'/f'table2b_{method}_13x7.csv'
            if not path.exists():
                self.skipTest('Archived inputs are not included in this checkout')
            with path.open(newline='') as stream:
                saved=next(r for r in csv.DictReader(stream) if r['# Escorts']=='16' and r['seed']=='1')
            row=direct.generated_row(method,13,7,[(6,0)],16,4,1,'leave','v4')
            self.assertEqual(lp.parse_instance(row),lp.parse_instance(saved))

    def test_v4_and_v5_use_different_coefficients_without_changing_initial_state(self):
        old=direct.generated_row('escortflow',3,2,[(0,0)],2,2,7,'leave','v4')
        new=direct.generated_row('escortflow',3,2,[(0,0)],2,2,7,'leave','v5')
        distance=sum(min(abs(x-a)+abs(y-b) for a,b in [(0,0)]) for x,y in ast.literal_eval(old['Target Loads']))
        self.assertEqual(old['flow_weight']-new['flow_weight'],distance)
        self.assertEqual(old['Target Loads'],new['Target Loads'])

    def prepare(self):
        with patch.object(direct.subprocess,'call',return_value=0) as replay,redirect_stdout(io.StringIO()):
            self.assertEqual(direct.run_standard(self.args,'escortflow'),0)
        return replay.call_args.args[0]

    def test_generation_needs_no_integer_csv_and_does_not_solve_an_integer_model(self):
        arguments=self.prepare()
        self.assertIn('--input',arguments)
        with Path(self.args.csv+'.instances.csv').open(newline='') as stream:
            rows=list(csv.DictReader(stream))
        self.assertEqual(len(rows),2)
        self.assertEqual({r['seed'] for r in rows},{'1','2'})
        self.assertTrue(Path(self.args.csv+'.generation.json').exists())
        self.assertFalse(Path(self.args.csv).exists())

    def test_resume_reuses_snapshot_without_recomputing_the_greedy_reference(self):
        self.prepare()
        self.args.resume=True
        with patch.object(direct,'generated_row',side_effect=AssertionError('regenerated')), \
             patch.object(direct.subprocess,'call',return_value=0) as replay,redirect_stdout(io.StringIO()):
            self.assertEqual(direct.run_standard(self.args,'escortflow'),0)
        self.assertIn('--resume',replay.call_args.args[0])

    def test_changed_seed_range_or_input_snapshot_is_rejected(self):
        self.prepare()
        self.args.resume=True
        self.args.reps_range='1-3'
        with self.assertRaisesRegex(ValueError,'settings changed'):
            direct.run_standard(self.args,'escortflow')
        self.args.reps_range='1-2'
        inputs=Path(self.args.csv+'.instances.csv')
        inputs.write_text(inputs.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'snapshot changed'):
            direct.run_standard(self.args,'escortflow')

    def test_existing_integer_output_is_not_overwritten(self):
        output=Path(self.args.csv)
        output.write_text('integer result\n')
        with self.assertRaisesRegex(ValueError,'output already exists'):
            direct.run_standard(self.args,'escortflow')
        self.assertEqual(output.read_text(),'integer result\n')

    def test_modified_frozen_source_is_rejected_before_replay(self):
        self.prepare()
        self.args.resume=True
        path=Path(self.args.csv+'.source')/'static_weighted_certification.py'
        path.write_text(path.read_text()+'\n# changed\n')
        with self.assertRaisesRegex(ValueError,'Frozen LP source changed'):
            direct.run_standard(self.args,'escortflow')

    def test_duplicate_ranges_and_integer_phase_limits_are_rejected(self):
        for key,value in (('reps_range','1,1'),('escorts_range','2,2'),('phase1_time_limit',0),('work_limit',0),('time_limit',0)):
            args=SimpleNamespace(**vars(self.args))
            setattr(args,key,value)
            with self.subTest(key=key),self.assertRaises(ValueError):
                direct.run_standard(args,'escortflow')


class DirectLPCLIIntegrationTests(unittest.TestCase):
    def test_regular_runners_and_safe_runner_solve_the_same_tiny_continuous_models(self):
        with tempfile.TemporaryDirectory() as temporary:
            for mode in ('leave','continue'):
                outputs=[]
                for method,script in (('escortflow','EscortFlowStatic.py'),('loadflow','LoadFlowStatic.py'),
                                      ('escortflow','RunSafeWeightedStatic.py'),('loadflow','RunSafeWeightedStatic.py')):
                    output=Path(temporary)/(mode+'_'+method+'_'+script+'.csv')
                    command=[sys.executable,str(direct.HERE/script),'-x','3','-y','2','-O','0','0','-l','2','-e','2','-r','1',
                             '-m',mode,'--lp','-f',str(output)]
                    if script=='RunSafeWeightedStatic.py':
                        command+=['--formulation',method,'--threads','1','--weighted-time-limit','20']
                    else:
                        command+=['--num_threads','1','-t','20']
                    completed=subprocess.run(command,capture_output=True,text=True,timeout=60)
                    self.assertEqual(completed.returncode,0,completed.stdout+'\n'+completed.stderr)
                    with output.open(newline='') as stream:
                        row=next(csv.DictReader(stream))
                    self.assertEqual(row['status'],'OPTIMAL')
                    self.assertEqual(row['loads'],'2')
                    self.assertEqual(row['retrieval_mode'],mode)
                    outputs.append(row)
                    resumed=subprocess.run(command+['--resume'],capture_output=True,text=True,timeout=60)
                    self.assertEqual(resumed.returncode,0,resumed.stdout+'\n'+resumed.stderr)
                    self.assertIn('0 to solve',resumed.stdout)
                self.assertAlmostEqual(float(outputs[0]['lp_objective']),float(outputs[2]['lp_objective']),places=8)
                self.assertAlmostEqual(float(outputs[1]['lp_objective']),float(outputs[3]['lp_objective']),places=8)
                self.assertLessEqual(float(outputs[1]['lp_objective']),float(outputs[0]['lp_objective'])+1e-6)


if __name__=='__main__':
    unittest.main()
