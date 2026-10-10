"""Regression checks for first-output retrieval against the actual Gurobi model."""
import unittest
from unittest.mock import patch
import gurobipy as gp
from load_flow_static_gurobi import LoadFlowStaticGurobiConfig, LoadFlowStaticGurobiSolver
from escort_flow_static_gurobi import StaticGurobiConfig, StaticEscortFlowGurobiSolver

class OutputRuleTests(unittest.TestCase):
    def solver(self, width, height, outputs, **overrides):
        config=dict(Lx=width,Ly=height,output_cells=outputs,move_method='BM',
                    alpha=0,beta=1,gamma=.01,time_limit=30,threads=1)
        config.update(overrides)
        solver=LoadFlowStaticGurobiSolver(LoadFlowStaticGurobiConfig(**config))
        self.addCleanup(solver.close)
        return solver

    def fix_assignment(self, solver, targets, escorts, horizon, values):
        keys=[('x',move,t,k) for move in solver.network['moves']
              for t in range(horizon+1) for k in (1,2)]
        keys += [('q',output,t) for output in solver.output_cells for t in range(horizon+1)]
        keys += [('z',)]
        optimize=gp.Model.optimize
        def fixed(model,*args,**kwargs):
            model.update()
            variables=model.getVars()
            self.assertEqual(len(variables),len(keys))
            for key,var in zip(keys,variables): var.LB=var.UB=values.get(key,0)
            return optimize(model,*args,**kwargs)
        with patch.object(gp.Model,'optimize',fixed): return solver.solve(targets,escorts,horizon)

    def test_single_target_delay_rejected_in_integer_and_fractional_models(self):
        output=(0,0)
        for fraction in (1,.5):
            with self.subTest(fraction=fraction):
                solver=self.solver(2,2,(output,),lp=fraction!=1)
                values={('x',(1,0,0,0),0,1):1,('x',(0,0,0,0),1,1):fraction,
                        ('q',output,1):1-fraction,('q',output,2):fraction}
                for loc in ((0,1),(1,1)):
                    for t in range(3): values['x',loc+loc,t,2]=1
                result=self.fix_assignment(solver,[(1,0)],[output],2,values)
                self.assertEqual(result['status_name'],'INFEASIBLE')

    def test_multi_target_counterexample_has_same_corrected_optimum_as_escort_flow(self):
        outputs=((0,0),(1,0));targets=[(2,0),(3,0),(4,0)]
        solver=self.solver(5,1,outputs)
        result=solver.solve(targets,outputs,5)
        self.assertEqual(result['status_name'],'OPTIMAL')
        self.assertAlmostEqual(result['flowtime'],9)
        self.assertAlmostEqual(result['movements'],6)
        ef=StaticEscortFlowGurobiSolver(StaticGurobiConfig(Lx=5,Ly=1,output_cells=outputs,
            retrieval_mode='leave',beta=1,gamma=.01,time_limit=30,threads=1))
        self.addCleanup(ef.close)
        escort_result=ef.solve(targets,outputs,4)
        self.assertEqual(escort_result['status_name'],'OPTIMAL')
        self.assertAlmostEqual(escort_result['flowtime'],result['flowtime'])
        self.assertAlmostEqual(escort_result['movements'],result['movements'])

    def test_target_passage_through_output_is_rejected(self):
        outputs=((0,0),(1,0));targets=[(2,0),(3,0),(4,0)]
        values={}
        for origin in (2,3,4): values['x',(origin,0,origin-1,0),0,1]=1
        for origin in (1,2,3): values['x',(origin,0,origin-1,0),1,1]=1
        values['x',(2,0,2,0),2,1]=1
        values['x',(2,0,1,0),3,1]=1
        for output,t in (((0,0),2),((1,0),2),((1,0),4)): values['q',output,t]=1
        solver=self.solver(5,1,outputs)
        self.assertEqual(self.fix_assignment(solver,targets,outputs,5,values)['status_name'],'INFEASIBLE')

    def test_target_occupancy_identity_holds_without_an_added_lp_row(self):
        solver=self.solver(3,2,((0,0),),lp=True)
        horizon=4
        optimize=gp.Model.optimize
        residuals=[]
        def inspect_solution(model,*args,**kwargs):
            model.update()
            self.assertIsNone(model.getConstrByName('flow_time_occupancy_identity'))
            keys=[('x',move,t,k) for move in solver.network['moves']
                  for t in range(horizon+1) for k in (1,2)]
            keys += [('q',output,t) for output in solver.output_cells for t in range(horizon+1)]
            keys += [('z',)]
            variables=dict(zip(keys,model.getVars()))
            optimize(model,*args,**kwargs)
            self.assertEqual(model.Status,gp.GRB.OPTIMAL)
            flow=sum(key[2]*var.X for key,var in variables.items() if key[0]=='q')
            occupancy=sum(var.X for key,var in variables.items()
                          if key[0]=='x' and key[3]==1 and key[1][:2] not in solver.output_set)
            residuals.append(flow-occupancy)
        with patch.object(gp.Model,'optimize',inspect_solution):
            result=solver.solve([(1,0),(2,0)],[(0,0),(0,1)],horizon)
        self.assertEqual(result['status_name'],'OPTIMAL')
        self.assertEqual(len(residuals),1)
        self.assertAlmostEqual(residuals[0],0)

    def test_terminal_blocker_cannot_move_beyond_physical_horizon(self):
        solver=self.solver(3,1,((0,0),))
        values={('x',(2,0,1,0),0,2):1}
        self.assertEqual(self.fix_assignment(solver,[],[(0,0),(1,0)],0,values)['status_name'],'INFEASIBLE')

if __name__=='__main__': unittest.main()
