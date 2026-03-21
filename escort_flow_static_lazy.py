from dataclasses import dataclass
import time

import gurobipy as gp
from gurobipy import GRB

from escort_flow_static_gurobi import StaticGurobiConfig as BaseStaticGurobiConfig, StaticEscortFlowGurobiSolver


@dataclass(frozen=True)
class StaticGurobiConfig(BaseStaticGurobiConfig):
    lazy_master_steps: int = 0


class LazyStaticEscortFlowGurobiSolver(StaticEscortFlowGurobiSolver):
    lazy_tol = 1e-6

    def _build_lazy_callback(self, x_a, x_e, lazy_tr):
        x_a_items = list(x_a.items())
        x_e_items = list(x_e.items())
        network = self.network
        tol = self.lazy_tol

        def callback(model, where):
            if where != GRB.Callback.MIPSOL:
                return

            x_a_values = model.cbGetSolution([var for _, var in x_a_items])
            x_e_values = model.cbGetSolution([var for _, var in x_e_items])
            x_a_sol = {key: value for (key, _), value in zip(x_a_items, x_a_values)}
            x_e_sol = {key: value for (key, _), value in zip(x_e_items, x_e_values)}

            # (8) in the paper - target load movements allowed and enforced.
            for loc in network["locations"]:
                cell_cover = network["cell_cover"][loc]
                stay_move = network["stay_move"][loc]
                nonstay_target_moves = network["nonstay_outgoing_a"][loc]

                for t in lazy_tr:
                    cover_value = sum(x_e_sol[(move, t)] for move in cell_cover)
                    stay_cover_value = x_a_sol[(stay_move, t)] + cover_value
                    if stay_cover_value > 1 + tol:
                        model.cbLazy(
                            x_a[(stay_move, t)] +
                            gp.quicksum(x_e[(move, t)] for move in cell_cover) <= 1
                        )

                    # Target load movements allowed only when induced by escort movements.
                    for move in nonstay_target_moves:
                        move_cover_value = sum(
                            x_e_sol[(escort_move, t)]
                            for escort_move in network["move_cover"][move]
                        )
                        if x_a_sol[(move, t)] > move_cover_value + tol:
                            model.cbLazy(
                                x_a[(move, t)] <=
                                gp.quicksum(
                                    x_e[(escort_move, t)]
                                    for escort_move in network["move_cover"][move]
                                )
                            )

        return callback

    def solve(self, target_positions, escort_positions, T, warmstart=None):
        if self.config.lp:
            raise ValueError("Lazy static Gurobi backend does not support --lp")

        target_set = set(target_positions)
        escort_set = set(escort_positions)
        loads_to_retrieve = len(target_set - self.output_set)
        master_tr = range(min(T + 1, self.config.lazy_master_steps))
        lazy_tr = range(self.config.lazy_master_steps, T + 1)

        solve_start = time.perf_counter()
        model = gp.Model("escort_flow_static_lazy", env=self.env)
        model.Params.OutputFlag = 1
        model.Params.StartNodeLimit = 100000
        model.Params.LazyConstraints = 1
        if self.config.time_limit is not None:
            model.Params.TimeLimit = self.config.time_limit
        if self.config.work_limit is not None:
            model.Params.WorkLimit = self.config.work_limit

        tr = range(T + 1)
        x_a = {
            (move, t): model.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY)
            for move in self.network["moves_a"]
            for t in tr
        }
        x_e = {
            (move, t): model.addVar(lb=0.0, ub=1.0, vtype=GRB.BINARY)
            for move in self.network["moves_e"]
            for t in tr
        }
        q = {
            output: model.addVar(lb=0.0, vtype=GRB.INTEGER)
            for output in self.output_cells
        }

        number_of_movements = gp.quicksum(
            self.network["move_cost_e"][move] * x_e[(move, t)]
            for move in self.network["moves_e"]
            for t in tr
        )
        model.setObjective(
            self.config.gamma * number_of_movements +
            self.config.beta * gp.quicksum(q[output] for output in self.output_cells),
            GRB.MINIMIZE,
        )

        if self.config.retrieval_mode == "stay":
            for loc in self.network["locations"]:
                for t in range(1, T + 1):
                    model.addConstr(
                        gp.quicksum(x_e[(move, t - 1)] for move in self.network["incoming_e"][loc]) ==
                        gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc])
                    )
                    model.addConstr(
                        gp.quicksum(x_a[(move, t - 1)] for move in self.network["incoming_a"][loc]) ==
                        gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc])
                    )
        elif self.config.retrieval_mode == "continue":
            for loc in self.network["locations"]:
                for t in range(1, T + 1):
                    model.addConstr(
                        gp.quicksum(x_e[(move, t - 1)] for move in self.network["incoming_e"][loc]) ==
                        gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc])
                    )
            for loc in self.network["not_outputs"]:
                for t in range(1, T + 1):
                    model.addConstr(
                        gp.quicksum(x_a[(move, t - 1)] for move in self.network["incoming_a"][loc]) ==
                        gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc])
                    )
        elif self.config.retrieval_mode == "leave":
            for output in self.output_cells:
                stay_move = self.network["stay_move"][output]
                incoming_output_moves = self.network["incoming_output_moves"][output]
                for t in range(1, T + 1):
                    model.addConstr(
                        x_a[(stay_move, t)] ==
                        gp.quicksum(x_a[(move, t - 1)] for move in incoming_output_moves)
                    )
                    model.addConstr(
                        x_a[(stay_move, t - 1)] +
                        gp.quicksum(x_e[(move, t - 1)] for move in self.network["incoming_e"][output]) ==
                        gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][output])
                    )
            for loc in self.network["not_outputs"]:
                for t in range(1, T + 1):
                    model.addConstr(
                        gp.quicksum(x_a[(move, t - 1)] for move in self.network["incoming_a"][loc]) ==
                        gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc])
                    )
                    model.addConstr(
                        gp.quicksum(x_e[(move, t - 1)] for move in self.network["incoming_e"][loc]) ==
                        gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc])
                    )
        else:
            raise ValueError(f"Unsupported retrieval_mode '{self.config.retrieval_mode}'")

        # (3) in the paper - target loads stay once they reach an output cell.
        for output in self.output_cells:
            nonstay_output_moves = self.network["nonstay_outgoing_a"][output]
            for t in tr:
                model.addConstr(
                    gp.quicksum(x_a[(move, t)] for move in nonstay_output_moves) == 0
                )

        # (4) and (5) in the paper - supply at the initial locations of target loads and escorts.
        for loc in self.network["locations"]:
            supply_a = 1 if loc in target_set else 0
            supply_e = 1 if loc in escort_set else 0
            model.addConstr(
                gp.quicksum(x_a[(move, 0)] for move in self.network["outgoing_a"][loc]) == supply_a
            )
            model.addConstr(
                gp.quicksum(x_e[(move, 0)] for move in self.network["outgoing_e"][loc]) == supply_e
            )

        # (7) in the paper - avoid conflicts for the initial master time steps.
        for loc in self.network["locations"]:
            for t in master_tr:
                model.addConstr(
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc]) <= 1
                )

        # (8) and (9) in the renumbered comments for the initial master time steps.
        for loc in self.network["locations"]:
            stay_move = self.network["stay_move"][loc]
            nonstay_target_moves = self.network["nonstay_outgoing_a"][loc]
            cell_cover = self.network["cell_cover"][loc]
            for t in master_tr:
                # Constraint (9) if target is crossed by an escort it can't stay still
                model.addConstr(
                    x_a[(stay_move, t)] +
                    gp.quicksum(x_e[(move, t)] for move in cell_cover) <= 1
                )

                # Constraint (8) A target can move at particular direction only uf crossed by an escort in the opposite direction
                for move in nonstay_target_moves:
                    model.addConstr(
                        x_a[(move, t)] <=
                        gp.quicksum(
                            x_e[(escort_move, t)]
                            for escort_move in self.network["move_cover"][move]
                        )
                    )

        # Constraint family (8)-(9) for later time steps is enforced lazily in the callback below.

        # (10) All target loads arrived
        model.addConstr(
            gp.quicksum(x_a[(move, t)] for move in self.network["arrival_moves"] for t in tr) ==
            loads_to_retrieve
        )

        # Integrality cut: q stores the summed arrival times at each output cell.
        for output in self.output_cells:
            model.addConstr(
                gp.quicksum(
                    (t + 1) * x_a[(move, t)]
                    for move in self.network["incoming_output_moves"][output]
                    for t in tr
                ) == q[output]
            )

        model.update()
        callback = self._build_lazy_callback(x_a, x_e, lazy_tr)
        if warmstart is not None:
            # First let Gurobi search without bias. If it fails to find any incumbent,
            # restart once and use the greedy solution as a fallback start.
            model.optimize(callback)
            status_name = self._status_name(model.Status)
            if model.SolCount == 0 and status_name not in {"INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED", "INTERRUPTED"}:
                model.reset()
                model.Params.SolutionLimit = 1
                model.Params.StartNodeLimit = -2
                self._apply_warmstart(x_a, x_e, q, warmstart)
                model.optimize(callback)
        else:
            model.optimize(callback)
        cpu_time = time.perf_counter() - solve_start

        status_name = self._status_name(model.Status)
        best_bound = getattr(model, "ObjBound", None)
        has_solution = model.SolCount > 0

        result = {
            "has_solution": has_solution,
            "status_name": status_name,
            "cpu_time": cpu_time,
            "work": getattr(model, "Work", None),
            "best_bound": best_bound,
            "makespan": None,
            "flowtime": None,
            "movements": None,
            "objective": None,
            "animation_moves": None,
        }

        try:
            if has_solution:
                arrival_values = [
                    (t + 1) * x_a[(move, t)].X
                    for move in self.network["arrival_moves"]
                    for t in tr
                ]
                result["makespan"] = max(arrival_values, default=0.0)
                result["flowtime"] = sum(
                    (t + 1) * x_a[(move, t)].X
                    for output in self.output_cells
                    for move in self.network["incoming_output_moves"][output]
                    for t in range(1, T + 1)
                )
                result["movements"] = sum(
                    self.network["move_cost_e"][move] * x_e[(move, t)].X
                    for move in self.network["moves_e"]
                    for t in tr
                )
                result["objective"] = model.ObjVal
                result["animation_moves"] = self._extract_animation_moves(
                    x_a,
                    x_e,
                    result["makespan"],
                )
        finally:
            model.dispose()

        return result
