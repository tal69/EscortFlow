from dataclasses import dataclass
import time

import gurobipy as gp
from gurobipy import GRB

from escort_flow_static_gurobi import StaticGurobiConfig as BaseStaticGurobiConfig, StaticEscortFlowGurobiSolver


@dataclass(frozen=True)
class StaticGurobiConfig(BaseStaticGurobiConfig):
    bnc_cut_limit: int | None = None


class BnCStaticEscortFlowGurobiSolver(StaticEscortFlowGurobiSolver):
    cut_tol = 1e-6
    lazy_cut_limit = 20

    def _effective_bnc_cut_limit(self, T):
        if self.config.bnc_cut_limit is not None:
            return self.config.bnc_cut_limit
        return 2 * T

    def _separate_movement_coupling(self, model, x_a_sol, x_e_sol, stay_specs_by_t, move_specs_by_t, tr, add_cut, max_cuts):
        tol = self.cut_tol
        added = 0

        # (8) Target load movements allowed.
        # (9) Target load movements enforced.
        for t in tr:
            for stay_key, cover_keys, stay_expr in stay_specs_by_t[t]:
                cover_value = sum(x_e_sol[key] for key in cover_keys)
                if x_a_sol[stay_key] + cover_value > 1 + tol:
                    add_cut(stay_expr <= 1)
                    added += 1
                    if added >= max_cuts:
                        return added

            for move_key, move_cover_keys, move_expr in move_specs_by_t[t]:
                if x_a_sol[move_key] <= tol:
                    continue
                move_cover_value = sum(x_e_sol[key] for key in move_cover_keys)
                if x_a_sol[move_key] > move_cover_value + tol:
                    add_cut(move_expr <= 0)
                    added += 1
                    if added >= max_cuts:
                        return added

        return added

    def _build_bnc_callback(self, x_a, x_e, stay_specs_by_t, move_specs_by_t, tr, callback_stats, root_cut_limit):
        x_a_items = list(x_a.items())
        x_e_items = list(x_e.items())

        def callback(model, where):
            if where == GRB.Callback.MIPNODE:
                node_status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
                if node_status != GRB.OPTIMAL:
                    return
                node_count = int(model.cbGet(GRB.Callback.MIPNODE_NODCNT))
                if node_count != 0:
                    return

                callback_start = time.perf_counter()
                x_a_values = model.cbGetNodeRel([var for _, var in x_a_items])
                x_e_values = model.cbGetNodeRel([var for _, var in x_e_items])
                x_a_sol = {key: value for (key, _), value in zip(x_a_items, x_a_values)}
                x_e_sol = {key: value for (key, _), value in zip(x_e_items, x_e_values)}

                self._separate_movement_coupling(
                    model,
                    x_a_sol,
                    x_e_sol,
                    stay_specs_by_t,
                    move_specs_by_t,
                    tr,
                    model.cbCut,
                    root_cut_limit,
                )
                callback_stats["time"] += time.perf_counter() - callback_start

            elif where == GRB.Callback.MIPSOL:
                callback_start = time.perf_counter()
                x_a_values = model.cbGetSolution([var for _, var in x_a_items])
                x_e_values = model.cbGetSolution([var for _, var in x_e_items])
                x_a_sol = {key: value for (key, _), value in zip(x_a_items, x_a_values)}
                x_e_sol = {key: value for (key, _), value in zip(x_e_items, x_e_values)}

                self._separate_movement_coupling(
                    model,
                    x_a_sol,
                    x_e_sol,
                    stay_specs_by_t,
                    move_specs_by_t,
                    tr,
                    model.cbLazy,
                    self.lazy_cut_limit,
                )
                callback_stats["time"] += time.perf_counter() - callback_start

        return callback

    def solve(self, target_positions, escort_positions, T, warmstart=None):
        if self.config.lp:
            raise ValueError("Branch-and-cut static Gurobi backend does not support --lp")

        target_set = set(target_positions)
        escort_set = set(escort_positions)
        loads_to_retrieve = len(target_set - self.output_set)

        solve_start = time.perf_counter()
        model = gp.Model("escort_flow_static_bnc", env=self.env)
        model.Params.OutputFlag = 1
        model.Params.StartNodeLimit = 100000
        model.Params.LazyConstraints = 1
        model.Params.PreCrush = 1
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

        # (6) Generalized capacity constraint.
        for loc in self.network["locations"]:
            for t in tr:
                model.addConstr(
                    gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc]) +
                    gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc]) <= 1
                )

        # (7) Avoid conflicts.
        for loc in self.network["locations"]:
            for t in tr:
                model.addConstr(
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc]) <= 1
                )

        # (10) All target loads arrive at output cells.
        model.addConstr(
            gp.quicksum(x_a[(move, t)] for move in self.network["arrival_moves"] for t in tr) ==
            loads_to_retrieve
        )

        # Only the movement-coupling constraints are separated in the callback below.
        stay_specs_by_t = {t: [] for t in tr}
        move_specs_by_t = {t: [] for t in tr}
        for loc in self.network["locations"]:
            stay_move = self.network["stay_move"][loc]
            nonstay_target_moves = self.network["nonstay_outgoing_a"][loc]
            for t in tr:
                cover_keys = tuple((move, t) for move in self.network["cell_cover"][loc])
                stay_specs_by_t[t].append(
                    (
                        (stay_move, t),
                        cover_keys,
                        x_a[(stay_move, t)] + gp.quicksum(x_e[key] for key in cover_keys),
                    )
                )
                for move in nonstay_target_moves:
                    move_cover_keys = tuple((escort_move, t) for escort_move in self.network["move_cover"][move])
                    move_specs_by_t[t].append(
                        (
                            (move, t),
                            move_cover_keys,
                            x_a[(move, t)] - gp.quicksum(x_e[key] for key in move_cover_keys),
                        )
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
        callback_stats = {"time": 0.0}
        root_cut_limit = self._effective_bnc_cut_limit(T)
        callback = self._build_bnc_callback(x_a, x_e, stay_specs_by_t, move_specs_by_t, tr, callback_stats, root_cut_limit)
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
            "user_cut_time": callback_stats["time"],
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
