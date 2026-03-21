from dataclasses import dataclass
import math
import time

import gurobipy as gp
from gurobipy import GRB


@dataclass(frozen=True)
class LoadFlowStaticGurobiConfig:
    Lx: int
    Ly: int
    output_cells: tuple
    move_method: str
    alpha: float
    beta: float
    gamma: float
    time_limit: int | None
    work_limit: float | None = None
    lp: bool = False
    threads: int = 0


class LoadFlowStaticGurobiSolver:
    def __init__(self, config):
        self.config = config
        self.output_cells = tuple(config.output_cells)
        self.output_set = set(self.output_cells)
        self.env = gp.Env(empty=True)
        self.env.setParam("OutputFlag", 1)
        self.env.setParam("Threads", self.config.threads)
        self.env.start()
        self.network = self._build_network()

    def close(self):
        if self.env is not None:
            self.env.dispose()
            self.env = None

    def _build_network(self):
        locations = [(x, y) for x in range(self.config.Lx) for y in range(self.config.Ly)]
        moves = []
        move_cost = {}
        outgoing = {loc: [] for loc in locations}
        incoming = {loc: [] for loc in locations}
        incoming_nonstay = {loc: [] for loc in locations}
        outgoing_nonstay = {loc: [] for loc in locations}
        incoming_horizontal = {loc: [] for loc in locations}
        outgoing_horizontal = {loc: [] for loc in locations}
        incoming_vertical = {loc: [] for loc in locations}
        outgoing_vertical = {loc: [] for loc in locations}
        reverse_move = {}

        for x, y in locations:
            for dest_x, dest_y in ((x + 1, y), (x, y + 1), (x - 1, y), (x, y - 1), (x, y)):
                if not (0 <= dest_x < self.config.Lx and 0 <= dest_y < self.config.Ly):
                    continue
                move = (x, y, dest_x, dest_y)
                moves.append(move)
                cost = 0.0 if (x, y) == (dest_x, dest_y) else self.config.gamma
                move_cost[move] = cost
                outgoing[(x, y)].append(move)
                incoming[(dest_x, dest_y)].append(move)
                if (x, y) != (dest_x, dest_y):
                    incoming_nonstay[(dest_x, dest_y)].append(move)
                    outgoing_nonstay[(x, y)].append(move)
                    if x == dest_x:
                        incoming_vertical[(dest_x, dest_y)].append(move)
                        outgoing_vertical[(x, y)].append(move)
                    else:
                        incoming_horizontal[(dest_x, dest_y)].append(move)
                        outgoing_horizontal[(x, y)].append(move)
                    reverse_move[move] = (dest_x, dest_y, x, y)

        nonstay_moves = [move for move in moves if move_cost[move] > 0.0]

        return {
            "locations": locations,
            "non_output_locations": [loc for loc in locations if loc not in self.output_set],
            "moves": moves,
            "nonstay_moves": nonstay_moves,
            "move_cost": move_cost,
            "outgoing": outgoing,
            "incoming": incoming,
            "incoming_nonstay": incoming_nonstay,
            "outgoing_nonstay": outgoing_nonstay,
            "incoming_horizontal": incoming_horizontal,
            "outgoing_horizontal": outgoing_horizontal,
            "incoming_vertical": incoming_vertical,
            "outgoing_vertical": outgoing_vertical,
            "reverse_move": reverse_move,
        }

    @staticmethod
    def _status_name(status_code):
        status_names = {
            GRB.LOADED: "LOADED",
            GRB.OPTIMAL: "OPTIMAL",
            GRB.INFEASIBLE: "INFEASIBLE",
            GRB.INF_OR_UNBD: "INF_OR_UNBD",
            GRB.UNBOUNDED: "UNBOUNDED",
            GRB.CUTOFF: "CUTOFF",
            GRB.ITERATION_LIMIT: "ITERATION_LIMIT",
            GRB.NODE_LIMIT: "NODE_LIMIT",
            GRB.TIME_LIMIT: "TIME_LIMIT",
            GRB.SOLUTION_LIMIT: "SOLUTION_LIMIT",
            GRB.INTERRUPTED: "INTERRUPTED",
            GRB.NUMERIC: "NUMERIC",
            GRB.SUBOPTIMAL: "SUBOPTIMAL",
            GRB.INPROGRESS: "INPROGRESS",
            GRB.USER_OBJ_LIMIT: "USER_OBJ_LIMIT",
            GRB.WORK_LIMIT: "WORK_LIMIT",
            GRB.MEM_LIMIT: "MEM_LIMIT",
        }
        return status_names.get(status_code, str(status_code))

    @staticmethod
    def _format_result_value(value):
        if value is None:
            return "-"
        if isinstance(value, float):
            if not math.isfinite(value):
                return "-"
            if abs(value - round(value)) < 1e-9:
                return str(int(round(value)))
        return str(value)

    def _extract_animation_moves(self, x, q, horizon):
        export_horizon = max(0, int(math.ceil(horizon - 1e-9)))
        moves = []
        for t in range(export_horizon + 1):
            one_step_moves = []
            for move in self.network["nonstay_moves"]:
                if sum(x[(move, t, commodity)].X for commodity in (1, 2)) > 0.99:
                    one_step_moves.append(((move[0], move[1]), (move[2], move[3])))
            for output in self.output_cells:
                if q[(output, t)].X > 0.99:
                    one_step_moves.append((output, (None, None)))
            moves.append(one_step_moves)
        return moves

    def solve(self, target_positions, escort_positions, T):
        target_set = set(target_positions)
        escort_set = set(escort_positions)
        blocking_set = set(self.network["locations"]) - target_set - escort_set

        solve_start = time.perf_counter()
        model = gp.Model("load_flow_static", env=self.env)
        model.Params.OutputFlag = 1
        if self.config.time_limit is not None:
            model.Params.TimeLimit = self.config.time_limit
        if self.config.work_limit is not None:
            model.Params.WorkLimit = self.config.work_limit

        flow_vtype = GRB.CONTINUOUS if self.config.lp else GRB.BINARY
        q_vtype = GRB.CONTINUOUS if self.config.lp else GRB.BINARY
        z_vtype = GRB.CONTINUOUS if self.config.lp else GRB.INTEGER

        tr = range(T + 1)
        x = {
            (move, t, commodity): model.addVar(lb=0.0, ub=1.0, vtype=flow_vtype)
            for move in self.network["moves"]
            for t in tr
            for commodity in (1, 2)
        }
        q = {
            (output, t): model.addVar(lb=0.0, ub=1.0, vtype=q_vtype)
            for output in self.output_cells
            for t in tr
        }
        z = model.addVar(lb=0.0, vtype=z_vtype)

        movement_expr = gp.quicksum(
            x[(move, t, commodity)]
            for move in self.network["nonstay_moves"]
            for t in tr
            for commodity in (1, 2)
        )
        flow_time_expr = gp.quicksum(
            t * q[(output, t)]
            for output in self.output_cells
            for t in tr
        )
        model.setObjective(
            self.config.alpha * z + self.config.gamma * movement_expr + self.config.beta * flow_time_expr,
            GRB.MINIMIZE,
        )

        for loc in self.network["locations"]:
            for t in range(1, T + 1):
                q_term = q[(loc, t)] if loc in self.output_set else 0.0
                model.addConstr(
                    gp.quicksum(x[(move, t - 1, 1)] for move in self.network["incoming"][loc]) ==
                    gp.quicksum(x[(move, t, 1)] for move in self.network["outgoing"][loc]) + q_term
                )
                model.addConstr(
                    gp.quicksum(x[(move, t - 1, 2)] for move in self.network["incoming"][loc]) ==
                    gp.quicksum(x[(move, t, 2)] for move in self.network["outgoing"][loc])
                )

        for output in self.output_cells:
            for t in tr:
                model.addConstr(
                    gp.quicksum(
                        x[(move, t, commodity)]
                        for commodity in (1, 2)
                        for move in self.network["incoming_nonstay"][output]
                    ) <= 1 - q[(output, t)]
                )

        for output in self.output_cells:
            supply_target = 1 if output in target_set else 0
            model.addConstr(
                gp.quicksum(x[(move, 0, 1)] for move in self.network["outgoing"][output]) + q[(output, 0)] ==
                supply_target
            )

        for loc in self.network["non_output_locations"]:
            supply_target = 1 if loc in target_set else 0
            model.addConstr(
                gp.quicksum(x[(move, 0, 1)] for move in self.network["outgoing"][loc]) == supply_target
            )

        for loc in self.network["locations"]:
            supply_blocking = 1 if loc in blocking_set else 0
            model.addConstr(
                gp.quicksum(x[(move, 0, 2)] for move in self.network["outgoing"][loc]) == supply_blocking
            )

        model.addConstr(gp.quicksum(q.values()) == len(target_set))

        for loc in self.network["locations"]:
            for t in range(1, T + 1):
                model.addConstr(
                    gp.quicksum(
                        x[(move, t - 1, commodity)]
                        for commodity in (1, 2)
                        for move in self.network["incoming"][loc]
                    ) <= 1
                )

        if self.config.move_method == "LM":
            for loc in self.network["locations"]:
                for t in tr:
                    q_term = q[(loc, t)] if loc in self.output_set else 0.0
                    model.addConstr(
                        gp.quicksum(
                            x[(move, t, commodity)]
                            for commodity in (1, 2)
                            for move in self.network["incoming_nonstay"][loc]
                        ) + q_term + gp.quicksum(
                            x[(move, t, commodity)]
                            for commodity in (1, 2)
                            for move in self.network["outgoing_nonstay"][loc]
                        ) <= 1
                    )
        else:
            for loc in self.network["locations"]:
                for t in tr:
                    model.addConstr(
                        gp.quicksum(
                            x[(move, t, commodity)]
                            for commodity in (1, 2)
                            for move in self.network["incoming_vertical"][loc]
                        ) + gp.quicksum(
                            x[(move, t, commodity)]
                            for commodity in (1, 2)
                            for move in self.network["outgoing_horizontal"][loc]
                        ) <= 1
                    )
                    model.addConstr(
                        gp.quicksum(
                            x[(move, t, commodity)]
                            for commodity in (1, 2)
                            for move in self.network["incoming_horizontal"][loc]
                        ) + gp.quicksum(
                            x[(move, t, commodity)]
                            for commodity in (1, 2)
                            for move in self.network["outgoing_vertical"][loc]
                        ) <= 1
                    )

            seen_pairs = set()
            for move in self.network["nonstay_moves"]:
                reverse_move = self.network["reverse_move"][move]
                pair_key = tuple(sorted((move, reverse_move)))
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)
                for t in range(T):
                    model.addConstr(
                        gp.quicksum(x[(arc, t, commodity)] for arc in pair_key for commodity in (1, 2)) <= 1
                    )

        if self.config.alpha > 0:
            for output in self.output_cells:
                for t in tr:
                    model.addConstr(t * q[(output, t)] <= z)

        model.optimize()
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
                actual_makespan = max(
                    (t for output in self.output_cells for t in tr if q[(output, t)].X > 1e-6),
                    default=0,
                )
                result["makespan"] = z.X if self.config.alpha > 0 else actual_makespan
                result["flowtime"] = sum(
                    t * q[(output, t)].X
                    for output in self.output_cells
                    for t in tr
                )
                result["movements"] = sum(
                    x[(move, t, commodity)].X
                    for move in self.network["nonstay_moves"]
                    for t in tr
                    for commodity in (1, 2)
                )
                result["objective"] = model.ObjVal
                result["animation_moves"] = self._extract_animation_moves(x, q, actual_makespan)
        finally:
            model.dispose()

        return result

    def build_csv_suffix(self, result):
        if result["has_solution"]:
            return (
                f",{self._format_result_value(result['makespan'])}, "
                f"{self._format_result_value(result['flowtime'])}, "
                f"{self._format_result_value(result['movements'])}, "
                f"{self._format_result_value(result['objective'])}, "
                f"{self._format_result_value(result['best_bound'])}, "
                f"{result['cpu_time']:.4f}, "
                f"{self._format_result_value(result.get('work'))}"
            )

        return (
            f",-,-,-,-,{self._format_result_value(result['best_bound'])}, "
            f"{result['cpu_time']:.4f}, "
            f"{self._format_result_value(result.get('work'))}"
        )
