from dataclasses import dataclass
import math
import time

import gurobipy as gp
from gurobipy import GRB


@dataclass(frozen=True)
class StaticGurobiConfig:
    Lx: int
    Ly: int
    output_cells: tuple
    retrieval_mode: str
    beta: float
    gamma: float
    time_limit: int
    lp: bool = False


class StaticEscortFlowGurobiSolver:
    def __init__(self, config):
        self.config = config
        self.output_cells = tuple(config.output_cells)
        self.output_set = set(self.output_cells)
        self.env = gp.Env(empty=True)
        self.env.setParam("OutputFlag", 0)
        self.env.start()
        self.network = self._build_network()

    def close(self):
        if self.env is not None:
            self.env.dispose()
            self.env = None

    def _build_network(self):
        locations = [(x, y) for x in range(self.config.Lx) for y in range(self.config.Ly)]

        na = {loc: [] for loc in locations}
        ne = {loc: [] for loc in locations}
        moves_a = []
        moves_e = []
        stay_move = {}
        outgoing_a = {loc: [] for loc in locations}
        incoming_a = {loc: [] for loc in locations}
        outgoing_e = {loc: [] for loc in locations}
        incoming_e = {loc: [] for loc in locations}
        move_cost_e = {}

        for x, y in locations:
            loc = (x, y)

            if x < self.config.Lx - 1:
                move = (x, y, x + 1, y)
                na[loc].append((x + 1, y))
                moves_a.append(move)
                outgoing_a[loc].append(move)
                incoming_a[(x + 1, y)].append(move)
            if y < self.config.Ly - 1:
                move = (x, y, x, y + 1)
                na[loc].append((x, y + 1))
                moves_a.append(move)
                outgoing_a[loc].append(move)
                incoming_a[(x, y + 1)].append(move)
            if x > 0:
                move = (x, y, x - 1, y)
                na[loc].append((x - 1, y))
                moves_a.append(move)
                outgoing_a[loc].append(move)
                incoming_a[(x - 1, y)].append(move)
            if y > 0:
                move = (x, y, x, y - 1)
                na[loc].append((x, y - 1))
                moves_a.append(move)
                outgoing_a[loc].append(move)
                incoming_a[(x, y - 1)].append(move)

            move = (x, y, x, y)
            na[loc].append(loc)
            moves_a.append(move)
            stay_move[loc] = move
            outgoing_a[loc].append(move)
            incoming_a[loc].append(move)

            for xx in range(self.config.Lx):
                move = (x, y, xx, y)
                ne[loc].append((xx, y))
                moves_e.append(move)
                outgoing_e[loc].append(move)
                incoming_e[(xx, y)].append(move)
                move_cost_e[move] = abs(x - xx)

            for yy in range(self.config.Ly):
                if yy == y:
                    continue
                move = (x, y, x, yy)
                ne[loc].append((x, yy))
                moves_e.append(move)
                outgoing_e[loc].append(move)
                incoming_e[(x, yy)].append(move)
                move_cost_e[move] = abs(y - yy)

        cell_cover = {loc: [] for loc in locations}
        move_cover = {move: [] for move in moves_a}
        for move in moves_e:
            orig_x, orig_y, dest_x, dest_y = move
            if orig_y == dest_y:
                for x in range(min(orig_x, dest_x), max(orig_x, dest_x) + 1):
                    cell_cover[(x, orig_y)].append(move)
                for x in range(orig_x, dest_x):
                    move_cover[(x + 1, orig_y, x, orig_y)].append(move)
                for x in range(orig_x, dest_x, -1):
                    move_cover[(x - 1, orig_y, x, orig_y)].append(move)
            else:
                for y in range(min(orig_y, dest_y), max(orig_y, dest_y) + 1):
                    cell_cover[(orig_x, y)].append(move)
                for y in range(orig_y, dest_y):
                    move_cover[(orig_x, y + 1, orig_x, y)].append(move)
                for y in range(orig_y, dest_y, -1):
                    move_cover[(orig_x, y - 1, orig_x, y)].append(move)

        moves_e_set = set(moves_e)
        conflicts = []

        for y in range(self.config.Ly):
            for x_end in range(self.config.Lx):
                for x_start in range(x_end):
                    conflict_moves = set()
                    for x1 in range(x_end + 1):
                        for x2 in range(x1 + 1, self.config.Lx):
                            conflict_moves.add((x1, y, x2, y))
                            conflict_moves.add((x2, y, x1, y))
                    for x in range(x_start, x_end + 1):
                        for y1 in range(y + 1):
                            for y2 in range(y, self.config.Ly):
                                if y1 != y2:
                                    conflict_moves.add((x, y1, x, y2))
                                    conflict_moves.add((x, y2, x, y1))
                    conflicts.append(sorted(conflict_moves & moves_e_set))

        for x in range(self.config.Lx):
            for y_end in range(self.config.Ly):
                for y_start in range(y_end):
                    conflict_moves = set()
                    for y1 in range(y_end + 1):
                        for y2 in range(y1 + 1, self.config.Ly):
                            conflict_moves.add((x, y1, x, y2))
                            conflict_moves.add((x, y2, x, y1))
                    for y in range(y_start, y_end + 1):
                        for x1 in range(x + 1):
                            for x2 in range(x, self.config.Lx):
                                if x1 != x2:
                                    conflict_moves.add((x1, y, x2, y))
                                    conflict_moves.add((x2, y, x1, y))
                    conflicts.append(sorted(conflict_moves & moves_e_set))

        incoming_output_moves = {
            output: [
                move for move in incoming_a[output]
                if (move[0], move[1]) != (move[2], move[3])
            ]
            for output in self.output_cells
        }
        nonstay_outgoing_a = {
            loc: [
                move for move in outgoing_a[loc]
                if (move[0], move[1]) != (move[2], move[3])
            ]
            for loc in locations
        }
        arrival_moves = [
            move
            for output in self.output_cells
            for move in incoming_output_moves[output]
        ]

        return {
            "locations": locations,
            "not_outputs": [loc for loc in locations if loc not in self.output_set],
            "na": na,
            "ne": ne,
            "moves_a": moves_a,
            "moves_e": moves_e,
            "stay_move": stay_move,
            "outgoing_a": outgoing_a,
            "incoming_a": incoming_a,
            "outgoing_e": outgoing_e,
            "incoming_e": incoming_e,
            "move_cost_e": move_cost_e,
            "cell_cover": cell_cover,
            "move_cover": move_cover,
            "conflicts": conflicts,
            "incoming_output_moves": incoming_output_moves,
            "nonstay_outgoing_a": nonstay_outgoing_a,
            "arrival_moves": arrival_moves,
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

    def _extract_animation_moves(self, x_a, x_e, calc_makespan):
        export_horizon = int(math.floor(calc_makespan + 1e-9))
        moves = []

        for t in range(export_horizon + 1):
            one_step_moves = []
            for move in self.network["moves_e"]:
                if x_e[(move, t)].X <= 0.99:
                    continue

                orig_x, orig_y, dest_x, dest_y = move
                if (orig_x, orig_y) == (dest_x, dest_y):
                    continue

                if dest_x < orig_x:
                    for x in range(dest_x, orig_x):
                        one_step_moves.append(((x, dest_y), (x + 1, orig_y)))
                elif dest_x > orig_x:
                    for x in range(orig_x, dest_x):
                        one_step_moves.append(((x + 1, dest_y), (x, orig_y)))
                elif dest_y < orig_y:
                    for y in range(dest_y, orig_y):
                        one_step_moves.append(((dest_x, y), (orig_x, y + 1)))
                elif dest_y > orig_y:
                    for y in range(orig_y, dest_y):
                        one_step_moves.append(((dest_x, y + 1), (orig_x, y)))

            if self.config.retrieval_mode == "leave" and t > 0:
                for move in self.network["moves_a"]:
                    if x_a[(move, t - 1)].X <= 0.99:
                        continue
                    if (move[0], move[1]) == (move[2], move[3]):
                        continue
                    if (move[2], move[3]) in self.output_set:
                        one_step_moves.append(((move[2], move[3]), (None, None)))

            moves.append(one_step_moves)

        return moves

    def solve(self, target_positions, escort_positions, T):
        target_set = set(target_positions)
        escort_set = set(escort_positions)
        loads_to_retrieve = len(target_set - self.output_set)

        solve_start = time.perf_counter()
        model = gp.Model("escort_flow_static", env=self.env)
        model.Params.OutputFlag = 0
        model.Params.TimeLimit = self.config.time_limit

        load_vtype = GRB.CONTINUOUS if self.config.lp else GRB.BINARY
        escort_vtype = GRB.CONTINUOUS if self.config.lp else GRB.BINARY
        q_vtype = GRB.CONTINUOUS if self.config.lp else GRB.INTEGER

        tr = range(T + 1)
        x_a = {
            (move, t): model.addVar(lb=0.0, ub=1.0, vtype=load_vtype)
            for move in self.network["moves_a"]
            for t in tr
        }
        x_e = {
            (move, t): model.addVar(lb=0.0, ub=1.0, vtype=escort_vtype)
            for move in self.network["moves_e"]
            for t in tr
        }
        q = {
            output: model.addVar(lb=0.0, vtype=q_vtype)
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

        for output in self.output_cells:
            nonstay_output_moves = self.network["nonstay_outgoing_a"][output]
            for t in tr:
                model.addConstr(
                    gp.quicksum(x_a[(move, t)] for move in nonstay_output_moves) == 0
                )

        for loc in self.network["locations"]:
            supply_a = 1 if loc in target_set else 0
            supply_e = 1 if loc in escort_set else 0
            model.addConstr(
                gp.quicksum(x_a[(move, 0)] for move in self.network["outgoing_a"][loc]) == supply_a
            )
            model.addConstr(
                gp.quicksum(x_e[(move, 0)] for move in self.network["outgoing_e"][loc]) == supply_e
            )

        for loc in self.network["locations"]:
            stay_move = self.network["stay_move"][loc]
            nonstay_target_moves = self.network["nonstay_outgoing_a"][loc]
            for t in tr:
                model.addConstr(
                    gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc]) +
                    gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc]) <= 1
                )
                model.addConstr(
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc]) <= 1
                )
                model.addConstr(
                    1 - x_a[(stay_move, t)] >=
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc])
                )
                for move in nonstay_target_moves:
                    model.addConstr(
                        x_a[(move, t)] <=
                        gp.quicksum(x_e[(escort_move, t)] for escort_move in self.network["move_cover"][move])
                    )

        model.addConstr(
            gp.quicksum(x_a[(move, t)] for move in self.network["arrival_moves"] for t in tr) ==
            loads_to_retrieve
        )

        for output in self.output_cells:
            model.addConstr(
                gp.quicksum(
                    (t + 1) * x_a[(move, t)]
                    for move in self.network["incoming_output_moves"][output]
                    for t in tr
                ) == q[output]
            )

        for conflict_set in self.network["conflicts"]:
            for t in tr:
                model.addConstr(
                    gp.quicksum(x_e[(move, t)] for move in conflict_set) <= 1
                )

        model.update()
        model.optimize()
        cpu_time = time.perf_counter() - solve_start

        status_name = self._status_name(model.Status)
        best_bound = getattr(model, "ObjBound", None)
        has_solution = model.SolCount > 0

        result = {
            "has_solution": has_solution,
            "status_name": status_name,
            "cpu_time": cpu_time,
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

    def build_csv_suffix(self, result):
        if result["has_solution"]:
            return (
                f",{self._format_result_value(result['makespan'])}, "
                f"{self._format_result_value(result['flowtime'])}, "
                f"{self._format_result_value(result['movements'])}, "
                f"{self._format_result_value(result['objective'])}, "
                f"{self._format_result_value(result['best_bound'])}, "
                f"{result['cpu_time']:.4f}, {result['status_name']}"
            )

        return (
            f",-,-,-,-,{self._format_result_value(result['best_bound'])}, "
            f"{result['cpu_time']:.4f}, {result['status_name']}"
        )
