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
    time_limit: int | None
    work_limit: float | None = None
    lp: bool = False
    omit_capacity_constraints: bool = False


class StaticEscortFlowGurobiSolver:
    def __init__(self, config):
        self.config = config
        self.output_cells = tuple(config.output_cells)
        self.output_set = set(self.output_cells)
        self.env = gp.Env(empty=True)
        self.env.setParam("OutputFlag", 1)
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
        export_horizon = max(0, int(math.ceil(calc_makespan - 1e-9)) - 1)
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

    def build_warmstart_from_trace(self, target_positions, escort_positions, T, target_move_history, escort_move_history):
        warmstart = {
            "x_a": {},
            "x_e": {},
            "q": {output: 0.0 for output in self.output_cells},
        }

        active_targets = {loc for loc in target_positions if loc not in self.output_set}
        current_output_stays = {loc for loc in target_positions if loc in self.output_set}
        active_escorts = set(escort_positions)

        for t in range(T + 1):
            step_target_moves = target_move_history[t] if t < len(target_move_history) else {}
            step_escort_moves = escort_move_history[t] if t < len(escort_move_history) else []

            load_dest_by_source = {
                source: dest
                for source, dest in step_target_moves.values()
            }
            escort_dest_by_source = {
                (orig_x, orig_y): (dest_x, dest_y)
                for orig_x, orig_y, dest_x, dest_y in step_escort_moves
            }

            for loc in current_output_stays:
                warmstart["x_a"][(self.network["stay_move"][loc], t)] = 1.0

            for loc in active_targets:
                dest = load_dest_by_source.get(loc, loc)
                warmstart["x_a"][((loc[0], loc[1], dest[0], dest[1]), t)] = 1.0

            for loc in active_escorts:
                dest = escort_dest_by_source.get(loc, loc)
                warmstart["x_e"][((loc[0], loc[1], dest[0], dest[1]), t)] = 1.0

            next_targets = set()
            next_output_stays = set()

            for loc in active_targets:
                dest = load_dest_by_source.get(loc, loc)
                if dest in self.output_set and dest != loc:
                    warmstart["q"][dest] += t + 1
                    if self.config.retrieval_mode == "stay":
                        next_targets.add(dest)
                    elif self.config.retrieval_mode == "leave":
                        next_output_stays.add(dest)
                else:
                    next_targets.add(dest)

            next_escorts = {
                escort_dest_by_source.get(loc, loc)
                for loc in active_escorts
            }
            if self.config.retrieval_mode == "leave":
                next_escorts.update(current_output_stays)
            elif self.config.retrieval_mode == "stay":
                next_output_stays.update(current_output_stays)

            active_targets = next_targets
            current_output_stays = next_output_stays
            active_escorts = next_escorts

        return self._densify_warmstart(warmstart, T)

    def build_feasible_leave_warmstart(self, target_positions, escort_positions):
        import OneStepHeuristic_v2

        warmstart = {
            "x_a": {},
            "x_e": {},
            "q": {output: 0.0 for output in self.output_cells},
        }

        ordered_targets = sorted(
            set(target_positions),
            key=lambda a: (min(abs(a[0] - o[0]) + abs(a[1] - o[1]) for o in self.output_cells), a[0], a[1]),
        )
        all_targets = {loc: idx + 1 for idx, loc in enumerate(ordered_targets)}
        active_targets = {loc: target_id for loc, target_id in all_targets.items() if loc not in self.output_set}
        current_output_stays = {loc: target_id for loc, target_id in all_targets.items() if loc in self.output_set}
        active_escorts = set(escort_positions)
        dist_map = OneStepHeuristic_v2.build_dist_map(self.config.Lx, self.config.Ly, self.output_cells)

        def move_touches_blocked_cells(move, blocked_cells):
            if move is None or not blocked_cells:
                return False
            orig_x, orig_y, dest_x, dest_y = move
            dir_x = int(math.copysign(1, dest_x - orig_x)) if dest_x != orig_x else 0
            dir_y = int(math.copysign(1, dest_y - orig_y)) if dest_y != orig_y else 0
            x, y = orig_x, orig_y
            while True:
                if (x, y) in blocked_cells:
                    return True
                if (x, y) == (dest_x, dest_y):
                    return False
                x += dir_x
                y += dir_y

        t = 0
        while active_targets or current_output_stays:
            moving_escort = None
            if active_targets:
                step_result = OneStepHeuristic_v2.OneStep(
                    self.config.Lx,
                    self.config.Ly,
                    set(self.output_cells),
                    active_targets,
                    active_escorts,
                    dist_map,
                    retrieval_mode="continue",
                    return_escort_moves=True,
                )
                _, _, _, escort_moves = step_result
                if not escort_moves and not current_output_stays:
                    raise RuntimeError("Greedy leave warmstart got stuck without an escort move")
                for escort_move in escort_moves:
                    if not move_touches_blocked_cells(escort_move, current_output_stays.keys()):
                        moving_escort = escort_move
                        break
                if moving_escort is None and not current_output_stays:
                    raise RuntimeError("Greedy leave warmstart could not find a feasible escort move")

            escort_dest_by_source = {}
            target_dest_by_source = {}

            if moving_escort is not None:
                orig_x, orig_y, dest_x, dest_y = moving_escort
                escort_dest_by_source[(orig_x, orig_y)] = (dest_x, dest_y)

                dir_x = int(math.copysign(1, dest_x - orig_x)) if dest_x != orig_x else 0
                dir_y = int(math.copysign(1, dest_y - orig_y)) if dest_y != orig_y else 0

                x, y = orig_x, orig_y
                while (x, y) != (dest_x, dest_y):
                    next_loc = (x + dir_x, y + dir_y)
                    if next_loc in active_targets:
                        target_dest_by_source[next_loc] = (x, y)
                    x, y = next_loc

            for loc in current_output_stays:
                warmstart["x_a"][(self.network["stay_move"][loc], t)] = 1.0

            for loc in active_targets:
                dest = target_dest_by_source.get(loc, loc)
                warmstart["x_a"][((loc[0], loc[1], dest[0], dest[1]), t)] = 1.0

            for loc in active_escorts:
                dest = escort_dest_by_source.get(loc, loc)
                warmstart["x_e"][((loc[0], loc[1], dest[0], dest[1]), t)] = 1.0

            next_targets = {}
            next_output_stays = {}
            for loc, target_id in active_targets.items():
                dest = target_dest_by_source.get(loc, loc)
                if dest in self.output_set and dest != loc:
                    warmstart["q"][dest] += t + 1
                    next_output_stays[dest] = target_id
                else:
                    next_targets[dest] = target_id

            next_escorts = {
                escort_dest_by_source.get(loc, loc)
                for loc in active_escorts
            }
            next_escorts.update(current_output_stays.keys())

            active_targets = next_targets
            current_output_stays = next_output_stays
            active_escorts = next_escorts
            t += 1

        T = max(0, t - 1)
        return T, self._densify_warmstart(warmstart, T)

    def _densify_warmstart(self, warmstart, T):
        dense_warmstart = {
            "x_a": {},
            "x_e": {},
            "q": {},
        }

        for t in range(T + 1):
            for move in self.network["moves_a"]:
                dense_warmstart["x_a"][(move, t)] = warmstart["x_a"].get((move, t), 0.0)
            for move in self.network["moves_e"]:
                dense_warmstart["x_e"][(move, t)] = warmstart["x_e"].get((move, t), 0.0)

        for output in self.output_cells:
            dense_warmstart["q"][output] = warmstart["q"].get(output, 0.0)

        return dense_warmstart

    @staticmethod
    def _apply_warmstart(x_a, x_e, q, warmstart):
        for key, var in x_a.items():
            var.Start = warmstart["x_a"].get(key, 0.0)
        for key, var in x_e.items():
            var.Start = warmstart["x_e"].get(key, 0.0)
        for output, var in q.items():
            var.Start = warmstart["q"].get(output, 0.0)

    def solve(self, target_positions, escort_positions, T, warmstart=None):
        target_set = set(target_positions)
        escort_set = set(escort_positions)
        loads_to_retrieve = len(target_set - self.output_set)

        solve_start = time.perf_counter()
        model = gp.Model("escort_flow_static", env=self.env)
        model.Params.OutputFlag = 1
        model.Params.StartNodeLimit = 100000
        if self.config.time_limit is not None:
            model.Params.TimeLimit = self.config.time_limit
        if self.config.work_limit is not None:
            model.Params.WorkLimit = self.config.work_limit

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
            # Flow conservation at nodes for escorts and target loads in stay mode.
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
            # Flow conservation at nodes for escorts in continue mode.
            for loc in self.network["locations"]:
                for t in range(1, T + 1):
                    model.addConstr(
                        gp.quicksum(x_e[(move, t - 1)] for move in self.network["incoming_e"][loc]) ==
                        gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc])
                    )
            # Flow conservation at non-output nodes for target loads in continue mode.
            for loc in self.network["not_outputs"]:
                for t in range(1, T + 1):
                    model.addConstr(
                        gp.quicksum(x_a[(move, t - 1)] for move in self.network["incoming_a"][loc]) ==
                        gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc])
                    )
        elif self.config.retrieval_mode == "leave":
            # Leave mode at outputs: arrivals become one-period load stays, then convert into escorts.
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
            # Leave mode away from outputs: standard flow conservation for target loads and escorts.
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

        # (3) in the paper: target loads stay once they reach an output cell.
        for output in self.output_cells:
            nonstay_output_moves = self.network["nonstay_outgoing_a"][output]
            for t in tr:
                model.addConstr(
                    gp.quicksum(x_a[(move, t)] for move in nonstay_output_moves) == 0
                )

        # (4) and (5) in the paper: supply at the initial locations of target loads and escorts.
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
                # Generalized capacity constraint, (6) in the paper
                if not self.config.omit_capacity_constraints:
                    model.addConstr(
                        gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc]) +
                        gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc]) <= 1
                    )

                # Escort conflict avoidance constraint, (7)
                model.addConstr(
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc]) <= 1
                )

                # A target load must move if crossed by an escort (9)
                model.addConstr(
                    1 - x_a[(stay_move, t)] >=
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc])
                )

                # A target load cannot move unless crossed by an escort (8)
                for move in nonstay_target_moves:
                    model.addConstr(
                        x_a[(move, t)] <=
                        gp.quicksum(x_e[(escort_move, t)] for escort_move in self.network["move_cover"][move])
                    )

        # (10) in the paper: every target load must eventually arrive at an output cell.
        model.addConstr(
            gp.quicksum(x_a[(move, t)] for move in self.network["arrival_moves"] for t in tr) ==
            loads_to_retrieve
        )

        # Integrality "cut": q stores the summed arrival times at each output cell.  (15)
        for output in self.output_cells:
            model.addConstr(
                gp.quicksum(
                    (t + 1) * x_a[(move, t)]
                    for move in self.network["incoming_output_moves"][output]
                    for t in tr
                ) == q[output]
            )

        # Explicit conflict cuts are disabled; rely on the cell-cover constraints above.

        model.update()
        if warmstart is not None and not self.config.lp:
            # First let Gurobi search without bias. If it fails to find any incumbent,
            # restart once and use the greedy solution as a fallback start.
            model.optimize()
            status_name = self._status_name(model.Status)
            if model.SolCount == 0 and status_name not in {"INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED", "INTERRUPTED"}:
                model.reset()
                model.Params.SolutionLimit = 1
                model.Params.StartNodeLimit = -2
                self._apply_warmstart(x_a, x_e, q, warmstart)
                model.optimize()
        else:
            model.optimize()
        cpu_time = time.perf_counter() - solve_start

        status_name = self._status_name(model.Status)
        best_bound = getattr(model, "ObjBound", None)
        has_solution = model.SolCount > 0

        result = {
            "has_solution": has_solution,
            "status_name": status_name,
            "cpu_time": cpu_time,
            "user_cut_time": 0.0,
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

    def build_csv_suffix(self, result):
        if result["has_solution"]:
            return (
                f",{self._format_result_value(result['makespan'])}, "
                f"{self._format_result_value(result['flowtime'])}, "
                f"{self._format_result_value(result['movements'])}, "
                f"{self._format_result_value(result['objective'])}, "
                f"{self._format_result_value(result['best_bound'])}, "
                f"{result['cpu_time']:.4f}, {self._format_result_value(result.get('work'))}, "
                f"{result.get('user_cut_time', 0.0):.4f}, {result['status_name']}"
            )

        return (
            f",-,-,-,-,{self._format_result_value(result['best_bound'])}, "
            f"{result['cpu_time']:.4f}, {self._format_result_value(result.get('work'))}, "
            f"{result.get('user_cut_time', 0.0):.4f}, {result['status_name']}"
        )
