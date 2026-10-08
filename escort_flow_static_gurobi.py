from dataclasses import dataclass
import math
import time

import gurobipy as gp
from gurobipy import GRB

OPTIMALITY_TOLERANCE = 1e-4  # 0.01%
OBJECTIVE_CUTOFF_TOLERANCE = 1e-6


@dataclass(frozen=True)
class StaticGurobiConfig:
    Lx: int
    Ly: int
    output_cells: tuple
    retrieval_mode: str
    beta: float
    gamma: float
    time_limit: float | None
    work_limit: float | None = None
    mip_focus: int = 0
    lp: bool = False
    threads: int = 0
    lexicographic: bool = False
    phase1_time_limit: float | None = None
    objective_mode: str = "legacy"
    certification_target: int | None = None
    weight_scale: int = 100
    flow_proof_extension_time_limit: float | None = None
    stop_on_flow_proof: bool = True


class StaticEscortFlowGurobiSolver:
    def __init__(self, config):
        from static_weighted_certification import validate_objective_mode
        validate_objective_mode(config)
        if config.lexicographic and config.lp:
            raise ValueError("Lexicographic integer stopping is not supported for LP relaxations")
        if config.phase1_time_limit is not None and not config.lexicographic:
            raise ValueError("phase1_time_limit requires lexicographic optimization")
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
        current_output_stays = ({loc for loc in target_positions if loc in self.output_set}
                                if self.config.retrieval_mode != "continue" else set())
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

    def build_warmstart_from_solution(self, result, T):
        """Copy a weighted incumbent into a horizon covering all its arrivals."""
        if self.config.lp or self.config.retrieval_mode != "leave":
            raise ValueError("Solution transfer requires an integer leave model")
        if not result.get("has_solution"):
            raise ValueError("Solution transfer requires a saved incumbent")
        old_T = result["solution_horizon"]
        if T < 0 or int(T) != T or T + 1 < result["makespan"]:
            raise ValueError("The receiving horizon must cover every saved target arrival")
        saved = result["solution_warmstart"]
        start = {name: {key: value for key, value in saved[name].items() if key[1] <= T}
                 for name in ("x_a", "x_e")}
        start["q"] = dict(saved["q"])
        start["removed_post_retrieval_movements"] = sum(
            self.network["move_cost_e"][move] for (move, t), value in saved["x_e"].items()
            if t > T and value > .5)
        if T <= old_T:
            return start

        escorts = {move[2:] for (move, t), value in saved["x_e"].items()
                   if t == old_T and value > .5}
        arrivals = set()
        for (move, t), value in saved["x_a"].items():
            if t != old_T or value <= .5:
                continue
            if move[:2] in self.output_set:
                escorts.add(move[2:])  # Completed one full period of output service.
            elif move[2:] in self.output_set:
                arrivals.add(move[2:])  # Needs its output-service period next.
            else:
                raise ValueError("Saved solution leaves a target unretrieved at the horizon")
        for t in range(old_T + 1, T + 1):
            for loc in escorts:
                start["x_e"][(self.network["stay_move"][loc], t)] = 1.0
            for loc in arrivals:
                start["x_a"][(self.network["stay_move"][loc], t)] = 1.0
            escorts.update(arrivals)
            arrivals = set()
        return start

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

    def summarize_warmstart(self, warmstart):
        arrival_times = [
            t + 1
            for (move, t), value in warmstart["x_a"].items()
            if value > 0.5 and move in self.network["arrival_moves"]
        ]
        flowtime = sum(warmstart["q"].values())
        movements = sum(
            self.network["move_cost_e"][move] * value
            for (move, _), value in warmstart["x_e"].items()
            if value > 0.5
        )
        return {
            "makespan": max(arrival_times, default=0),
            "flowtime": flowtime,
            "movements": movements,
            "objective": self.config.beta * flowtime + self.config.gamma * movements,
        }

    @staticmethod
    def _apply_warmstart(x_a, x_e, q, warmstart):
        for key, var in x_a.items():
            var.Start = warmstart["x_a"].get(key, 0.0)
        for key, var in x_e.items():
            var.Start = warmstart["x_e"].get(key, 0.0)
        for output, var in q.items():
            var.Start = warmstart["q"].get(output, 0.0)

    def solve(self, target_positions, escort_positions, T, warmstart=None, objective_cutoff=None):
        if self.config.objective_mode != "legacy" and objective_cutoff is not None:
            raise ValueError("Weighted/certification modes cannot use an objective cutoff")
        if self.config.lexicographic and objective_cutoff is not None:
            raise ValueError("A weighted objective cutoff is not valid for lexicographic optimization")
        target_set = set(target_positions)
        escort_set = set(escort_positions)
        loads_to_retrieve = len(target_set - self.output_set)

        solve_start = time.perf_counter()
        model = gp.Model("escort_flow_static", env=self.env)
        model.Params.OutputFlag = 1
        model.Params.StartNodeLimit = 100000
        if not self.config.lp:
            model.Params.MIPFocus = self.config.mip_focus
        model.Params.MIPGap = OPTIMALITY_TOLERANCE
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
        objective_expr = (
            self.config.gamma * number_of_movements +
            self.config.beta * gp.quicksum(q[output] for output in self.output_cells)
        )
        model.setObjective(objective_expr, GRB.MINIMIZE)

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
            supply_a = int(loc in target_set and not (
                self.config.retrieval_mode == "continue" and loc in self.output_set))
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
                # Escort conflict avoidance constraint, (6)
                model.addConstr(
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc]) <= 1
                )

                # A target load must move if crossed by an escort (8)
                model.addConstr(
                    1 - x_a[(stay_move, t)] >=
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc])
                )

                # A target load cannot move unless crossed by an escort (7)
                for move in nonstay_target_moves:
                    model.addConstr(
                        x_a[(move, t)] <=
                        gp.quicksum(x_e[(escort_move, t)] for escort_move in self.network["move_cover"][move])
                    )

        # (9) in the paper: every target load must eventually arrive at an output cell.
        model.addConstr(
            gp.quicksum(x_a[(move, t)] for move in self.network["arrival_moves"] for t in tr) ==
            loads_to_retrieve
        )

        # Integrality "cut": q stores the summed arrival times at each output cell.  (14)
        for output in self.output_cells:
            model.addConstr(
                gp.quicksum(
                    (t + 1) * x_a[(move, t)]
                    for move in self.network["incoming_output_moves"][output]
                    for t in tr
                ) == q[output]
            )

        if objective_cutoff is not None:
            cutoff_value = objective_cutoff + OBJECTIVE_CUTOFF_TOLERANCE
            if self.config.lp:
                model.addConstr(objective_expr <= cutoff_value)
            else:
                model.Params.Cutoff = cutoff_value

        # Explicit conflict cuts are disabled; rely on the cell-cover constraints above.

        model.update()
        if warmstart is not None and not self.config.lp:
            self._apply_warmstart(x_a, x_e, q, warmstart)
        if self.config.lexicographic or self.config.objective_mode != "legacy":
            flow_proof_context = None
            if self.config.flow_proof_extension_time_limit is not None:
                flow_proof_context = dict(
                    targets=tuple(target_set), outputs=self.output_cells,
                    cell_count=self.config.Lx * self.config.Ly, escort_count=len(escort_set),
                    physical_horizon=T + 1, weighted_time_limit=self.config.time_limit,
                    extension_time_limit=self.config.flow_proof_extension_time_limit,
                    stop_on_flow_proof=self.config.stop_on_flow_proof)
            return self._solve_lexicographic_model(
                model, x_a, x_e, q, T, solve_start, flow_proof_context=flow_proof_context)
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
                    for t in tr
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

    def _solve_lexicographic_model(self, model, x_a, x_e, q, T, solve_start, callback=None,
                                  flow_proof_context=None):
        from static_lexicographic import solve_lexicographic

        flowtime_expr = gp.quicksum(q.values())
        movement_expr = gp.quicksum(
            self.network["move_cost_e"][move] * x_e[(move, t)]
            for move in self.network["moves_e"] for t in range(T + 1)
        )

        def extract_solution():
            makespan = max(
                (t + 1 for move in self.network["arrival_moves"]
                 for t in range(T + 1) if x_a[(move, t)].X > 0.5),
                default=0,
            )
            result = {
                "makespan": makespan,
                "animation_moves": self._extract_animation_moves(x_a, x_e, makespan),
            }
            if self.config.objective_mode == "weighted_integer":
                # Save the incumbent before the shared solve helper disposes its model.
                result["solution_horizon"] = T
                result["solution_warmstart"] = {
                    "x_a": {key: 1.0 for key, var in x_a.items() if var.X > .5},
                    "x_e": {key: 1.0 for key, var in x_e.items() if var.X > .5},
                    "q": {key: float(round(var.X)) for key, var in q.items()},
                }
            return result

        if self.config.objective_mode != "legacy":
            from static_weighted_certification import solve_weighted_or_certificate

            def extract_callback_metrics(callback_model):
                arrivals = [(t + 1, x_a[(move, t)]) for move in self.network["arrival_moves"]
                            for t in range(T + 1)]
                values = callback_model.cbGetSolution([var for _, var in arrivals])
                return dict(makespan=max((arrival for (arrival, _), value in zip(arrivals, values)
                                          if value > .5), default=0))

            return solve_weighted_or_certificate(
                model, flowtime_expr, movement_expr, extract_solution,
                status_name=self._status_name, solve_start=solve_start,
                mode=self.config.objective_mode, target=self.config.certification_target,
                weight_scale=self.config.weight_scale,
                flow_proof_context=flow_proof_context,
                extract_callback_metrics=extract_callback_metrics,
            )
        return solve_lexicographic(
            model, flowtime_expr, movement_expr, extract_solution,
            status_name=self._status_name, solve_start=solve_start,
            time_limit=self.config.time_limit,
            phase1_time_limit=self.config.phase1_time_limit,
            work_limit=self.config.work_limit, callback=callback,
        )

    def build_csv_suffix(self, result):
        if self.config.lexicographic:
            from static_lexicographic import build_lex_csv_suffix
            return build_lex_csv_suffix(result)
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
