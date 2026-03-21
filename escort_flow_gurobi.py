from dataclasses import dataclass
import time

import gurobipy as gp
from gurobipy import GRB

import OneStepHeuristic_v2


OPTIMALITY_TOLERANCE = 1e-4  # 0.01%

WARMSTART_SUCCESS_MARKERS = (
    "User MIP start produced solution with objective",
    "Loaded user MIP start with objective",
    "MIP start from previous solve produced solution with objective",
    "Loaded MIP start from previous solve with objective",
)
WARMSTART_FAILURE_MARKERS = (
    "User MIP start did not produce a new incumbent solution",
    "MIP start from previous solve did not produce a new incumbent solution",
)
WARMSTART_REPAIR_MARKERS = (
    "Warning: Completing partial solution",
    "Another try with MIP start",
    "User MIP start violates ",
    "MIP start from previous solve violates ",
)


@dataclass(frozen=True)
class SolverConfig:
    Lx: int
    Ly: int
    output_cells: tuple
    distance_penalty: float
    time_penalty: float
    gamma: float
    full: bool
    acyclic: bool
    fractional_horizon: int
    integer_horizon: int
    epoch: int
    time_limit: int
    num_threads: int
    mip_focus: int
    warmstart_mode: str
    warmstart_zero_fill: bool
    optimality_tolerance: float = OPTIMALITY_TOLERANCE


class RollingHorizonGurobiSolver:
    def __init__(self, config):
        self.config = config
        self.output_cells = tuple(config.output_cells)
        self.output_set = set(self.output_cells)
        self.env = gp.Env(empty=True)
        self.env.setParam("OutputFlag", 0)
        self.env.start()
        self.dist_map = OneStepHeuristic_v2.build_dist_map(
            config.Lx, config.Ly, self.output_cells
        )
        self.network = self._build_gurobi_network()
        self.model = None
        self.x_a = {}
        self.x_e = {}
        self.q_vars = None
        self.supply_a_constraints = {}
        self.supply_e_constraints = {}
        self.retrieve_total_constraint = None
        self._pending_initial_model_build_time = self._build_persistent_model()

    def _configure_model_parameters(self, model):
        model.Params.OutputFlag = 0
        model.Params.TimeLimit = self.config.time_limit
        model.Params.MIPGap = self.config.optimality_tolerance
        model.Params.MIPFocus = self.config.mip_focus
        if self.config.num_threads > 0:
            model.Params.Threads = self.config.num_threads

    def _build_gurobi_network(self):
        locations = [
            (x, y)
            for x in range(self.config.Lx)
            for y in range(self.config.Ly)
        ]
        not_outputs = [loc for loc in locations if loc not in self.output_set]

        na = {loc: [loc] for loc in locations}
        ne = {loc: [] for loc in locations}
        moves_a = []
        moves_e = []
        move_cost_e = {}
        outgoing_a = {loc: [] for loc in locations}
        incoming_a = {loc: [] for loc in locations}
        outgoing_e = {loc: [] for loc in locations}
        incoming_e = {loc: [] for loc in locations}
        stay_move = {}

        location_penalty = {
            loc: self.config.time_penalty
            + (self.config.Lx + self.config.Ly + 1) * self.config.distance_penalty
            for loc in locations
        }
        for loc in locations:
            for output in self.output_cells:
                distance = abs(output[0] - loc[0]) + abs(output[1] - loc[1])
                location_penalty[loc] = min(
                    location_penalty[loc],
                    self.config.time_penalty + self.config.distance_penalty * distance,
                )
            if loc in self.output_set:
                location_penalty[loc] = 0

        for x, y in locations:
            move = (x, y, x, y)
            moves_a.append(move)
            outgoing_a[(x, y)].append(move)
            incoming_a[(x, y)].append(move)
            stay_move[(x, y)] = move

            if x < self.config.Lx - 1:
                move = (x, y, x + 1, y)
                na[(x, y)].append((x + 1, y))
                moves_a.append(move)
                outgoing_a[(x, y)].append(move)
                incoming_a[(x + 1, y)].append(move)
            if y < self.config.Ly - 1:
                move = (x, y, x, y + 1)
                na[(x, y)].append((x, y + 1))
                moves_a.append(move)
                outgoing_a[(x, y)].append(move)
                incoming_a[(x, y + 1)].append(move)
            if x > 0:
                move = (x, y, x - 1, y)
                na[(x, y)].append((x - 1, y))
                moves_a.append(move)
                outgoing_a[(x, y)].append(move)
                incoming_a[(x - 1, y)].append(move)
            if y > 0:
                move = (x, y, x, y - 1)
                na[(x, y)].append((x, y - 1))
                moves_a.append(move)
                outgoing_a[(x, y)].append(move)
                incoming_a[(x, y - 1)].append(move)

            for x1 in range(self.config.Lx):
                move = (x, y, x1, y)
                ne[(x, y)].append((x1, y))
                moves_e.append(move)
                move_cost_e[move] = abs(x - x1)
                outgoing_e[(x, y)].append(move)
                incoming_e[(x1, y)].append(move)
            for y1 in range(self.config.Ly):
                if y1 == y:
                    continue
                move = (x, y, x, y1)
                ne[(x, y)].append((x, y1))
                moves_e.append(move)
                move_cost_e[move] = abs(y - y1)
                outgoing_e[(x, y)].append(move)
                incoming_e[(x, y1)].append(move)

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

        return {
            "locations": locations,
            "output_set": self.output_set,
            "not_outputs": not_outputs,
            "na": na,
            "ne": ne,
            "moves_a": moves_a,
            "moves_e": moves_e,
            "move_cost_e": move_cost_e,
            "outgoing_a": outgoing_a,
            "incoming_a": incoming_a,
            "outgoing_e": outgoing_e,
            "incoming_e": incoming_e,
            "cell_cover": cell_cover,
            "move_cover": move_cover,
            "stay_move": stay_move,
            "location_penalty": location_penalty,
        }

    def _build_persistent_model(self):
        model_build_start = time.perf_counter()
        model = gp.Model("escort_flow_rh_v7", env=self.env)
        self._configure_model_parameters(model)

        T = self.config.fractional_horizon
        x_a = {}
        x_e = {}
        for t in range(T + 1):
            a_vtype = GRB.BINARY if self.config.full or t <= self.config.integer_horizon else GRB.CONTINUOUS
            e_vtype = GRB.BINARY if self.config.full or t <= self.config.integer_horizon else GRB.CONTINUOUS
            for move in self.network["moves_a"]:
                x_a[(move, t)] = model.addVar(lb=0.0, ub=1.0, vtype=a_vtype)
            for move in self.network["moves_e"]:
                x_e[(move, t)] = model.addVar(lb=0.0, ub=1.0, vtype=e_vtype)

        q_vars = None
        if self.config.full:
            q_vars = {output: model.addVar(lb=0.0, vtype=GRB.INTEGER) for output in self.output_cells}

        movement_term = gp.quicksum(
            self.network["move_cost_e"][move] * x_e[(move, t)]
            for move in self.network["moves_e"]
            for t in range(T + 1)
        )
        if self.config.full:
            model.setObjective(
                gp.quicksum(q_vars[output] for output in self.output_cells)
                + self.config.gamma * movement_term,
                GRB.MINIMIZE,
            )
        else:
            surrogate_term = gp.quicksum(
                self.network["location_penalty"][(move[2], move[3])] * x_a[(move, t)]
                for move in self.network["moves_a"]
                for t in range(T + 1)
            )
            model.setObjective(
                self.config.gamma * movement_term + surrogate_term,
                GRB.MINIMIZE,
            )

        for t in range(1, T + 1):
            for loc in self.network["locations"]:
                model.addConstr(
                    gp.quicksum(x_e[(move, t - 1)] for move in self.network["incoming_e"][loc]) ==
                    gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc]),
                    name=f"escort_flow_{loc[0]}_{loc[1]}_t{t}",
                )
            for loc in self.network["not_outputs"]:
                model.addConstr(
                    gp.quicksum(x_a[(move, t - 1)] for move in self.network["incoming_a"][loc]) ==
                    gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc]),
                    name=f"target_flow_{loc[0]}_{loc[1]}_t{t}",
                )

        for output in self.output_cells:
            nonstay_output_moves = [
                move for move in self.network["outgoing_a"][output]
                if (move[0], move[1]) != (move[2], move[3])
            ]
            for t in range(T + 1):
                model.addConstr(
                    gp.quicksum(x_a[(move, t)] for move in nonstay_output_moves) == 0,
                    name=f"output_nonstay_{output[0]}_{output[1]}_t{t}",
                )

        supply_a_constraints = {}
        supply_e_constraints = {}
        for loc in self.network["locations"]:
            supply_a_constraints[loc] = model.addConstr(
                gp.quicksum(x_a[(move, 0)] for move in self.network["outgoing_a"][loc]) == 0,
                name=f"supply_a_{loc[0]}_{loc[1]}",
            )
            supply_e_constraints[loc] = model.addConstr(
                gp.quicksum(x_e[(move, 0)] for move in self.network["outgoing_e"][loc]) == 0,
                name=f"supply_e_{loc[0]}_{loc[1]}",
            )

        for loc in self.network["locations"]:
            stay_move = self.network["stay_move"][loc]
            nonstay_target_moves = [
                move for move in self.network["outgoing_a"][loc]
                if move != stay_move
            ]
            for t in range(T + 1):
                model.addConstr(
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc]) <= 1,
                    name=f"cell_cover_{loc[0]}_{loc[1]}_t{t}",
                )
                model.addConstr(
                    1 - x_a[(stay_move, t)] >= gp.quicksum(
                        x_e[(move, t)] for move in self.network["cell_cover"][loc]
                    ),
                    name=f"stay_cover_{loc[0]}_{loc[1]}_t{t}",
                )
                for move in nonstay_target_moves:
                    model.addConstr(
                        x_a[(move, t)] <= gp.quicksum(
                            x_e[(escort_move, t)]
                            for escort_move in self.network["move_cover"][move]
                        ),
                        name=f"move_cover_{move[0]}_{move[1]}_{move[2]}_{move[3]}_t{t}",
                    )

        retrieve_total_constraint = None
        if self.config.full:
            arrival_moves_to_outputs = [
                move for move in self.network["moves_a"]
                if (move[2], move[3]) in self.output_set and (move[0], move[1]) != (move[2], move[3])
            ]
            retrieve_total_constraint = model.addConstr(
                gp.quicksum(
                    x_a[(move, t)] for move in arrival_moves_to_outputs for t in range(T + 1)
                ) == 0
            )
            for output in self.output_cells:
                incoming_output_moves = [
                    move for move in self.network["incoming_a"][output]
                    if (move[0], move[1]) != (move[2], move[3])
                ]
                model.addConstr(
                    gp.quicksum(
                        (t + 1) * x_a[(move, t)]
                        for move in incoming_output_moves
                        for t in range(T + 1)
                    ) == q_vars[output]
                )

        model.update()
        self.model = model
        self.x_a = x_a
        self.x_e = x_e
        self.q_vars = q_vars
        self.supply_a_constraints = supply_a_constraints
        self.supply_e_constraints = supply_e_constraints
        self.retrieve_total_constraint = retrieve_total_constraint
        return time.perf_counter() - model_build_start

    def _update_dynamic_model_state(self, target_positions, escort_positions):
        target_set = set(target_positions)
        escort_set = set(escort_positions)

        for loc in self.network["locations"]:
            self.supply_a_constraints[loc].RHS = 1 if loc in target_set else 0
            self.supply_e_constraints[loc].RHS = 1 if loc in escort_set else 0

        if self.retrieve_total_constraint is not None:
            self.retrieve_total_constraint.RHS = len(target_set - self.output_set)

    def _gurobi_status_name(self, status_code):
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

    def _extract_planned_moves(self, x_e_vars):
        planned_moves = []
        for t in range(self.config.epoch):
            one_step_moves = []
            for move in self.network["moves_e"]:
                if x_e_vars[(move, t)].X <= 0.99:
                    continue
                orig_x, orig_y, dest_x, dest_y = move
                if orig_x == dest_x and orig_y == dest_y:
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

            planned_moves.append(one_step_moves)
        return planned_moves

    def _infer_escort_moves_from_step_moves(self, step_moves):
        if not step_moves:
            return []

        escort_moves = []
        buckets = {}
        for src, dst in step_moves:
            if src[1] == dst[1]:
                direction = dst[0] - src[0]
                buckets.setdefault(("h", src[1], direction), []).append((src, dst))
            else:
                direction = dst[1] - src[1]
                buckets.setdefault(("v", src[0], direction), []).append((src, dst))

        for (axis, line, direction), bucket_moves in buckets.items():
            if axis == "h" and direction == -1:
                ordered = sorted(bucket_moves, key=lambda mv: mv[1][0])
                chains = []
                chain = [ordered[0]]
                for mv in ordered[1:]:
                    if mv[1][0] == chain[-1][1][0] + 1:
                        chain.append(mv)
                    else:
                        chains.append(chain)
                        chain = [mv]
                chains.append(chain)
                for chain in chains:
                    escort_moves.append((chain[0][1][0], line, chain[-1][0][0], line))

            elif axis == "h" and direction == 1:
                ordered = sorted(bucket_moves, key=lambda mv: mv[0][0])
                chains = []
                chain = [ordered[0]]
                for mv in ordered[1:]:
                    if mv[0][0] == chain[-1][0][0] + 1:
                        chain.append(mv)
                    else:
                        chains.append(chain)
                        chain = [mv]
                chains.append(chain)
                for chain in chains:
                    escort_moves.append((chain[-1][1][0], line, chain[0][0][0], line))

            elif axis == "v" and direction == -1:
                ordered = sorted(bucket_moves, key=lambda mv: mv[1][1])
                chains = []
                chain = [ordered[0]]
                for mv in ordered[1:]:
                    if mv[1][1] == chain[-1][1][1] + 1:
                        chain.append(mv)
                    else:
                        chains.append(chain)
                        chain = [mv]
                chains.append(chain)
                for chain in chains:
                    escort_moves.append((line, chain[0][1][1], line, chain[-1][0][1]))

            elif axis == "v" and direction == 1:
                ordered = sorted(bucket_moves, key=lambda mv: mv[0][1])
                chains = []
                chain = [ordered[0]]
                for mv in ordered[1:]:
                    if mv[0][1] == chain[-1][0][1] + 1:
                        chain.append(mv)
                    else:
                        chains.append(chain)
                        chain = [mv]
                chains.append(chain)
                for chain in chains:
                    escort_moves.append((line, chain[-1][1][1], line, chain[0][0][1]))

        return escort_moves

    def _build_greedy_warmstart_segment(self, target_positions, escort_positions, start_t, end_t):
        start_values = {}
        A_state = {loc: idx for idx, loc in enumerate(sorted(target_positions), start=1)}
        E_state = set(escort_positions)

        for t in range(start_t, end_t + 1):
            if A_state:
                A_next, E_next, step_moves = OneStepHeuristic_v2.OneStep(
                    self.config.Lx,
                    self.config.Ly,
                    self.output_set,
                    A_state,
                    set(E_state),
                    self.dist_map,
                    acyclic=self.config.acyclic,
                )
            else:
                A_next, E_next, step_moves = {}, set(E_state), []

            escort_moves = self._infer_escort_moves_from_step_moves(step_moves)
            moved_escort_origins = {(move[0], move[1]) for move in escort_moves}
            for move in escort_moves:
                start_values[("e", move, t)] = 1.0
            for escort in E_state:
                if escort not in moved_escort_origins:
                    start_values[("e", (escort[0], escort[1], escort[0], escort[1]), t)] = 1.0

            next_loc_by_target_id = {target_id: loc for loc, target_id in A_next.items()}
            for orig_loc, target_id in A_state.items():
                next_loc = next_loc_by_target_id.get(target_id)
                if next_loc is None:
                    if orig_loc in self.output_set:
                        continue
                    next_loc = orig_loc
                start_values[("a", (orig_loc[0], orig_loc[1], next_loc[0], next_loc[1]), t)] = 1.0

            A_state = A_next
            E_state = E_next

        return start_values

    def _build_greedy_warmstart(self, target_positions, escort_positions):
        return self._build_greedy_warmstart_segment(
            target_positions, escort_positions, 0, self.config.fractional_horizon
        )

    def _build_shifted_previous_start(self, previous_solution):
        start_values = {}
        for (move, old_t), value in previous_solution["x_a"].items():
            if old_t < self.config.epoch:
                continue
            if value >= 0.99:
                start_values[("a", move, old_t - self.config.epoch)] = 1.0
        for (move, old_t), value in previous_solution["x_e"].items():
            if old_t < self.config.epoch:
                continue
            if value >= 0.99:
                start_values[("e", move, old_t - self.config.epoch)] = 1.0
        return start_values

    def _extract_positions_from_start_values(self, start_values, kind, t):
        return {
            (move[0], move[1])
            for (entry_kind, move, entry_t), value in start_values.items()
            if entry_kind == kind and entry_t == t and value >= 0.99
        }

    def _extract_state_after_warmstart_step(self, start_values, t):
        target_positions = {
            (move[2], move[3])
            for (kind, move, layer_t), value in start_values.items()
            if kind == "a" and layer_t == t and value >= 0.99 and (move[2], move[3]) not in self.output_set
        }
        escort_positions = {
            (move[2], move[3])
            for (kind, move, layer_t), value in start_values.items()
            if kind == "e" and layer_t == t and value >= 0.99
        }
        return target_positions, escort_positions

    def _merge_shifted_with_greedy_tail(self, shifted_start):
        tail_start = max(0, self.config.fractional_horizon - self.config.epoch + 1)
        merged = dict(shifted_start)

        if tail_start == 0:
            tail_target_positions = self._extract_positions_from_start_values(shifted_start, "a", 0)
            tail_escort_positions = self._extract_positions_from_start_values(shifted_start, "e", 0)
        else:
            tail_target_positions, tail_escort_positions = self._extract_state_after_warmstart_step(
                shifted_start, tail_start - 1
            )

        greedy_tail = self._build_greedy_warmstart_segment(
            tail_target_positions, tail_escort_positions,
            tail_start, self.config.fractional_horizon
        )
        merged.update(greedy_tail)
        return merged

    def _fill_unset_warmstart_arcs_with_zero(self, start_values):
        dense_start = dict(start_values)
        for t in range(self.config.fractional_horizon + 1):
            for move in self.network["moves_a"]:
                dense_start.setdefault(("a", move, t), 0.0)
            for move in self.network["moves_e"]:
                dense_start.setdefault(("e", move, t), 0.0)
        return dense_start

    def _build_full_model_q_warmstart(self, start_values):
        if not self.config.full:
            return {}

        q_values = {}
        for output in self.output_cells:
            incoming_output_moves = [
                move for move in self.network["incoming_a"][output]
                if (move[0], move[1]) != (move[2], move[3])
            ]
            q_values[output] = int(round(sum(
                (t + 1) * start_values.get(("a", move, t), 0.0)
                for move in incoming_output_moves
                for t in range(self.config.fractional_horizon + 1)
            )))
        return q_values

    def _finalize_warmstart_vector(self, start_values, source):
        if not start_values:
            return None

        if self.config.warmstart_zero_fill:
            start_values = self._fill_unset_warmstart_arcs_with_zero(start_values)
            q_values = self._build_full_model_q_warmstart(start_values)
            is_complete = True
        else:
            q_values = {}
            is_complete = False

        return {
            "values": start_values,
            "q_values": q_values,
            "source": source,
            "is_complete": is_complete,
        }

    def _shifted_ilp_warmstart_is_state_compatible(self, shifted_start, target_positions, escort_positions):
        return (
            self._extract_positions_from_start_values(shifted_start, "a", 0) == set(target_positions)
            and self._extract_positions_from_start_values(shifted_start, "e", 0) == set(escort_positions)
        )

    def _can_use_shifted_ilp_warmstart(
        self,
        current_target_loads,
        previous_target_loads,
        previous_solution,
        shifted_start,
        target_positions,
        escort_positions,
    ):
        if (
            previous_solution is None
            or previous_target_loads is None
            or not set(current_target_loads).issubset(previous_target_loads)
        ):
            return False

        return self._shifted_ilp_warmstart_is_state_compatible(
            shifted_start, target_positions, escort_positions
        )

    def select_warmstart_vector(
        self,
        target_positions,
        escort_positions,
        current_target_loads,
        previous_solution,
        previous_target_loads,
    ):
        if self.config.warmstart_mode == "none":
            return None

        if self.config.warmstart_mode == "greedy":
            return self._finalize_warmstart_vector(
                self._build_greedy_warmstart(target_positions, escort_positions),
                "greedy",
            )

        shifted_start = None
        if previous_solution is not None:
            shifted_start = self._build_shifted_previous_start(previous_solution)
        can_shift_previous = self._can_use_shifted_ilp_warmstart(
            current_target_loads,
            previous_target_loads,
            previous_solution,
            shifted_start,
            target_positions,
            escort_positions,
        )

        if can_shift_previous:
            return self._finalize_warmstart_vector(
                self._merge_shifted_with_greedy_tail(shifted_start),
                "ilp",
            )

        if self.config.warmstart_mode == "ilp-greedy":
            return self._finalize_warmstart_vector(
                self._build_greedy_warmstart(target_positions, escort_positions),
                "greedy",
            )

        return None

    def _apply_warmstart_vector(self, model, x_a, x_e, q_vars, warmstart_vector):
        if warmstart_vector is None:
            return False

        model.NumStart = 1
        model.Params.StartNumber = 0
        for (kind, move, t), value in warmstart_vector["values"].items():
            if kind == "a":
                x_a[(move, t)].Start = value
            else:
                x_e[(move, t)].Start = value
        if q_vars is not None:
            for output, value in warmstart_vector["q_values"].items():
                q_vars[output].Start = value
        return True

    def _classify_warmstart_outcome(self, gurobi_messages, warmstart_vector):
        if warmstart_vector is None:
            return "none"

        joined_messages = "".join(gurobi_messages)
        success = any(marker in joined_messages for marker in WARMSTART_SUCCESS_MARKERS)
        repaired = any(marker in joined_messages for marker in WARMSTART_REPAIR_MARKERS)
        failed = any(marker in joined_messages for marker in WARMSTART_FAILURE_MARKERS)

        if success:
            if warmstart_vector["is_complete"] and not repaired:
                return "feasible"
            return "repaired"
        if failed or warmstart_vector is not None:
            return "failed"
        return "none"

    def _warmstart_log_excerpt(self, gurobi_messages):
        relevant_lines = [
            line.strip()
            for line in gurobi_messages
            if "MIP start" in line or "partial solution" in line
        ]
        return " | ".join(line for line in relevant_lines if line)

    def run_model(self, target_positions, escort_positions, warmstart_vector=None):
        model_build_start = time.perf_counter()
        model = self.model
        model.reset(1)
        self._update_dynamic_model_state(target_positions, escort_positions)
        T = self.config.fractional_horizon
        T_exec = self.config.epoch
        x_a = self.x_a
        x_e = self.x_e
        q = self.q_vars

        warmstart_applied = self._apply_warmstart_vector(model, x_a, x_e, q, warmstart_vector)
        model.update()

        gurobi_messages = []

        def message_callback(cb_model, where):
            if where == GRB.Callback.MESSAGE:
                gurobi_messages.append(cb_model.cbGet(GRB.Callback.MSG_STRING))

        model_construction_time = (
            time.perf_counter() - model_build_start + self._pending_initial_model_build_time
        )
        self._pending_initial_model_build_time = 0.0
        model.optimize(message_callback if warmstart_applied else None)

        status_code = model.Status
        result = {
            "A": [],
            "E": [],
            "solver_status": status_code,
            "solver_status_name": self._gurobi_status_name(status_code),
            "sol_count": model.SolCount,
            "solver_time_iter": model.Runtime,
            "model_construction_time_iter": model_construction_time,
            "makespan_rh": 0,
            "flowtime_rh": 0,
            "NumberOfMovements_rh": 0,
            "obj_rh": 0.0,
            "lb_rh": 0.0,
            "planned_moves": [],
            "ilp_gap": 1.0,
            "is_optimal_with_tolerance": False,
            "solution_snapshot": None,
            "warmstart_applied": warmstart_applied,
            "warmstart_source": "none" if warmstart_vector is None else warmstart_vector["source"],
            "warmstart_outcome": self._classify_warmstart_outcome(gurobi_messages, warmstart_vector),
            "warmstart_log_excerpt": self._warmstart_log_excerpt(gurobi_messages),
        }

        if model.SolCount <= 0:
            return result

        result["A"] = [
            (move[0], move[1]) for move in self.network["moves_a"] if x_a[(move, T_exec)].X > 0.99
        ]
        result["E"] = [
            (move[0], move[1]) for move in self.network["moves_e"] if x_e[(move, T_exec)].X > 0.99
        ]
        result["planned_moves"] = self._extract_planned_moves(x_e)
        result["NumberOfMovements_rh"] = round(sum(
            self.network["move_cost_e"][move] * x_e[(move, t)].X
            for move in self.network["moves_e"]
            for t in range(T_exec)
        ))
        result["makespan_rh"] = round(max(
            [
                (t + 1) * x_a[(move, t)].X
                for move in self.network["moves_a"]
                for t in range(T_exec)
                if (move[0], move[1]) not in self.output_set
            ] or [0]
        ))
        result["flowtime_rh"] = round(sum(
            (t + 1) * x_a[(move, t)].X
            for move in self.network["moves_a"]
            for t in range(T_exec)
            if (move[2], move[3]) in self.output_set and (move[0], move[1]) != (move[2], move[3])
        ))
        result["obj_rh"] = model.ObjVal
        result["lb_rh"] = model.ObjBound
        result["ilp_gap"] = model.MIPGap if model.IsMIP else 0.0
        result["is_optimal_with_tolerance"] = (
            result["solver_status"] == GRB.OPTIMAL
            or result["ilp_gap"] <= self.config.optimality_tolerance
        )
        if self.config.warmstart_mode != "none":
            result["solution_snapshot"] = {
                "x_a": {
                    (move, t): x_a[(move, t)].X
                    for move in self.network["moves_a"]
                    for t in range(T + 1)
                    if x_a[(move, t)].X > 1e-9
                },
                "x_e": {
                    (move, t): x_e[(move, t)].X
                    for move in self.network["moves_e"]
                    for t in range(T + 1)
                    if x_e[(move, t)].X > 1e-9
                },
            }
        return result
