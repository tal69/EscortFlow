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
    bnc: bool = False
    full_all_integer_horizon: bool = False
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
        self.current_fractional_horizon = config.fractional_horizon
        self.current_integer_horizon = min(config.integer_horizon, config.fractional_horizon)
        self.supply_a_constraints = {}
        self.supply_e_constraints = {}
        self.retrieve_total_constraint = None
        self.bnc_stay_specs_by_t = None
        self.bnc_move_specs_by_t = None
        if self.config.full:
            self._pending_initial_model_build_time = 0.0
        else:
            self._pending_initial_model_build_time = self._build_persistent_model()

    def _configure_model_parameters(self, model):
        model.Params.OutputFlag = 0
        model.Params.TimeLimit = self.config.time_limit
        model.Params.MIPGap = self.config.optimality_tolerance
        model.Params.MIPFocus = self.config.mip_focus
        if self.config.bnc:
            model.Params.LazyConstraints = 1
            model.Params.PreCrush = 1
        if self.config.num_threads > 0:
            model.Params.Threads = self.config.num_threads

    def _dispose_model(self):
        if self.model is not None:
            self.model.dispose()
            self.model = None

    def _estimate_full_horizon(self, target_positions, escort_positions):
        if not target_positions:
            return self.config.epoch

        max_steps = max(
            1,
            (self.config.Lx + self.config.Ly) * len(target_positions) * 20 // max(len(escort_positions), 1),
        )
        makespan, _, _ = OneStepHeuristic_v2.SolveGreedy(
            self.config.Lx,
            self.config.Ly,
            self.output_set,
            set(target_positions),
            set(escort_positions),
            verbal=False,
            max_steps=max_steps,
            acyclic=self.config.acyclic,
        )
        return max(self.config.epoch, makespan)

    def _ensure_model_for_state(self, target_positions, escort_positions, count_build_time_as_pending=False):
        if self.config.full:
            fractional_horizon = self._estimate_full_horizon(target_positions, escort_positions)
            if self.config.full_all_integer_horizon:
                integer_horizon = fractional_horizon
            else:
                integer_horizon = min(self.config.integer_horizon, fractional_horizon)
        else:
            fractional_horizon = self.config.fractional_horizon
            integer_horizon = self.config.integer_horizon

        if (
            self.model is not None
            and self.current_fractional_horizon == fractional_horizon
            and self.current_integer_horizon == integer_horizon
        ):
            return

        self._dispose_model()
        self.current_fractional_horizon = fractional_horizon
        self.current_integer_horizon = integer_horizon
        build_time = self._build_persistent_model()
        if count_build_time_as_pending:
            self._pending_initial_model_build_time += build_time

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
        model = gp.Model("escort_flow_rh_v8", env=self.env)
        self._configure_model_parameters(model)

        T = self.current_fractional_horizon
        integer_horizon = self.current_integer_horizon
        x_a = {}
        x_e = {}
        for t in range(T + 1):
            a_vtype = GRB.BINARY if t <= integer_horizon else GRB.CONTINUOUS
            e_vtype = GRB.BINARY if t <= integer_horizon else GRB.CONTINUOUS
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

        # (1) Escort flow conservation.
        # (2) Target flow conservation away from outputs.
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

        # (3) Target loads stay once they reach an output cell.
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

        # (4) and (5) Supply at the initial locations of target loads and escorts.
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
                # (6) Generalized capacity constraint. - found to be redundant
                # model.addConstr(
                #     gp.quicksum(x_a[(move, t)] for move in self.network["outgoing_a"][loc]) +
                #     gp.quicksum(x_e[(move, t)] for move in self.network["outgoing_e"][loc]) <= 1,
                #     name=f"cap_{loc[0]}_{loc[1]}_t{t}",
                # )
                # (7) Escort conflict avoidance.
                model.addConstr(
                    gp.quicksum(x_e[(move, t)] for move in self.network["cell_cover"][loc]) <= 1,
                    name=f"cell_cover_{loc[0]}_{loc[1]}_t{t}",
                )
                if not self.config.bnc:
                    # (9) A target load must move if crossed by an escort.
                    model.addConstr(
                        1 - x_a[(stay_move, t)] >= gp.quicksum(
                            x_e[(move, t)] for move in self.network["cell_cover"][loc]
                        ),
                        name=f"stay_cover_{loc[0]}_{loc[1]}_t{t}",
                    )
                    # (8) A target load cannot move unless crossed by an escort.
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
            # (10) Every target load must eventually arrive at an output cell.
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
                # (15) q stores the summed arrival times at each output cell.
                model.addConstr(
                    gp.quicksum(
                        (t + 1) * x_a[(move, t)]
                        for move in incoming_output_moves
                        for t in range(T + 1)
                    ) == q_vars[output]
                )

        bnc_stay_specs_by_t = None
        bnc_move_specs_by_t = None
        if self.config.bnc:
            # Prepare separated representations of (8) and (9) for the callback.
            bnc_stay_specs_by_t = {t: [] for t in range(T + 1)}
            bnc_move_specs_by_t = {t: [] for t in range(T + 1)}
            for loc in self.network["locations"]:
                stay_move = self.network["stay_move"][loc]
                nonstay_target_moves = [
                    move for move in self.network["outgoing_a"][loc]
                    if move != stay_move
                ]
                for t in range(T + 1):
                    cover_keys = tuple((move, t) for move in self.network["cell_cover"][loc])
                    bnc_stay_specs_by_t[t].append(
                        (
                            (stay_move, t),
                            cover_keys,
                            x_a[(stay_move, t)] + gp.quicksum(x_e[key] for key in cover_keys),
                        )
                    )
                    for move in nonstay_target_moves:
                        move_cover_keys = tuple(
                            (escort_move, t)
                            for escort_move in self.network["move_cover"][move]
                        )
                        bnc_move_specs_by_t[t].append(
                            (
                                (move, t),
                                move_cover_keys,
                                x_a[(move, t)] - gp.quicksum(x_e[key] for key in move_cover_keys),
                            )
                        )

        model.update()
        self.model = model
        self.x_a = x_a
        self.x_e = x_e
        self.q_vars = q_vars
        self.supply_a_constraints = supply_a_constraints
        self.supply_e_constraints = supply_e_constraints
        self.retrieve_total_constraint = retrieve_total_constraint
        self.bnc_stay_specs_by_t = bnc_stay_specs_by_t
        self.bnc_move_specs_by_t = bnc_move_specs_by_t
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
            target_positions, escort_positions, 0, self.current_fractional_horizon
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
        tail_start = max(0, self.current_fractional_horizon - self.config.epoch + 1)
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
            tail_start, self.current_fractional_horizon
        )
        merged.update(greedy_tail)
        return merged

    def _fill_unset_warmstart_arcs_with_zero(self, start_values):
        dense_start = dict(start_values)
        for t in range(self.current_fractional_horizon + 1):
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
                for t in range(self.current_fractional_horizon + 1)
            )))
        return q_values

    def _finalize_warmstart_vector(self, start_values, source):
        if not start_values:
            return None

        start_values = {
            key: value
            for key, value in start_values.items()
            if key[2] <= self.current_fractional_horizon
        }

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

    def _cut_time_order(self, include_fractional):
        integer_limit = min(self.current_integer_horizon, self.current_fractional_horizon)
        if not include_fractional:
            return list(range(integer_limit + 1))
        return list(range(integer_limit + 1)) + list(range(integer_limit + 1, self.current_fractional_horizon + 1))

    def _effective_bnc_cut_limit(self):
        return 2 * self.current_fractional_horizon

    def _separate_movement_coupling(self, x_a_sol, x_e_sol, add_cut, max_cuts, include_fractional):
        if max_cuts <= 0:
            return 0

        tol = 1e-6
        added = 0

        for t in self._cut_time_order(include_fractional):
            for stay_key, cover_keys, stay_expr in self.bnc_stay_specs_by_t[t]:
                cover_value = sum(x_e_sol[key] for key in cover_keys)
                if x_a_sol[stay_key] + cover_value > 1 + tol:
                    add_cut(stay_expr <= 1)
                    added += 1
                    if added >= max_cuts:
                        return added

            for move_key, move_cover_keys, move_expr in self.bnc_move_specs_by_t[t]:
                if x_a_sol[move_key] <= tol:
                    continue
                move_cover_value = sum(x_e_sol[key] for key in move_cover_keys)
                if x_a_sol[move_key] > move_cover_value + tol:
                    add_cut(move_expr <= 0)
                    added += 1
                    if added >= max_cuts:
                        return added

        return added

    def _build_runtime_callback(self, x_a, x_e, gurobi_messages, callback_stats, capture_messages):
        use_bnc = self.config.bnc
        x_a_items = list(x_a.items()) if use_bnc else None
        x_e_items = list(x_e.items()) if use_bnc else None
        cut_limit = self._effective_bnc_cut_limit() if use_bnc else 0

        def callback(model, where):
            if capture_messages and where == GRB.Callback.MESSAGE:
                gurobi_messages.append(model.cbGet(GRB.Callback.MSG_STRING))
                return

            if not use_bnc:
                return

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
                    x_a_sol,
                    x_e_sol,
                    model.cbCut,
                    cut_limit,
                    include_fractional=False,
                )
                callback_stats["time"] += time.perf_counter() - callback_start

            elif where == GRB.Callback.MIPSOL:
                callback_start = time.perf_counter()
                x_a_values = model.cbGetSolution([var for _, var in x_a_items])
                x_e_values = model.cbGetSolution([var for _, var in x_e_items])
                x_a_sol = {key: value for (key, _), value in zip(x_a_items, x_a_values)}
                x_e_sol = {key: value for (key, _), value in zip(x_e_items, x_e_values)}
                self._separate_movement_coupling(
                    x_a_sol,
                    x_e_sol,
                    model.cbLazy,
                    cut_limit,
                    include_fractional=True,
                )
                callback_stats["time"] += time.perf_counter() - callback_start

        return callback

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

        self._ensure_model_for_state(
            target_positions,
            escort_positions,
            count_build_time_as_pending=True,
        )

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
        self._ensure_model_for_state(target_positions, escort_positions)
        model = self.model
        model.reset(1)
        self._update_dynamic_model_state(target_positions, escort_positions)
        T = self.current_fractional_horizon
        T_exec = self.config.epoch
        x_a = self.x_a
        x_e = self.x_e
        q = self.q_vars

        warmstart_applied = self._apply_warmstart_vector(model, x_a, x_e, q, warmstart_vector)
        model.update()

        gurobi_messages = []
        callback_stats = {"time": 0.0}
        callback = self._build_runtime_callback(
            x_a,
            x_e,
            gurobi_messages,
            callback_stats,
            capture_messages=warmstart_applied,
        )

        model_construction_time = (
            time.perf_counter() - model_build_start + self._pending_initial_model_build_time
        )
        self._pending_initial_model_build_time = 0.0
        if self.config.bnc or warmstart_applied:
            model.optimize(callback)
        else:
            model.optimize()

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
            "user_cut_time": callback_stats["time"],
            "fractional_horizon_used": self.current_fractional_horizon,
            "integer_horizon_used": self.current_integer_horizon,
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
