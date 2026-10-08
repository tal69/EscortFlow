from dataclasses import dataclass
import math
import time

import gurobipy as gp
from gurobipy import GRB

from static_lexicographic import build_lex_csv_suffix, solve_lexicographic

OPTIMALITY_TOLERANCE = 1e-4  # 0.01%
OBJECTIVE_CUTOFF_TOLERANCE = 1e-6


@dataclass(frozen=True)
class LoadFlowStaticGurobiConfig:
    Lx: int
    Ly: int
    output_cells: tuple
    move_method: str
    alpha: float
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


class LoadFlowStaticGurobiSolver:
    def __init__(self, config):
        from static_weighted_certification import validate_objective_mode
        validate_objective_mode(config)
        if config.lexicographic and config.lp:
            raise ValueError("Lexicographic optimization requires an integer model")
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

        nonstay_moves = [move for move in moves if move[:2] != move[2:]]

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

    def build_warmstart_from_trace(self, target_positions, escort_positions, T,
                                   target_move_history, escort_move_history):
        """Encode the same BM leave trace used by the escort-flow solver.

        A straight escort move shifts each load on its path by one cell in
        the opposite direction. Commodity 1 tracks target loads; commodity 2
        tracks all remaining loads. A target arriving after transition t is
        retrieved through q at t + 1. During that retrieval period the output
        stays unavailable, as in the greedy trace and escort-flow model.

        Store nonzero values sparsely; _apply_warmstart explicitly supplies
        zero for every other variable, so Gurobi receives a complete start.
        """
        if self.config.lp or self.config.move_method != "BM":
            raise ValueError("Heuristic warm starts require an integer BM load-flow model")
        if len(target_move_history) != len(escort_move_history):
            raise ValueError("Target and escort traces must have the same length")
        if T < len(escort_move_history):
            raise ValueError(
                f"Warm-start trace needs load-flow horizon at least {len(escort_move_history)}, got {T}"
            )
        targets, escorts = set(target_positions), set(escort_positions)
        locations = set(self.network["locations"])
        if targets & escorts or not (targets | escorts) <= locations:
            raise ValueError("Target and escort locations must be disjoint and inside the grid")
        occupied = {loc: 1 if loc in targets else 2 for loc in locations - escorts}
        warmstart = {"x": {}, "q": {}, "z": 0.0}

        for t in range(T + 1):
            retrieved = {loc for loc, commodity in occupied.items()
                         if commodity == 1 and loc in self.output_set}
            for loc in retrieved:
                warmstart["q"][(loc, t)] = 1.0
                warmstart["z"] = float(t)
                del occupied[loc]

            escort_moves = escort_move_history[t] if t < len(escort_move_history) else []
            target_moves = target_move_history[t] if t < len(target_move_history) else {}
            load_destinations = {}
            used_cells = set()
            for orig_x, orig_y, dest_x, dest_y in escort_moves:
                origin, destination = (orig_x, orig_y), (dest_x, dest_y)
                if origin in occupied or origin in retrieved:
                    raise ValueError(f"Trace escort is unavailable at {origin}, period {t}")
                if origin == destination:
                    continue
                if (orig_x != dest_x and orig_y != dest_y) or destination not in locations:
                    raise ValueError("Trace escort moves must be straight and inside the grid")
                # OneStep can return numpy integer coordinates, whose boolean
                # comparisons do not support direct subtraction.
                dx = int(dest_x > orig_x) - int(dest_x < orig_x)
                dy = int(dest_y > orig_y) - int(dest_y < orig_y)
                current = origin
                path = {current}
                while current != destination:
                    source = (current[0] + dx, current[1] + dy)
                    if source not in occupied:
                        raise ValueError(f"Trace crosses an empty cell at {source}, period {t}")
                    load_destinations[source] = current
                    path.add(source)
                    current = source
                if path & (used_cells | retrieved):
                    raise ValueError(f"Trace has conflicting movements or output service at period {t}")
                used_cells.update(path)

            expected_target_moves = {
                (source, dest) for source, dest in load_destinations.items()
                if occupied[source] == 1
            }
            if expected_target_moves != set(target_moves.values()):
                raise ValueError(f"Target and escort traces disagree at period {t}")

            next_occupied = {}
            for source, commodity in occupied.items():
                dest = load_destinations.get(source, source)
                if dest in next_occupied:
                    raise ValueError(f"Trace moves two loads to {dest}, period {t}")
                move = source + dest
                warmstart["x"][(move, t, commodity)] = 1.0
                next_occupied[dest] = commodity
            occupied = next_occupied

        if any(commodity == 1 for commodity in occupied.values()):
            raise ValueError("Warm-start trace does not retrieve every target within the horizon")
        return warmstart

    def build_warmstart_from_solution(self, result, T):
        """Retarget a weighted incumbent's horizon without changing retrievals."""
        if self.config.lp or self.config.move_method != "BM":
            raise ValueError("Solution transfer requires an integer BM load-flow model")
        if not result.get("has_solution"):
            raise ValueError("Solution transfer requires a saved incumbent")
        old_T = result["solution_horizon"]
        if T < 0 or int(T) != T or T < result["makespan"]:
            raise ValueError("The receiving horizon must cover every saved target retrieval")
        saved = result["solution_warmstart"]
        start = {"x": {key: value for key, value in saved["x"].items() if key[1] <= T},
                 "q": {key: value for key, value in saved["q"].items() if key[1] <= T},
                 "z": saved["z"],
                 "removed_post_retrieval_movements": sum(
                     1 for (move, t, _), value in saved["x"].items()
                     if t > T and value > .5 and move[:2] != move[2:])}
        if T == old_T:
            return start

        boundary = min(T, old_T)
        # All targets have left by boundary. Final-period blocker moves are
        # unnecessary and their endpoints lack capacity/anti-swap constraints
        # in the shorter model. Replace them with stays before extending.
        blockers = set()
        for (move, t, commodity), value in saved["x"].items():
            if t != boundary or value <= .5:
                continue
            if commodity == 1:
                raise ValueError("Saved solution leaves a target unretrieved at the horizon")
            blockers.add(move[:2])
            if move[:2] != move[2:]:
                start["removed_post_retrieval_movements"] += 1
            start["x"].pop((move, t, commodity), None)
        for t in range(boundary, T + 1):
            for loc in blockers:
                start["x"][(loc + loc, t, 2)] = 1.0
        return start

    @staticmethod
    def _apply_warmstart(x, q, z, warmstart):
        for key, var in x.items():
            var.Start = warmstart["x"].get(key, 0.0)
        for key, var in q.items():
            var.Start = warmstart["q"].get(key, 0.0)
        z.Start = warmstart["z"]

    def solve(self, target_positions, escort_positions, T, objective_cutoff=None, warmstart=None):
        if self.config.objective_mode != "legacy" and objective_cutoff is not None:
            raise ValueError("Weighted/certification modes cannot use an objective cutoff")
        if self.config.lexicographic and objective_cutoff is not None:
            raise ValueError("A weighted objective cutoff cannot be used with lexicographic optimization")
        if warmstart is not None and (self.config.lp or self.config.move_method != "BM"):
            raise ValueError("Heuristic warm starts require an integer BM load-flow model")
        target_set = set(target_positions)
        escort_set = set(escort_positions)
        blocking_set = set(self.network["locations"]) - target_set - escort_set

        solve_start = time.perf_counter()
        model = gp.Model("load_flow_static", env=self.env)
        model.Params.OutputFlag = 1
        if not self.config.lp:
            model.Params.MIPFocus = self.config.mip_focus
        model.Params.MIPGap = OPTIMALITY_TOLERANCE
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
        objective_expr = self.config.alpha * z + self.config.gamma * movement_expr + self.config.beta * flow_time_expr
        model.setObjective(objective_expr, GRB.MINIMIZE)

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

        if warmstart is not None:
            self._apply_warmstart(x, q, z, warmstart)

        if self.config.lexicographic or self.config.objective_mode != "legacy":
            def extract_solution():
                actual_makespan = max(
                    (t for output in self.output_cells for t in tr if q[(output, t)].X > 1e-6),
                    default=0,
                )
                result = {
                    "makespan": actual_makespan,
                    "animation_moves": self._extract_animation_moves(x, q, actual_makespan),
                }
                if self.config.objective_mode == "weighted_integer":
                    result["solution_horizon"] = T
                    result["solution_warmstart"] = {
                        "x": {key: 1.0 for key, var in x.items() if var.X > .5},
                        "q": {key: 1.0 for key, var in q.items() if var.X > .5},
                        "z": float(round(z.X)),
                    }
                return result

            if self.config.objective_mode != "legacy":
                from static_weighted_certification import solve_weighted_or_certificate
                flow_proof_context = None
                if self.config.flow_proof_extension_time_limit is not None:
                    flow_proof_context = dict(
                        targets=tuple(target_set), outputs=self.output_cells,
                        cell_count=self.config.Lx * self.config.Ly, escort_count=len(escort_set),
                        physical_horizon=T, weighted_time_limit=self.config.time_limit,
                        extension_time_limit=self.config.flow_proof_extension_time_limit)

                def extract_callback_metrics(callback_model):
                    retrievals = [(t, q[(output, t)]) for output in self.output_cells for t in tr]
                    values = callback_model.cbGetSolution([var for _, var in retrievals])
                    return dict(makespan=max((retrieval for (retrieval, _), value in zip(retrievals, values)
                                              if value > .5), default=0))

                return solve_weighted_or_certificate(
                    model, flow_time_expr, movement_expr, extract_solution,
                    status_name=self._status_name, solve_start=solve_start,
                    mode=self.config.objective_mode, target=self.config.certification_target,
                    weight_scale=self.config.weight_scale,
                    flow_proof_context=flow_proof_context,
                    extract_callback_metrics=extract_callback_metrics,
                )
            return solve_lexicographic(
                model,
                flow_time_expr,
                movement_expr,
                extract_solution,
                status_name=self._status_name,
                solve_start=solve_start,
                time_limit=self.config.time_limit,
                phase1_time_limit=self.config.phase1_time_limit,
                work_limit=self.config.work_limit,
            )

        if objective_cutoff is not None:
            cutoff_value = objective_cutoff + OBJECTIVE_CUTOFF_TOLERANCE
            if self.config.lp:
                model.addConstr(objective_expr <= cutoff_value)
            else:
                model.Params.Cutoff = cutoff_value

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
        if self.config.lexicographic:
            return build_lex_csv_suffix(result)
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
