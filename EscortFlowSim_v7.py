"""Simulate dynamic PBS retrieval under a rolling-horizon control policy.

The script generates a stream of retrieval requests, runs either the MILP-based
controller or the greedy fallback heuristic over successive execution horizons,
and records both aggregate performance statistics and a raw pickle trace.

Timing convention:
- Requests with arrival time `t` enter the simulation at time `t`.
- In offline mode, such requests are immediately eligible for the next
  planning call. In online mode, they become eligible only after one full
  execution-horizon delay, so a request may wait through one horizon before it
  is first considered by the optimizer.
- A planning decision taken at time `t` produces moves for the next
  `epoch` simulated time steps, so a newly eligible request may still
  wait until the current solve finishes before its load starts moving.
- The move list for the current time step is then applied to the PBS state.
- Requests whose loads occupy an output cell after those moves depart at time
  `t + 1`.

Outputs:
- CSV summary line with configuration and performance statistics.
- Raw pickle containing the simulation trace needed by `CI_Calculation.py` and
  `PBSAnimation.py`.

Author:      Tal Raviv,  talraviv@tau.ac.il
Copyright:   (c) Tal Raviv 2023, 2024, 2026
Licence:     Free but please let me know that you are using it
Depends on   PBSCom.py, OneStepHeuristic_v2.py, gurobipy

"""

import copy
import random
import itertools
import pickle
import time
import os
import sys
import socket
import numpy as np
import argparse
import gurobipy as gp
from gurobipy import GRB

import OneStepHeuristic_v2
from PBSCom import *
import CI_Calculation

sim_log_file = ""
OPTIMALITY_TOLERANCE = 1e-4  # 0.01%


def log_message(*messages, echo_to_screen=False):
    """Write messages to the log file and optionally echo them to stdout."""
    if echo_to_screen:
        for message in messages:
            print(message, end="")

    if not args.log or not sim_log_file:
        return

    with open(sim_log_file, "a") as log_file:
        for message in messages:
            log_file.write(message)


def print_progress(current, total):
    """Print a simple one-line progress bar."""
    bar_width =100
    filled = bar_width if total == 0 else int(bar_width * current / total)
    bar = "#" * filled + "-" * (bar_width - filled)
    print(f"\rProgress [{bar}] {current}/{total}", end="", flush=True)
    if current >= total:
        print("")

parser = argparse.ArgumentParser()
default_num_threads = 0 if sys.platform == "darwin" else 0

parser.add_argument("-x", "--Lx", type=int, help="Horizontal dimension of the PBS unit", required=True)
parser.add_argument("-y", "--Ly", type=int, help="Vertical dimension of the PBS unit", required=True)
parser.add_argument("-O", "--output_cells", nargs='+',
                    help="List of output locations, can be built from ranges in PBSCom format. e.g., 0-4-2 can be used to express cells (0,0),(0,2),(0,4)",
                    default=['0', '0'])
parser.add_argument("-e", "--escorts_num", help="Number of escorts", type=int, default=8)

parser.add_argument("-f", "--csv", help="File name of the result csv file", default="sim_escort_flow.csv")
parser.add_argument("--gamma", type=float, help="Weight of the movements in the objective function (default 0.01)",
                    default=0.01)
parser.add_argument("--distance_penalty", type=float,
                    help="penalty for each unit of distance at each time unit in the MILP model (default, 1)",
                    default=1)
parser.add_argument("--time_penalty", type=float,
                    help="fixed penalty for each unit of time outside output cell in the MILP model (default, 1)",
                    default=1)
parser.add_argument("--num_threads", type=int,
                    help=f"Number of threads to be used by the MIP solver. 0 means use machine default (default {default_num_threads})", default=default_num_threads)
# On macOS, 8 threads is a safer default than the machine-wide solver default.

parser.add_argument("-S", "--number_of_requests", type=int,
                    help="Total number of requests in the simulation (default 4000)", default=1000)
parser.add_argument("-T", "--fractional_horizon", type=int,
                    help="Number of fractional periods in model (default: max(--epoch + 4, --integer_horizon))",
                    default=None)
parser.add_argument("-I", "--integer_horizon", type=int,
                    help="Number of periods in the planning horizon represented by boolean variables (default: --epoch)",
                    default=None)
parser.add_argument("-E", "--epoch", type=int,
                    help="Number of periods in the execution epoch (default 1)",
                    default=1)
parser.add_argument("-t", "--time_limit", type=int,
                    help="Time limit for Gurobi calls in seconds (default: --epoch)",
                    default=None)
parser.add_argument("-m", "--offline", action="store_true",
                    help="Offline rolling horizon: use target loads visible at the current decision time; cannot be combined with --greedy or --hybrid")
parser.add_argument("--full", action="store_true",
                    help="Use the full MILP instead of the surrogate model, in either realtime or offline mode")

parser.add_argument("--greedy", action="store_true",
                    help="Use the greedy heuristic only; ignores Gurobi and currently requires --epoch 1")
parser.add_argument("--hybrid", action="store_true",
                    help="Use rolling horizon as usual, but switch to a greedy epoch whenever the number of new open requests is at least the number of old open requests times --hybrid_ratio")
parser.add_argument("--hybrid_ratio", type=float, default=1.0,
                    help="Ratio for --hybrid: use greedy when number_of_new_open_requests >= number_of_old_open_requests * --hybrid_ratio (default 1.0). Must be positive.")
parser.add_argument("--acyclic", action="store_true",
                    help="Use greedy target order based on target seniority instead of the default dynamic distance-based order. This is inline with our acyclicty proof but it provde poorer result than the distance-based priority rule.")

parser.add_argument("-L", "--log", action="store_true",
                    help="Write report after each iteration into log file (sim_log<time>.txt), overwrite previous report")
parser.add_argument("-R", "--request_rate", type=float, help="Request arrival rate (default 0.1 request/time step)",
                    default=0.1)
parser.add_argument("-o", "--max_opt_gap", type=float, help="Maximum optimality gap to use ILP solution, otherwise use heuristic solution (default 0.4)",
                    default=1)
parser.add_argument("--seed", type=int, help="random seed (default 0)", default=0)
parser.add_argument("-H", "--header_line", action="store_true",
                    help="Print header line in result csv file; a new or empty file gets a header automatically")
parser.add_argument("-a", "--save_raw", action="store_true",
                    help="Save the raw pickle trace used by CI_Calculation.py and PBSAnimation.py")
parser.add_argument("-M", "--max_balls_in_air", type=int,
                    help="Maximum number of target loads that are considered concurrently; if omitted, all cells are eligible",
                    default=None)
parser.add_argument('-q', '--queue_management', choices=['fifo', 'spt'],
                    help='Select queue management strategy when there are to many open requests to hande (more than '
                         'max_balls_in__air spt (default) or spt - closest load with most requests'
                         , default='spt')


args = parser.parse_args()

if args.integer_horizon is None:
    args.integer_horizon = args.epoch
if args.time_limit is None:
    args.time_limit = args.epoch
if args.fractional_horizon is None:
    args.fractional_horizon = max(args.epoch + 4, args.integer_horizon)

if args.integer_horizon < args.epoch:
    parser.error("--integer_horizon must be at least --epoch")
if args.fractional_horizon < args.integer_horizon:
    parser.error("--fractional_horizon must be at least --integer_horizon")
if args.greedy and args.epoch != 1:
    parser.error("--greedy requires --epoch == 1")
if args.greedy and args.full:
    parser.error("--greedy cannot be combined with --full")
if args.greedy and args.hybrid:
    parser.error("--greedy cannot be combined with --hybrid")
if args.offline and args.greedy:
    parser.error("--offline cannot be combined with --greedy")
if args.offline and args.hybrid:
    parser.error("--offline cannot be combined with --hybrid")
if args.hybrid_ratio <= 0:
    parser.error("--hybrid_ratio must be positive. Use --greedy to run the heuristic at every step")

result_csv_file = args.csv
csv_name_prefix = os.path.splitext(os.path.basename(result_csv_file))[0] or "sim"

gamma = args.gamma  # movement weight
time_limit = args.time_limit  # Time limit for gurobi run (seconds)
Lx = args.Lx
Ly = args.Ly
max_balls_in_air_csv = "n/a" if args.max_balls_in_air is None else args.max_balls_in_air
if args.max_balls_in_air is None:
    args.max_balls_in_air = Lx * Ly

if args.greedy:
    args.max_balls_in_air = Lx*Ly

if len(args.output_cells) % 2:
    print(f"Error: output cell coordinate list must be of even length - {args.output_cells}")
    exit(1)

O = []
for i in range(len(args.output_cells) // 2):
    for x in str2range(args.output_cells[i * 2]):
        for y in str2range(args.output_cells[i * 2 + 1]):
            O.append((x, y))
            if x not in range(Lx) or y not in range(Ly):
                print(f"Error: escort ({x},{y}) not in range (0-{Lx - 1}) x (0-{Ly - 1})")
                exit(1)

if len(set(O)) < len(O):
    print(f"Error: All output cells location must be unique {sorted(O)}")
    exit(1)

Locations = sorted(set(itertools.product(range(Lx), range(Ly))))

""" Create distance table for queue priority """

dist = np.ones((Lx, Ly), dtype=np.int16) * (Lx + Ly)
dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)]

for (x, y) in O:
    dist[x, y] = 0

B_new = set(O)

while B_new:
    B = copy.copy(B_new)
    B_new = set([])
    for (x, y) in B:
        for d in range(4):
            x1, y1 = x + dirs[d][0], y + dirs[d][1]
            if x1 in range(Lx) and y1 in range(Ly):
                if dist[x1, y1] > dist[x, y] + 1:
                    dist[x1, y1] = dist[x, y] + 1
                    B_new.add((x1, y1))

if args.log:
    sim_log_file = f"{csv_name_prefix}_log{time.strftime('%Y-%m-%d_%H%M%S', time.localtime())}.txt"
    log_message(f"{args}\n", "Start simulation logging ...\n")

pickle_file_name = ""
if args.save_raw:
    pickle_file_name = f"{csv_name_prefix}_raw{time.strftime('%Y-%m-%d_%H%M%S', time.localtime())}.p"
curr_t = 0
np.random.seed(args.seed)
random.seed(args.seed + 1)

E = random.sample(Locations, args.escorts_num)
E_orig = copy.copy(E)
number_of_requests = args.number_of_requests
max_simulation_length = (int(number_of_requests/args.request_rate)+ 2500)  # +2500 for a good measure to allow cool down period after the arrival of the last request
moves = [[] for _ in range((int(number_of_requests/args.request_rate)+ 2500))]
cpu_time = 0
total_lead_time = 0
non_optimal = 0
heuristic_sol = 0
fallback_heuristic_sol = 0
hybrid_heuristic_sol = 0
sim_iter = 0
NumberOfMovements = 0
idle_takt = 0
max_gap = 0
actual_max_balls = 0

requests_on_load = dict([])
load_loc = dict([])
pbs = np.zeros((Lx, Ly), dtype=int)
number_of_loads = 0
for x in range(Lx):
    for y in range(Ly):
        if (x, y) not in E:
            number_of_loads += 1
            pbs[x, y] = number_of_loads
            load_loc[number_of_loads] = (x, y)
            requests_on_load[number_of_loads] = []

req_count = 0

arrivals = []
req2load = []
enter_via_cell = []  # request i sees its load enter via this cell; used only to reconstruct animation coloring
load_queue = []  # loads waiting to be handled when there is enough capacity

# create all arrivals in advance (for variability reduction)
count = 0
t = 0
while count < number_of_requests:
    load_num_curr_t = np.random.poisson(args.request_rate)
    if count + load_num_curr_t > number_of_requests:
        load_num_curr_t = number_of_requests - count
    for i in range(load_num_curr_t):
        arrivals.append(t)
        req2load.append(random.randint(1, number_of_loads))
    count += load_num_curr_t
    t += 1


departures = np.zeros(number_of_requests, dtype=np.int32)
# First time a request enters the control logic. This is used to compute
# waiting time as "arrival until first active handling attempt", not merely
# until physical movement starts.
start_move = np.full(number_of_requests, np.iinfo(np.int32).max, dtype=np.int32)
orig_distance = np.zeros(number_of_requests, dtype=np.int32)  # the distance of the request from the output cell upon arrival


header_cols = [
    "Machine Name", "Time Stamp", "version", "cpu_time", "Non optimal", "Greedy runs",
    "Fallback Greedy Runs", "Hybrid-Ratio Greedy Runs",
    "Algorithm Name", "Queue Management", "Seed", "Request Rate", "PBS Dimensions", "Output Cells",
    "Number of outputs", "Number of Escorts", "Number of Requests", "Simulation End Time",
    "Fractional Horizon", "Integer Horizon", "Epoch", "Time Limit", "Max Balls In Air",
    "Actual Max Balls", "Max actual opt gap", "Max allowed opt gap",
    "Lead Time Mean", "Lead Time CI Half Width 95%",
    "Waiting Time Mean", "Waiting Time CI Half Width 95%",
    "Flow Time Mean", "Flow Time CI Half Width 95%",
    "Excess time Mean", "Excess Time CI Half Width 95%",
    "Lead Time Deleted", "Lead Time Used", "Lead Time Batch Size", "Lead Time Number of Batches",
    "Lead Time Lag1 Autocorr", "Lead Time Lag1 Threshold", "Lead Time Batch Means Variance",
    "Waiting Time Deleted", "Waiting Time Used", "Waiting Time Batch Size", "Waiting Time Number of Batches",
    "Waiting Time Lag1 Autocorr", "Waiting Time Lag1 Threshold", "Waiting Time Batch Means Variance",
    "Flow Time Deleted", "Flow Time Used", "Flow Time Batch Size", "Flow Time Number of Batches",
    "Flow Time Lag1 Autocorr", "Flow Time Lag1 Threshold", "Flow Time Batch Means Variance",
    "Excess Time Deleted", "Excess Time Used", "Excess Time Batch Size", "Excess Time Number of Batches",
    "Excess Time Lag1 Autocorr", "Excess Time Lag1 Threshold", "Excess Time Batch Means Variance",
    "Log File Name", "Raw Pickle File Name"
]
expected_header = ",".join(header_cols)
write_header = args.header_line or not os.path.exists(result_csv_file) or os.path.getsize(result_csv_file) == 0
if not write_header:
    with open(result_csv_file, "r") as existing_csv:
        first_line = existing_csv.readline().strip()
    write_header = first_line != expected_header

f = open(result_csv_file, 'a')
if write_header:
    f.write(",".join(header_cols) + "\n")


f.close()  # we want to open and close the file anyway just to make sure that the file is available for writing.

req_count = 0
last_start_sol = 0
last_progress_served = -1

print(f"Running {os.path.basename(__file__)}")
for arg_name, arg_value in sorted(vars(args).items()):
    print(f"  {arg_name} = {arg_value}")

log_message(f"\n{time.ctime()} :Instance {Lx}x{Ly}, O={O}, E={E}\n")

open_requests = []

dist_map = OneStepHeuristic_v2.build_dist_map(Lx, Ly, O)  # for the greedy heuristic


def collect_target_loads(cutoff_time):
    """Return unique requested loads visible by `cutoff_time`, in queue priority order."""
    target_loads = []
    target_set = set()
    for r in open_requests:
        if arrivals[r] <= cutoff_time and req2load[r] not in target_set:
            target_set.add(req2load[r])
            target_loads.append(req2load[r])
    if args.queue_management == 'fifo':
        return target_loads
    return sorted(
        target_loads,
        key=lambda q: dist[load_loc[q][0], load_loc[q][1]] / len(requests_on_load[q])
    )


def collect_visible_requests(cutoff_time):
    """Return open requests visible by `cutoff_time`."""
    return [r for r in open_requests if arrivals[r] <= cutoff_time]


def build_gurobi_network():
    locations = list(Locations)
    output_set = set(O)
    not_outputs = [loc for loc in locations if loc not in output_set]

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
        loc: args.time_penalty + (Lx + Ly + 1) * args.distance_penalty
        for loc in locations
    }
    for loc in locations:
        for output in O:
            distance = abs(output[0] - loc[0]) + abs(output[1] - loc[1])
            location_penalty[loc] = min(
                location_penalty[loc],
                args.time_penalty + args.distance_penalty * distance,
            )
        if loc in output_set:
            location_penalty[loc] = 0

    for x, y in locations:
        move = (x, y, x, y)
        moves_a.append(move)
        outgoing_a[(x, y)].append(move)
        incoming_a[(x, y)].append(move)
        stay_move[(x, y)] = move

        if x < Lx - 1:
            move = (x, y, x + 1, y)
            na[(x, y)].append((x + 1, y))
            moves_a.append(move)
            outgoing_a[(x, y)].append(move)
            incoming_a[(x + 1, y)].append(move)
        if y < Ly - 1:
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

        for x1 in range(Lx):
            move = (x, y, x1, y)
            ne[(x, y)].append((x1, y))
            moves_e.append(move)
            move_cost_e[move] = abs(x - x1)
            outgoing_e[(x, y)].append(move)
            incoming_e[(x1, y)].append(move)
        for y1 in range(Ly):
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
        "output_set": output_set,
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


NETWORK = build_gurobi_network()


def schedule_greedy_epoch(start_time, candidate_loads, escort_positions):
    """Plan greedy moves for the next epoch from the current state."""
    global cpu_time, NumberOfMovements, heuristic_sol, actual_max_balls

    actual_max_balls = max(actual_max_balls, len(candidate_loads))
    for load_id in candidate_loads:
        for req in requests_on_load[load_id]:
            start_move[req] = min(start_move[req], start_time)

    A = {
        load_loc[load_id]: min(requests_on_load[load_id])
        for load_id in candidate_loads
    }
    E = set(escort_positions)
    heuristic_start_time = time.perf_counter()
    for i in range(args.epoch):
        A, E, one_step_move = OneStepHeuristic_v2.OneStep(
            Lx, Ly, set(O), A, set(E), dist_map, acyclic=args.acyclic
        )
        moves[start_time + i] = one_step_move
        NumberOfMovements += len(one_step_move)
    cpu_time += time.perf_counter() - heuristic_start_time
    heuristic_sol += 1
    return list(A.keys()), list(E)


def gurobi_status_name(status_code):
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


def extract_planned_moves(x_e_vars):
    planned_moves = []
    for t in range(args.epoch):
        one_step_moves = []
        for move in NETWORK["moves_e"]:
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


def run_gurobi_model(target_positions, escort_positions):
    """Solve the current rolling-horizon model directly in Gurobi."""
    model = gp.Model("escort_flow_rh_v7")
    model.Params.OutputFlag = 0
    model.Params.TimeLimit = time_limit
    model.Params.MIPGap = OPTIMALITY_TOLERANCE
    model.Params.MIPFocus = 1
    if args.num_threads > 0:
        model.Params.Threads = args.num_threads

    T = args.fractional_horizon
    T_exec = args.epoch
    output_set = NETWORK["output_set"]
    target_set = set(target_positions)
    escort_set = set(escort_positions)
    loads_to_retrieve = len(target_set - output_set)

    x_a = {}
    x_e = {}
    for t in range(T + 1):
        a_vtype = GRB.BINARY if args.full or t <= args.integer_horizon else GRB.CONTINUOUS
        e_vtype = GRB.BINARY if args.full or t <= args.integer_horizon else GRB.CONTINUOUS
        for move in NETWORK["moves_a"]:
            x_a[(move, t)] = model.addVar(lb=0.0, ub=1.0, vtype=a_vtype)
        for move in NETWORK["moves_e"]:
            x_e[(move, t)] = model.addVar(lb=0.0, ub=1.0, vtype=e_vtype)

    q = None
    if args.full:
        q = {output: model.addVar(lb=0.0, vtype=GRB.INTEGER) for output in O}

    model.update()

    movement_term = gp.quicksum(
        NETWORK["move_cost_e"][move] * x_e[(move, t)]
        for move in NETWORK["moves_e"]
        for t in range(T + 1)
    )
    if args.full:
        model.setObjective(gp.quicksum(q[output] for output in O) + gamma * movement_term, GRB.MINIMIZE)
    else:
        surrogate_term = gp.quicksum(
            NETWORK["location_penalty"][(move[2], move[3])] * x_a[(move, t)]
            for move in NETWORK["moves_a"]
            for t in range(T + 1)
        )
        model.setObjective(gamma * movement_term + surrogate_term, GRB.MINIMIZE)

    for t in range(1, T + 1):
        for loc in NETWORK["locations"]:
            model.addConstr(
                gp.quicksum(x_e[(move, t - 1)] for move in NETWORK["incoming_e"][loc]) ==
                gp.quicksum(x_e[(move, t)] for move in NETWORK["outgoing_e"][loc])
            )
        for loc in NETWORK["not_outputs"]:
            model.addConstr(
                gp.quicksum(x_a[(move, t - 1)] for move in NETWORK["incoming_a"][loc]) ==
                gp.quicksum(x_a[(move, t)] for move in NETWORK["outgoing_a"][loc])
            )

    for output in O:
        nonstay_output_moves = [
            move for move in NETWORK["outgoing_a"][output]
            if (move[0], move[1]) != (move[2], move[3])
        ]
        for t in range(T + 1):
            model.addConstr(gp.quicksum(x_a[(move, t)] for move in nonstay_output_moves) == 0)

    for loc in NETWORK["locations"]:
        supply_a = 1 if loc in target_set else 0
        supply_e = 1 if loc in escort_set else 0
        model.addConstr(gp.quicksum(x_a[(move, 0)] for move in NETWORK["outgoing_a"][loc]) == supply_a)
        model.addConstr(gp.quicksum(x_e[(move, 0)] for move in NETWORK["outgoing_e"][loc]) == supply_e)

    for loc in NETWORK["locations"]:
        stay_move = NETWORK["stay_move"][loc]
        nonstay_target_moves = [
            move for move in NETWORK["outgoing_a"][loc]
            if move != stay_move
        ]
        for t in range(T + 1):
            model.addConstr(
                gp.quicksum(x_a[(move, t)] for move in NETWORK["outgoing_a"][loc]) +
                gp.quicksum(x_e[(move, t)] for move in NETWORK["outgoing_e"][loc]) <= 1
            )
            model.addConstr(gp.quicksum(x_e[(move, t)] for move in NETWORK["cell_cover"][loc]) <= 1)
            model.addConstr(
                1 - x_a[(stay_move, t)] >= gp.quicksum(x_e[(move, t)] for move in NETWORK["cell_cover"][loc])
            )
            for move in nonstay_target_moves:
                model.addConstr(
                    x_a[(move, t)] <= gp.quicksum(x_e[(escort_move, t)] for escort_move in NETWORK["move_cover"][move])
                )

    if args.full:
        arrival_moves_to_outputs = [
            move for move in NETWORK["moves_a"]
            if (move[2], move[3]) in output_set and (move[0], move[1]) != (move[2], move[3])
        ]
        model.addConstr(
            gp.quicksum(x_a[(move, t)] for move in arrival_moves_to_outputs for t in range(T + 1)) == loads_to_retrieve
        )
        for output in O:
            incoming_output_moves = [
                move for move in NETWORK["incoming_a"][output]
                if (move[0], move[1]) != (move[2], move[3])
            ]
            model.addConstr(
                gp.quicksum((t + 1) * x_a[(move, t)] for move in incoming_output_moves for t in range(T + 1)) == q[output]
            )

    model.optimize()

    status_code = model.Status
    result = {
        "A": [],
        "E": [],
        "solver_status": status_code,
        "solver_status_name": gurobi_status_name(status_code),
        "sol_count": model.SolCount,
        "cpu_time_iter": model.Runtime,
        "makespan_rh": 0,
        "flowtime_rh": 0,
        "NumberOfMovements_rh": 0,
        "obj_rh": 0.0,
        "lb_rh": 0.0,
        "planned_moves": [],
        "ilp_gap": 1.0,
        "is_optimal_with_tolerance": False,
    }

    if model.SolCount <= 0:
        return result

    result["A"] = [
        (move[0], move[1]) for move in NETWORK["moves_a"] if x_a[(move, T_exec)].X > 0.99
    ]
    result["E"] = [
        (move[0], move[1]) for move in NETWORK["moves_e"] if x_e[(move, T_exec)].X > 0.99
    ]
    result["planned_moves"] = extract_planned_moves(x_e)
    result["NumberOfMovements_rh"] = round(sum(
        NETWORK["move_cost_e"][move] * x_e[(move, t)].X
        for move in NETWORK["moves_e"]
        for t in range(T_exec)
    ))
    result["makespan_rh"] = round(max(
        [(t + 1) * x_a[(move, t)].X for move in NETWORK["moves_a"] for t in range(T_exec) if (move[0], move[1]) not in output_set] or [0]
    ))
    result["flowtime_rh"] = round(sum(
        (t + 1) * x_a[(move, t)].X
        for move in NETWORK["moves_a"]
        for t in range(T_exec)
        if (move[2], move[3]) in output_set and (move[0], move[1]) != (move[2], move[3])
    ))
    result["obj_rh"] = model.ObjVal
    result["lb_rh"] = model.ObjBound
    result["ilp_gap"] = model.MIPGap if model.IsMIP else 0.0
    result["is_optimal_with_tolerance"] = (
        result["solver_status"] == GRB.OPTIMAL or result["ilp_gap"] <= OPTIMALITY_TOLERANCE
    )
    return result


def is_usable_model_solution(model_result, old_A, old_E):
    """Return whether a model result should be executed, and why if not."""
    if model_result["sol_count"] <= 0:
        return False, f"Gurobi status {model_result['solver_status_name']} produced no feasible solution"

    if model_result["ilp_gap"] > args.max_opt_gap:
        return False, (
            f"ILP gap {model_result['ilp_gap']:.4f} exceeds the allowed optimality gap threshold {args.max_opt_gap:.4f}"
        )

    if old_A != [] and set(model_result["A"]) == set(old_A) and set(model_result["E"]) == set(old_E):
        return False, "the model returned the same A/E state as the input state"

    return True, None

while True:
    served_requests = req_count - len(open_requests)
    if not args.log and served_requests != last_progress_served:
        print_progress(served_requests, number_of_requests)
        last_progress_served = served_requests

    if args.log:
        log_message(f"Time: {curr_t}: \n")

        if curr_t > max_simulation_length:
            log_message(f"Panic: stop because get stuck for more than {max_simulation_length} taktst\n")
            exit(1)

    # extract all requests that arrive at this takt
    while req_count < number_of_requests and arrivals[req_count] <= curr_t:
        orig_distance[req_count] = dist[load_loc[req2load[req_count]][0], load_loc[req2load[req_count]][1]]
        if orig_distance[req_count] > 0:  # request is not for a load currently located on an output
            requests_on_load[req2load[req_count]].append(req_count)
            open_requests.append(req_count)
            log_message(
                f"\trequest #{req_count} arrive, load:{req2load[req_count]}, currently @ {load_loc[req2load[req_count]]}\n")
            enter_via_cell.append(load_loc[req2load[req_count]])
            # if ll not in O:
            #     moves[curr_t].insert(0,(ll,ll))  # change the color of the load
        else:  # request happened to be for a load located on an output cell
            enter_via_cell.append(load_loc[req2load[req_count]])
            log_message(
                f"\trequest #{req_count} arrive, load:{req2load[req_count]}, currently @ {load_loc[req2load[req_count]]} - will be removed immediately\n")
            departures[req_count] = curr_t
            start_move[req_count] = curr_t
        req_count += 1

    # idle takt
    if not open_requests:
        log_message(f"    Skipping an idle takt @ {curr_t}\n")
        idle_takt += 1

    elif curr_t % args.epoch == 0:
        E = set(Locations) - set(load_loc.values())
        old_visible_requests = collect_visible_requests(curr_t - args.epoch)
        current_visible_requests = collect_visible_requests(curr_t)
        old_target_loads = collect_target_loads(curr_t - args.epoch)
        current_time_target_loads = collect_target_loads(curr_t)
        old_visible_request_set = set(old_visible_requests)
        new_visible_requests = [r for r in current_visible_requests if r not in old_visible_request_set]
        hybrid_ratio_active = args.hybrid and len(new_visible_requests) >= len(old_visible_requests) * args.hybrid_ratio

        if args.greedy:
            if current_time_target_loads:
                schedule_greedy_epoch(curr_t, current_time_target_loads, E)
        elif hybrid_ratio_active:
            if current_time_target_loads:
                log_message(
                    f"\tHybrid mode: using greedy because new open requests={len(new_visible_requests)} >= old open requests={len(old_visible_requests)} * ratio {args.hybrid_ratio}\n"
                )
                greedy_A, greedy_E = schedule_greedy_epoch(curr_t, current_time_target_loads, E)
                hybrid_heuristic_sol += 1
                log_message(f"Greedy forecast after epoch: E = {greedy_E}, A = {greedy_A}\n")

        else:  # run the ILP model
            eligible_target_loads = current_time_target_loads if args.offline else old_target_loads
            target_loads = eligible_target_loads[:args.max_balls_in_air]

            for l in target_loads:
                for r in requests_on_load[l]:
                    start_move[r] = min(start_move[r], curr_t)

            A = set(load_loc[ll] for ll in target_loads)
            if set(A) - set(O):
                old_A, old_E = copy.copy(A), copy.copy(E)
                actual_max_balls = max(actual_max_balls, len(A))
                sim_iter += 1
                model_result = run_gurobi_model(A, E)
                A, E = model_result["A"], model_result["E"]
                cpu_time_iter = model_result["cpu_time_iter"]
                ilp_gap = model_result["ilp_gap"]
                cpu_time += cpu_time_iter

                log_message(
                    f"\tfinish running gurobi, iteration {sim_iter}, status:{model_result['solver_status_name']} LB={model_result['lb_rh']}, ob={model_result['obj_rh']} "
                    f"flowtime={len(A) * args.epoch + model_result['flowtime_rh']}, open requests: {open_requests}\n"
                    f"{'****** ' if model_result['solver_status'] == GRB.TIME_LIMIT and not model_result['is_optimal_with_tolerance'] else ''} cpu_time={cpu_time_iter:.2f} "
                    f"Gap={100 * model_result['ilp_gap']:.4f}%\n")

                use_model_solution, unusable_reason = is_usable_model_solution(model_result, old_A, old_E)
                if use_model_solution:
                    if not model_result["is_optimal_with_tolerance"]:
                        non_optimal += 1
                    NumberOfMovements += model_result["NumberOfMovements_rh"]
                    max_gap = max(max_gap, ilp_gap)
                    log_message(f"State after: E = {E}, A = {A}\n")
                    if model_result["planned_moves"]:
                        for i, one_step_move in enumerate(model_result["planned_moves"][:args.epoch]):
                            moves[curr_t + i] = one_step_move
                else:
                    log_message(
                        f"Did not use ILP model solution because {unusable_reason}\n"
                        "Running *** greedy *** heuristic on all target loads visible at the current decision time\n",
                        echo_to_screen=True,
                    )
                    greedy_A, greedy_E = schedule_greedy_epoch(curr_t, current_time_target_loads, old_E)
                    fallback_heuristic_sol += 1
                    log_message(f"Greedy forecast after epoch: E = {greedy_E}, A = {greedy_A}\n")

    # Apply the move already planned for the current takt, even if the next solve is scheduled for later.
    if moves[curr_t]:
        for (loc1, loc2) in moves[curr_t]:
            if loc1 != loc2:  # can happen that they are equal because of changing color for the animation
                load_loc[pbs[loc1[0], loc1[1]]] = loc2
                log_message(
                    f"\tload {pbs[loc1[0], loc1[1]]} moves {loc1}->{loc2} \n")

        pbs = np.zeros((Lx, Ly), dtype=int)
        for l, loc in load_loc.items():
            pbs[loc[0], loc[1]] = l

        for loc in O:
            if pbs[loc[0], loc[1]] != 0:
                leaving_load = pbs[loc[0], loc[1]]
                for req in requests_on_load[leaving_load]:
                    departures[req] = curr_t + 1
                    total_lead_time += (departures[req] - arrivals[req])
                    start_move[req] = min(curr_t+1,start_move[req])  # just in case the request arrive by chance before starting to be handled by the optimization

                    open_requests.remove(req)
                    log_message(
                        f"\trequest #{req} departs via ({loc[0]}, {loc[1]}) with load {leaving_load} at the end of the takt (time  {departures[req]})\n")
                    #  For QA
                    if departures[req] - arrivals[req] < orig_distance[req]:
                        log_message("Panic: lead time is smaller than distance\n", echo_to_screen=True)
                        exit(1)

                requests_on_load[leaving_load] = []

    if curr_t > arrivals[-1] and not open_requests:
        if not args.log:
            print_progress(number_of_requests, number_of_requests)
        log_message(
            "Simulation end\n",
            f"Total lead time: {total_lead_time}\n",
            f"Mean lead time: {total_lead_time / req_count:.2f}\n",
            echo_to_screen=True,
        )
        break

    if idle_takt > max_simulation_length:
        log_message(
            "Error: simulation got stuck on idle state. Stopping and reporting results for debugging only\n",
            echo_to_screen=True,
        )
        break

    curr_t += 1  # advance the clock one takt

arrivals = np.array(arrivals, dtype=int) # convert to np vector
simulation_end_time = int(np.max(departures)) if departures.size else 0
late_departure_indices = np.flatnonzero(departures > arrivals[-1])
if late_departure_indices.size > 0:
    cooldown_cutoff_arrival_time = arrivals[late_departure_indices[0]]
    cooldown_trim_mask = arrivals <= cooldown_cutoff_arrival_time
else:
    cooldown_trim_mask = np.ones_like(arrivals, dtype=bool)

lead_times = (departures - arrivals)[cooldown_trim_mask]
waiting_times = (start_move - arrivals)[cooldown_trim_mask]
flow_times = (departures - start_move)[cooldown_trim_mask]
excess_times = (departures - arrivals - orig_distance)[cooldown_trim_mask]

lead_stats = CI_Calculation.summarize_metric(lead_times)
waiting_stats = CI_Calculation.summarize_metric(waiting_times)
flow_stats = CI_Calculation.summarize_metric(flow_times)
excess_stats = CI_Calculation.summarize_metric(excess_times)

for metric_name, stats in [
    ("Lead Time", lead_stats),
    ("Waiting Time", waiting_stats),
    ("Flow Time", flow_stats),
    ("Excess Time", excess_stats),
]:
    if stats["res"] is None:
        log_message(
            f"{metric_name} steady-state summary skipped: {stats.get('error', 'unknown error')}\n",
            echo_to_screen=True,
        )


# just for debugging
if min(departures - arrivals) < 0 and args.log:
    log_message(
        "Error: some arrivals occurs after departure\n",
        f"{np.where(departures - arrivals < 0)}\n",
        echo_to_screen=True,
    )

if max(start_move) == np.iinfo(np.int32).max:
    print("Error: the start moves time value for some request was not recorede")
    print(f"{np.where(start_move==np.iinfo(np.int32).max)}")

if args.greedy:
    alg_name = "greedy"
else:
    horizon_suffix = "-PLPR" if args.fractional_horizon > args.integer_horizon else ""
    alg_prefix = f"hybrid{args.hybrid_ratio}-" if args.hybrid else ""
    if not alg_prefix:
        alg_prefix = "offline-" if args.offline else ""
    alg_name = alg_prefix
    if not args.full:
        alg_name += "surrogate"

    if max_balls_in_air_csv != "n/a":
        alg_name += f"-Q{args.max_balls_in_air}"
    alg_name += horizon_suffix

if alg_name == "":
    alg_name = "-"

f = open(result_csv_file, 'a')
script_version = f"{os.path.basename(__file__)} ({time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(__file__)))})"
machine_name = f"{socket.gethostname()}-T{args.num_threads}"
row_vals = [
    machine_name, time.ctime(), script_version, f"{cpu_time:.2f}", non_optimal, heuristic_sol,
    fallback_heuristic_sol, hybrid_heuristic_sol,
    alg_name, args.queue_management, args.seed, args.request_rate, f"{Lx}x{Ly}", tuple_opl(O),
    len(O), len(E_orig), args.number_of_requests, simulation_end_time,
    args.fractional_horizon, args.integer_horizon, args.epoch, time_limit,
    max_balls_in_air_csv, actual_max_balls, f"{max_gap:.4f}", args.max_opt_gap,
    lead_stats["mean"], lead_stats["half_width_95"],
    waiting_stats["mean"], waiting_stats["half_width_95"],
    flow_stats["mean"], flow_stats["half_width_95"],
    excess_stats["mean"], excess_stats["half_width_95"],
    lead_stats["obs_deleted"], lead_stats["obs_used"], lead_stats["batch_size"], lead_stats["num_batches"],
    lead_stats["lag1_autocorr"], lead_stats["lag1_threshold"], lead_stats["batch_means_variance"],
    waiting_stats["obs_deleted"], waiting_stats["obs_used"], waiting_stats["batch_size"], waiting_stats["num_batches"],
    waiting_stats["lag1_autocorr"], waiting_stats["lag1_threshold"], waiting_stats["batch_means_variance"],
    flow_stats["obs_deleted"], flow_stats["obs_used"], flow_stats["batch_size"], flow_stats["num_batches"],
    flow_stats["lag1_autocorr"], flow_stats["lag1_threshold"], flow_stats["batch_means_variance"],
    excess_stats["obs_deleted"], excess_stats["obs_used"], excess_stats["batch_size"], excess_stats["num_batches"],
    excess_stats["lag1_autocorr"], excess_stats["lag1_threshold"], excess_stats["batch_means_variance"],
    sim_log_file, pickle_file_name,
]
f.write(",".join(CI_Calculation.csv_cell(v) for v in row_vals) + "\n")
f.close()


if args.save_raw:
    moves = moves[:curr_t+1]  # Trim trailing empty entries of the realized move list

    for i in range(number_of_requests):
        if enter_via_cell[i] not in O:
            # A self-move is inserted at the request arrival time so animation can
            # recolor that load as a target even when the raw pickle, rather than a
            # dedicated animation file, is used as input.
            moves[arrivals[i]].insert(0, (enter_via_cell[i], enter_via_cell[i]))

    # This is consumed by CI_Calculation.py and can also be read directly by
    # PBSAnimation.py.
    pickle.dump((alg_name, args.queue_management, args.seed, args.request_rate, args.fractional_horizon, args.integer_horizon,
                 args.epoch,time_limit, args.max_balls_in_air, args.max_opt_gap, Lx, Ly, O, E_orig,arrivals,
                 departures, start_move, moves, actual_max_balls,non_optimal, max_gap,heuristic_sol), open(pickle_file_name,"wb"))
