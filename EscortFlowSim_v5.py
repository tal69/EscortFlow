# -------------------------------------------------------------------------------
# Name:        EscortFlowSim_v5
#
# Purpose:     Simulate dynamic PBS system with parallel retrival using our rolling horizon framework for
#              the escort flow ILP model, SBM/SLM  (formerly BM/LM-NBM),  multi-loads,
#              continue retrival mode
#              Poisson arrival of requests
#              Works only with SBM regime. Adding SLM should be simple but I never test it
#
# Author:      Tal Raviv,  talraviv@tau.ac.il
#
# Created:     3/12/2023, first version adapted from the load flow formulation implementation, warm start for single load
#              5/12/2023 Added animation export file and argparse parameters
#              8/12/2023 Added block movement support (running pbs_escorts_bm.mod)
#              9/12/2023 added using DP heuristic to create upper bound for BM with a single load (without warm-start for now)
#              11/12/2023 Implementation of the RH framework
#              13/12/2023 Block movement support for RH (running  escort_flow_bm_rh.mod)
#              15/12/2023 Dynamic simulation
#              16/12/2023 Add max_balls_in_air and queuing, collecting statistics in blocks and simulation warmup
#              17/12/2023 Real time version that adjust to the solution time
#              21/12/2023 Implementation of the simple heuristic for SBM
#              31/12/2023 Record lead time at each block,   apply heuristic if ilp gap is above some threshold (default 0.3)
#              13/1/2023 (v2), accept the requests while performing the solution, so requests that are
#                         arrive at the output cells accidentally are removed immediately.
#                         for this end we have to create a forecast of the state of the system at the beginning
#                         of each execution horizon.
#              18/1/2023 Keep record of the distance of the requests from the output cells upon their arrival
#              15/1/2026 (v3) back to a version without warm start (cleaner). Also limits number of threads to eight by default to fix mac sudio bug
#              26/1/2026 Option to run full static model at offline rolling horizon, report about actual_max_balls
#              2/3/2026  Stop if get stuck
#              4/3/2026  Collect real waiting times data
#              6/3/2026  v4 - applying a new heuristic from the file OneStepHeuristic.py
#              10/3/2026 applying v2 of one step heuristic
#              11/3/2026 save detailed results in "raw<time>.p" file
#              11/3/2026 v5 - remove block/warmup/CI estimation and use number_of_requests directly
#
# Copyright:   (c) Tal Raviv 2023, 2024, 2026
# Licence:     Free but please let me know that you are using it
# Depends on   PBSCom.py, escort_flow_lm_rh.mod, escort_flow_bm_rh.mod, SimpleHeuristic.py
#              Assumes oplrun is installed and on the path
#              Creates input for FlowAnimation (2023 version)
# -------------------------------------------------------------------------------
import copy
import random
import itertools
import subprocess
import pickle
import time
import os
import socket
import numpy as np
import argparse

import OneStepHeuristic_v2
from PBSCom import *
import CI_Calculation

parser = argparse.ArgumentParser()

parser.add_argument("-x", "--Lx", type=int, help="Horizontal dimension of the PBS unit", required=True)
parser.add_argument("-y", "--Ly", type=int, help="Vertical dimension of the PBS unit", required=True)
parser.add_argument("-O", "--output_cells", nargs='+',
                    help="List of output locations, can be built from ranges in PBSCom format. e.g., 0-4-2 0  can be used to express cells (0,0),(0,2),(0,4)",
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
                    help="Number of threads to be used by the MIP solver - larger is not always better (default 8)", default=8)
# parser.add_argument("-S", "--simulation_length", type=int,
#                     help="Number of periods in the not including time it takes to empty the system (default, 100 time steps)",
#                     default=100)
parser.add_argument("-S", "--number_of_requests", type=int,
                    help="Total number of requests in the simulation (default 4000)", default=1000)
parser.add_argument("-T", "--fractional_horizon", type=int,
                    help="Number of fractional periods in model (default: --exec_horizon + 4)",
                    default=None)
parser.add_argument("-I", "--integer_horizon", type=int,
                    help="Number of periods in the planning horizon represented by boolean variables (default: --exec_horizon)",
                    default=None)
parser.add_argument("-E", "--exec_horizon", type=int,
                    help="Number of periods in the execution horizon (default 5)",
                    default=5)
parser.add_argument("-t", "--time_limit", type=int,
                    help="Time limit for CPLEX calls in seconds (default: --exec_horizon)",
                    default=None)
parser.add_argument("-a", "--export_animation", action="store_true",
                    help="Export animation files, one for each instance")
parser.add_argument("-m", "--offline", action="store_true",
                    help="The actual computation time is "
                         "completely ignored. With -E 1 this should give the best lead time. But it is not realistic ")
parser.add_argument( "--static", action="store_true",
                    help="Use the full static model instead of the surogate one. Should be used with --offline and long --time_limit ")

parser.add_argument( "--greedy", action="store_true",
                    help="Use greedy heuristic only (ignore the --static flag and force --offline)")

parser.add_argument("-L", "--log", action="store_true",
                    help="Write report after each iteration into log file (sim_log<time>.txt), overwrite previous report")
parser.add_argument("-R", "--request_rate", type=float, help="Request arrival rate (default 0.1 request/time step)",
                    default=0.1)
parser.add_argument("-o", "--max_opt_gap", type=float, help="Maximum optimality gap to use ILP solution, otherwise use heuristic solution (default 0.4)",
                    default=1)
parser.add_argument("--seed", type=int, help="random seed (default 0)", default=0)
parser.add_argument("-H", "--header_line", action="store_true", help="Print header line in result csv file")
parser.add_argument("-F", "--full_model", action="store_true", help="run the full model instead of the surrogate model)")
parser.add_argument("-M", "--max_balls_in_air", type=int,
                    help="Maximum number of target loads that are considered concurrently, the rest are queued (default 5)",
                    default=5)
parser.add_argument('-q', '--queue_management', choices=['fifo', 'spt'],
                    help='Select queue management strategy when there are to many open requests to hande (more than '
                         'max_balls_in__air spt (default) or spt - closest load with most requests'
                         , default='spt')


args = parser.parse_args()

if args.integer_horizon is None:
    args.integer_horizon = args.exec_horizon
if args.time_limit is None:
    args.time_limit = args.exec_horizon
if args.fractional_horizon is None:
    args.fractional_horizon = args.exec_horizon + 4

if args.integer_horizon < args.exec_horizon:
    parser.error("--integer_horizon must be at least --exec_horizon")
if args.fractional_horizon < args.integer_horizon:
    parser.error("--fractional_horizon must be at least --integer_horizon")

result_csv_file = args.csv
csv_name_prefix = os.path.splitext(os.path.basename(result_csv_file))[0] or "sim"
file_export = "out.txt"
dat_file = "escort_flow_sim.dat"

gamma = args.gamma  # movement weight
time_limit = args.time_limit  # Time limit for cplex run (seconds)
Lx = args.Lx
Ly = args.Ly

if args.greedy:
    args.offline = True
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

sim_log_file = ""
if args.log:
    sim_log_file = f"{csv_name_prefix}_log{time.strftime('%Y-%m-%d_%H%M%S', time.localtime())}.txt"
    f = open(sim_log_file, "w")
    f.write(f"{args}\n")
    f.write("Start simulation logging ...\n")
    f.close()

pickle_file_name = f"{csv_name_prefix}_raw{time.strftime('%Y-%m-%d_%H%M%S', time.localtime())}.p"
curr_t = 0
np.random.seed(args.seed)
random.seed(args.seed + 1)

E = random.sample(Locations, args.escorts_num)
A_orig, E_orig = [], copy.copy(E)  # save for script file
number_of_requests = args.number_of_requests
max_simulation_length = (int(number_of_requests/args.request_rate)+ 2500)  # +2500 for a good measure to allow cool down period after the arrival of the last request
moves = [[] for _ in range((int(number_of_requests/args.request_rate)+ 2500))]
cpu_time = 0
total_lead_time = 0
non_optimal = 0
heuristic_sol = 0
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
enter_via_cell = []  # for the animation
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
start_move = np.full(number_of_requests, np.iinfo(np.int32).max, dtype=np.int32) # Time the request start to be considered to allow calculation of waiting times
orig_distance = np.zeros(number_of_requests, dtype=np.int32)  # the distance of the request from the output cell upon arrival


f = open(result_csv_file, 'a')
if args.header_line:
    header_cols = [
        "Machine Name", "Time Stamp", "version", "cpu_time", "Non optimal", "Greedy runs", "Max gap",
        "Algorithm Name", "Queue Management", "Seed", "Request Rate", "PBS Dimensions", "Output Cells",
        "Number of outputs", "Number of Escorts",
        "Fractional Horizon", "Integer Horizon", "Execution Horizon", "Time Limit", "Max Balls In Air",
        "Max Opt Gap", "Actual Max Balls",
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

    f.write(",".join(header_cols) + "\n")


f.close()  # we want to open and close the file anyway just to make sure that the file is available for writing.

req_count = 0
last_start_sol = 0

if args.log:
    f = open(sim_log_file, "a")
    f.write(f"\n{time.ctime()} :Instance {Lx}x{Ly}, O={O}, E={E}\n")
    f.close()

open_requests = []
next_execution_time = args.exec_horizon * 2  # start moving loads on after two execution horizons

dist_map = OneStepHeuristic_v2.build_dist_map(Lx, Ly, O)  # for the greedy heuristic

while True:
    if args.log:
        f = open(sim_log_file, "a")
        f.write(f"Time: {curr_t}: \n")

        if curr_t > max_simulation_length:
            f.write(f"Panic: stop because get stuck for more than {max_simulation_length} taktst\n")
            exit(1)
        f.close()

    # extract all requests that arrive at this takt
    while req_count < number_of_requests and arrivals[req_count] <= curr_t:
        orig_distance[req_count] = dist[load_loc[req2load[req_count]][0], load_loc[req2load[req_count]][1]]
        if orig_distance[req_count] > 0:  # request is not for a load currently located on an output
            requests_on_load[req2load[req_count]].append(req_count)
            open_requests.append(req_count)
            if args.log:
                f = open(sim_log_file, "a")
                f.write(
                    f"\trequest #{req_count} arrive, load:{req2load[req_count]}, currently @ {load_loc[req2load[req_count]]}\n")
                f.close()
            if args.export_animation:
                enter_via_cell.append(load_loc[req2load[req_count]])
                # if ll not in O:
                #     moves[curr_t].insert(0,(ll,ll))  # change the color of the load
        else:  # request happened to be for a load located on an output cell
            if args.export_animation:
                enter_via_cell.append(load_loc[req2load[req_count]])
            if args.log:
                f = open(sim_log_file, "a")
                f.write(
                    f"\trequest #{req_count} arrive, load:{req2load[req_count]}, currently @ {load_loc[req2load[req_count]]} - will be removed immediately\n")
                f.close()
            departures[req_count] = curr_t
            start_move[req_count] = curr_t
        req_count += 1

    # idle takt
    if not open_requests:
        if args.log:
            f = open(sim_log_file, "a")
            f.write(f"    Skipping an idle takt @ {curr_t}\n")
            f.close()
        #curr_t += 1
        idle_takt += 1

    # solve next execution horizon based on requests that where available at curt_t- execution_horizon
    # or in the case of offline, available now
    elif curr_t >= next_execution_time:
        E = set(Locations) - set(load_loc.values())
        # locations of target loads
        target_loads = set([])  # target loads
        if args.queue_management == 'fifo':
            for r in open_requests:
                if arrivals[r] <= curr_t - (1 - args.offline) * args.exec_horizon:
                    target_loads.add(req2load[r])
                if len(target_loads) >= args.max_balls_in_air:
                    break
        else: # spt
            for r in open_requests:
                if arrivals[r] <= curr_t - (1 - args.offline) * args.exec_horizon:
                    target_loads.add(req2load[r])
            target_loads = sorted(target_loads, key=lambda q: dist[load_loc[q][0], load_loc[q][1]] / len(requests_on_load[q]))[:args.max_balls_in_air]

        for l in target_loads:
            for r in requests_on_load[l]:
                start_move[r] = min(start_move[r], curr_t)  # update when item is considered by the MILP first

        A = set([load_loc[ll] for ll in target_loads])

        if set(A) - set(O):
            # run optimization only if we have target loads outside output cells

            old_A, old_E = copy.copy(A), copy.copy(E)
            actual_max_balls = max(actual_max_balls, len(A))

            if not args.greedy:
                next_execution_time = curr_t + args.exec_horizon
                sim_iter += 1
                f = open(dat_file, "w")
                f.write('file_export = "%s";\n' % file_export)
                f.write('time_limit = %d;\n' % time_limit)
                f.write('num_threads = %d;\n' % args.num_threads)
                f.write(f'distance_penalty={args.distance_penalty};\n')
                f.write(f'time_penalty={args.time_penalty};\n')
                f.write('gamma=%f;\n' % gamma)
                f.write('Lx=%d;\n' % Lx)
                f.write('Ly=%d;\n' % Ly)
                f.write(f'T={args.fractional_horizon};\n')
                f.write(f'T_exec={args.exec_horizon};\n')
                if not args.static:
                    f.write(f'T_int={args.integer_horizon};\n')
                f.write('E=%s;\n' % tuple_opl(E))
                f.write('A=%s;\n' % tuple_opl(A))
                f.write('O=%s;\n' % tuple_opl(O))
                f.write(f'retrieval_mode = "continue";\n')
                f.close()

                try:
                    if args.static:   # using the same model as v3 so I didn't change the names
                        subprocess.run(["oplrun", "escort_flow_bm_rh_static_v3.mod", dat_file], check=True)
                    else:
                        subprocess.run(["oplrun", "escort_flow_bm_rh_v3.mod", dat_file], check=True)
                except:
                    print("Panic: Could not solve the model")
                    if args.log:
                        f = open(sim_log_file, "a")
                        f.write("Panic: Could not solve the model")
                        f.close()
                    exit(1)
                else:
                    f = open("end_of_exec_horizon.txt", "r")
                    s = f.readlines()
                    f.close()
                    A, E = eval(s[0]), eval(s[1])
                    cplex_status, cpu_time_iter = int(s[2]), float(s[3])
                    obj_rh, lb_rh = float(s[7]), float(s[8])
                    if obj_rh != 0:
                        ilp_gap = (obj_rh - lb_rh) / obj_rh
                    else:
                        ilp_gap = 1
                    cpu_time += cpu_time_iter

                    makespan_rh, flowtime_rh, NumberOfMovements_rh = round(float(s[4])), round(float(s[5])), round(
                        float(s[6]))

                    if args.log:
                        f = open(sim_log_file, "a")

                        f.write(
                            f"\tfinish running cplex, iteration {sim_iter}, cplex status:{cplex_status} LB={lb_rh}, ob={obj_rh} "
                            f"flowtime={len(A) * args.exec_horizon + flowtime_rh}, open requests: {open_requests}\n"
                            f"{'****** ' if cpu_time_iter >= time_limit else ''} cpu_time={cpu_time_iter:.2f} "
                            f"Gap={100 * (obj_rh - lb_rh) / max(obj_rh, 1):.2f}%\n")
                        f.close()


            if args.greedy or cplex_status not in [1, 11, 101, 102, 127] or ilp_gap > args.max_opt_gap or (
                    old_A != [] and set(A) == set(old_A) and set(E) == set(old_E)):

                print("Did not solve ILP model because greedy= true or solution obtained for the model is infeasible or trivial\n running *** greedy *** heuristic")
                if args.log:
                    f = open(sim_log_file, "a")
                    f.write("Did not solve ILP model because greedy= true or solution obtained for the model is infeasible or trivial\n Running *** greedy *** heuristic\n")

                A, E = old_A, old_E  # no need to copy
                for i in range(args.exec_horizon):
                    A, E, one_step_move = OneStepHeuristic_v2.OneStep(Lx, Ly, set(O), set(A), set(E), dist_map)
                    moves[curr_t+i] = one_step_move
                    NumberOfMovements += len(one_step_move)

                A,E = list(A), list(E)
                heuristic_sol += 1

                if args.log:
                    f = open(sim_log_file, "a")
                    f.write(
                        f"      iteration {sim_iter}")
                    f.write(f"      old_E = {old_E}, old_A = {old_A}, open_requests = {open_requests}\n")
                    f.write(f"      E = {E}, A = {A}\n")
                    f.close()

                # print(f"obf_rh={obj_rh}    lb_rh={lb_rh}, oldA={old_A}    A={A}   old_E={old_E}  E={E}")
                # exit(1)
            else:

                # we use round(float(s[4])) and not int(s[4]) because in some rare cases it is non integer but very close, e.g., 4.999999999998
                # In such a case int(s[4]) will produce an error

                if cpu_time_iter >= time_limit:
                    non_optimal += 1

                # total_lead_time += (len(A) * args.exec_horizon + flowtime_rh)
                NumberOfMovements += NumberOfMovements_rh
                max_gap = max(max_gap, ilp_gap)
                if args.log:
                    f = open(sim_log_file, "a")
                    f.write(f"State after: E = {E}, A = {A}\n")

                    f.close()
                f = open(file_export)
                mvs = f.readlines()
                f.close()
                moves[curr_t:(curr_t+args.exec_horizon)] = eval(mvs[-1])

    # Apply the move already planned for the current takt, even if the next solve is scheduled for later.
    if moves[curr_t]:
        for (loc1, loc2) in moves[curr_t]:
            if loc1 != loc2:  # can happen that they are equal because of changing color for the animation
                load_loc[pbs[loc1[0], loc1[1]]] = loc2
                if args.log:
                    f = open(sim_log_file, "a")
                    f.write(
                        f"\tload {pbs[loc1[0], loc1[1]]} moves {loc1}->{loc2} \n")
                    f.close()

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
                    if args.log:
                        f = open(sim_log_file, "a")
                        f.write(
                            f"\trequest #{req} departs via ({loc[0]}, {loc[1]}) with load {leaving_load} at the end of the takt (time  {departures[req]})\n")
                        f.close()
                        #  For QA
                        if departures[req] - arrivals[req] < orig_distance[req]:
                            print("Panic: lead time is smaller than distance")
                            if args.log:
                                f = open(sim_log_file, "a")
                                f.write("Panic: lead time is smaller than distance")
                                f.close()
                            exit(1)

                requests_on_load[leaving_load] = []

    if curr_t > arrivals[-1] and not open_requests:
        print(f"Simulation end")
        print(f"Total lead time: {total_lead_time}")
        print(f"Mean lead time: {total_lead_time / req_count:.2f}")
        break

    if idle_takt > max_simulation_length:
        print("Error: simulation got stuck on idle state. Stopping and reporting results for debugging only")
        if args.log:
            f = open(sim_log_file, "a")
            f.write(
                "Error: simulation got stuck on idle state. Stopping and reporting results for debugging only\n")
            f.close()
        break

    curr_t += 1  # advance the clock one takt

arrivals = np.array(arrivals, dtype=int) # convert to np vector
lead_times = departures - arrivals
waiting_times = start_move - arrivals
flow_times = departures - start_move
excess_times = departures - arrivals - orig_distance

lead_stats = CI_Calculation.summarize_metric(lead_times)
waiting_stats = CI_Calculation.summarize_metric(waiting_times)
flow_stats = CI_Calculation.summarize_metric(flow_times)
excess_stats = CI_Calculation.summarize_metric(excess_times)



# just for debugging
if min(departures - arrivals) < 0 and args.log:
    print("Error: some arrivals occurs after departure")
    print(f"{np.where(departures - arrivals < 0)}")
    f = open(sim_log_file, "a")
    f.write("Error: some arrivals occurs after departure\n")
    f.write(f"{np.where(departures - arrivals < 0)}\n")
    f.close()

if max(start_move) == np.iinfo(np.int32).max:
    print("Error: the start moves time value for some request was not recorede")
    print(f"{np.where(start_move==np.iinfo(np.int32).max)}")

if args.greedy:
    alg_name = "Greedy"
else:
    if args.offline:
        alg_name = "offline-"
    else:
        alg_name = "realtime-"
    if args.static:
        alg_name += "static"
    else:
        alg_name += "surrogate"

    alg_name += f"Q{args.max_balls_in_air}"

f = open(result_csv_file, 'a')
script_version = f"{os.path.basename(__file__)} ({time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(__file__)))})"
machine_name = socket.gethostname()
row_vals = [
    machine_name, time.ctime(), script_version, f"{cpu_time:.2f}", non_optimal, heuristic_sol, f"{max_gap:.4f}",
    alg_name, args.queue_management, args.seed, args.request_rate, f"{Lx}x{Ly}", tuple_opl(O),
    len(O), len(E_orig), args.fractional_horizon, args.integer_horizon, args.exec_horizon, time_limit,
    args.max_balls_in_air, args.max_opt_gap, actual_max_balls,
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
f.write(",".join(map(str, row_vals)) + "\n")
f.close()


moves = moves[:curr_t+1]  # Trim trailing empty entries of the list

pickle.dump((alg_name, args.queue_management, args.seed, args.request_rate, args.fractional_horizon, args.integer_horizon,
             args.exec_horizon,time_limit, args.max_balls_in_air, args.max_opt_gap, Lx, Ly, O, E_orig,arrivals,
             departures, start_move, moves, actual_max_balls,non_optimal, max_gap,heuristic_sol), open(pickle_file_name,"wb"))

if args.export_animation:

    for i in range(number_of_requests):
        if enter_via_cell[i] not in O:
            moves[arrivals[i]].insert(0, (enter_via_cell[i], enter_via_cell[i]))

    pickle.dump((Lx, Ly, O, E_orig, A_orig, moves),
                # remove empty moves  periods at the end
                open(
                    f"script_sim_v5_{Lx}_{Ly}_{args.escorts_num}_{args.request_rate}_{number_of_requests}.p",
                    "wb"))
