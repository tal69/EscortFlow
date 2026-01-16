# -------------------------------------------------------------------------------
# Name:        EscortFlowSim_v2, second try, request are added and handled immediately
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
#
# Copyright:   (c) Tal Raviv 2023, 2024
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
import numpy as np
import argparse

import SimpleHeuristic
from PBSCom import *

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
parser.add_argument("-S", "--simulation_length", type=int,
                    help="Number of periods in the not including time it takes to empty the system (default, 100 time steps)",
                    default=100)
parser.add_argument("-T", "--fractional_horizon", type=int,
                    help="Number of fractional periods in model (default 10)",
                    default=10)
parser.add_argument("-I", "--integer_horizon", type=int,
                    help="Number of periods in the planning horizon represented by boolean variables (default 5)",
                    default=5)
parser.add_argument("-E", "--exec_horizon", type=int,
                    help="Number of periods in the execution horizon (default 5)",
                    default=5)
parser.add_argument("-t", "--time_limit", type=int, help="Time limit for CPLEX calls (default 10 sec.)", default=10)
parser.add_argument("-a", "--export_animation", action="store_true",
                    help="Export animation files, one for each instance")
parser.add_argument("-m", "--offline", action="store_true",
                    help="The actual computation time is "
                         "completely ignored. With -E 1 this should give the best lead time. But it is not realistic ")
parser.add_argument("-L", "--log", action="store_true",
                    help="Write report after each iteration into log file (sim_log<time>.txt), overwrite previous report")
parser.add_argument("-R", "--request_rate", type=float, help="Request arrival rate (default 0.1 request/time step)",
                    default=0.1)
parser.add_argument("-o", "--max_opt_gap", type=float, help="Maximum optimality gap to use ILP solution, otherwise use heuristic solution (default 0.3)",
                    default=0.4)
parser.add_argument("--seed", type=int, help="random seed (default 0)", default=0)
parser.add_argument("-H", "--header_line", action="store_true", help="Print header line in result csv file")
parser.add_argument("-M", "--max_balls_in_air", type=int,
                    help="Maximum number of loads that are retried concurrently, the rest are queued (default 5)",
                    default=5)
parser.add_argument("-B", "--block_length", type=int, help="Number of requests in each simulation block (default 50)", default=50)
parser.add_argument('-q', '--queue_management', choices=['fifo', 'spt'],
                    help='Select queue management strategy when there are to many open requests to hande (more than '
                         'max_balls_in__air fifo (default) or spt - closest load with most requests'
                         , default='fifo')
parser.add_argument("-w", "--warmup_arrivals", type=int,
                    help="Number of first loads to ignore when estimating mean lead time (default 10)", default=10)

args = parser.parse_args()
result_csv_file = args.csv
file_export = "out.txt"
dat_file = "escort_flow_sim.dat"

gamma = args.gamma  # movement weight
time_limit = args.time_limit  # Time limit for cplex run (seconds)
Lx = args.Lx
Ly = args.Ly

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
    sim_log_file = f"sim_log{time.strftime('%H%M%S', time.localtime())}.txt"
    f = open(sim_log_file, "w")
    f.write(f"{args}\n")
    f.write("Start simulation logging ...\n")
    f.close()

curr_t = 0
np.random.seed(args.seed)
random.seed(args.seed + 1)

E = random.sample(Locations, args.escorts_num)
A_orig, E_orig = [], copy.copy(E)  # save for script file
moves = [[] for _ in range((args.simulation_length + 1000))]
cpu_time = 0
total_lead_time = 0
non_optimal = 0
heuristic_sol = 0
sim_iter = 0
NumberOfMovements = 0
idle_takt = 0
max_gap = 0

requests_on_load = dict([])
load_loc = dict([])
pbs = np.zeros((Lx, Ly), dtype=int)
count = 0
for x in range(Lx):
    for y in range(Ly):
        if (x, y) not in E:
            count += 1
            pbs[x, y] = count
            load_loc[count] = (x, y)
            requests_on_load[count] = []

req_count = 0

arrivals = []
req2load = []
enter_via_cell = []  # for the animation
load_queue = []  # loads waiting to be handled when there is enough capacity

# create all arrivals in advance (for variability reduction)
for t in range(args.simulation_length):
    load_num = np.random.poisson(args.request_rate)
    for i in range(load_num):
        arrivals.append(t)
        req2load.append(random.randint(1, count))
number_of_requests = len(arrivals)

departures = np.zeros(number_of_requests, dtype=int)
orig_distance = np.zeros(number_of_requests, dtype=int)  # the distance of the request from the output cell upon arrival


f = open(result_csv_file, 'a')
if args.header_line:
    f.write(
        f"\ndate, version, Solution mode, Queue Management, Lx x Ly, #IOs, # Escorts, IOs, Request rate, gamma, "
        f"distance_penalty, time_penalty, seed , T fractional, T integer, T execution, cpu time_limit, max balls in air, max opt gap, "
        f"T sim, Warmup loads, T actual, Total lead time, # loads entered, mean lead time,  95% C.I,mean excess time, 95% C.I, #load movements, total CPU time, "
        f"Sim iterations, Idle takts, Non optimal iter, max gap, Heuristic solution, Steady state block, blocks,..\n")
f.close()  # we want to open and close the file anyway just to make sure that the file is available for writing.

req_count = 0
last_start_sol = 0

if args.log:
    f = open(sim_log_file, "a")
    f.write(f"\n{time.ctime()} :Instance {Lx}x{Ly}, O={O}, E={E}\n")
    f.close()

open_requests = []
next_execution_time = args.exec_horizon * 2  # start moving loads on after two execution horizons

while True:
    if args.log:
        f = open(sim_log_file, "a")
        f.write(
            f"Time: {curr_t}: \n")
        f.close()

    # extract all requests that arrive at this takt
    while req_count < number_of_requests and arrivals[req_count] <= curr_t:
        open_requests.append(req_count)
        requests_on_load[req2load[req_count]].append(req_count)
        orig_distance[req_count] = dist[load_loc[req2load[req_count]][0], load_loc[req2load[req_count]][1]]
        if args.log:
            f = open(sim_log_file, "a")
            f.write(
                f"\trequest #{req_count} arrive, load:{req2load[req_count]}, currently @ {load_loc[req2load[req_count]]}\n")
            f.close()
        if args.export_animation:
            if req_count == 10:
                print("hi")
            enter_via_cell.append(load_loc[req2load[req_count]])
            # if ll not in O:
            #     moves[curr_t].insert(0,(ll,ll))  # change the color of the load
        req_count += 1

    # idle takt
    if not open_requests:
        if args.log:
            f = open(sim_log_file, "a")
            f.write(f"    Skipping an idle takt @ {curr_t}\n")
            f.close()
        curr_t += 1
        idle_takt += 1
        # if curr_t < args.simulation_length:
        #     continue
    # solve next execution horizon based on requests that where available at curt_t- execution_horizon
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

        A = [load_loc[ll] for ll in target_loads]

        if set(A) - set(O):  # run optimization only if we have target loads outside output cells
            next_execution_time = curr_t + args.exec_horizon
            sim_iter += 1
            f = open(dat_file, "w")
            f.write('file_export = "%s";\n' % file_export)
            f.write('time_limit = %d;\n' % time_limit)
            f.write(f'distance_penalty={args.distance_penalty};\n')
            f.write(f'time_penalty={args.time_penalty};\n')
            f.write('gamma=%f;\n' % gamma)
            f.write('Lx=%d;\n' % Lx)
            f.write('Ly=%d;\n' % Ly)
            f.write(f'T={args.fractional_horizon};\n')
            f.write(f'T_int={args.integer_horizon};\n')
            f.write(f'T_exec={args.exec_horizon};\n')
            f.write('E=%s;\n' % tuple_opl(E))
            f.write('A=%s;\n' % tuple_opl(A))
            f.write('O=%s;\n' % tuple_opl(O))
            f.write(f'retrieval_mode = "continue";\n')
            f.write(f'warm_start = 0;\n')
            f.write('init_solA = {};\n')
            f.write('init_solE = {};\n')
            f.close()

            try:
                subprocess.run(["oplrun", "escort_flow_bm_rh.mod", dat_file], check=True)
            except:
                print("Panic: Could not solve the model")
                exit(1)
            else:
                old_A, old_E = copy.copy(A), copy.copy(E)

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

                if cplex_status not in [1, 11, 102] or ilp_gap > args.max_opt_gap or (
                        old_A != [] and set(A) == set(old_A) and set(E) == set(old_E)):
                    A_arr, E_arr, O_arr = np.array(old_A), np.array(list(old_E)), np.array(O)

                    for i in range(args.exec_horizon):
                        old_Aa, old_Ea = copy.copy(A_arr), copy.copy(E_arr)  # # just for debugging and QA
                        A_arr, E_arr, one_step_move = SimpleHeuristic.OneStep(Lx, Ly, O_arr, A_arr, E_arr, True)
                        moves[curr_t+i] = one_step_move
                        NumberOfMovements += len(one_step_move)

                        # Just for debugging and QA
                        if not SimpleHeuristic.test_step(Lx, Ly, old_Aa, old_Ea, A_arr, E_arr, one_step_move):
                            print("Panic")
                            exit(1)

                    A = list(map(tuple, A_arr))
                    E = list(map(tuple, E_arr))
                    heuristic_sol += 1

                    if args.log:
                        f = open(sim_log_file, "a")
                        f.write(
                            f"      iteration {sim_iter}, status:{cplex_status} Cplex failed - ====== applying simple heuristic"
                            f",  makespan={curr_t} =======\n")
                        f.write(f"      E = {old_E}, A = {old_A}, open_requests = {open_requests}\n")
                        f.write(f"      E = {E}, A = {A}\n")
                        f.close()
                else:
                    makespan_rh, flowtime_rh, NumberOfMovements_rh = round(float(s[4])), round(float(s[5])), round(float(s[6]))
                    # we use round(float(s[4])) and not int(s[4]) because in some rare cases it is non integer but very close, e.g., 4.999999999998
                    # In such a case int(s[4]) will produce an error

                    if cpu_time_iter > time_limit:
                        non_optimal += 1

                    # total_lead_time += (len(A) * args.exec_horizon + flowtime_rh)
                    NumberOfMovements += NumberOfMovements_rh
                    max_gap = max(max_gap, ilp_gap)
                    if args.log:
                        f = open(sim_log_file, "a")
                        f.write(
                            f"\tfinish running cplex, iteration {sim_iter}, status:{cplex_status} LB={lb_rh}, ob={obj_rh} "
                            f"flowtime={len(A) * args.exec_horizon + flowtime_rh},  makespan={curr_t}  "
                            f"{'****** ' if cpu_time_iter >= time_limit else ''} cpu_time={cpu_time_iter:.2f} "
                            f"Gap={100 * (obj_rh - lb_rh) / max(obj_rh, 1):.2f}%\n")
                        f.write(f"\t\tE = {E}, A = {A}\n")
                        f.close()
                    f = open(file_export)
                    mvs = f.readlines()
                    f.close()
                    moves[curr_t:(curr_t+args.exec_horizon)] = eval(mvs[-1])

    if moves[curr_t]:
        for (loc1, loc2) in moves[curr_t]:
            if loc1 != loc2:  # can happen because of changing color in animation
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

                    open_requests.remove(req)
                    if args.log:
                        f = open(sim_log_file, "a")
                        f.write(
                            f"\trequest #{req} departs via ({loc[0]}, {loc[1]}) with load {leaving_load} at the end of the takt (time  {departures[req]})\n")
                        f.close()
                        #  For QA
                        if departures[req] - arrivals[req] < orig_distance[req]:
                            print("Panic: lead time is smaller than distance")
                            exit(1)


                requests_on_load[leaving_load] = []

    if curr_t > args.simulation_length and not open_requests:
        print(f"Simulation end")
        print(f"Total lead time: {total_lead_time}")
        print(f"Mean lead time: {total_lead_time / req_count:.2f}  (including warmup period)")
        break

    if idle_takt > args.simulation_length:
        print("Error: simulation got stuck on idle state. Stopping and reporting results for debugging only")
        if args.log:
            f = open(sim_log_file, "a")
            f.write(
                "Error: simulation got stuck on idle state. Stopping and reporting results for debugging only\n")
            f.close()
        break

    curr_t += 1  # advance the clock one takt

if args.export_animation:
    moves = moves[:curr_t+1]
    for i in range(number_of_requests):
        if enter_via_cell[i] not in O:
            moves[arrivals[i]].insert(0, (enter_via_cell[i], enter_via_cell[i]))

    pickle.dump((Lx, Ly, O, E_orig, A_orig, moves),
                # remove empty moves  periods at the end
                open(
                    f"script_sim_v2_{Lx}_{Ly}_{args.escorts_num}_{args.request_rate}_{args.simulation_length}.p",
                    "wb"))

arrivals = np.array(arrivals, dtype=int)
if args.warmup_arrivals < number_of_requests-1:
    lead_time_std = np.std(departures[args.warmup_arrivals:] - arrivals[args.warmup_arrivals:], ddof=1)
    lead_time_mean = np.mean(departures[args.warmup_arrivals:] - arrivals[args.warmup_arrivals:])
    half_ci = 1.96 * lead_time_std / (number_of_requests - args.warmup_arrivals) ** 0.5
else:
    total_lead_time = 0  # just to indicate that we don't have statistics
    lead_time_mean = 0
    half_ci = 0
    lead_time_std = 0
# evaluate steady state -  by comparing lead time of last block of requests with the previous blocks
steady_state_block = 0
if number_of_requests > 2*args.block_length:
    num_of_blocks = number_of_requests // args.block_length
    last_block_lead_time = np.mean(departures[-args.block_length:] - arrivals[-args.block_length:])
    for i in range(1, num_of_blocks):
        curr_block_leadtime = np.mean(departures[(-(i+1)*args.block_length):-i*args.block_length] - arrivals[(-(i+1)*args.block_length):-i*args.block_length])
        if last_block_lead_time - curr_block_leadtime < (1.96 * args.block_length**0.5 )*lead_time_std:
            steady_state_block = num_of_blocks - i
else:
    num_of_blocks = 0
    last_block_lead_time = 0

# Calculate mean excess time
excess_time = departures[args.warmup_arrivals:] - arrivals[args.warmup_arrivals:] - orig_distance[args.warmup_arrivals:]
mean_excess_time = np.mean(excess_time)
std_excess_time = np.std(excess_time, ddof=1)
half_ci_excess_time = 1.96 * std_excess_time / (number_of_requests-args.warmup_arrivals) ** 0.5

# just for debugging
if min(departures - arrivals) < 0 and args.log:
    print("Error: some arrivals occurs after departure")
    print(f"{np.where(departures - arrivals < 0)}")
    f = open(sim_log_file, "a")
    f.write("Error: some arrivals occurs after departure\n")
    f.write(f"{np.where(departures - arrivals < 0)}\n")
    f.close()

f = open(result_csv_file, 'a')
f.write(
    f"{time.ctime()},v2, {'offline' if args.offline else 'real-time'},{args.queue_management},{Lx}x{Ly}, {len(O)}, {len(E_orig)}, "
    f"{tuple_opl(O)},{args.request_rate}, {gamma}, {args.distance_penalty}, {args.time_penalty},{args.seed},"
    f" {args.fractional_horizon}, {args.integer_horizon}, {args.exec_horizon},{time_limit}, {args.max_balls_in_air}, {args.max_opt_gap}, "
    f"{args.simulation_length}, {args.warmup_arrivals}, {curr_t}, {total_lead_time},{number_of_requests}, "
    f"{lead_time_mean:.3f}, {half_ci:.3f}, {mean_excess_time:.3f}, {half_ci_excess_time:.3f}, {NumberOfMovements},{cpu_time:.2f}, {sim_iter}, {idle_takt}, {non_optimal}, "
    f"{max_gap:.4f}, {heuristic_sol}, {steady_state_block}")

if number_of_requests > 2*args.block_length:
    for i in range(num_of_blocks-1, 0, -1):
        f.write(f" ,{np.mean(departures[(-(i+1)*args.block_length):-i*args.block_length] - arrivals[(-(i+1)*args.block_length):-i*args.block_length])}")
    f.write(f" ,{last_block_lead_time}")
f.write("\n")
f.close()