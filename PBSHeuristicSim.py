# -------------------------------------------------------------------------------
# Name:        EscortFlowSimHeuristic
# Purpose:     Simulate dynamic PBS system with parallel retrival using our simple heuristic , SBM/SLM  (formerly BM/LM-NBM),  multi-loads,
#              continue retrival mode. Works only on edge configurations
#              Poisson arrival of requests
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
#              22/12/2023  Creating a version that exploit the one step heuristic in real time
#
# Copyright:   (c) Tal Raviv 2020, 2023
# Licence:     Free but please let me know that you are using it
# Depends on   PBSCom.py, SimpleHeuristic.py
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
parser.add_argument("-O", "--output_cells", nargs='+', help="List of output locations, can be built from ranges in PBSCom format. e.g., 0-4-2 0  can be used to express cells (0,0),(0,2),(0,4)", default=[0, 0])
parser.add_argument("-e", "--escorts_num", help="Number of escorts", type=int, default=8)
parser.add_argument('-m', '--retrieval_mode', choices=['continue'],
                    help='Select a retrieval mode. Default is stay. Also note that for stay mode the number '
                         'of output cells must be greater or equal the number target loads', default='continue')
parser.add_argument("-f", "--csv", help="File name of the result csv file", default="sim_escort_flow_heuristic.csv")
parser.add_argument("--gamma", type=float, help="Weight of the movements in the objective function (default 0.01)",
                    default=0.01)

parser.add_argument("-S", "--simulation_length", type=int,
                    help="Number of periods in the not including time it takes to empty the system (default, 100 time steps)",
                    default=100)
parser.add_argument("-a", "--export_animation", action="store_true",
                    help="Export animation files, one for each instance")
parser.add_argument("-b", "--bm", action="store_true", help="Simultaneous block movement regime (otherwise SLM, aka LM, NBM)")
parser.add_argument("-L", "--log", action="store_true",
                    help="Write report after each iteration into log file (sim_log<time>.txt), overwrite previous report")
parser.add_argument("-R", "--request_rate", type=float, help="Request arrival rate (default 0.1 request/time step)",
                    default=0.1)
parser.add_argument("--seed", type=int, help="random seed (default 0)", default=0)
parser.add_argument("-H", "--header_line", action="store_true", help="Print header line in result csv file")
parser.add_argument("-M", "--max_balls_in_air", type=int,
                    help="Maximum number of loads that are retried concurrently, the rest are queued (default 5)  - DONT SURE WE NEED IT HERE",
                    default=5)
parser.add_argument("-B", "--block_length", type=int, help="Number of time step in each simulation block (default 50)",
                    default=50)
parser.add_argument("-w", "--warmup_arrivals", type=int,
                    help="Number of first loads to ignore when estimating mean lead time (default 10)", default=10)
parser.add_argument("-P", "--print_blocks", action="store_true", help="Print block data")

args = parser.parse_args()
result_csv_file = args.csv
file_export = "out.txt"
dat_file = "escort_flow_sim.dat"

Lx = args.Lx
Ly = args.Ly


if len( args.output_cells) % 2:
    print(f"Error: output cell coordinate list must be of even length - {args.output_cells}")
    exit(1)

O = []
for i in range(len(args.output_cells) // 2):
    for x in str2range( args.output_cells[i * 2]):
        for y in str2range( args.output_cells[i * 2+1]):
            O.append((x, y))
            if x not in range(Lx) or y not in range(Ly):
                print(f"Error: escort ({x},{y}) not in range (0-{Lx-1}) x (0-{Ly-1})")
                exit(1)

if len(set(O)) < len(O):
    print(f"Error: All output cells location must be unique {sorted(O)}")
    exit(1)

if args.simulation_length % args.block_length != 0:
    print(
        f"Error: Simulation length ({args.simulation_length}) must be an integer multiplication of block length ({args.block_length})")
    exit(1)

num_of_blocks = args.simulation_length // args.block_length

block_max_queue = np.zeros(num_of_blocks, dtype=int)

Locations = sorted(set(itertools.product(range(Lx), range(Ly))))

if args.log:
    sim_log_file = f"sim_heuristic_log{time.strftime('%H%M%S', time.localtime())}.txt"
    f = open(sim_log_file, "w")
    f.write("Start simulation logging ...\n")
    f.close()

curr_t = 0
np.random.seed(args.seed)
random.seed(args.seed + 1)

E = random.sample(Locations, args.escorts_num)
A = []
A_orig, E_orig = copy.copy(A), copy.copy(E)  # save for script file
moves = []
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
load_queue = []  # loads waiting to be handled when there is enough capacity

# create all arrivals in advance (for variability reduction)
for t in range(args.simulation_length):
    load_num = np.random.poisson(args.request_rate)
    for i in range(load_num):
        arrivals.append(t)
        req2load.append(random.randint(1, count))
number_of_requests = len(arrivals)

departures = np.zeros(number_of_requests, dtype=int)

f = open(result_csv_file, 'a')
if args.header_line:
    f.write(
        f"\ndate, SBM/SLM, Retrieval Mode, Lx x Ly, #IOs, # Escorts, IOs, Request rate, "
        f" seed , max balls in air, T sim, Warmup loads,  "
        f"T actual, Total lead time, # loads entered, mean lead time,  95% C.I, #load movements, total CPU time, "
        f"Sim iterations, Idle takts")
    if args.print_blocks:
        for i in range(num_of_blocks):
            f.write(f", B{i} - max queue")
    f.write("\n")
f.close()  # we want to open and close the file anyway just to make sure that the file is available for writing.

req_count = 0
last_start_sol = 0

if args.log:
    f = open(sim_log_file, "a")
    f.write(f"\n{time.ctime()} :Instance {Lx}x{Ly}, O={O}, E={E}, A={A}\n")
    f.close()

while True:
    if curr_t % 100 == 0:
        print(f"time: {curr_t}")
    new_target_loads_xy = set([])
    while load_queue and len(A) < args.max_balls_in_air:
        load_from_queue = load_queue.pop(0)
        new_target_loads_xy.add(load_loc[load_from_queue])  # just for the animation
        A.append(load_loc[load_from_queue])
        # CAN BE IMPROVED BY SERVING THE LOAD WITH MOST REQUESTS, CONSIDER USING A PRIORITY HEAP

    while req_count < number_of_requests:
        if arrivals[req_count] == curr_t:

            if load_loc[req2load[req_count]] in O:
                departures[req_count] = curr_t+1  # items is actually not introduced to the problem
                if args.log:
                    f = open(sim_log_file, "a")
                    f.write(f"Time: {arrivals[req_count]}: request #{req_count} arrive, load:{req2load[req_count]}, currently @ {load_loc[req2load[req_count]]} - lucky me, the load is on an output cell - departure time is {curr_t+1}\n")
                    f.close()
            else:
                requests_on_load[req2load[req_count]].append(req_count)
                if args.log:
                    f = open(sim_log_file, "a")
                    f.write(f"Time: {arrivals[req_count]}: request #{req_count} arrive, load:{req2load[req_count]}, currently @ {load_loc[req2load[req_count]]}\n")
                    f.close()
                if load_loc[req2load[req_count]] not in A:
                    if len(A) < args.max_balls_in_air:
                        new_target_loads_xy.add(load_loc[req2load[req_count]])  # just for the animation
                        A.append(load_loc[req2load[req_count]])
                    else:  # load is not already on the way
                        if req2load[req_count] not in load_queue:
                            load_queue.append(req2load[req_count])
                            block_max_queue[arrivals[req_count] // args.block_length] = max(
                                block_max_queue[arrivals[req_count] // args.block_length], len(load_queue))
            req_count += 1
        else:
            break

    if not A:
        moves.append([])  # idle cycle
        curr_t += 1
        idle_takt += 1
        if args.log:
            f = open(sim_log_file, "a")
            f.write(f"    Skipping an idle takt @ {curr_t}\n")
            f.close()
        if curr_t < args.simulation_length:
            continue
        else:
            print(f"Simulation end")
            print(f"Total lead time: {total_lead_time}")
            print(f"Mean lead time: {total_lead_time / req_count:.2f}  (including warmup period)")
            break

    sim_iter += 1
    curr_t += 1
    #We can do a bit better by advancing the simulation clock until next arrival when there is nothing to do

    old_A, old_E = copy.copy(A), copy.copy(E)

    A_arr, E_arr, O_arr = np.array(old_A), np.array(old_E), np.array(O)

    old_Aa, old_Ea = copy.copy(A_arr), copy.copy(E_arr)  # # just for debugging and QA
    A_arr, E_arr, one_step_move = SimpleHeuristic.OneStep(Lx, Ly, O_arr, A_arr, E_arr, args.bm)
    NumberOfMovements += len(one_step_move)

    # Just for debugging and QA
    if not SimpleHeuristic.test_step(Lx, Ly, old_Aa, old_Ea,A_arr, E_arr, one_step_move):
        print("Panic")
        exit(1)

    A = list(map(tuple, A_arr))
    E = list(map(tuple, E_arr))

    if args.log:
        f = open(sim_log_file, "a")
        f.write(
            f"Simulation Time: {curr_t}, \n")
        f.write(f"      Before operation: E = {E}, A = {A}, Q = {load_queue}\n")

    A = list(map(tuple, A_arr))
    E = list(map(tuple, E_arr))

    if args.log:
        f.write(f"      After operation: E = {E}, A = {A}, Q = {load_queue}\n")
        f.close()


    for (loc1, loc2) in one_step_move:
        load_loc[pbs[loc1[0], loc1[1]]] = loc2

    pbs = np.zeros((Lx, Ly), dtype=int)
    for l, loc in load_loc.items():
        pbs[loc[0], loc[1]] = l

    for loc in O:
        if pbs[loc[0], loc[1]] != 0:
            leaving_load = pbs[loc[0], loc[1]]
            for req in requests_on_load[leaving_load]:
                total_lead_time += (curr_t - arrivals[req])
                departures[req] = curr_t
                if args.log:
                    f = open(sim_log_file, "a")
                    f.write(
                        f"Time: {curr_t } request #{req} departs via ({loc[0]}, {loc[1]}) with load {leaving_load}\n")
                    #f.write(f"      E = {E}, A = {A}, Q = {load_queue}\n")
                    f.close()
            requests_on_load[leaving_load] = []
            if leaving_load in load_queue:
                if args.log:
                    f = open(sim_log_file, "a")
                    f.write(f"Time: {curr_t} load {leaving_load} accidentally arrived at an output cell and was removed from the queue\n")
                    f.close()
                load_queue.pop(load_queue.index(leaving_load))

    """  For QA only - check the integrity of the data """
    for l, loc in load_loc.items():
        if requests_on_load[l]: # there are still requests on the load
            if loc not in A and l not in load_queue:
                print(f"************ Panic: Time: {curr_t}   load {l} at location {loc} with requests {requests_on_load[l]} got lost")
                if args.log:
                    f = open(sim_log_file, "a")
                    f.write(f"************ Panic: load  Panic: load {l} at location {loc} with requests {requests_on_load[l]} got lost\n")
                    f.close()


    for l in new_target_loads_xy:
        one_step_move.insert(0, (l, l))  # to make the target load in the animation turn red
    moves.append(one_step_move)  # read moves in the current horizon

    if len(A) == 0 and len(load_queue) == 0 and curr_t > args.simulation_length and req_count >= number_of_requests:
        print(f"Simulation end")
        print(f"Total lead time: {total_lead_time}")
        print(f"Mean lead time: {total_lead_time / req_count:.2f}  (including warmup period)")
        break

    if curr_t > args.simulation_length*5:
        print(f"Something went wrong, exit at time {curr_t}")
        break

if args.export_animation:
    pickle.dump((Lx, Ly, O, E_orig, A_orig, moves),
                # remove empty moves  periods at the end
                open(
                    f"script_sim_heur_{'SBM' if args.bm else 'SLM'}_{Lx}_{Ly}_{args.escorts_num}_{args.request_rate}_{args.simulation_length}.p",
                    "wb"))

arrivals = np.array(arrivals, dtype=int)
lead_time_std = np.std(departures[args.warmup_arrivals:] - arrivals[args.warmup_arrivals:], ddof=1)
lead_time_mean = np.mean(departures[args.warmup_arrivals:] - arrivals[args.warmup_arrivals:])
half_ci = 1.96 * lead_time_std / (number_of_requests - args.warmup_arrivals) ** 0.5

# just for debugging and QA
if min(departures - arrivals)< 0 and args.log:
    print("Error: some arrivals occurs after departure")
    print(f"{np.where(departures - arrivals< 0)}")

    f = open(sim_log_file, "a")
    f.write("Error: some arrivals occurs after departure\n")
    f.write(f"{np.where(departures - arrivals < 0)}\n")
    f.close()


f = open(result_csv_file, 'a')
f.write(
    f"{time.ctime()},{'SBM' if args.bm else 'SLM'}, continue, {Lx}x{Ly}, {len(O)}, {len(E_orig)}, "
    f"{tuple_opl(O)},{args.request_rate},"
    f"{args.seed},{args.max_balls_in_air}, {args.simulation_length}, "
    f"{args.warmup_arrivals}, {curr_t}, {total_lead_time},{number_of_requests}, "
    f"{lead_time_mean:.3f}, {half_ci:.3f}, {NumberOfMovements},{cpu_time:.2f}, {idle_takt}, {max_gap:.4f}")

if args.print_blocks:
    for i in range(num_of_blocks):
        f.write(f",{block_max_queue[i]}")
f.write("\n")
f.close()