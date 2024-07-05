# -------------------------------------------------------------------------------
# Name:        load_flow_multy
# Purpose:     Run the experiment with the full model for multi retrival and under LM and BM
#
# Author:      Tal Raviv,  talraviv@tau.ac.il
#
# Created:     23/08/2020, updated 18/9/2020, 15/3/2023, 10/12/2023, 18/2/2024
# Copyright:   (c) Tal Raviv 2020, 2023
# Licence:     Free but please let me know that you are using it
# -------------------------------------------------------------------------------

# import os
import random
import sys
import itertools
import subprocess
import pickle
import time
import copy
from PBSCom import *
import numpy as np
import argparse
import PBS_DPHeuristic_lm
import PBS_DPHeuristic_bm

parser = argparse.ArgumentParser()

parser.add_argument("-x", "--Lx", type=int, help="Horizontal dimension of the PBS unit", required=True)
parser.add_argument("-y", "--Ly", type=int, help="Vertical dimension of the PBS unit", required=True)
parser.add_argument("-O", "--output_cells", nargs='+', type=int, help="List of output locations", required=True)
parser.add_argument("-e", "--escorts_range", help="Range of number of escorts  3, 3-10, or 1-30-3", default="5")
parser.add_argument("-r", "--reps_range", help="Replication range (seed range) 1, 1-10, or 1-30-3", default="1")
parser.add_argument("-l", "--load_num", type=int, help="Number of target loads (default 1)", default=1)
parser.add_argument('-m', '--retrieval_mode', choices=['stay', 'leave', 'continue'],
                    help='Select a retrieval mode. Default is stay. Also note that for stay mode the number '
                         'of output cells must be greater or equal the number target loads  (ONLY LEAVE IS SUPPORTED FOR NOW)',
                    default='leave')
parser.add_argument("-f","--csv", help="File name of the result csv file", default="res_load_flow.csv")
parser.add_argument("--alpha", type=float, help="Weight of the makespan in the objective function (default 0.0)",
                    default=0)
parser.add_argument("--beta", type=float, help="Weight of the flowtime in the objective function (default 1.0)",
                    default=1.0)
parser.add_argument("--gamma", type=float, help="Weight of the movements in the objective function (default 0.01)",
                    default=0.01)
parser.add_argument("-T", "--T_factor", type=float,
                    help="multiplier of the planning horizon length when no heuristic is used to create an upper bound (default 2.0)",
                    default=2.0)
parser.add_argument("-t", "--time_limit", type=int, help="Time limit for CPLEX (default 300)", default=300)

parser.add_argument("-a", "--export_animation", action="store_true",
                    help="Export animation files, one for each instance")

parser.add_argument("-b", "--bm", action="store_true", help="Block movement regime (otherwise NBM, aka LM)")
parser.add_argument("--dp_file",
                    help="Dynamic programming file for upper bound and warm start (only relevant for single target load for now)", default="")
parser.add_argument("-k", "--k_prime", type=int,
                    help="k' value for the dp based heuristic (provide only if you provided dp_file)", default=0)



args = parser.parse_args()
result_csv_file = args.csv
file_export = "out.txt" if args.export_animation else ""

dat_file = "escorts_flow3.dat"

""" consider getting this input from an external source """
alpha = args.alpha  # C_max weight
beta = args.beta  # flowtime weight
gamma = args.gamma  # movement weight
time_limit = args.time_limit  # Time limit for cplex run (seconds)
Lx = args.Lx
Ly = args.Ly





reps_range = str2range(args.reps_range)
load_num = args.load_num
escorts_range = str2range(args.escorts_range)

ez = [int(a) for a in args.output_cells]
if len(ez) % 2:
    print(f"Warning: output cell coordinate list must be of even length - {ez}")
    exit(1)

O = []
for i in range(len(ez) // 2):
    O.append((ez[i * 2], ez[i * 2 + 1]))

if len(O) < load_num and args.retrieval_mode == "stay":
    print(
        f"Panic: In 'stay' settings number of output cells ({len(O)}) can not be smaller than the number of target loads ({load_num})")
    exit(1)


if args.k_prime > 0:
    if load_num > 1:
        print(
            "Panic: DP-k' heuristic works now only for retrieval of one load and only in 'stay' mode")
        exit(1)

if args.dp_file:
    print(f"Loading DP file {args.dp_file}...", flush=True)
    S = pickle.load(open(args.dp_file, "rb"))
    print("Done", flush=True)

Locations = sorted(set(itertools.product(range(Lx), range(Ly))))

f = open(result_csv_file, 'a')
f.write(
    "\ndate, Lx x Ly, #IOs, # Escorts, #Loads, IOs, Escorts, Target Loads, alpha, beta, gamma, k', seed, makespan, flowtime, #load movements, obj, LB, CPU timee\n")
f.close()

for escort_num in escorts_range:
    for rep in reps_range:
        """ Generate instance and solve the multi retrival problem """
        random.seed(rep)

        R, E = GeneretaeRandomInstance(rep, Locations, escort_num, load_num)
        if args.dp_file:
            if args.bm:
                moves = PBS_DPHeuristic_bm.DOHueristicBM(S, R[0], E, Lx, Ly, O, args.k_prime, False)
            else:
                moves = PBS_DPHeuristic_lm.DOHueristicLM(S, R[0], E, Lx, Ly, O, args.k_prime, False)
            T = len(moves)
        else:  # just guess T
            moves = []  # so it prints 0
            T = int((Lx + Ly + len(R) ** 0.7 - len(O) - escort_num ** 0.5) * args.T_factor)
            if args.bm:
                T = int(T * .8)


        f = open("pbs_load_flow.dat", "w")

        f.write('file_export = "%s";\n' % file_export)
        f.write(f'file_res = "{result_csv_file}";\n')
        f.write('time_limit = %d;\n' % time_limit)
        f.write(f'MoveMethod = {"BM" if args.bm else "LM"};\n')
        f.write('alpha=%f;\n' % alpha)
        f.write('beta=%f;\n' % beta)
        f.write('gamma=%f;\n' % gamma)
        f.write('Lx=%d;\n' % Lx)
        f.write('Ly=%d;\n' % Ly)
        f.write('T=%d;\n' % T)
        f.write('E=%s;\n' % tuple_opl(E))
        f.write('R=%s;\n' % tuple_opl(R))
        f.write('O=%s;\n' % tuple_opl(O))
        f.close()

        f = open(result_csv_file, 'a')
        f.write(
            f"{time.ctime()},{Lx}x{Ly}, {len(O)}, {len(E)}, {len(R)}, {tuple_opl(O)}, {tuple_opl(E)}, {tuple_opl(R)}, {alpha},{beta},{gamma},{args.k_prime},{rep}")
        f.close()

        try:
            subprocess.run(["oplrun", "pbs_load_flow_multi.mod", "pbs_load_flow.dat"], check=True)
            f = open(result_csv_file, 'a')
            f.write("\n")
            f.close()
        except:
            print("Could not solve the model")
            f = open(result_csv_file, 'a')
            f.write(f", No solution, -,-,-,-,{time_limit + 0.01}\n")
            f.close()
        else:
            # Create script file
            if file_export != "":
                f = open(file_export)
                s = f.readlines()
                moves = eval(s[-1])  # read moves in the current horizon
                f.close()
                pickle.dump((Lx, Ly, O, E, R, moves),
                            open(f"script_load_flow_{'BM' if args.bm else 'LM'}_{args.retrieval_mode}_{Lx}_{Ly}_"
                                 f"{escort_num}_{load_num}_{rep}.p", "wb"))
