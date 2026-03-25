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
import OneStepHeuristic_v2

parser = argparse.ArgumentParser()


def solver_thread_count():
    if sys.platform.startswith("linux"):
        return 12
    if sys.platform == "darwin":
        return 8
    return 0


default_num_threads = solver_thread_count()

parser.add_argument("-x", "--Lx", type=int, help="Horizontal dimension of the PBS unit", required=True)
parser.add_argument("-y", "--Ly", type=int, help="Vertical dimension of the PBS unit", required=True)
parser.add_argument("-O", "--output_cells", nargs='+', type=int, help="List of output locations", required=True)
parser.add_argument("-e", "--escorts_range", help="Range of number of escorts  3, 3-10, or 1-30-3", default="5")
parser.add_argument("-r", "--reps_range", help="Replication range (seed range) 1, 1-10, or 1-30-3", default="1")
parser.add_argument(
    "-l",
    "--load_num",
    dest="load_num",
    type=int,
    help="Number of target loads (default 1)",
    default=1,
)
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
parser.add_argument("-t", "--time_limit", type=int,
                    help="Time limit for CPLEX/Gurobi in seconds (default 300 unless omitted while using only --work_limit)",
                    default=None)
parser.add_argument(
    "--num_threads",
    type=int,
    help=(
        "Number of solver threads for the static backends. "
        f"Default {default_num_threads}; 0 means use the solver default."
    ),
    default=default_num_threads,
)
parser.add_argument("--work_limit", type=float,
                    help="Work limit for Gurobi solves in work units (default none)", default=None)
parser.add_argument(
    "--mip_emphasis",
    choices=["balanced", "feasibility", "optimality", "bound"],
    default="balanced",
    help="Gurobi MIP emphasis for static MILP solves (default balanced)",
)

parser.add_argument("-a", "--export_animation", action="store_true",
                    help="Export animation files, one for each instance")

parser.add_argument("--lm", "--LM", dest="lm", action="store_true",
                    help="Run the static load-flow model in LM mode instead of the default BM mode")
parser.add_argument("--dp_file",
                    help="Dynamic programming file for upper bound and warm start (only relevant for single target load for now)", default="")
parser.add_argument("-k", "--k_prime", type=int,
                    help="k' value for the dp based heuristic (provide only if you provided dp_file)", default=0)
parser.add_argument("--lp", action="store_true",
                    help="Run lp relaxation work only with bm movement regime (default False)")
parser.add_argument("--gurobi", action="store_true",
                    help="Solve the static load-flow model with the Gurobi Python API (default)")
parser.add_argument("--opl", dest="opl", action="store_true",
                    help="Solve the static load-flow model with oplrun/CPLEX instead of the default Gurobi backend")
parser.add_argument("--cplex", dest="opl", action="store_true", help=argparse.SUPPRESS)


args = parser.parse_args()
if args.gurobi and args.opl:
    print("Panic: --gurobi and --opl cannot be combined")
    exit(1)
args.gurobi = not args.opl
result_csv_file = args.csv
file_export = "out.txt" if args.export_animation else ""
is_bm = not args.lm

dat_file = "escorts_flow3.dat"

""" consider getting this input from an external source """
alpha = args.alpha  # C_max weight
beta = args.beta  # flowtime weight
gamma = args.gamma  # movement weight
time_limit = args.time_limit  # Time limit for cplex run (seconds)
if time_limit is None and args.work_limit is None:
    time_limit = 300
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

if args.work_limit is not None and args.work_limit <= 0:
    print("Panic: --work_limit must be positive")
    exit(1)

if args.num_threads < 0:
    print("Panic: --num_threads must be nonnegative")
    exit(1)

if args.work_limit is not None and not args.gurobi:
    print("Panic: --work_limit is supported only with --gurobi")
    exit(1)

if args.dp_file:
    print(f"Loading DP file {args.dp_file}...", flush=True)
    S = pickle.load(open(args.dp_file, "rb"))
    print("Done", flush=True)

Locations = sorted(set(itertools.product(range(Lx), range(Ly))))
solver_threads = args.num_threads


def mip_emphasis_focus():
    return {
        "balanced": 0,
        "feasibility": 1,
        "optimality": 2,
        "bound": 3,
    }[args.mip_emphasis]


def mip_emphasis_for_csv():
    if args.gurobi and not args.lp:
        return args.mip_emphasis
    return "-"


def greedy_upper_bound(target_positions, escort_positions):
    max_steps = max(1, (Lx + Ly) * len(target_positions) * 20 // max(len(escort_positions), 1))
    makespan, _, _ = OneStepHeuristic_v2.SolveGreedy(
        Lx,
        Ly,
        set(O),
        set(target_positions),
        set(escort_positions),
        verbal=False,
        max_steps=max_steps,
        retrieval_mode=args.retrieval_mode,
    )
    return makespan

f = open(result_csv_file, 'a')
f.write(
    "\ndate, Moves, Model, MIP Emphasis, Retrieval Mode, Lx x Ly, #IOs, # Escorts, #Loads, IOs, Escorts, Target Loads, alpha, beta, gamma, k', seed, makespan, flowtime, #load movements, obj, LB, Wall Clock Time, Work\n")
f.close()

model_name = "LP-Gurobi" if args.gurobi and args.lp else (
    "ILP-Gurobi" if args.gurobi else ("LP" if args.lp else "ILP")
)

gurobi_solver = None
if args.gurobi:
    try:
        from load_flow_static_gurobi import LoadFlowStaticGurobiConfig, LoadFlowStaticGurobiSolver
    except ModuleNotFoundError as exc:
        if exc.name == "gurobipy":
            print("Panic: --gurobi requires the gurobipy package. Install it or run without --gurobi.")
            exit(1)
        raise

    gurobi_solver = LoadFlowStaticGurobiSolver(
        LoadFlowStaticGurobiConfig(
            Lx=Lx,
            Ly=Ly,
            output_cells=tuple(O),
            move_method="BM" if is_bm else "LM",
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            time_limit=time_limit,
            work_limit=args.work_limit,
            mip_focus=mip_emphasis_focus(),
            lp=args.lp,
            threads=solver_threads,
        )
    )

try:
    for escort_num in escorts_range:
        for rep in reps_range:
            """ Generate instance and solve the multi retrival problem """
            random.seed(rep)

            R, E = GeneretaeRandomInstance(rep, Locations, escort_num, load_num)
            if args.dp_file:
                if is_bm:
                    moves = PBS_DPHeuristic_bm.DOHueristicBM(S, R[0], E, Lx, Ly, O, args.k_prime, False)
                else:
                    moves = PBS_DPHeuristic_lm.DOHueristicLM(S, R[0], E, Lx, Ly, O, args.k_prime, False)
                T = len(moves)
            elif is_bm and args.retrieval_mode in ["continue", "leave"]:
                moves = []
                T = greedy_upper_bound(R, E)
            else:  # just guess T
                moves = []  # so it prints 0
                T = int((Lx + Ly + len(R) ** 0.7 - len(O) - escort_num ** 0.5) * args.T_factor)
                if is_bm:
                    T = int(T * .8)

            if not args.gurobi:
                f = open("pbs_load_flow.dat", "w")
                f.write('file_export = "%s";\n' % file_export)
                f.write(f'file_res = "{result_csv_file}";\n')
                f.write('time_limit = %d;\n' % time_limit)
                f.write('threads = %d;\n' % solver_threads)
                f.write(f'MoveMethod = {"BM" if is_bm else "LM"};\n')
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
                f"{time.ctime()},{'BM' if is_bm else 'LM'}, {model_name}, {mip_emphasis_for_csv()},"
                f"{args.retrieval_mode},{Lx}x{Ly}, {len(O)}, {len(E)}, {len(R)}, {tuple_opl(O)}, "
                f"{tuple_opl(E)}, {tuple_opl(R)}, {alpha},{beta},{gamma},{args.k_prime},{rep}")
            f.close()

            if args.gurobi:
                try:
                    result = gurobi_solver.solve(R, E, T)
                except Exception as exc:
                    print(f"Could not solve the model with Gurobi: {exc}")
                    f = open(result_csv_file, 'a')
                    f.write(",-,-,-,-,-,-,-\n")
                    f.close()
                else:
                    if result["status_name"] == "INTERRUPTED":
                        sys.exit(130)
                    if file_export != "" and result["has_solution"]:
                        pickle.dump(
                            (Lx, Ly, O, E, R, result["animation_moves"]),
                            open(
                                f"script_load_flow_{'BM' if is_bm else 'LM'}_{args.retrieval_mode}_{Lx}_{Ly}_"
                                f"{escort_num}_{load_num}_{rep}.p",
                                "wb",
                            ),
                        )
                    f = open(result_csv_file, 'a')
                    f.write(gurobi_solver.build_csv_suffix(result) + "\n")
                    f.close()
            else:
                try:
                    if args.lp:
                        subprocess.run(["oplrun", "pbs_load_flow_multi_lp.mod", "pbs_load_flow.dat"], check=True)
                    else:
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
                                    open(f"script_load_flow_{'BM' if is_bm else 'LM'}_{args.retrieval_mode}_{Lx}_{Ly}_"
                                         f"{escort_num}_{load_num}_{rep}.p", "wb"))
finally:
    if gurobi_solver is not None:
        gurobi_solver.close()
