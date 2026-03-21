# -------------------------------------------------------------------------------
# Name:        EscortFlowStatic
# Purpose:     Run the escort flow ILP model, SBM (formerly LM)/ BM, multi-loads, stay/continue/leave
#              Use DP to create upper bound
#              For now works with one target load only
#
# Author:      Tal Raviv,  talraviv@tau.ac.il
#
# Created:     3/12/2023, first version adapted from the load flow formulation implementation, warm start for single load
#              5/12/2023 Added animation export file and argparse parameters
#              8/12/2023 Added block movement support (running pbs_escorts_bm.mod)
#              9/12/2023 added using DP heuristic to create upper bound for BM with a single load (without warm-start for now)
#              14/12/2023 adding lp relaxation with "--lp" switch
#              20/12/2025 (v3) adding cuts as described in the revised paper, removing support for external network, SLM, removing makespan objective (alpha)
#              17/3/2026 Rename to EscortFlowStatic.py and add support for the GreedyHeurisitc
#              19/3/2026 Add optional Gurobi API backend with --gurobi
#
# Copyright:   (c) Tal Raviv 2020, 2023, 2025, 2026
# Licence:     Free but please let me know that you are using it
# Depends on   PBSCom.py, PBS_DPHeuristic_lm.py, pbs_escorts3.mod, , pbs_escorts_bm.mod
#              Assumes oplrun is installed and on the path
# -------------------------------------------------------------------------------
import copy
import random
import sys
import itertools
import subprocess
import pickle
import time
import os
import socket

import numpy as np

from PBSCom import *
import PBS_DPHeuristic_lm
import PBS_DPHeuristic_bm
import OneStepHeuristic_v2
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-x", "--Lx", type=int, help="Horizontal dimension of the PBS unit", required=True)
parser.add_argument("-y", "--Ly", type=int, help="Vertical dimension of the PBS unit", required=True)
parser.add_argument("-O", "--output_cells", nargs='+', type=int, help="List of output locations", required=True)
parser.add_argument("-e", "--escorts_range", help="Range of number of escorts  3, 3-10, or 1-30-3", default="5")
parser.add_argument("-r", "--reps_range", help="Replication range (seed range) 1, 1-10, or 1-30-3", default="1")
parser.add_argument("-l", "--load_num", type=int, help="Number of target loads (default 1)", default=1)
parser.add_argument('-m', '--retrieval_mode', choices=['stay', 'leave', 'continue'],
                    help='Select a retrieval mode. Default is stay. Also note that for stay mode the number '
                         'of output cells must be greater or equal the number target loads', default='leave')
parser.add_argument("-f", "--csv", help="File name of the result csv file", default="res_escort_flow.csv")
parser.add_argument("--beta", type=float, help="Weight of the flowtime in the objective function (default 1.0)",
                    default=1.0)
parser.add_argument("--gamma", type=float, help="Weight of the movements in the objective function (default 0.01)",
                    default=0.01)
parser.add_argument("-T", "--T_factor", type=float,
                    help="multiplier of the planning horizon length when no heuristic is used to create an upper bound (default 2.0)",
                    default=1.6)
parser.add_argument("-t", "--time_limit", type=int,
                    help="Time limit for CPLEX/Gurobi in seconds (default 300 unless omitted while using only --work_limit)",
                    default=None)
parser.add_argument("--work_limit", type=float,
                    help="Work limit for Gurobi solves in work units (default none)", default=None)
parser.add_argument("--dp_file",
                    help="Dynamic programming file for upper bound and warm start (only relevant for single target load for now)",
                    default="")
parser.add_argument("-k", "--k_prime", type=int,
                    help="k' value for the dp based heuristic (provide only if you provided dp_file)", default=0)
parser.add_argument("-a", "--export_animation", action="store_true",
                    help="Export animation files, one for each instance")
parser.add_argument("--lp", action="store_true",
                    help="Run lp relaxation work only with bm movement regime (default False)")
parser.add_argument("--greedy", action="store_true",
                    help="Solve the static instance with OneStepHeuristic_v2.SolveGreedy instead of OPL")
parser.add_argument("--gurobi", action="store_true",
                    help="Solve the static model with the Gurobi Python API instead of oplrun/CPLEX")
parser.add_argument("--warmstart", action="store_true",
                    help="Use a heuristic MIP start with the Gurobi static backend (default False)")
parser.add_argument("--lazy", nargs="?", const=0, default=None, type=int,
                    help="Use the lazy-constraint Gurobi static backend; optional value is the number of initial time steps kept in the master problem, default 0; implies --gurobi")
parser.add_argument("--bnc", nargs="?", const=-1, default=None, type=int,
                    help="Use the branch-and-cut Gurobi static backend; optional value is the maximum number of user cuts added per separated node, default 2*T when omitted; implies --gurobi")

args = parser.parse_args()
if args.lazy is not None or args.bnc:
    args.gurobi = True

result_csv_file = args.csv
file_export = "out.txt" if args.export_animation else ""

dat_file = "escorts_flow3.dat"
network_file = "escort_network.dat"  # fixed file common for all number of escorts and load number""

""" consider getting this input from an external source """
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
dp_file = args.dp_file
k_prime = args.k_prime

if load_num == 1 and min(escorts_range) < k_prime and k_prime > 0:
    print(f"Warning: minimal number of escorts {min(escorts_range)} must be greater then k'={k_prime}")
    exit(1)

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

if k_prime > 0:
    if load_num > 1:
        print(
            "Panic: DP-k' heuristic works now only for retrieval of one load and only in 'stay' mode")
        exit(1)

    if args.retrieval_mode != "stay":
        print(
            "Warning: warm start for works with 'stay' mode only but we will use to DP-k' heuristic to create upper bound")

if len(set(O)) < len(O):
    print(f"Panic: All output cells location must be unique {sorted(O)}")
    exit(1)

if args.greedy and args.lp:
    print("Panic: --greedy cannot be combined with --lp")
    exit(1)

if args.greedy and args.retrieval_mode not in ["continue", "leave"]:
    print("Panic: --greedy currently supports only --retrieval_mode continue or leave")
    exit(1)

if args.work_limit is not None and args.work_limit <= 0:
    print("Panic: --work_limit must be positive")
    exit(1)

if args.work_limit is not None and not args.gurobi and args.lazy is None and not args.bnc:
    print("Panic: --work_limit is supported only with the Gurobi static backends")
    exit(1)

if args.lazy is not None and args.lazy < 0:
    print("Panic: --lazy must be a nonnegative integer")
    exit(1)

if args.bnc is not None and args.bnc < -1:
    print("Panic: --bnc must be a nonnegative integer")
    exit(1)

if args.lazy is not None and args.bnc:
    print("Panic: --lazy and --bnc cannot be combined")
    exit(1)

if args.lazy is not None and args.greedy:
    print("Panic: --lazy cannot be combined with --greedy")
    exit(1)

if args.lazy is not None and args.lp:
    print("Panic: --lazy is supported only for integer Gurobi solves, not with --lp")
    exit(1)

if args.bnc and args.greedy:
    print("Panic: --bnc cannot be combined with --greedy")
    exit(1)

if args.bnc and args.lp:
    print("Panic: --bnc is supported only for integer Gurobi solves, not with --lp")
    exit(1)

if args.warmstart and not args.gurobi:
    print("Panic: --warmstart is supported only together with --gurobi")
    exit(1)

if args.warmstart and args.greedy:
    print("Panic: --warmstart cannot be combined with --greedy")
    exit(1)

if args.warmstart and args.lp:
    print("Panic: --warmstart is supported only for integer Gurobi solves, not with --lp")
    exit(1)

if args.warmstart and args.retrieval_mode not in ["continue", "leave"]:
    print("Panic: --warmstart currently supports only --retrieval_mode continue or leave")
    exit(1)

if args.warmstart and dp_file:
    print("Panic: --warmstart is currently supported only on the greedy Gurobi upper-bound path, not with --dp_file")
    exit(1)

Locations = sorted(set(itertools.product(range(Lx), range(Ly))))
f = open(network_file, "w")

f.write('locations = {')
f.write(' '.join([f"<{x},{y}>" for x in range(Lx) for y in range(Ly)]))
f.write("};\n")

f.write('NA = [')
for x in range(Lx):
    for y in range(Ly):
        f.write("{")
        if x < Lx - 1:
            f.write(f"<{x + 1}, {y}> ")  # right
        if y < Ly - 1:
            f.write(f"<{x}, {y+1}> ")  # up
        if x > 0:
            f.write(f"<{x - 1}, {y}> ")  # left
        if y > 0:
            f.write(f"<{x}, {y-1}> ")  # down
        f.write(f"<{x}, {y}> ")
        f.write("}\n")
f.write("];\n")

f.write('NE = [')
for x in range(Lx):
    for y in range(Ly):
        f.write("{")
        for xx in range(Lx):
            f.write(f"<{xx}, {y}> ")
        for yy in range(Ly):
            if yy != y:
                f.write(f"<{x}, {yy}> ")
        f.write(f"<{x}, {y} > ")
        f.write("}\n")
f.write("];\n")

f.write('movesE = {')
for x in range(Lx):
    for y in range(Ly):
        for xx in range(Lx):
            f.write(f"<{x} {y} {xx} {y}> ")
        for yy in range(Ly):
            if yy != y:
                f.write(f"<{x} {y} {x} {yy}> ")
        f.write("\n")
f.write("};\n")





f.write('movesA = {')
for x in range(Lx):
    for y in range(Ly):
        if x < Lx - 1:
            f.write(f"<{x} {y} {x + 1} {y}> ")  # right
        if y < Ly - 1:
            f.write(f"<{x} {y} {x} {y + 1}> ")  # up
        if x > 0:
            f.write(f"<{x} {y} {x - 1} {y}> ")  # left
        if y > 0:
            f.write(f"<{x} {y} {x} {y - 1}> ")  # down
        f.write(f"<{x} {y} {x} {y}> ")
f.write("};\n")

f.write('CellCover = [')
for x in range(Lx):
    for y in range(Ly):
        f.write("{")
        for x1 in range(x+1):
            for x2 in range(x, Lx):
                if x1 != x2:
                    f.write(f"<{x1} {y} {x2} {y}> ")
                    f.write(f"<{x2} {y} {x1} {y}> ")
        for y1 in range(y+1):
            for y2 in range(y, Ly):
                if y1 != y2:
                    f.write(f"<{x} {y1} {x} {y2}> ")
                    f.write(f"<{x} {y2} {x} {y1}> ")
        f.write(f"<{x} {y} {x} {y}> ")
        f.write("}\n")
f.write("];\n")


# Conflict-pair generation is disabled; the static models rely on CellCover/avoid_conflicts.


f.write('MoveCover = [')
for x in range(Lx):
    for y in range(Ly):
        # right movement to (x+1,y)
        if x < Lx-1:
            f.write("{")
            for x1 in range(x+1):
                for x2 in range(x, Lx):
                    if x1 != x2:
                        f.write(f"<{x2} {y} {x1} {y}> ")
            f.write("}\n")

        # up movement to (x,y+1)
        if y < Ly-1:
            f.write("{")
            for y1 in range(y + 1):
                for y2 in range(y, Ly):
                    if y1 != y2:
                        f.write(f"<{x} {y2} {x} {y1}> ")
            f.write("}\n")

        # left movement to (x-1,y)
        if x> 0:
            f.write("{")
            for x1 in range(x+1):
                for x2 in range(x, Lx):
                    if x1 != x2:
                        f.write(f"<{x1} {y} {x2} {y}> ")
            f.write("}\n")

        # down movement to (x,y-1)
        if y > 0:
            f.write("{")
            for y1 in range(y+1):
                for y2 in range(y, Ly):
                    if y1 != y2:
                        f.write(f"<{x} {y1} {x} {y2}> ")
            f.write("}\n")
        f.write("{}\n")

f.write("];\n")
f.close()

header_line = (
    "Machine Name, Time Stamp, version, Moves, Model, Retrieval Mode, Lx x Ly, #IOs, # Escorts, #Loads, IOs, Escorts, "
    "Target Loads, beta, gamma, seed, T, Heuristic UB makespan, Heuristic UB movements , "
    "ILP makespan, ILP flowtime, #load movements, obj, LB, Wall Clock Time, Work, User Cut Time, Solver Status"
)
if not csv_file_contains_row(result_csv_file, header_line):
    append_csv_text(result_csv_file, header_line + "\n", ensure_record_start=True)

script_version = f"{os.path.basename(__file__)} ({time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(__file__)))})"
machine_name = socket.gethostname()

if dp_file:
    print(f"Loading DP file {dp_file}...", flush=True)
    S = pickle.load(open(dp_file, "rb"))
    print("Done", flush=True)

def model_name_for_instance(T):
    if args.greedy:
        return "Greedy"
    if args.gurobi:
        if args.bnc is not None:
            bnc_limit = 2 * T if args.bnc == -1 else args.bnc
            return f"ILP-Gurobi-BnC-{bnc_limit}"
        if args.lazy is not None:
            return "ILP-Gurobi-Lazy"
        return "LP-Gurobi" if args.lp else "ILP-Gurobi"
    return "LP" if args.lp else "ILP"

gurobi_solver = None
if args.gurobi and not args.greedy:
    try:
        if args.bnc:
            from escort_flow_static_bnc import StaticGurobiConfig, BnCStaticEscortFlowGurobiSolver
        elif args.lazy is not None:
            from escort_flow_static_lazy import StaticGurobiConfig, LazyStaticEscortFlowGurobiSolver
        else:
            from escort_flow_static_gurobi import StaticGurobiConfig, StaticEscortFlowGurobiSolver
    except ModuleNotFoundError as exc:
        if exc.name == "gurobipy":
            print("Panic: --gurobi requires the gurobipy package. Install it or run without --gurobi.")
            exit(1)
        raise

    if args.bnc is not None:
        solver_class = BnCStaticEscortFlowGurobiSolver
    elif args.lazy is not None:
        solver_class = LazyStaticEscortFlowGurobiSolver
    else:
        solver_class = StaticEscortFlowGurobiSolver
    config_kwargs = dict(
        Lx=Lx,
        Ly=Ly,
        output_cells=tuple(O),
        retrieval_mode=args.retrieval_mode,
        beta=beta,
        gamma=gamma,
        time_limit=time_limit,
        work_limit=args.work_limit,
        lp=args.lp,
    )
    if args.lazy is not None:
        config_kwargs["lazy_master_steps"] = args.lazy
    if args.bnc is not None:
        config_kwargs["bnc_cut_limit"] = None if args.bnc == -1 else args.bnc
    gurobi_solver = solver_class(StaticGurobiConfig(**config_kwargs))


def greedy_upper_bound_and_warmstart(target_positions, escort_positions, build_warmstart=False):
    max_steps = max(1, (Lx + Ly) * len(target_positions) * 20 // max(len(escort_positions), 1))
    need_trace = build_warmstart and args.retrieval_mode != "leave"
    greedy_result = OneStepHeuristic_v2.SolveGreedy(
        Lx,
        Ly,
        set(O),
        set(target_positions),
        set(escort_positions),
        verbal=False,
        max_steps=max_steps,
        return_trace=need_trace,
        retrieval_mode=args.retrieval_mode,
    )
    if need_trace:
        makespan, _, movements, _, escort_moves, target_moves = greedy_result
    else:
        makespan, _, movements = greedy_result
        escort_moves = None
        target_moves = None
    warmstart = None
    if build_warmstart and gurobi_solver is not None and not args.lp:
        if args.retrieval_mode == "leave":
            makespan, warmstart = gurobi_solver.build_feasible_leave_warmstart(
                target_positions,
                escort_positions,
            )
        else:
            warmstart = gurobi_solver.build_warmstart_from_trace(
                target_positions,
                escort_positions,
                makespan,
                target_moves,
                escort_moves,
            )
    return makespan, movements, warmstart

try:
    for escort_num in escorts_range:
        for rep in reps_range:
            random.seed(rep)
            R, E = GeneretaeRandomInstance(rep, Locations, escort_num, load_num)
            heuristic_ub_makespan = 0
            heuristic_ub_movements = 0
            warmstart = None

            if dp_file:
                moves = PBS_DPHeuristic_bm.DOHueristicBM(S, R[0], E, Lx, Ly, O, k_prime, False)
                T = len(moves)
                heuristic_ub_makespan = len(moves)
                heuristic_ub_movements = sum(len(step) for step in moves)
            elif args.gurobi and args.retrieval_mode in ["continue", "leave"]:
                moves = []
                heuristic_ub_makespan, heuristic_ub_movements, warmstart = greedy_upper_bound_and_warmstart(
                    R,
                    E,
                    build_warmstart=args.warmstart,
                )
                T = heuristic_ub_makespan
            else:  # just guess T
                moves = []  # so it prints 0
                T = int((Lx + Ly + len(R) ** 0.7 - len(O) - escort_num ** 0.5) * args.T_factor)

            model_name = model_name_for_instance(T)

            if not args.greedy and not args.gurobi:
                f = open(dat_file, "w")
                f.write('file_export = "%s";\n' % file_export)
                f.write(f'file_res = "{result_csv_file}";\n')
                f.write('time_limit = %d;\n' % time_limit)
                f.write('beta=%f;\n' % beta)
                f.write('gamma=%f;\n' % gamma)
                f.write('Lx=%d;\n' % Lx)
                f.write('Ly=%d;\n' % Ly)
                f.write(f'T={T};\n')
                f.write('E=%s;\n' % tuple_opl(E))
                f.write('A=%s;\n' % tuple_opl(R))
                f.write('O=%s;\n' % tuple_opl(O))
                f.write(f'retrieval_mode = "{args.retrieval_mode}";\n')
                f.close()

            append_csv_text(
                result_csv_file,
                f"{machine_name}, {time.ctime()}, {script_version}, BM, {model_name},"
                f"{args.retrieval_mode}, {Lx}x{Ly}, {len(O)}, {len(E)}, {len(R)}, {tuple_opl(O)}, "
                f"{tuple_opl(E)}, {tuple_opl(R)},{beta},{gamma},{rep}, {T}, {heuristic_ub_makespan}, "
                f"{heuristic_ub_movements}",
                ensure_record_start=True,
            )

            if args.greedy:
                max_steps = max(T, (Lx + Ly) * len(R) * 20 // max(len(E), 1))
                greedy_start_time = time.perf_counter()
                makespan, flowtime, movements, greedy_moves = OneStepHeuristic_v2.SolveGreedy(
                    Lx, Ly, set(O), set(R), set(E), verbal=False, max_steps=max_steps, return_moves=True,
                    retrieval_mode=args.retrieval_mode
                )
                cpu_time = time.perf_counter() - greedy_start_time

                if file_export != "":
                    pickle.dump((Lx, Ly, O, E, R, greedy_moves),
                                open(f"script_BM_{args.retrieval_mode}_{Lx}_{Ly}_{escort_num}_{load_num}_{rep}.p", "wb"))

                append_csv_text(
                    result_csv_file,
                    f",{makespan}, {flowtime}, {movements}, {beta * flowtime + gamma * movements}, , {cpu_time:.4f}, , 0.0000, Greedy"
                )
            elif args.gurobi:
                try:
                    result = gurobi_solver.solve(R, E, T, warmstart=warmstart)
                except Exception as exc:
                    print(f"Could not solve the model with Gurobi: {exc}")
                    append_csv_text(result_csv_file, ",-,-,-,-,-,-,-,0.0000, ERROR")
                else:
                    if result["status_name"] == "INTERRUPTED":
                        sys.exit(130)
                    if file_export != "" and result["has_solution"]:
                        pickle.dump(
                            (Lx, Ly, O, E, R, result["animation_moves"]),
                            open(f"script_BM_{args.retrieval_mode}_{Lx}_{Ly}_{escort_num}_{load_num}_{rep}.p", "wb")
                        )
                    append_csv_text(result_csv_file, gurobi_solver.build_csv_suffix(result))
            else:
                try:
                    if args.lp:
                        subprocess.run(["oplrun", "pbs_escorts_bm_lp_v3.mod", dat_file, network_file], check=True)
                    else:
                        subprocess.run(["oplrun", "pbs_escorts_bm_v3.mod", dat_file, network_file], check=True)

                except:
                    print("Could not solve the model")
                else:
                    # Create script file
                    if file_export != "":
                        f = open(file_export)
                        s = f.readlines()
                        moves = eval(s[-1])  # read moves in the current horizon
                        f.close()
                        pickle.dump((Lx, Ly, O, E, R, moves),
                                    open(f"script_BM_{args.retrieval_mode}_{Lx}_{Ly}_{escort_num}_{load_num}_{rep}.p", "wb"))
            append_csv_text(result_csv_file, "\n")
finally:
    if gurobi_solver is not None:
        gurobi_solver.close()
