#-------------------------------------------------------------------------------
# Name:        PBS_DPHeuristic_bm.py
# Purpose:     load data, generate instance and apply the DP heurstic  for the LM case
#              PBSk|sim,lm,mIO| *
#              The objective function * determined by the DP table
#              Block movement version
#
# Author:      Tal Raviv talraviv@au.ac.il
#
# Created:     25-06-2020
# Copyright:   (c) TAL 2020
# Licence:     Free to use but please contact me
#-------------------------------------------------------------------------------

import argparse
import pickle
from PBSCom import *
import itertools
import random
import time
import copy
from OneStepHeuristic_v2 import SolveGreedy

scripts_path = "./"
reps = 20  # number of replications


def print_progress(current, total):
    bar_width = 30
    filled = bar_width if total == 0 else int(bar_width * current / total)
    bar = "#" * filled + "-" * (bar_width - filled)
    print(f"\rProgress [{bar}] {current}/{total}", end="", flush=True)
    if current == total:
        print("")



# assume two 2-tuples as input return a directional vector
def direction(a,b):

    if a[0]<b[0]:
        return ((1,0))
    elif a[0]>b[0]:
        return ((-1,0))
    elif a[1]<b[1]:
        return ((0,1))
    elif a[1]>b[1]:
        return ((0,-1))
    else:
        return ((0,0))


def movement_cells_and_load_moves(origin, destination):
    cells = {origin}
    load_moves = []
    curr_loc = origin
    d = direction(origin, destination)

    while curr_loc != destination:
        next_loc = (curr_loc[0] + d[0], curr_loc[1] + d[1])
        cells.add(next_loc)
        load_moves.append((next_loc, curr_loc))
        curr_loc = next_loc

    return cells, load_moves


def rank_candidate_moves(S, I, E_local, Lx, Ly, k):
    candidates = []
    for Et in itertools.combinations(E_local, k):
        escorts = sorted(list(Et))
        state_idx = listTuple2Int([I] + escorts, Lx, Ly)
        score, successor = S[state_idx]
        candidates.append((score, escorts, successor))

    candidates.sort(key=lambda item: (item[0], item[1]))
    return candidates


def normalize_loads(loads):
    if isinstance(loads, tuple) and len(loads) == 2 and all(isinstance(v, int) for v in loads):
        return [loads]
    return list(loads)


def distance_to_closest_output(load, terminals):
    return min(abs(load[0] - out[0]) + abs(load[1] - out[1]) for out in terminals)


def apply_non_conflicting_candidate(
    candidate_escorts,
    candidate_successor,
    occupied_escorts,
    used_cells,
    used_escorts,
    moved_load_cells,
    period_load_cells,
    period_load_ranks,
    current_priority_rank,
    protected_load,
    terminals,
    Lx,
    Ly,
):
    next_state = int2ListTuple(candidate_successor, Lx, Ly)
    next_escorts = next_state[1:]
    curr_move = []

    for origin, destination in zip(candidate_escorts, next_escorts):
        if origin == destination or origin in used_escorts or origin not in occupied_escorts:
            continue

        movement_cells, load_moves = movement_cells_and_load_moves(origin, destination)
        candidate_moved_loads = (movement_cells - {origin}) & period_load_cells
        violates_priority = False
        for old_load, new_load in load_moves:
            if old_load != protected_load:
                if current_priority_rank == 0:
                    continue
                old_rank = period_load_ranks.get(old_load)
                if old_rank is None or old_rank <= current_priority_rank:
                    continue
            if distance_to_closest_output(new_load, terminals) > distance_to_closest_output(old_load, terminals):
                violates_priority = True
                break

        if movement_cells & used_cells:
            continue
        if candidate_moved_loads & moved_load_cells:
            continue
        if violates_priority:
            continue
        if (movement_cells - {origin}) & occupied_escorts:
            continue

        used_escorts.add(origin)
        used_cells.update(movement_cells)
        moved_load_cells.update(candidate_moved_loads)
        curr_move.extend(load_moves)
        occupied_escorts.remove(origin)
        occupied_escorts.add(destination)

    return curr_move, occupied_escorts


# Prints the optimal path from state S
#  This is improted by  the path_print_dp.py program
def DOHueristicBM(S, I, E, Lx , Ly, Terminals,k, chat=True ):
    def print_state_snapshot(step, loads, escorts):
        print("Period ", step, " - ", loads, escorts)

        for y in range(Ly - 1, -1, -1):
            for x in range(Lx):
                if (x, y) in loads:
                    print("$", end="")
                elif (x, y) in escorts:
                    print(".", end="")
                elif (x, y) in Terminals:
                    print("T", end="")
                else:
                    print("*", end="")
            print("")
        print("")

    t = 0

    active_loads = list(enumerate(normalize_loads(I)))
    E_local = copy.copy(E)
    max_chat_steps = 4 * (Lx + Ly)
    seen_states = set()

    moves = []
    flowtime = 0

    while True:
        active_loads = [(load_id, load) for load_id, load in active_loads if load not in Terminals]
        if not active_loads:
            break

        current_state = (
            tuple(load for _, load in active_loads),
            tuple(sorted(E_local)),
        )
        if current_state in seen_states:
            if chat:
                print("Cycle detected. Switching to OneStepHeuristic_v2 fallback.")

            remaining_loads = [load for _, load in active_loads]
            greedy_max_steps = max(1, (Lx + Ly) * max(1, len(remaining_loads)) * 20 // max(1, len(E_local)))
            _, greedy_flowtime, _, greedy_moves = SolveGreedy(
                Lx,
                Ly,
                set(Terminals),
                set(remaining_loads),
                set(E_local),
                verbal=False,
                max_steps=greedy_max_steps,
                acyclic=False,
                return_moves=True,
                retrieval_mode="continue",
            )
            moves.extend(greedy_moves)
            flowtime += greedy_flowtime
            _, E_local = play_moves(remaining_loads, sorted(E_local), greedy_moves)
            active_loads = []
            break
        seen_states.add(current_state)

        flowtime += len(active_loads)

        if chat:
            print_state_snapshot(t, [load for _, load in active_loads], E_local)


        # ez = ShortestPath(I,E_local,Terminals)  # check if we can reuse this in the future
        # if ez:
        #     print("*****  Finish with shortest path shortcut ******")
        #     return moves+[ [a] for a in ez]

        period_loads = copy.copy(active_loads)
        period_load_positions = [load for _, load in period_loads]
        period_escorts = sorted(copy.copy(E_local))
        period_load_cells = set(period_load_positions)
        period_load_ranks = {load: idx for idx, load in enumerate(period_load_positions)}
        protected_load = period_load_positions[0]
        candidate_moves_by_load = [
            rank_candidate_moves(S, load, period_escorts, Lx, Ly, k)
            for load in period_load_positions
        ]
        max_rank = max(len(candidate_moves) for candidate_moves in candidate_moves_by_load)
        used_cells = set()
        used_escorts = set()
        moved_load_cells = set()
        occupied_escorts = set(period_escorts)
        curr_move = []
        t += 1

        for rank_idx in range(max_rank):
            if len(used_escorts) == len(period_escorts):
                break
            for load_rank, candidate_moves in enumerate(candidate_moves_by_load):
                if len(used_escorts) == len(period_escorts):
                    break
                if rank_idx >= len(candidate_moves):
                    continue

                _, candidate_escorts, candidate_successor = candidate_moves[rank_idx]
                if candidate_successor == "Sink":
                    continue

                candidate_move, occupied_escorts = apply_non_conflicting_candidate(
                    candidate_escorts,
                    candidate_successor,
                    occupied_escorts,
                    used_cells,
                    used_escorts,
                    moved_load_cells,
                    period_load_cells,
                    period_load_ranks,
                    load_rank,
                    protected_load,
                    Terminals,
                    Lx,
                    Ly,
                )
                curr_move.extend(candidate_move)

        moves.append(curr_move)
        new_load_positions, E_local = play_moves(period_load_positions, period_escorts, [curr_move])
        active_loads = [
            (period_loads[idx][0], new_load_positions[idx])
            for idx in range(len(period_loads))
        ]
        E_local.sort()

        if chat and len(moves) >= max_chat_steps:
            print(f"Stopping chat trace after {max_chat_steps} steps.")
            break

        if not curr_move:
            break

    if chat:
        final_loads = [load for _, load in active_loads if load not in Terminals]
        print_state_snapshot(t, final_loads, E_local)

    return moves, flowtime


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply the block-movement DP heuristic to generated PBS instances.",
    )
    parser.add_argument(
        "-x",
        "--Lx",
        type=int,
        required=True,
        help="Horizontal dimension of the PBS unit.",
    )
    parser.add_argument(
        "-y",
        "--Ly",
        type=int,
        required=True,
        help="Vertical dimension of the PBS unit.",
    )
    parser.add_argument(
        "-O",
        "--output_cells",
        nargs="+",
        type=int,
        required=True,
        help="Flat list of output-cell coordinates, e.g. '-O 2 0 4 0'.",
    )
    parser.add_argument(
        "-e",
        "--num_escorts",
        required=True,
        type=int,
        help="Number of actual escorts in the generated instance.",
    )
    parser.add_argument(
        "--dp_file",
        required=True,
        help="Pickle file containing the precomputed DP table.",
    )
    parser.add_argument(
        "-k",
        "--k_prime",
        type=int,
        required=True,
        help="k' value used by the DP heuristic.",
    )
    parser.add_argument(
        "-r",
        "--seed_range",
        help="Range of random seeds to test, parsed by str2range.",
        default="1",
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=0.01,
        help="Weight of the movements in the objective value (default 0.01).",
    )
    parser.add_argument(
        "-l",
        "--load_num",
        type=int,
        default=1,
        help="Number of target loads in the generated instance (default 1).",
    )
    parser.add_argument(
        "-c",
        "--chat",
        action="store_true",
        help="Print the period-by-period heuristic trace.",
    )
    parser.add_argument(
        "-a",
        "--export_animation",
        action="store_true",
        help="Export the solution to a pickle file for animation.",
    )
    return parser.parse_args()



if __name__ == '__main__':
    args = parse_args()

    if len(args.output_cells) % 2 != 0:
        raise ValueError("--output_cells must contain an even number of integers")

    Lx = args.Lx
    Ly = args.Ly
    Terminals = [
        (args.output_cells[i], args.output_cells[i + 1])
        for i in range(0, len(args.output_cells), 2)
    ]
    k = args.k_prime
    name = args.dp_file

    if k - 1 == args.num_escorts:
        print(f"Panic: cannot solve using k={k} instances with --num_escorts={args.num_escorts}")
        exit(1)

    rep_range = list(str2range(args.seed_range))

    Locations =  sorted(set(itertools.product(range(Lx), range(Ly))))
    print("Loading DP Table...")
    S = pickle.load(open(args.dp_file, "rb"))
    print("Finished loading DP Table.")
    f = open("res_dph_bm.csv", "a")
    f.write("Date, Name, Model, Lx x Ly, seed, IOs, #escorts, k', Load, Escorts, makespan, flowtime, #moves, objective, cpu time\n")
    f.close()
    K = args.num_escorts
    total_makespan = 0.0
    total_flowtime = 0.0
    total_num_moves = 0.0
    num_instances = 0
    for idx, i in enumerate(rep_range, start=1):
        A, E = GeneretaeRandomInstance(i, Locations, K, args.load_num)

        startTime = time.time()
        moves, flowtime = DOHueristicBM(S, A, E, Lx, Ly, Terminals, k, args.chat)
        makespan = len(moves)
        num_moves = sum(len(a) for a in moves)
        objective_value = makespan + args.gamma * num_moves
        total_makespan += makespan
        total_flowtime += flowtime
        total_num_moves += num_moves
        num_instances += 1

        f = open("res_dph_bm.csv", "a")
        f.write(
            f"{time.ctime()}, {name}, BM, {Lx}x{Ly}, {i}, "
            f"{str(Terminals).replace(',', '')},{K},  {k}, "
            f"{tuple_opl(A)}, {tuple_opl(E)}, {makespan}, {flowtime}, "
            f"{num_moves}, {objective_value:1.3f}, {time.time() - startTime:1.3f}\n"
        )
        f.close()

        if args.export_animation:
            f = open(f"{scripts_path}/DPH_BM_{Lx}_{Ly}_{len(E)}_{k}_{i}.p", "wb")
            pickle.dump((Lx, Ly, Terminals, E, A, moves), f)
            f.close()

        if not args.chat:
            print_progress(idx, len(rep_range))

    if num_instances:
        print(f"Average makespan: {total_makespan / num_instances:1.3f}")
        print(f"Average flowtime: {total_flowtime / num_instances:1.3f}")
        print(f"Average #moves: {total_num_moves / num_instances:1.3f}")
