"""
Test script file, 14/12/2023 Tal Raviv, talraviv@tau.ac.il
NOT WORKING FOR BM YET
"""
import copy
import pickle
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('filename')
parser.add_argument("-b", "--bm", action="store_true", help="Block movement regime (otherwise NBM, aka LM)")

args = parser.parse_args()

Lx, Ly, O, E, A, moves = pickle.load(open(args.filename, "rb"))

pbs_now = np.ones((Lx, Ly), dtype=np.int8)
target_load_grid = np.zeros((Lx, Ly), dtype=np.int8)


for e in E:
    if e[0] not in range(Lx) or e[1] not in range(Ly):
        print(f"Error: Initial location of escort {e} outside the {Lx}x{Ly} grid")
    pbs_now[e[0], e[1]] = 0

for a in A:
    if a[0] not in range(Lx) or a[1] not in range(Ly):
        print(f"Error: Initial location of target load {a} outside the {Lx}x{Ly} grid")
    target_load_grid[a[0], a[1]] = 1

number_of_movements = 0
flow_time = 0
remaining_loads = len(A)

print(f"Initial A: {A}")
print(f"Initial E: {E}")

for t in range(len(moves)):
    pbs_next = copy.copy(pbs_now)
    number_of_movements += len(moves[t])
    for (l1,l2) in moves[t]:

        if l1[0] not in range(Lx) or l1[1] not in range(Ly):
            print(f"Error: At time {t}  movement {l1} to {l2}: origin outside the {Lx}x{Ly} grid")
        else:
            if pbs_now[l1[0], l1[1]] == 0:
                print(f"Error: At time {t}  movement {l1} to {l2}: illegal movement from an empty location ")



        if l2 == (None, None):
            if l1 not in O:
                print(f"Error: At time {t}  movement {l1} to {l2}: a load leaves not via output cell")
            #elif pbs_now[l1[0], l1[1]] != 2:  not working now
            #    print(f"Error: At time {t}  movement {l1} to {l2}: a non target load leaves the PBS")

            print(f"At time {t} to {t+1} a load leaves the system via {l1}")
            remaining_loads -= 1
            flow_time += (t+1)
            pbs_next[l1[0], l1[1]] = 0
        else:
            if abs(l1[0] - l2[0]) + abs(l1[1] - l2[1]) > 1:
                print(f"Error: At time {t}  illegal move from {l1} to {l2}: an L1 distance of more than 1")

            if l2[0] not in range(Lx) or l2[1] not in range(Ly):
                print(f"Error: At time {t}  movement {l1} to {l2}: destination outside the {Lx}x{Ly} grid")
            else:
                if pbs_now[l2[0], l2[1]] != 0 and not args.bm:  # still allow perpendicular movements
                    print(f"Error: At time {t}  movement {l1} to {l2}: illegal movement onto a non-empty location at the current period - load type {pbs_now[l2[0], l2[1]]} in NBM regime")
                # if pbs_next[l2[0], l2[1]] != 0 and not args.bm:  # still allow perpendicular movements
                #     print(f"Error: At time {t}  movement {l1} to {l2}: illegal movement onto a non-empty location at the next period - load type  {pbs_now[l2[0], l2[1]]} a conflict")
                pbs_next[l2[0], l2[1]] += 1
                pbs_next[l1[0], l1[1]] -= 1

    over_crowded = list(zip(*np.where(pbs_now > 1)))
    under_crowded = list(zip(*np.where(pbs_now < 0)))
    if over_crowded:
        print(f"Time {t}: Over crowded cell {over_crowded}")
    if under_crowded:
        print(f"Time {t}: Under crowded cell {under_crowded}")
        print(f"{moves[t]}\n")

    pbs_now = copy.copy(pbs_next)


print("OK. Done")
print(f"Makespan: {t+1} (including additional step for the last load to leave the output cell)")
print(f"Flow time: {flow_time} (including additional step for each load to leave the output cell)")
print(f"Movements: {number_of_movements} (Including movements outside)")

if remaining_loads:
    print(f"Error: {remaining_loads} target loads left in the system")
