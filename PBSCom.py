# -------------------------------------------------------------------------------
# Name:        PBSCom
# Purpose:     Common funciton used by both block and laod movement
#
# Author:      Tal Raviv
#
# Created:     31/08/2019
# Updated      4/12/2023  add play_moves
# Copyright:   (c) Tal Raviv 2019
# Licence:     Free but please contact me before use talraviv@tauex.tau.ac.il
# -------------------------------------------------------------------------------
import copy
import math
import itertools
import random


# The general function in this file are used also by other programs
# So check dependency issue before modifying or removing these functions
def distance(a, b=(0, 0)):
    return sum([abs(a[i] - b[i]) for i in range(2)])


def sgn(a):
    if a > 0:
        return 1
    elif a < 0:
        return -1
    else:
        return 0


def ShortestPath(I, E, Terminals):
    # first find closeset IO
    ez = min([(distance(I, b), b) for b in Terminals])
    # L1dist = ez[0]
    IO = ez[1]

    if I == IO:
        return []  # just in case

    OpenNodes = [I]
    Nodes = dict()

    while OpenNodes:
        curr_loc = OpenNodes.pop()
        dir_x, dir_y = sgn(IO[0] - curr_loc[0]), sgn(IO[1] - curr_loc[1])

        if dir_x != 0:
            next_loc = (curr_loc[0] + dir_x, curr_loc[1])
            if next_loc in E and not next_loc in Nodes:
                Nodes[next_loc] = curr_loc
                OpenNodes.append(next_loc)
                if next_loc == IO:
                    break

        if dir_y != 0:
            next_loc = (curr_loc[0], curr_loc[1] + dir_y)
            if next_loc in E and not next_loc in Nodes:
                Nodes[next_loc] = curr_loc
                OpenNodes.append(next_loc)
                if next_loc == IO:
                    break

    if next_loc == IO and IO in Nodes:
        ret_path = []
        while True:
            prev_loc = Nodes[next_loc]
            ret_path.insert(0, (prev_loc, next_loc))
            if prev_loc == I:
                break
            next_loc = prev_loc
        return ret_path
    else:
        return []


def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


# print schematic of the satte
def print_state(s, Lx, Ly):
    ez = pz.int2ListTuple(s, Lx, Ly)
    k = len(ez) - 1
    Ix, Iy = ez[0]
    E = set(ez[1:])

    for y in range(Ly - 1, -1, -1):
        for x in range(Lx):

            if (Ix, Iy) == (x, y):
                print("$", end="")
            elif (x, y) in E:
                print(".", end="")
            else:
                print("*", end="")
        print("")


# Convert a list of tuples in the range 0..Lx-1, 0..ly-1 into an integer number that is used as the key of the state
# L is the list of tuple where the first is the location of the item and the reset are
# the location of the escrots
# The structure of the compact int is number_of_escorts+ sum(i in 0..1) (L[i][0]+Lx*L[i][1])*B^(i+1)

# The number of escorts is greater than the number of the loads the state is encoded by the locations of
# the loads and as a negative number
# MAYBE WE CAN REMOVE THIS FEATURE BECAUSE WE ARE ONLY INTERESTED IN SMALL NUMBER OF ESCORTS


def listTuple2Int(L, Lx, Ly):
    NumOfCells = Lx * Ly;
    b = NumOfCells
    n = len(L) - 1

    for i in L:
        n += (i[0] * Ly + i[1]) * b
        b *= NumOfCells

    return n


# Inverese of  listTuple2Int.
# n - the state encoding
# Lx, Ly - the dimension of the waerhouse
# The function returns a list of tuples including the location of the item and the locations of all the escorts

def int2ListTuple(n, Lx, Ly):
    l = []
    NumOfCells = Lx * Ly;
    n, k = divmod(n, NumOfCells)  # get number of escorts/loads
    for j in range(k + 1):
        n, r = divmod(n, NumOfCells)
        x, y = divmod(r, Ly)
        l.append((x, y))

    return l


# Input:
#    S - pointer to the DP table
#    old_s_int - state key (as integer)
#    Lx, Ly - diminsions of the warehouse
# Output:
#   list of 2-tuple pair that describes the movement of each escort from its origin
#   location to its destination location
def state2Act(S, old_s_int, Lx, Ly):
    new_s_int = S[old_s_int][1]
    old_s = int2ListTuple(old_s_int, Lx, Ly)
    new_s = int2ListTuple(new_s_int, Lx, Ly)
    return list(zip(old_s[1:], new_s[1:]))


# find if s is already a key of one the dicts S[d0],...,S[d]
def AlreadyInDicts(S, d0, d1, s):
    for i in range(max(0, d0), d1 + 1):
        if s in S[i]:
            return True
    return False


# obtain an integer representation of state that was created fromm unsorted
# escorts list and return the equvalent sorted one
# This function is used when the decision a state is decoded in order to find
# the next state in the sequence
# We keep the decision sorted by the order of the original state
# so we can know the exact movement of each escort
def sortIntState(cs, Lx, Ly):
    # first disabke cs to its integer components (not tuples)
    l = []
    NumOfCells = Lx * Ly;
    n, k = divmod(cs, NumOfCells)  # get number of escorts/loads
    for j in range(k + 1):
        n, r = divmod(n, NumOfCells)
        l.append((r))

    L = [k, l[0]] + sorted(l[1:])
    return (sum([L[i] * (NumOfCells ** i) for i in range(k + 2)]))


def checkComp(s1, s2, Lx, Ly):
    t1 = int2ListTuple(s1, Lx, Ly)
    t2 = int2ListTuple(s2, Lx, Ly)

    for i in range(len(t1)):
        if t1[i][0] != t2[i][0] and t1[i][1] != t2[i][1]:
            return False
    return True


"""
   [(1,1),(2,3)] -> {<1 1>  <2 3>}
"""


def tuple_opl(L):
    s = "{"

    for t in L:
        s += "<" + str(t[0]) + " " + str(t[1]) + "> "

    return s + "}"


# Sample locations of retrived load(s) and escorts for a given location list
# use
#     Locations =  sorted(set(itertools.product(range(Lx), range(Ly))))
# to create the list of locations first

def GeneretaeRandomInstance(seed, Locations, num_escorts, num_load=1):
    # We take the entire permutation to assure consistancey.
    # I.e., that larger sample with the same seed will have smaller one as subsets
    # This should have beem the default for random.sample but there is a bug
    random.seed(seed)

    # Not working correctly
    # ez = random.sample(Locations,num_escorts+num_load)
    # return  sorted(ez[:num_load]),  sorted(ez[num_load:])

    ez = random.sample(Locations, len(Locations))
    return sorted(ez[:num_load]), sorted(ez[num_load:(num_load + num_escorts)])


def str2range(s):
    a = s.split('-')
    if len(a) == 1:
        return range(int(a[0]), int(a[0]) + 1)
    elif len(a) == 2:
        return range(int(a[0]), int(a[1]) + 1)
    elif len(a) == 3:
        return range(int(a[0]), int(a[1]) + 1, int(a[2]))


"""  obtained  initial location of target loads+escorts and list of load moves 
and return these location after applying the moves
init_R, init_E  - initial locations
mv - list of load moves each element correspond to a single time 
     step and include all the moves at the time step  as ((x1,y1), (x2,y2))
     note that the moves are NOT escort moves 
"""


def play_moves(init_R, init_E, mv):
    R = copy.copy(init_R)
    E = copy.copy(init_E)
    for i in range(len(mv)):
        for j in range(len(mv[i])):
            e_index = E.index(mv[i][j][1])
            E[e_index] = mv[i][j][0]
            if mv[i][j][0] in R:
                r_index = R.index(mv[i][j][0])
                R[r_index] = mv[i][j][1]
    return R, E


""" Return all the target loads induced by escort movements """


# def target_load_path(init_R, init_E, mv):
#     A = copy.copy(init_R)
#     E = copy.copy(init_E)
#     moves_tup = []
#
#     for t in range(len(mv)):
#         for j in range(len(mv[t])):
#             e_index = E.index(mv[t][j][1])
#             E[e_index] = mv[t][j][0]
#             moves_tup.append((mv[t][j][1][0], mv[t][j][1][1], mv[t][j][0][0], mv[t][j][0][1], t, 2))
#             if mv[t][j][0] in A:
#                 r_index = A.index(mv[t][j][0])
#                 A[r_index] = mv[t][j][1]
#                 moves_tup.append((mv[t][j][0][0], mv[t][j][0][1], mv[t][j][1][0], mv[t][j][1][1], t, 1))
#     return moves_tup


def warm_start_movements(init_R, init_E, mv, T=0):
    R = copy.copy(init_R)
    E = copy.copy(init_E)
    moves_tup = []

    for t in range(len(mv)):
        for i in range(len(E)):
            e = E[i]
            em = [x for x in mv[t] if x[1] == e]
            if em:
                moves_tup.append((em[0][1][0], em[0][1][1], em[0][0][0], em[0][0][1], t, 2))
                E[i] = em[0][0]
            else:
                moves_tup.append((e[0], e[1], e[0], e[1], t, 2))

        for i in range(len(R)):
            r = R[i]
            rm = [x for x in mv[t] if x[0] == r]
            if rm:
                moves_tup.append((rm[0][0][0], rm[0][0][1], rm[0][1][0], rm[0][1][1], t, 1))
                R[i] = rm[0][1]
            else:
                moves_tup.append((r[0], r[1], r[0], r[1], t, 1))

    # sometimes we want T > len(mv) to give the search process some slack
    while t < T:
        t += 1
        for i in range(len(E)):
            e = E[i]
            moves_tup.append((e[0], e[1], e[0], e[1], t, 2))
        for i in range(len(R)):
            r = R[i]
            moves_tup.append((r[0], r[1], r[0], r[1], t, 1))

    return moves_tup


