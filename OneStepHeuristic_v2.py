"""
Tal Raviv, talraviv@tau.ac.il,   14/12/2023
           revised as OneStepHeuristic_v2.py 10/3/2026

"""
import copy
from operator import itemgetter

import numpy as np
import math
from collections import defaultdict
import random

from tornado.escape import xhtml_escape

""" 
    run one step of the heuristic and can apply it repeatedly to obtain a complete solution
    At each step the algorithm consider the a single load (the one closest to to its output)
    and greedily moves escorts to the best location
    First escort that moves the selected load closer to the output are preformed, and next escort are moved arround the
    load from zone B->A, the from C->B and D->C
    (see zoning examples and scheme below inside the function  
    
    10/3/2026 update to handle movement of several targets at each step
    
"""




""" Check that a move is feasible at the current state and given the loads that already moved in the current step """
def check_move( E, cell_used, x0, y0, x1,y1):
    if (x0,y0) not in E:
        return False  # only an escort can move
    dir_x, dir_y = int(np.sign(x1 - x0)), int(np.sign(y1 - y0))

    if dir_x != 0 and dir_y != 0:
        return False  # diagonal movements are not allowed

    x, y = x0, y0

    while (x, y) != (x1, y1):
        if (x,y) in cell_used:
            return False
        x += dir_x
        y += dir_y
        if (x, y) in E:
            return False

    if (x, y) in cell_used:
        return False
    else:
        return True




def build_dist_map(Lx, Ly, O):

    dist_map = defaultdict(lambda: (Lx+Ly, (-1-1), (-1-1)))
    for x in range(Lx):
        for y in range(Ly):
            for o in O:
                dist = abs(x - o[0]) + abs(y - o[1]) + x / Lx + y / (Lx*Ly)  # to break symetry
                if dist < dist_map[(x,y)][0]:
                    dist_map[(x,y)] = (dist, (int(np.sign(o[0] -x)), int(np.sign(o[1] - y))), o)
    return dist_map




"""   OneStep(Lx, Ly, O, _A, _E)
      This function accept a PBS configuration and state and return new locations of escorts, target loads and moves
      That can be done in one time step and are represented presumably better state. It is used to handle 
      situation when we are unable to proceed with the rolling horizon heuristic when the MILP can't be solved
      or return the same state as the initial state
      Lx, Ly : dimension of the PBS units
      O   location of the output cells - as set of tuples
      _A, _E initial locations of the target loads and escorts - as sets of tuples
      All locations are zero based index (botton left corner is (0,0) and tpp right is (Lx-1,Ly-1)
      
      :returns
            A,E - the new target and escort locations as set of tuples
            moves - List of movement, each movement is an ((x0,y0),(x1,y1)) tuple of tuples
      
      It is not clear if it leads to a deadlock free solution if applied repeatedly, but I don't have a counterexample yet
"""

def OneStep(Lx, Ly, O, _A, _E, dist_map):
    """ move escort from (x0,y0) to (x1,y0)  and update moves, and cell_used
      the input assumed to be valid movement (no checking is done)

      The function return nothing it just update A, A_new, E, E_new, moves, and cell_used"""

    def move_escort(x0, y0, x1, y1):

        E.remove((x0, y0))
        E_new.add((x1, y1))
        dir_x, dir_y = int(np.sign(x1 - x0)), int(np.sign(y1 - y0))

        x, y = x0, y0
        while (x, y) != (x1, y1):
            moves.append(((x + dir_x, y + dir_y), (x, y)))
            cell_used.add((x, y))

            if (x + dir_x, y + dir_y) in A:
                A_new.add((x, y))
                A.remove((x + dir_x, y + dir_y))
                target_moved[(x + dir_x, y + dir_y)] = (x, y)
            x += dir_x
            y += dir_y

        cell_used.add((x, y))  # last cell not included in the loop

    """" Return the location of the target loads with sorted by rect distance to its nearest output cell + their directions """
    def find_closest_targets():

        As = []
        for a in A:
            min_dist = math.inf
            for o in O:
                dist = abs(a[0] - o[0]) + abs(a[1] - o[1]) + a[0] * 0.01 + a[1] * 0.0001  # to break symetry
                if dist < min_dist:
                    min_dist = dist
                    a_dir = (int(np.sign(o[0] - a[0])), int(np.sign(o[1] - a[1])))
            As += [(min_dist, a, a_dir)]

        return sorted(As, key=lambda x: x[0])

    """  Find the closest escort that can be used to move a target load at location x_target, y_target in the right direction 
         (i.e., an escort in Zone A) and return its location, or (-1,-1) if no escort is found  """

    def find_closest_escort(x_target, y_target):
        min_dist = math.inf
        best_escort = (-1, -1)

        dist_x, dist_y = dist_map[(x_target, y_target)][2]

        for (x_escort, y_escort) in E:
            cell_dist, A_dir, ez  = dist_map[(x_target, y_target)]
            if (x_escort == x_target and np.sign(y_escort - y_target) == A_dir[1] or
                    y_escort == y_target and np.sign(x_escort - x_target) == A_dir[0]):
                dist = abs(x_escort - x_target) / (Lx**6) + abs(y_escort - y_target) / (Ly**6)
                dist = abs(x_escort - x_target) / (0.1+dist_x) + abs(y_escort - y_target) / (0.1+dist_y)
                #  Both seems ok


                if dist < min_dist and check_move(E, cell_used, x_escort, y_escort, x_target, y_target):
                    min_dist = dist
                    best_escort = (x_escort, y_escort)

        return best_escort

    A_new, E_new = set([]), set([])
    A = copy.deepcopy(_A)  # to make sure we are not changing the argument outside
    E = copy.deepcopy(_E)
    moves = []

    for a in _A:
        if a in O:
            A.remove(a)

    if len(A) == 0:
        return A, E, moves


    dist_map = build_dist_map(Lx, Ly,O)

    cell_used = set([])

    target_moved = {}

    for x in range(Lx):
        cell_used.add( (x,-1))
        cell_used.add( (x,Ly))

    for y in range(Ly):
        cell_used.add( (-1,y))
        cell_used.add( (Lx,y))


    # sort targets by distance from their output
    A_sorted = find_closest_targets()

    for i in range(len(A_sorted)):
        x0, y0 = A_sorted[i][1]

        # First check that it didn't move already in the current step
        if (x0,y0) not in target_moved.keys():
            # Try to move the item immediately in the right direction. An escort moves from Zone A->C/D
            A_dir =A_sorted[i][2]
            x_escort,y_escort = find_closest_escort(x0, y0)

            if (x_escort,y_escort) != (-1,-1): # found an escort
                if x_escort == x0:  # vertical movement
                    # Look for another load in the same column to move
                    yy, y1 = y0, y0
                    while yy - A_dir[1] in range(Ly) and (x_escort, yy - A_dir[1]) not in E and (x_escort, yy - A_dir[1]) not in cell_used:
                        yy -= A_dir[1]
                        if (x0, yy) in A:
                            y1 = yy
                    move_escort( x_escort,y_escort, x_escort, y1)
                    y0 += A_dir[1]
                else:
                    xx, x1 = x0, x0
                    while xx - A_dir[0] in range(Lx) and (xx - A_dir[0], y_escort) not in E and (xx - A_dir[0],y_escort) not in cell_used:
                        xx -= A_dir[0]
                        if (xx, y0) in A:
                            x1 = xx
                    move_escort(x_escort, y_escort, x1, y_escort)
                    x0 += A_dir[0]

            # check if loads arrive at outputs  - we dont need to remove them just yet this will happe in the next iteration
            # for (x,y) in O:
            #     if (x,y) in A_new:
            #         A_new.remove((x,y))
            #         target_moved[(x, y)] = (-1, -1)


            # make sure that no further escort movement in the current step moves high priority loads in the wrong direction
            cell_used.add((x0,y0))
        else:  # if the item was moved before while handling previous targets update its location
            continue  #  For some reason it seems better to skip an item that moved altogether also from considering B->A and C->B movments

            if target_moved[(x0,y0)] == (-1, -1):
                continue  # if it was ejected from the system nothing to do with this load
            else:
                x0, y0 = target_moved[(x0,y0)]



        """
        Zoning - output at (0,0)
    
        BBBBCCCCC
        BBBBCCCCC
        AAAAXCCCC
        BBBBABBBB
        BBBBABBBB
    
        DCCCC
        XCCCC
        ABBBB
        ABBBB
    
        BBBCCC
        BBBCCC 
        AAAXDD    
        
        
        BBBBBC
        BBBBBC
        AAAAAX
        
        Zoning - output at (3,0)
        
        CCBBBBB   Target at a general position
        CCBBBBB
        CXAAAAA
        BABBBBB
        BABBBBB
           *
        
        CCCDCCC    Target above output
        CCCDCCC
        CCCXCCC
        BBBABBB
        BBBABBB
           * 
        
        
        BBBBCCC   Target at the same row as output
        BBBBCCC
        BBBBCCC   
        BBBBCCC   
        AAAAXDD
           *
           
        """

        # Next try to move escorts from zone B -> A
        moved = False
        if A_dir[1] != 0:
            x, y = x0, y0 + A_dir[1]
            lst = []
            while y in range(Ly):
                if (x0, y) not in E | E_new:
                    lst += [(xx, y) for xx in range(Lx) if (xx, y) in E]
                y += A_dir[1]

            if lst:
                for e in lst:
                    if check_move(E, cell_used, e[0], e[1], x0, e[1]):
                        move_escort( e[0], e[1], x0, e[1])
                        moved = True
                        break  # move only one load B-> A in for each target load
                        # could be made better by moving the closest possible escort from b->a

        if not moved and A_dir[0] != 0:
            x, y = x0 + A_dir[0], y0
            lst = []
            while x in range(Lx):
                if (x, y0) not in E | E_new:
                    lst += [(x, yy) for yy in range(Ly) if (x, yy) in E]
                x += A_dir[0]

            if lst:
                for e in lst:
                    if check_move(E, cell_used, e[0], e[1], e[0], y0):
                        move_escort( e[0], e[1], e[0], y0)

        # Next try to move escorts from zone C -> B and  D -> C

        # C -> B
        moved = False
        if A_dir[0]!= 0 and A_dir[1] != 0:
            x= x0
            while x in range(Lx):
                y = y0
                while y in range(Ly):
                    if (x,y) in E:
                        # if the escort is "above" it can move vertically only
                        if x == x0 and x0+A_dir[0] in range(Lx) and check_move(E,cell_used,x,y,x0+A_dir[0],y):
                            x1 = x0+A_dir[0]
                            while x1+A_dir[0] in range(Lx) and check_move(E,cell_used,x,y,x1+A_dir[0],y):
                                x1 += A_dir[0]
                            move_escort(x, y, x1, y)
                            moved = True
                        # if the escort is "to the right" it can move horizontally only
                        elif y == y0 and y0+A_dir[1] in range(Ly) and check_move(E,cell_used,x,y,x,y0+A_dir[1]):
                            y1 = y0+A_dir[1]
                            while y1+A_dir[1] in range(Ly) and  check_move(E,cell_used,x,y,x,y1+A_dir[1]):
                                y1 += A_dir[1]
                            move_escort(x, y, x, y1)
                            moved = True

                        if x!=x0 and y!=y0:
                            # check if it is possible to move vertically beyond the item
                            if y0 + A_dir[1] in range(Ly) and check_move(E,cell_used,x,y,x,y0+A_dir[1]):
                                y1 = y0 + A_dir[1]
                                # if y1 + A_dir[1] in range(Ly)  and check_move(E, cell_used, x, y, x, y1 + A_dir[1]):
                                #      y1 += A_dir[1]   # Not helping
                                move_escort( x, y, x, y1)
                                moved = True
                            elif x0+A_dir[0] in range(Lx) and check_move(E,cell_used,x,y,x0+A_dir[0],y):
                                x1 = x0 + A_dir[0]
                                # if x1 + A_dir[0] in range(Lx) and check_move(E, cell_used, x, y, x1 + A_dir[0], y):
                                #        x1 += A_dir[0] # Not helping
                                move_escort(x, y, x1, y)
                                moved = True

                    y -= A_dir[1]
                    if moved:
                        break
                x -= A_dir[0]
                if moved:
                    break

        elif not moved and A_dir[1] == 0:   # target directly left or right to the output on the same row
            x = x0
            while not moved and x in range(Lx):
                for y in range(Ly):
                    if y!= y0 and (x, y) in E:
                        if check_move(E, cell_used, x, y, x0 + A_dir[0], y):
                            move_escort( x, y, x0 + A_dir[0], y)
                            moved = True
                            break
                x -= A_dir[0]


            #  Now check D->C
            x = x0 - A_dir[0]
            shift = 1 if (y0 < Ly - 1) else  -1  # move up or down  (maybe add randomization here)
            while not moved and x in range(Lx):
                if (x, y0) in E:
                    if check_move(E, cell_used, x, y0, x, y0+shift):
                        move_escort( x, y0, x, y0 + shift)
                        moved = True
                x -= A_dir[0]

        elif A_dir[0] == 0: # target "above" output
            for x in range(Lx):
                if x != x0:
                    y = y0
                    while not moved and y in range(Ly):
                        if (x, y) in E:
                            if check_move(E, cell_used, x, y, x, y0 + A_dir[1]):
                                move_escort(x, y, x, y0 + A_dir[1])
                                moved = True
                        y -= A_dir[1]

            #  Now check D->C
            y = y0 - A_dir[1]
            shift = 1 if (x0 < Lx - 1) else -1  # move right or left  (maybe add randomization here)
            while not moved and y in range(Ly):
                if (x0, y) in E:
                    if check_move(E, cell_used, x0, y, x0 + shift, y):
                        move_escort(x0, y, x0 + shift, y)
                        moved = True
                y -= A_dir[1]

    return A | A_new, E | E_new, moves


def SolveGreedy(Lx, Ly, O, _A, _E, verbal=False, max_steps=100):
    A = copy.deepcopy(_A)  # to make sure we are not changing the argument outside
    E = copy.deepcopy(_E)

    dist_map = build_dist_map(Lx, Ly,O)

    if verbal:
        print(f"initial A = {A}")
        print(f"initial E = {E}")

    flow_time = 0
    makespan = 0
    movements= 0
    while A:
        flow_time += len(A)

        A, E, mv = OneStep(Lx, Ly, O, A, E, dist_map)
        if verbal:
            print(f"After step: {makespan}")
            print(f"A {A}")
            print(f"E {E}")
            print(f"Moves {mv}")
        makespan += 1
        movements += len(mv)


        if makespan >= max_steps:
            print(f"Panic: could not solve in {max_steps} steps")
            break
    return makespan, flow_time, movements


if __name__ == '__main__':

    from PBSCom import *


    O = set([(4,0), (13,0)])
    Lx, Ly = 18,5

    Locations = sorted(set(itertools.product(range(Lx), range(Ly))))
    total_makespan = 0
    total_flow_time = 0
    total_movements= 0
    num_of_loads = 2
    n = 100
    for seed in range(n):
        if seed % 50 == 0:
          print(f"{seed} / {n}")
        A,E = GeneretaeRandomInstance(seed, Locations, 20, num_load=num_of_loads)
        makepan, flow_time, movements = SolveGreedy(Lx, Ly, O, set(A), set(E), verbal=False, max_steps=(Lx + Ly) * len(A) * 2)
        total_makespan += makepan
        total_flow_time += flow_time
        total_movements += movements

    print("")
    print(f"makespan {total_makespan/n}")
    print(f"total flow_time {total_flow_time/n}")
    print(f"mean flowtime {total_flow_time/n/num_of_loads}")
    print(f"movements {total_movements/n}")
    dist_map = build_dist_map(Lx, Ly, O)
    print(dist_map)