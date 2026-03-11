"""
Tal Raviv, talraviv@tau.ac.il,   14/12/2023
           revised as OneStepHeuristic.py 6/3/2026

"""
import copy
import numpy as np
import math
from collections import defaultdict
import random

""" 
    run one step of the heuristic and can apply it repeatedly to obtain a complete solution
    At each step the algorithm consider the a single load (the one closest to to its output)
    and greedily moves escorts to the best location
    First escort that moves the selected load closer to the output are preformed, and next escort are moved arround the
    load from zone B->A, the from C->B and D->C
    (see zoning examples and scheme below inside the function  
    
"""


""" move escort from (x0,y0) to (x1,y0)  and update moves, and cell_used
  the input assumed to be valid movement (no checking is done)
  
  The function return nothing it just update A, A_new, E, E_new, moves, and cell_used"""

def move_escort( E, E_new, A, A_new, cell_used, x0, y0, x1,y1, moves):

    E.remove((x0,y0))
    E_new.add((x1,y1))
    dir_x, dir_y = int(np.sign(x1-x0)),int(np.sign(y1-y0))

    x,y = x0,y0
    while (x,y) != (x1,y1):
        moves.append( ((x+dir_x,y+dir_y),(x,y)))
        cell_used.add((x, y))

        if (x + dir_x, y + dir_y) in A:
            A_new.add((x,y))
            A.remove((x + dir_x, y + dir_y))
        x += dir_x
        y += dir_y

    cell_used.add((x, y)) # last cell not included in the loop

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


""" Return the location of the target load with the minimal rect distance to its nearest output cell"""
def find_closest_target(A,O):
    A_dist = math.inf
    for a in A:
        for o in O:
            dist = abs(a[0] - o[0]) + abs(a[1]-o[1])
            if dist < A_dist:
                A_dist = dist
                As = a
                A_dir = (int(np.sign(o[0] - a[0])), int(np.sign(o[1] - a[1])))

    return As, A_dir



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

def OneStep(Lx, Ly, O, _A, _E):

    A_new, E_new = set([]), set([])
    A = copy.deepcopy(_A)  # to make sure we are not changing the argument outside
    E = copy.deepcopy(_E)

    moves = []

    for a in _A:
        if a in O:
            A.remove(a)

    if len(A) == 0:
        return A, E, moves

    cell_used = set([])

    # find closet load to its output cell - to avoid cycling we focus on a single load at each step
    (x0,y0), A_dir = find_closest_target(A,O)

    # first try to move the item immediately in the right direction  Zone A->C
    moved = False

    if A_dir[0] != 0:  # the load is not at the same column of its closest output
        x, y = x0 + A_dir[0], y0
        while x in range(Lx) and (x,y) not in E:
            x += A_dir[0]

        # Look for another load in the same row to move
        xx,x1 = x0 - A_dir[0], x0


        while xx in range(Lx) and (xx,y) not in E:
            if (xx,y0) not in A:
                x1 = xx
            xx -= A_dir[0]

        if x in range(Lx):  # we have found an escort
            move_escort(E, E_new,A, A_new, cell_used, x, y, x1, y,  moves)
            x0 += A_dir[0]
            moved = True


    if A_dir[1] != 0 and not moved:  # the load is not at the row column of its closest output
        x, y = x0 , y0+ A_dir[1]
        while y in range(Ly) and (x, y) not in E:
            y += A_dir[1]

        yy, y1 = y0 - A_dir[1], y0

        while yy in range(Ly) and (x, yy) not in E:
            if (x0, yy) not in A:
                y1 = yy
            yy -= A_dir[1]

        if y in range(Ly):  # we have found an escort
            move_escort(E, E_new, A, A_new, cell_used, x,y, x, y1,  moves)
            y0 += A_dir[1]

    if (x0,y0) in O: # if current item arrive at an output
        A_new.remove((x0,y0))
        return A | A_new, E | E_new, moves

    # Next try to move escorts from zone B -> A
    if A_dir[1] != 0:
        x, y = x0 , y0+ A_dir[1]
        lst =[]
        while y in range(Ly):
            if (x0,y) not in E | E_new:
                lst += [(xx,y) for xx in range(Lx) if (xx,y)in E]
            y += A_dir[1]

        if lst:
            for e in lst:
                if check_move(E,cell_used,e[0],e[1],x0,e[1]):
                    move_escort(E, E_new, A, A_new, cell_used, e[0],e[1],x0,e[1], moves)


    if A_dir[0] != 0:
        x, y = x0 +A_dir[0] , y0
        lst =[]
        while x in range(Lx):
            if (x,y0) not in E | E_new:
                lst += [(x,yy) for yy in range(Ly) if (x,yy)in E]
            x += A_dir[0]

        if lst:
            for e in lst:
                if check_move(E,cell_used,e[0],e[1],e[0],y0):
                    move_escort(E, E_new, A, A_new, cell_used, e[0],e[1],e[0],y0, moves)

    # Next try to move escorts from zone C -> B and  D -> C
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

    # C -> B
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
                        move_escort(E, E_new, A, A_new, cell_used, x, y, x1, y, moves)
                    # if the escort is "to the right" it can move horizontally only
                    elif y == y0 and y0+A_dir[1] in range(Ly) and check_move(E,cell_used,x,y,x,y0+A_dir[1]):
                        y1 = y0+A_dir[1]
                        while y1+A_dir[1] in range(Ly) and  check_move(E,cell_used,x,y,x,y1+A_dir[1]):
                            y1 += A_dir[1]
                        move_escort(E, E_new, A, A_new, cell_used, x, y, x, y1, moves)

                    if x!=x0 and y!=y0:
                        # check if it is possible to move vertically beyond the item
                        if y0 + A_dir[1] in range(Ly) and check_move(E,cell_used,x,y,x,y0+A_dir[1]):
                            y1 = y0 + A_dir[1]
                            # if y1 + A_dir[1] in range(Ly)  and check_move(E, cell_used, x, y, x, y1 + A_dir[1]):
                            #      y1 += A_dir[1]   # Not helping
                            move_escort(E, E_new, A, A_new, cell_used, x, y, x, y1, moves)
                        elif x0+A_dir[0] in range(Lx) and check_move(E,cell_used,x,y,x0+A_dir[0],y):
                            x1 = x0 + A_dir[0]
                            # if x1 + A_dir[0] in range(Lx) and check_move(E, cell_used, x, y, x1 + A_dir[0], y):
                            #        x1 += A_dir[0] # Not helping
                            move_escort(E, E_new, A, A_new, cell_used, x, y, x1, y, moves)

                y -= A_dir[1]
            x -= A_dir[0]

    elif A_dir[1] == 0:   # target directly left or right to the output on the same row
        x = x0
        while x in range(Lx):
            for y in range(Ly):
                if y!= y0 and (x, y) in E:
                    if check_move(E, cell_used, x, y, x0 + A_dir[0], y):
                        move_escort(E, E_new, A, A_new, cell_used, x, y, x0 + A_dir[0], y, moves)
            x -= A_dir[0]

        #  Now check D->C
        x = x0 - A_dir[0]
        shift = 1 if (y0 < Ly - 1) else  -1  # move up or down  (maybe add randomization here)
        while x in range(Lx):
            if (x, y0) in E:
                if check_move(E, cell_used, x, y0, x, y0+shift):
                    move_escort(E, E_new, A, A_new, cell_used, x, y0, x, y0 + shift, moves)
            x -= A_dir[0]

    elif A_dir[0] == 0: # target "above" output
        for x in range(Lx):
            if x != x0:
                y = y0
                while y in range(Ly):
                    if (x, y) in E:
                        if check_move(E, cell_used, x, y, x, y0 + A_dir[1]):
                            move_escort(E, E_new, A, A_new, cell_used, x, y, x, y0 + A_dir[1], moves)
                    y -= A_dir[1]

        #  Now check D->C
        y = y0 - A_dir[1]
        shift = 1 if (x0 < Lx - 1) else -1  # move right or left  (maybe add randomization here)
        while y in range(Ly):
            if (x0, y) in E:
                if check_move(E, cell_used, x0, y, x0 + shift, y):
                    move_escort(E, E_new, A, A_new, cell_used, x0, y, x0 + shift, y, moves)
            y -= A_dir[1]

    return A | A_new, E | E_new, moves


def SolveGreedy(Lx, Ly, O, _A, _E, verbal=False, max_steps=100):
    A = copy.deepcopy(_A)  # to make sure we are not changing the argument outside
    E = copy.deepcopy(_E)

    if verbal:
        print(f"initial A = {A}")
        print(f"initial E = {E}")

    count = 0
    while A:

        A, E, mv = OneStep(Lx, Ly, O, A, E)
        if verbal:
            print(f"After step: {count}")
            print(f"A {A}")
            print(f"E {E}")
            print(f"Moves {mv}")
        count += 1

        if count >= max_steps:
            print(f"Panic: could not solve in {max_steps} steps")
            break
    return count


if __name__ == '__main__':

    from PBSCom import *


    O = set([(9,0)])
    Lx, Ly = 20,10
    # E = set([(2, 2),(0,0)])
    # A = set([(0,1), (0,2)])
    # OneStep(Lx, Ly, O, A, E)

    Locations = sorted(set(itertools.product(range(Lx), range(Ly))))
    total = 0
    n = 100
    for seed in range(n):
        A,E = GeneretaeRandomInstance(seed, Locations, 12, num_load=9)
        total += SolveGreedy(Lx, Ly, O, set(A), set(E), verbal=False, max_steps=(Lx + Ly) * len(A) * 2)
        if seed % 40 == 0:
          print(f"{seed} / {n}")


    print("")
    print(total/n)