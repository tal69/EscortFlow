"""
Tal Raviv, talraviv@tau.ac.il,   14/12/2023
"""
import copy
import numpy as np

""" 
    run one step of the heuristic, avoid step that are opposite to the one in the tabu list
    (used by the full heuristic to avoid cycles) 
    
    Note that O, A, E are accepted as numpy arrays
    Works only for NBM (LM) at 'leave' mode for now 
    
    For now (14/12/23) this heuristic works nicely for single load retrieval problems but tends to deadlock when 
    applied on multiple parallel retrieval 
    
"""


def FindAttractiveCells(Lx, Ly, A, A2O, E, cell_used):
    A_dir = np.sign(A2O - A)  # recalculate desired directions of each target load after movement
    num_loads = len(A)
    """ Find set of attractive cells that can facilitate a move in the next takt """
    attractive_cells = set([])
    for i in range(num_loads):
        if A_dir[i, 0] == 1:  # right
            x = A[i, 0] + 1
            y = A[i, 1]

            while x < Lx and not cell_used[x, y] and len(np.where((E == np.array([x, y])).all(axis=1))[0]) == 0 and len(
                    np.where((A == np.array([x, y])).all(axis=1))[0]) == 0:
                attractive_cells.add((x, y))
                x += 1

        if A_dir[i, 0] == -1:  # left
            x = A[i, 0] - 1
            y = A[i, 1]
            while x >= 0 and not cell_used[x, y] and len(np.where((E == np.array([x, y])).all(axis=1))[0]) == 0 and len(
                    np.where((A == np.array([x, y])).all(axis=1))[0]) == 0:
                attractive_cells.add((x, y))
                x -= 1

        if A_dir[i, 1] == 1:  # up
            x, y = A[i, 0], A[i, 1] + 1
            while y < Ly and not cell_used[x, y] and len(np.where((E == np.array([x, y])).all(axis=1))[0]) == 0 and len(
                    np.where((A == np.array([x, y])).all(axis=1))[0]) == 0:
                attractive_cells.add((x, y))
                y += 1

        if A_dir[i, 1] == -1:  # down
            x, y = A[i, 0], A[i, 1] - 1
            while y >= 0 and not cell_used[x, y] and len(np.where((E == np.array([x, y])).all(axis=1))[0]) == 0 and len(
                    np.where((A == np.array([x, y])).all(axis=1))[0]) == 0:
                attractive_cells.add((x, y))
                y -= 1
    return attractive_cells


""" Find location from which attractive cell are reachable in the next iteration """


def FindAttractiveCells2(Lx, Ly, attractive_cells, E, A):
    attractive_cells2 = set()
    for (x0, y0) in attractive_cells:
        x = x0 + 1  # cells to the right
        c = np.array([x, y0])
        while x < Lx and (x, y0) not in attractive_cells and len(np.where((E == c).all(axis=1))[0]) == 0 and len(
                np.where((A == c).all(axis=1))[0]) == 0:
            attractive_cells2.add((x, y0))
            x += 1

        x = x0 - 1  # cells to the left
        c = np.array([x, y0])
        while x >= 0 and (x,y0) not in attractive_cells and  len(np.where((E == c).all(axis=1))[0]) == 0 and len(
                np.where((A == c).all(axis=1))[0]) == 0:
            attractive_cells2.add((x, y0))
            x -= 1

        y = y0 + 1  # cells upward
        c = np.array([x0, y])
        while y < Ly and (x0,y) not in attractive_cells and  len(np.where((E == c).all(axis=1))[0]) == 0 and len(
                np.where((A == c).all(axis=1))[0]) == 0:
            attractive_cells2.add((x0, y))
            y += 1

        y = y0 - 1  # cells upward
        c = np.array([x0, y])
        while y >= 0 and (x0,y) not in attractive_cells and len(np.where((E == c).all(axis=1))[0]) == 0 and len(
                np.where((A == c).all(axis=1))[0]) == 0:
            attractive_cells2.add((x0, y))
            y -= 1
    return attractive_cells2


""" Move  escorts to first attractive cells if possible
    The tie breaking here is tuned toward output cells at the bottom and top
    For a system with output cells at the side check first left and right movements 
"""
def Move2Attractive(Lx, Ly, attractive_cells, E, A, cell_used, moves):

    for i in range(len(E)):

        # try moving escort down
        if not cell_used[E[i, 0], E[i, 1]]:
            x, y0, y = E[i, 0], E[i, 1], E[i, 1] - 1
            while y >= 0 and not cell_used[x, y] and len(np.where((E == np.array([x, y])).all(axis=1))[0]) == 0 and len(
                    np.where((A == np.array([x, y])).all(axis=1))[0]) == 0:
                if (x, y) in attractive_cells:
                    for yy in range(y, y0 + 1):
                        cell_used[x, yy] = 1
                    for yy in range(y, y0):
                        moves.append(((x, yy), (x, yy + 1)))
                    E[i, 1] = y
                    break
                y -= 1

        # try moving escort up
        if not cell_used[E[i, 0], E[i, 1]]:
            x, y0, y = E[i, 0], E[i, 1], E[i, 1] + 1
            while y < Ly and not cell_used[x, y] and len(
                    np.where((E == np.array([x, y])).all(axis=1))[0]) == 0 and len(
                    np.where((A == np.array([x, y])).all(axis=1))[0]) == 0:
                if (x, y) in attractive_cells:
                    for yy in range(y0, y + 1):
                        cell_used[x, yy] = 1
                    for yy in range(y0, y):
                        moves.append(((x, yy + 1), (x, yy)))
                    E[i, 1] = y
                    break
                y += 1

        if not cell_used[E[i, 0], E[i, 1]]:
            # try moving escort right
            x0, x, y = E[i, 0], E[i, 0] + 1, E[i, 1]
            while x < Lx and not cell_used[x, y] and len(np.where((E == np.array([x, y])).all(axis=1))[0]) == 0 and len(
                    np.where((A == np.array([x, y])).all(axis=1))[0]) == 0:
                if (x, y) in attractive_cells:
                    for xx in range(x0, x + 1):
                        cell_used[xx, y] = 1
                    for xx in range(x0, x):
                        moves.append(((xx + 1, y), (xx, y)))
                    E[i, 0] = x
                    break
                x += 1

        if not cell_used[E[i, 0], E[i, 1]]:
            # try moving escort left
            x0, x, y = E[i, 0], E[i, 0] - 1, E[i, 1]
            while x >= 0 and not cell_used[x, y] and len(np.where((E == np.array([x, y])).all(axis=1))[0]) == 0 and len(
                    np.where((A == np.array([x, y])).all(axis=1))[0]) == 0:
                if (x, y) in attractive_cells:
                    for xx in range(x, x0 + 1):
                        cell_used[xx, y] = 1
                    for xx in range(x, x0):
                        moves.append(((xx, y), (xx + 1, y)))
                    E[i, 0] = x
                    break
                x -= 1





""" Move escorts randomly"""
def MoveRandom(Lx, Ly, E, A, cell_used, moves):
    for i in range(len(E)):
        if not cell_used[E[i, 0], E[i, 1]]:
            # try moving escort right
            x0, y0 = E[i, 0],  E[i, 1]
            (dx, dy) = random.sample([(0,1), (0,-1), (1,0), (-1,0)],1)[0]
            x1, y1 = x0+dx, y0+dy
            if x1 in range(Lx) and y1 in range(Ly):
                c = np.array([x1, y1])
                if not cell_used[x1, y1] and len(np.where((E == c).all(axis=1))[0]) == 0 and len(
                    np.where((A == c).all(axis=1))[0]) == 0:
                    E[i, 0], E[i,1] = x1, y1
                    moves.append(((x1, y1),(x0,y0)))
                    cell_used[x0, y0] = 1
                    cell_used[x1, y1] = 1
                    targets_here = np.where((A == c).all(axis=1))[0]
                    if len(targets_here) > 0:
                        A[targets_here[0]] = (x0, y0)




"""   This function accept a PBS configuration and state and return new locations of escorts, target loads and moves
      That can be done in one time step and are represented presumably better state. It is used to handle 
      situation when we are unable to proceed with the rolling horizon heuristic when the MILP can't be solved
      or return the same state as the initial state
      Lx, Ly : dimension of the PBS units
      O   location of the output cells - as numpy array
      A, E initial locations of the target loads and escorts - as numpy array
      bm = block movement (if true)   NOT WORKING YET
"""

"""  For debugging and QA purposes """
""" Works for 'continue' mode for now"""
def test_step(Lx, Ly, _A,_E, A1, E1, mv):


    A, E = copy.copy(_A), copy.copy(_E)
    pbs_now = np.ones((Lx, Ly), dtype=np.int8)  # number of loads on each cell
    target_load_grid = np.zeros((Lx, Ly), dtype=np.int8)

    for e in E:
        if e[0] not in range(Lx) or e[1] not in range(Ly):
            print(f"Error: Initial location of escort {e} outside the {Lx}x{Ly} grid")
        pbs_now[e[0], e[1]] = 0

    for a in A:
        if a[0] not in range(Lx) or a[1] not in range(Ly):
            print(f"Error: Initial location of target load {a} outside the {Lx}x{Ly} grid")
        target_load_grid[a[0], a[1]] = 1

    pbs_next = copy.copy(pbs_now)

    for (l1, l2) in mv:

        if l1[0] not in range(Lx) or l1[1] not in range(Ly):
            print(f"Error:  movement {l1} to {l2}: origin outside the {Lx}x{Ly} grid")
            return False
        # else:
        #     if pbs_now[l1[0], l1[1]] == 0:
        #         print(f"Error: movement {l1} to {l2}: illegal movement from an empty location ")
        #         return False

        if abs(l1[0] - l2[0]) + abs(l1[1] - l2[1]) > 1:
            print(f"Error: illegal move from {l1} to {l2}: an L1 distance of more than 1")
            return False

        if l2[0] not in range(Lx) or l2[1] not in range(Ly):
            print(f"Error: movement {l1} to {l2}: destination outside the {Lx}x{Ly} grid")
            return False

        pbs_next[l2[0], l2[1]] += 1
        pbs_next[l1[0], l1[1]] -= 1

    over_crowded = list(zip(*np.where(pbs_now > 1)))
    under_crowded = list(zip(*np.where(pbs_now < 0)))
    if over_crowded:
        print(f"Error: Over crowded cell {over_crowded}")
        return False
    if under_crowded:
        print(f"Error: Under crowded cell {under_crowded}")
        return False

    for e in E1:
        if pbs_next[e[0], e[1]] != 0:
            print(f"Error: there should be an escort at {e}")
            return False

    for x in range(Lx):
        for y in range(Ly):
            if pbs_next[x,y] == 0 and (x,y) not in E1:
                print(f"Error: there is an escort in ({x}, {y}) but there is not")
                return False

    # for i in range(len(A)):
    #     if (A[i] != A1[i]).any() and ((A[i][0], A[i][1]), (A1[i][0], A1[i][1])) not in mv:
    #         print(f"Target load moved from {A[i]} to {A1[i]} but there is no such movement")
    #         return False


    return True
def OneStep(Lx, Ly, O, _A, _E, BM, retrieval_method='continue', tabu_set=set([])):
    A, E = copy.copy(_A), copy.copy(_E)
    L = np.array([Lx, Ly])
    cell_used = np.zeros((Lx, Ly), dtype=np.int8)
    if not BM:
        for a in A:
            cell_used[a[0], a[1]] = 1
    # first see if there are target loads on output cell and make them escorts (leave method for now

    retain = []
    new_escorts = 0
    moves = []
    for i in range(A.shape[0]):
        arrived = np.where((O == A[i]).all(axis=1))[0]
        if len(arrived) == 0:
            retain.append(i)
        else:
            if retrieval_method == 'leave':
                E = np.concatenate((E, A[i].reshape(1, 2)), axis=0)
                new_escorts += 1
                moves.append((tuple(A[i]), (None, None)))
            elif retrieval_method == 'continue':
                # moves.append((tuple(A[i]), A[i]))
                pass
    A = A[retain]

    num_loads = A.shape[0]
    num_output = O.shape[0]
    num_escorts = E.shape[0]

    E_moved = np.zeros(num_escorts, dtype=np.int8)
    if new_escorts:
        E_moved[-new_escorts:] = 1

    if num_loads == 0:
        return A, E, moves

    # find closet output cell for each target load
    A2O = np.array([O[0]] * len(A))
    A_dist = np.array([abs(a[0] - O[0][0]) + abs(a[1] - O[0][1]) for a in A])
    for j in range(1, num_output):
        for i in range(num_loads):
            dist = np.sum(np.abs(O[j] - A[i]))
            if dist < A_dist[i]:
                A2O[i] = O[j]
                A_dist[i] = dist

    A_dir = np.sign(A2O - A)  # find desired directions of each target load

    if BM:
        # targets_not_yet_moved = set([(A[i,0], A[i,1]) for i in range(A.shape[0])])
        """  first look for escorts that can promote target loads in one move  """
        for moving_escort in range(num_escorts):
            x0, y0 = tuple(E[moving_escort])
            targets_move = []

            # search escort moves up
            move_to = y0
            y1 = y0 + 1
            while y1 < Ly and not cell_used[x0, y1]:
                c = np.array([x0, y1])
                escorts_here = np.where((E == c).all(axis=1))[0]
                if len(escorts_here) > 0:
                    break
                targets_here = np.where((A == c).all(axis=1))[0]
                if len(targets_here) > 0:
                    if A_dir[targets_here[0], 1] == -1:
                        move_to = y1
                        targets_move.append(targets_here[0])
                    else:
                        break
                y1 += 1
            if move_to > y0:
                for y in range(y0, move_to + 1):
                    cell_used[x0, y] = 1  # mark used cells
                for target in targets_move:  # move affected target loads down
                    A[target, 1] -= 1
                    # targets_not_yet_moved.remove(tuple(A[target]))
                E[moving_escort] = np.array([x0, move_to])
                for y in range(y0, move_to):  # moving loads
                    moves.append(((x0, y + 1), (x0, y)))
                continue

            # search escort moves down
            y1 = y0 - 1
            while y1 >= 0 and not cell_used[x0, y1]:
                c = np.array([x0, y1])
                escorts_here = np.where((E == c).all(axis=1))[0]
                if len(escorts_here) > 0:
                    break
                targets_here = np.where((A == c).all(axis=1))[0]
                if len(targets_here) > 0:
                    if A_dir[targets_here[0], 1] == +1:
                        move_to = y1
                        targets_move.append(targets_here[0])
                    else:
                        break
                y1 -= 1
            if move_to < y0:
                for y in range(move_to, y0 + 1):
                    cell_used[x0, y] = 1  # mark used cells
                for target in targets_move:  # move affected target loads up
                    A[target, 1] += 1
                    # targets_not_yet_moved.remove(tuple(A[target]))
                E[moving_escort] = np.array([x0, move_to])
                for y in range(move_to, y0):  # moving loads
                    moves.append(((x0, y), (x0, y + 1)))
                continue

            # escort move right
            move_to = x0
            x1 = x0 + 1
            while x1 < Lx and not cell_used[x1, y0]:
                c = np.array([x1, y0])
                escorts_here = np.where((E == c).all(axis=1))[0]
                if len(escorts_here) > 0:
                    break
                targets_here = np.where((A == c).all(axis=1))[0]
                if len(targets_here) > 0:
                    if A_dir[targets_here[0], 0] == -1:
                        move_to = x1
                        targets_move.append(targets_here[0])
                    else:
                        break
                x1 += 1
            if move_to > x0:
                for x in range(x0, move_to + 1):
                    cell_used[x, y0] = 1  # mark used cells
                for target in targets_move:  # move affected target loads left
                    A[target, 0] -= 1
                    # targets_not_yet_moved.remove(tuple(A[target]))
                E[moving_escort] = np.array([move_to, y0])
                for x in range(x0, move_to):  # moving loads
                    moves.append(((x + 1, y0), (x, y0)))
                continue

            # escort move left
            x1 = x0 - 1
            while x1 >= 0 and not cell_used[x1, y0]:
                c = np.array([x1, y0])
                escorts_here = np.where((E == c).all(axis=1))[0]
                if len(escorts_here) > 0:
                    break
                targets_here = np.where((A == c).all(axis=1))[0]
                if len(targets_here) > 0:
                    if A_dir[targets_here[0], 0] == 1:
                        move_to = x1
                        targets_move.append(targets_here[0])
                    else:
                        break
                x1 -= 1
            if move_to < x0:
                for x in range(move_to, x0 + 1):
                    cell_used[x, y0] = 1  # mark used cells
                for target in targets_move:  # move affected target loads right
                    A[target, 0] += 1
                    # targets_not_yet_moved.remove(tuple(A[target]))
                E[moving_escort] = np.array([move_to, y0])
                for x in range(move_to, x0):  # moving loads
                    moves.append(((x, y0), (x + 1, y0)))
                continue

        attractive_cells = FindAttractiveCells(Lx, Ly, A, A2O, E, cell_used)

        """ Move remaining escorts to attractive cells if possible"""
        Move2Attractive(Lx, Ly, attractive_cells, E, A, cell_used, moves)

        """ Recalculate attractive cells with respect to new locations """
        attractive_cells = FindAttractiveCells(Lx, Ly, A, A2O, E, cell_used)

        """ Now find 2nd order attractive locations  - places from which attractive cells are reachable"""
        attractive_cells2 = FindAttractiveCells2(Lx, Ly, attractive_cells, E, A)

        """ And move the rest of the escort toward these cells"""
        Move2Attractive(Lx, Ly, attractive_cells2, E, A, cell_used, moves)

        attractive_cells3 = FindAttractiveCells2(Lx, Ly, attractive_cells2, E, A)

        Move2Attractive(Lx, Ly, attractive_cells3, E, A, cell_used, moves)

        """ if all else went wrong move each escort in some random feasible direction if any """
        # while np.sum(cell_used) == 0: # nothing moved
        #     MoveRandom(Lx, Ly, E, A, cell_used, moves)

    else:  # LM
        # See first if there are escorts that can move loads immediately at their desired directions
        for i in range(num_loads):
            for d in range(2):
                if A_dir[i, d] != 0:
                    new_loc = copy.copy(A[i])
                    new_loc[d] += A_dir[i, d]
                    if not cell_used[new_loc[0], new_loc[1]]:
                        moving_escort = np.where((E == new_loc).all(axis=1))[0]
                        if len(moving_escort) > 0:  # we found an escort that can do the job
                            moves.append((tuple(A[i]), tuple(E[moving_escort[0]])))
                            E[moving_escort[0]], A[i] = copy.copy(A[i]), copy.copy(E[moving_escort[0]])
                            cell_used[moves[-1][0][0], moves[-1][0][1]] = 1
                            cell_used[moves[-1][1][0], moves[-1][1][1]] = 1
                            E_moved[moving_escort[0]] = 1
                            break
        # next see if we can move the remaining escorts to better positions in one step

        for e in E:
            cell_used[e[0], e[1]] = 1  # We can do it only after moving target loads (that can be moved into escorts)

        E_dir = np.zeros((num_escorts, 2))
        for j in range(num_escorts):
            if not E_moved[j]:
                # find nearest load  (measure distance to next cell of the load in the desired direction
                best_dist = Lx + Ly + 1
                for i in range(num_loads):
                    for d in range(2):
                        if A_dir[i, d] != 0:
                            new_loc = copy.copy(A[i])
                            new_loc[d] += A_dir[i, d]
                            if new_loc[d] in range(L[d]):  # check that the target location is within the gird
                                if not cell_used[new_loc[0], new_loc[1]]:
                                    dist = np.sum(np.abs(E[j] - new_loc))
                                    if dist < best_dist:
                                        best_dist = dist
                                        nearest_loc = copy.copy(new_loc)

                if best_dist < Lx + Ly + 1:  # we found somewhere to move
                    E_dir[j] = np.sign(nearest_loc - E[j])
                    for d in range(2):
                        if E_dir[j, d] != 0:
                            new_loc = copy.copy(E[j])
                            new_loc[d] += E_dir[j, d]

                            if not cell_used[new_loc[0], new_loc[1]] and (tuple(E[j]), tuple(new_loc)) not in tabu_set:
                                moves.append((tuple(new_loc), tuple(E[j])))
                                E[j] = copy.copy(new_loc)
                                E_moved[j] = 1
                                cell_used[new_loc[0], new_loc[1]] = 1
                                break
        #  Now do U movements with the remaining escorts
        for j in range(num_escorts):
            if not E_moved[j]:
                for d in range(2):
                    if E_dir[j, d] != 0:  # but we could not go in this direction
                        alt_d = 1 - d  # so check perpendicular direction
                        new_loc = copy.copy(E[j])
                        new_loc[alt_d] += 1
                        if new_loc[alt_d] < L[alt_d] and not cell_used[new_loc[0], new_loc[1]]:
                            moves.append((tuple(new_loc), tuple(E[j])))
                            E[j] = copy.copy(new_loc)
                            E_moved[
                                j] = 1  # not sure that we need just in case we try other movements in the future (random?)
                            cell_used[new_loc[0], new_loc[1]] = 1
                            break  # cannot move in two directions
                        else:
                            new_loc[alt_d] -= 2  # the other side
                            if new_loc[alt_d] > 0 and not cell_used[new_loc[0], new_loc[1]]:
                                moves.append((tuple(new_loc), tuple(E[j])))
                                E[j] = copy.copy(new_loc)
                                E_moved[
                                    j] = 1  # not sure that we need just in case we try other movements in the future (random?)
                                cell_used[new_loc[0], new_loc[1]] = 1
                                break  # cannot move in two directions
    return A, E, moves


"""  This function applies OneStep repeatedly to create a solution for a static problem """


def SimpleHeuristic(Lx, Ly, O, A, E, BM, verbose=False):
    moves_current = []
    moves = []
    iter_num = 0
    loads_num = A.shape[0]

    while A.shape[0]:
        iter_num += 1
        if verbose:
            print(f"Step {iter_num}")
            print(f"A0= {A}")
            print(f"E0= {E}")

        A0, E0 = copy.copy(A), copy.copy(E)
        A, E, moves_current = OneStep(Lx, Ly, O, A, E, BM, args.retrieval_mode, set(moves_current))

        if not test_step(Lx, Ly, A0 ,E0 ,A, E,moves_current):
            print(f"Stopped at iteration {iter_num}")
            break
        if verbose:
            print(f"moves = {moves_current}")
            print("***************************")
        if iter_num > (Lx + Ly) * int(loads_num ** 0.5):
            print("It seems that I got into a deadlock")
            break

        """  For debugging """
        origs = [a[0] for a in moves_current]
        dests = [a[0] for a in moves_current]

        if len(dests) != len(set(dests)):
            print(f"Panic: iteration {iter_num} non unique origs {moves_current}")

        if len(origs) != len(set(origs)):
            print(f"Panic: iteration {iter_num} non unique dests{moves_current}")

        moves.append(moves_current)

    return moves


if __name__ == '__main__':
    # A = np.array([[3, 3]])
    # E = np.array([[4, 4], [0,3]])
    # print(f"A={A}, \n E={E}")
    # A, E, moves = OneStep(5, 5, np.array([[0, 0]]), A, E, True)
    # print('----------------')
    # print(f"A={A}, \n E={E}")
    # print(f"moves={moves}")
    #
    # exit(1)
    import pickle
    from PBSCom import *
    import argparse

    E = np.array([[3,0]])
    A = np.array([[2,1]])
    O = np.array([[2,0]])
    Lx, Ly = 5,5

    A1, E1, mv = OneStep(Lx,Ly,O,A,E, True)
    exit(0)


    parser.add_argument("-x", "--Lx", type=int, help="Horizontal dimension of the PBS unit", required=True)
    parser.add_argument("-y", "--Ly", type=int, help="Vertical dimension of the PBS unit", required=True)
    parser.add_argument("-O", "--output_cells", nargs='+', type=int, help="List of output locations", default=[0, 0])
    parser.add_argument("-e", "--escorts_range", help="Range of number of escorts  3, 3-10, or 1-30-3", default="5")
    parser.add_argument("-r", "--reps_range", help="Replication range (seed range) 1, 1-10, or 1-30-3", default="1")
    parser.add_argument("-l", "--load_num", type=int, help="Number of target loads (default 1)", default=1)
    parser.add_argument('-m', '--retrieval_mode', choices=['stay', 'leave', 'continue'],
                        help='Select a retrieval mode. Default is stay. Also note that for stay mode the number '
                             'of output cells must be greater or equal the number target loads', default='leave')
    parser.add_argument("-f", "--csv", help="File name of the result csv file", default="simple_pbs.csv")
    parser.add_argument("-a", "--export_animation", action="store_true",
                        help="Export animation files, one for each instance")
    parser.add_argument("-b", "--bm", action="store_true", help="Block movement regime (otherwise NBM, aka LM)")

    args = parser.parse_args()
    result_csv_file = args.csv
    reps_range = str2range(args.reps_range)
    escorts_range = str2range(args.escorts_range)

    if len(args.output_cells) % 2:
        print(f"Error: output cell coordinate list must be of even length - {args.output_cells}")
        exit(1)

    O = []
    for i in range(len(args.output_cells) // 2):
        for x in str2range(args.output_cells[i * 2]):
            for y in str2range(args.output_cells[i * 2 + 1]):
                O.append((x, y))
                if x not in range(args.Lx) or y not in range(args.Ly):
                    print(f"Error: escort ({x},{y}) not in range (0-{args.Lx - 1}) x (0-{args.Ly - 1})")
                    exit(1)

    if len(O) < args.load_num and args.retrieval_mode == "stay":
        print(
            f"Panic: In 'stay' settings number of output cells ({len(O)}) can not be smaller than the number of target loads ({load_num})")
        exit(1)

    if len(set(O)) < len(O):
        print(f"Panic: All output cells location must be unique {sorted(O)}")
        exit(1)

    Locations = sorted(set(itertools.product(range(args.Lx), range(args.Ly))))

    for escort_num in escorts_range:
        for rep in reps_range:
            random.seed(rep)
            A, E = GeneretaeRandomInstance(rep, Locations, escort_num, args.load_num)
            moves = SimpleHeuristic(args.Lx, args.Ly, np.array(O), np.array(A), np.array(E), args.bm, True)
            if args.export_animation:
                pickle.dump((args.Lx, args.Ly, O, E, A, moves), open(
                    f"script_simple_{'BM' if args.bm else 'LM'}_{args.retrieval_mode}_{args.Lx}_{args.Ly}_"
                    f"{escort_num}_{args.load_num}_{rep}.p", "wb"))