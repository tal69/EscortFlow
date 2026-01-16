"""
Tal Raviv, talraviv@tau.ac.il,   22/12/2023
Second try with the LM heuristic

Basic idea
-   At a preprocessing step, find the closest output to each target cell and the directions (1 or 2) from which one should
    proceed from it toward the escort. In case of ties there may be even 3 such cells.
    A cell is called "used" with respect to the current takt if it was involved in a movement either as an origin or
    destination.

- At each takt
  1. Greedily move target loads to escort toward their output.
  2. Find all unused cells attractive_1 in the sense that they can facilitate q targets move in the next takt
     and that such a move is not possible otherwise (i.e) there is no other escort next to this target in a desired direction.
  3. Identify unused escorts that can move in the current takt to attractive_1 cells. Move
     these escorts to their adjacent attractive_1 cells  whenever this is possible collision wise. For cells that have
     more than one such option, prioritize movements in the direction of the nearest output cell.
  4. Recalculate attractive_1 after all the changes of 3. Remove from this set currently occupied cells and cells adjacent
     to targets that can already move in a desired direction in the next takt
  5. Repeat the following  for n=2,3,4,...
        - Indentify attractive_n cells which are unused cells that are adjacent to attractive_{n-1} cells.
        - Whenever possible, move all unused escorts to attractive_n cells, break ties in favor of movements toward the
          closest output cell
    Until no nore attratcive cell can be found are all escorts moved.
"""
import copy
import numpy as np

dirs = [(1,0), (0,1), (-1,0),(0,-1) ]
def CreateRoutingTable(Lx, Ly, O):

    dist = np.ones((Lx, Ly), dtype=np.int16) * (Lx + Ly)
    best_dir = np.empty((Lx, Ly), dtype=set)

    for (x, y) in O:
        dist[x, y] = 0
        best_dir[x, y] = set([])
    B_new = set(O)

    while B_new:
        B = copy.copy(B_new)
        B_new  = set([])
        for (x,y) in B:
            for d in range(4):
                x1, y1 = x + dirs[d][0], y + dirs[d][1]
                if x1 in range(Lx) and y1 in range(Ly):
                    if dist[x1,y1] > dist[x,y] + 1:
                        dist[x1, y1] = dist[x, y] + 1
                        best_dir[x1,y1] = set([(-dirs[d][0], -dirs[d][1])])
                        B_new.add((x1,y1))
                    elif dist[x1,y1] == dist[x,y] + 1:
                        best_dir[x1, y1].add((-dirs[d][0], -dirs[d][1]))

    return best_dir



def GreedyMoves(_A, _E, _used_cells, best_dir):
    moves = []
    for i in range(len(_A)):
        x, y = _A[i][0], _A[i][1]
        for dir in best_dir[x, y]:
            if (x, y) not in _used_cells:
                x1, y1 = x + dir[0], y + dir[1]
                if (x1, y1) not in _used_cells and (x1, y1) in _E:
                    _A[i] = (x1, y1)
                    _E[_E.index((x1, y1))] = (x, y)
                    _used_cells.update([(x, y), (x1, y1)])
                    moves.append(((x, y), (x1, y1)))
    return A,E, moves

def FindAttractive(_A, _E, _used_cells, best_dir):
    attractive = set([])
    for i in range(len(_A)):
        x, y = _A[i][0], _A[i][1]
        for dir in best_dir[x, y]:
            x1, y1 = x + dir[0], y + dir[1]
            if (x1, y1) not in _used_cells and (x1, y1) not in _E and (x1,y1) not in _A:
                attractive.add((x1,y1))

    return attractive

def Move2Attractive(_A, _E, _used_cells, best_dir, attractive):


best_dir = CreateRoutingTable(5, 5, [(2,0)])
used_cells = set([])


A = [(1,3), (2,4), (1,1)]
E = [(1,2), (3,2), (2,1),(1,0)]

print(f"A={A},  E={E}")

A,E, moves = GreedyMoves(A, E, used_cells, best_dir)
print(f"A={A},  E={E}")
print(moves)
print(used_cells)
print(FindAttractive(A, E, used_cells, best_dir))



