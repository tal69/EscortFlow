"""Greedy one-step escort heuristic for the PBS retrieval problem.

Originally written by Tal Raviv on 14/12/2023 and revised as
``OneStepHeuristic_v2.py`` on 10/3/2026.

The main entry point is ``OneStep()``, which advances the system by one takt.
``SolveGreedy()`` repeatedly applies ``OneStep()`` until all target loads
reach output cells.
"""
import itertools
from collections import defaultdict

import numpy as np


def build_dist_map(Lx, Ly, O):
    """Map each cell to its closest output and preferred movement direction."""
    dist_map = defaultdict(lambda: (-1, -1))
    best_dist = defaultdict(lambda: Lx + Ly)
    for x in range(Lx):
        for y in range(Ly):
            for o in O:
                dist = abs(x - o[0]) + abs(y - o[1]) + x / Lx + y / (Lx * Ly)  # tie-break symmetry
                if dist < best_dist[(x, y)]:
                    best_dist[(x, y)] = dist
                    dist_map[(x, y)] = (int(np.sign(o[0] - x)), int(np.sign(o[1] - y)))
    return dist_map


def OneStep(Lx, Ly, O, _A, _E, dist_map, acyclic=False, retrieval_mode="continue"):
    """Advance the heuristic by one takt.

    Args:
        Lx, Ly: PBS dimensions.
        O: Output cell locations as a set of ``(x, y)`` tuples.
        _A: Dictionary mapping target-load locations to target ids.
        _E: Escort locations as a set of tuples.
        dist_map: Optional precomputed distance map from ``build_dist_map``.
        acyclic: If true, scan targets by increasing target id. Otherwise use
            dynamic priority by distance to the closest output.
        retrieval_mode: Target-handling mode after arrival at an output. In
            ``continue`` mode, loads at outputs stop being targets but do not
            create new escorts. In ``leave`` mode, they become escorts at the
            beginning of the next step.

    Returns:
        Tuple ``(A, E, moves)`` where ``A`` maps updated target locations to
        target ids, ``E`` is the updated escort set, and ``moves`` is the list
        of movements executed during the current time step.
    """

    def check_move(x0, y0, x1, y1):
        """Return whether an escort move is feasible in the current one-step plan."""
        if (x0, y0) not in E:
            return False  # only an escort can move
        dir_x, dir_y = int(np.sign(x1 - x0)), int(np.sign(y1 - y0))

        if dir_x != 0 and dir_y != 0:
            return False  # diagonal movements are not allowed

        x, y = x0, y0

        while (x, y) != (x1, y1):
            if (x, y) in cell_used:
                return False
            x += dir_x
            y += dir_y
            if (x, y) in E:
                return False

        return (x, y) not in cell_used

    def move_escort(x0, y0, x1, y1):
        E.remove((x0, y0))
        E_new.add((x1, y1))
        dir_x, dir_y = int(np.sign(x1 - x0)), int(np.sign(y1 - y0))

        x, y = x0, y0
        while (x, y) != (x1, y1):
            moves.append(((x + dir_x, y + dir_y), (x, y)))
            cell_used.add((x, y))

            if (x + dir_x, y + dir_y) in A:
                target_id = A.pop((x + dir_x, y + dir_y))
                A_new[(x, y)] = target_id
                moved_target_ids.add(target_id)
            x += dir_x
            y += dir_y

        cell_used.add((x, y))  # last cell not included in the loop

    def get_targets_in_priority_order():
        """Return target loads in the selected heuristic priority order."""
        if acyclic:
            return sorted(A.items(), key=lambda item: item[1])

        def priority_key(item):
            loc, _ = item
            closest_output = min(
                O,
                key=lambda o: (abs(loc[0] - o[0]) + abs(loc[1] - o[1]), o[0], o[1]),
            )
            return (
                abs(loc[0] - closest_output[0]) + abs(loc[1] - closest_output[1]),
                closest_output[0],
                closest_output[1],
                loc[0],
                loc[1],
            )

        return sorted(A.items(), key=priority_key)


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

    def find_zone_escorts(x_target, y_target):
        """Return escort sets for zones A, B, C, and D for the given target.

        The zones follow the paper definitions:
        A: same row/column as the target, in the output direction
        B: not in A, but sharing a row/column with an A-cell, and not on the
           target row/column
        C: not in A or B, but sharing a row/column with a B-cell
        D: not in A, B, or C, but sharing a row/column with a C-cell
        """
        A_dir = dist_map[(x_target, y_target)]
        has_horizontal_a = A_dir[0] != 0
        has_vertical_a = A_dir[1] != 0

        def in_horizontal_a(x):
            return has_horizontal_a and np.sign(x - x_target) == A_dir[0]

        def in_vertical_a(y):
            return has_vertical_a and np.sign(y - y_target) == A_dir[1]

        def is_zone_a(x, y):
            return (y == y_target and in_horizontal_a(x)) or (x == x_target and in_vertical_a(y))

        def is_zone_b(x, y):
            return not is_zone_a(x, y) and x != x_target and y != y_target and (in_horizontal_a(x) or in_vertical_a(y))

        def row_has_zone_b(y):
            return y != y_target and (has_horizontal_a or in_vertical_a(y))

        def col_has_zone_b(x):
            return x != x_target and (has_vertical_a or in_horizontal_a(x))

        def is_zone_c(x, y):
            return not is_zone_a(x, y) and not is_zone_b(x, y) and (row_has_zone_b(y) or col_has_zone_b(x))

        def is_zone_d(x, y):
            if is_zone_a(x, y) or is_zone_b(x, y) or is_zone_c(x, y):
                return False
            if has_horizontal_a and not has_vertical_a:
                return y == y_target and np.sign(x - x_target) == -A_dir[0]
            if has_vertical_a and not has_horizontal_a:
                return x == x_target and np.sign(y - y_target) == -A_dir[1]
            return False

        zone_a_escorts = set()
        zone_b_escorts = set()
        zone_c_escorts = set()
        zone_d_escorts = set()

        for x_escort, y_escort in E:
            escort = (x_escort, y_escort)
            if is_zone_a(x_escort, y_escort):
                zone_a_escorts.add(escort)
            elif is_zone_b(x_escort, y_escort):
                zone_b_escorts.add(escort)
            elif is_zone_c(x_escort, y_escort):
                zone_c_escorts.add(escort)
            elif is_zone_d(x_escort, y_escort):
                zone_d_escorts.add(escort)

        return zone_a_escorts, zone_b_escorts, zone_c_escorts, zone_d_escorts

    def sort_escorts_by_distance(x, y, EE):
        """Return escorts in ``EE`` sorted by distance to ``(x, y)`` and then lexicographically."""
        return sorted(EE, key=lambda escort: (abs(escort[0] - x) + abs(escort[1] - y), escort[0], escort[1]))

    def extend_move_for_lower_priority_targets(x0, y0, x1, y1):
        """Extend a zone-to-zone escort move if it can also promote lower-priority targets."""
        dir_x, dir_y = int(np.sign(x1 - x0)), int(np.sign(y1 - y0))
        if dir_x == 0 and dir_y == 0:
            return x1, y1

        best_x, best_y = x1, y1
        xx, yy = x1, y1

        while (xx + dir_x, yy + dir_y) not in E and (xx + dir_x, yy + dir_y) not in cell_used:
            xx += dir_x
            yy += dir_y
            if (xx, yy) in A:
                target_dir = dist_map[(xx, yy)]
                if (dir_x != 0 and target_dir[0] == dir_x) or (dir_y != 0 and target_dir[1] == dir_y):
                    best_x, best_y = xx, yy

        return best_x, best_y

    A_new, E_new = dict(), set()
    A = dict(_A)  # to make sure we are not changing the argument outside
    E = set(_E)
    moves = []

    for a in list(_A):
        if a in O:
            A.pop(a, None)
            if retrieval_mode == "leave":
                E.add(a)

    if len(A) == 0:
        return A, E, moves

    if dist_map is None:
        dist_map = build_dist_map(Lx, Ly, O)

    cell_used = set()

    moved_target_ids = set()

    for x in range(Lx):
        cell_used.add((x, -1))
        cell_used.add((x, Ly))

    for y in range(Ly):
        cell_used.add((-1, y))
        cell_used.add((Lx, y))


    A_sorted = get_targets_in_priority_order()

    for i in range(len(A_sorted)):
        (x0, y0), target_id = A_sorted[i]

        # First check that it didn't move already in the current step due to the movement of a load with higher priority
        if target_id not in moved_target_ids:
            # Try to move the item immediately in the right direction. An escort moves from Zone A->C/D
            #A_dir = A_sorted[i][2]  # we can also take it from dist_map
            zone_A_escorts, zone_B_escorts, zone_C_escorts, zone_D_escorts = find_zone_escorts(x0, y0)
            lst = sort_escorts_by_distance(x0, y0, zone_A_escorts)
            for (x_escort, y_escort ) in lst:
                if check_move(x_escort, y_escort, x0, y0):
                    dx, dy =np.sign(x0-x_escort), np.sign(y0-y_escort)  # only one of them is non zero
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x0, y0)
                    move_escort(x_escort, y_escort, x1, y1)
                    x0, y0 = x0-dx, y0-dy  # update the location of the current target load
                    zone_A_escorts, zone_B_escorts, zone_C_escorts, zone_D_escorts = find_zone_escorts(x0, y0)
                    break

            # make sure that no further escort movement in the current step moves high priority loads in the wrong direction
            cell_used.add((x0, y0))

        else:
            # If a higher-priority target already moved this load, skip further re-positioning attempts in the current takt.
            continue


        # try to move an escort from zone B -> A
        lst = sort_escorts_by_distance(x0, y0, zone_B_escorts)
        for (x_escort, y_escort) in lst:
            if np.sign(y_escort - y0) == dist_map[(x0, y0)][1]: # escort is "below" the load
                if check_move(x_escort, y_escort, x0, y_escort):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x0, y_escort)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
            else:  # escort is "above" the load
                if check_move(x_escort, y_escort, x_escort, y0):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort, y0)
                    move_escort(x_escort, y_escort, x1, y1)
                    break

        # try to move an escort from zone C to zone B
        lst = sort_escorts_by_distance(x0, y0, zone_C_escorts)
        for (x_escort, y_escort) in lst:
            dir_x, dir_y = dist_map[(x0, y0)]
            if dir_x == 0:
                if check_move(x_escort, y_escort, x_escort, y0 + dir_y):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort, y0 + dir_y)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
            elif dir_y == 0:
                if check_move(x_escort, y_escort, x0 + dir_x, y_escort):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x0 + dir_x, y_escort)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
            elif np.sign(y0 - y_escort) == dir_y:  # escort is "above" the load
                if check_move(x_escort, y_escort, x0 + dir_x, y_escort):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x0 + dir_x, y_escort)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
            else:  # escort is "below" the load
                if check_move(x_escort, y_escort, x_escort, y0 + dir_y):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort, y0 + dir_y)
                    move_escort(x_escort, y_escort, x1, y1)
                    break

        # Finally, check D -> C
        lst = sort_escorts_by_distance(x0, y0, zone_D_escorts)
        for (x_escort, y_escort) in lst:
            dir_x, dir_y = dist_map[(x0, y0)]
            if dir_y == 0:
                if check_move(x_escort, y_escort, x_escort, y_escort + 1):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort, y_escort + 1)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
                elif check_move(x_escort, y_escort, x_escort, y_escort - 1):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort, y_escort - 1)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
            elif dir_x == 0:
                if check_move(x_escort, y_escort, x_escort + 1, y_escort):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort + 1, y_escort)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
                elif check_move(x_escort, y_escort, x_escort - 1, y_escort):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort - 1, y_escort)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
            elif np.sign(y0 - y_escort) == dir_y:  # escort is "above" the load
                if check_move(x_escort, y_escort, x_escort + 1, y_escort):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort + 1, y_escort)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
                elif check_move(x_escort, y_escort, x_escort - 1, y_escort):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort - 1, y_escort)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
            else: # escort to the left or to the right of the load
                if check_move(x_escort, y_escort, x_escort, y_escort + 1):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort, y_escort + 1)
                    move_escort(x_escort, y_escort, x1, y1)
                    break
                elif check_move(x_escort, y_escort, x_escort, y_escort - 1):
                    x1, y1 = extend_move_for_lower_priority_targets(x_escort, y_escort, x_escort, y_escort - 1)
                    move_escort(x_escort, y_escort, x1, y1)
                    break

    A.update(A_new)
    return A, E | E_new, moves


def SolveGreedy(
    Lx, Ly, O, _A, _E, verbal=False, max_steps=500, acyclic=False, return_moves=False,
    retrieval_mode="continue"
):
    """Repeatedly apply ``OneStep`` until all target loads are retrieved."""
    ordered_targets = sorted(
        set(_A),
        key=lambda a: (min(abs(a[0] - o[0]) + abs(a[1] - o[1]) for o in O), a[0], a[1]),
    )
    A = {loc: idx + 1 for idx, loc in enumerate(ordered_targets)}
    E = set(_E)

    dist_map = build_dist_map(Lx, Ly, O)

    if verbal:
        print(f"initial A = {A}")
        print(f"initial E = {E}")

    flow_time = 0
    makespan = 0
    movements = 0
    move_history = []
    while A:
        flow_time += len(A)

        A, E, mv = OneStep(
            Lx, Ly, O, A, E, dist_map, acyclic=acyclic, retrieval_mode=retrieval_mode
        )
        if return_moves:
            move_history.append(mv)
        if verbal:
            print(f"After step: {makespan}")
            print(f"A {A}")
            print(f"E {E}")
            print(f"Moves {mv}")
        makespan += 1
        movements += len(mv)


        if makespan >= max_steps:
            print(f"Panic: could not solve in {max_steps} steps")
            exit(1)

    if return_moves:
        return makespan, flow_time, movements, move_history

    return makespan, flow_time, movements


if __name__ == '__main__':

    from PBSCom import *


    O = set([(0,0)])
    Lx, Ly = 10,10

    Locations = sorted(set(itertools.product(range(Lx), range(Ly))))
    total_makespan = 0
    total_flow_time = 0
    total_movements = 0
    num_of_loads = 4
    n = 1000
    for seed in range(n):
        if seed % 1 == 0:
          print(f"{seed} / {n}")
        A,E = GeneretaeRandomInstance(seed, Locations, 12, num_load=num_of_loads)
        makepan, flow_time, movements = SolveGreedy(Lx, Ly, O, set(A), set(E), verbal=False, max_steps=(Lx + Ly) * len(A) * 20 // len(E))
        total_makespan += makepan
        total_flow_time += flow_time
        total_movements += movements

    print("")
    print(f"makespan {total_makespan/n}")
    print(f"total flow_time {total_flow_time/n}")
    print(f"mean flowtime {total_flow_time/n/num_of_loads}")
    print(f"movements {total_movements/n}")
