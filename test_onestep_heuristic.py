"""Regression test for OneStepHeuristic_v2.

Checks, for a battery of random instances:

1. Termination: all target loads are retrieved, no panic exit.
2. Theoretical bound: in acyclic (fixed-priority) mode, the makespan respects
   the bound 4*n*d_max of the termination theorem in the papers.
3. Move legality per time step: axis-parallel moves only, every move origin is
   an escort, the sets of involved cells of the moves of a step are mutually
   disjoint, and no path crosses a stationary escort.
4. Optionally, a comparison against a baseline copy of the module:

       python3 test_onestep_heuristic.py /path/to/baseline_OneStepHeuristic_v2.py

Run without arguments to test the current module only:

       python3 test_onestep_heuristic.py
"""
import importlib.util
import os
import random
import sys

CODE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "OneStepHeuristic_v2.py")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def gen_instance(seed, Lx, Ly, n_esc, n_tgt, n_out):
    rnd = random.Random(seed)
    cells = [(x, y) for x in range(Lx) for y in range(Ly)]
    boundary = [c for c in cells if c[0] in (0, Lx - 1) or c[1] in (0, Ly - 1)]
    O = set(rnd.sample(boundary, n_out))
    rest = [c for c in cells if c not in O]
    picks = rnd.sample(rest, n_esc + n_tgt)
    return O, set(picks[n_esc:]), set(picks[:n_esc])  # O, A, E


def dmax(Lx, Ly, O):
    return max(
        min(abs(x - ox) + abs(y - oy) for ox, oy in O)
        for x in range(Lx) for y in range(Ly)
    )


def involved_cells(x0, y0, x1, y1):
    if x0 == x1:
        lo, hi = sorted((y0, y1))
        return {(x0, y) for y in range(lo, hi + 1)}
    lo, hi = sorted((x0, x1))
    return {(x, y0) for x in range(lo, hi + 1)}


def validate_run(mod, Lx, Ly, O, A0, E0, cap):
    """Run OneStep in a loop with per-step legality checks."""
    A = {loc: i + 1 for i, loc in enumerate(sorted(A0))}
    E = set(E0)
    for step in range(cap):
        A2, E2, mv, emv = mod.OneStep(
            Lx, Ly, set(O), A, set(E), None, acyclic=True,
            return_escort_moves=True,
        )
        used = set()
        origins = set()
        for (x0, y0, x1, y1) in emv:
            assert x0 == x1 or y0 == y1, "diagonal move"
            assert (x0, y0) in E, "origin not an escort at step start"
            cells = involved_cells(x0, y0, x1, y1)
            assert not (cells & used), "conflicting moves share a cell"
            used |= cells
            origins.add((x0, y0))
        for (x0, y0, x1, y1) in emv:
            cells = involved_cells(x0, y0, x1, y1)
            others = (E - origins) & (cells - {(x0, y0)})
            assert not others, "path crosses a stationary escort"
        assert len(E2) == len(E), "escort count changed"
        assert not (set(A2) & set(E2)), "target and escort share a cell"
        A, E = A2, set(E2)
        if not A or all(loc in O for loc in A):
            return step + 1
    raise AssertionError("did not finish within cap")


def bench(mod, acyclic, seeds=25):
    grids = [(5, 5), (9, 5), (7, 7), (13, 7)]
    tot = dict(ft=0, ms=0, mv=0, n=0)
    fails, viol = [], []
    per_instance = {}
    for (Lx, Ly) in grids:
        for seed in range(seeds):
            n_esc = 2 + seed % 5
            n_tgt = 1 + seed % 5
            n_out = 1 + seed % 2
            O, A, E = gen_instance(seed * 7 + Lx * 131 + Ly, Lx, Ly, n_esc, n_tgt, n_out)
            bound = 4 * len(A) * dmax(Lx, Ly, O)
            cap = bound + 50 if acyclic else 20000
            try:
                ms, ft, mv = mod.SolveGreedy(
                    Lx, Ly, O, set(A), set(E), max_steps=cap, acyclic=acyclic
                )
            except SystemExit:
                fails.append((Lx, Ly, seed))
                continue
            tot["ft"] += ft; tot["ms"] += ms; tot["mv"] += mv; tot["n"] += 1
            per_instance[(Lx, Ly, seed)] = (ms, ft, mv)
            if acyclic and ms > bound:
                viol.append((Lx, Ly, seed, ms, bound))
    return tot, fails, viol, per_instance


def main():
    baseline_path = sys.argv[1] if len(sys.argv) > 1 else None
    fixed = load("fixed_h", CODE)
    base = load("base_h", baseline_path) if baseline_path else None

    print("== Legality validation (current code, acyclic, 60 instances) ==")
    ok = 0
    for (Lx, Ly) in [(5, 5), (9, 5), (7, 7)]:
        for seed in range(20):
            O, A, E = gen_instance(9000 + seed * 13 + Lx, Lx, Ly, 2 + seed % 4, 1 + seed % 4, 1)
            cap = 4 * len(A) * dmax(Lx, Ly, O) + 50
            validate_run(fixed, Lx, Ly, O, A, E, cap)
            ok += 1
    print(f"  all {ok} runs legal and terminated")

    for acyclic in (True, False):
        mode = "acyclic (fixed priorities)" if acyclic else "dynamic ordering"
        tf, ff, vf, pf = bench(fixed, acyclic)
        print(f"\n== Mode: {mode} ==")
        print(f"  current : solved {tf['n']}/100  mean_flowtime={tf['ft']/tf['n']:.3f} "
              f"mean_makespan={tf['ms']/tf['n']:.3f} mean_moves={tf['mv']/tf['n']:.3f} "
              f"panics={len(ff)} bound_violations={len(vf)}")
        if base is not None:
            tb, fb, vb, pb = bench(base, acyclic)
            print(f"  baseline: solved {tb['n']}/100  mean_flowtime={tb['ft']/tb['n']:.3f} "
                  f"mean_makespan={tb['ms']/tb['n']:.3f} mean_moves={tb['mv']/tb['n']:.3f} "
                  f"panics={len(fb)} bound_violations={len(vb)}")
            diff = sum(1 for k in pb if k in pf and pb[k] != pf[k])
            print(f"  instances with different trajectories: {diff}/{min(len(pb), len(pf))}")
        assert not vf, f"bound violations: {vf[:3]}"
        assert not ff, f"panics: {ff[:3]}"

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
