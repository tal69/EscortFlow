"""Validated, solver-independent aggregation of revised PBS paper experiments.

Called automatically by ReproducePaper.py. All percentages are calculated per
instance before averaging. Initial-cutoff bounds/times remain method-specific;
best-known solutions include both models and their conditional extensions.
"""

import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

from RunStaticLP import MODEL_FILES, fingerprint, parse_instance, reference_key, validate_value
from RunStaticCampaign import LAYOUTS

METHODS = ("loadflow", "escortflow")
EPSILON = 0.001


def read_csv(path):
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Malformed or truncated CSV: {path}")
    return rows


def write_csv(path, rows):
    if not rows:
        raise ValueError("Cannot write an empty table: " + str(path))
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def instance_key(row):
    return (row.get("retrieval_mode", "leave"), row["Lx x Ly"], int(row["#Loads"]),
            int(row["# Escorts"]), int(row["seed"]))


def coordinates(row):
    return tuple(tuple(sorted(map(tuple, ast.literal_eval(row[field]))))
                 for field in ("IOs", "Escorts", "Target Loads"))


def integer(row, field):
    value = float(row[field])
    if not math.isfinite(value) or int(value) != value or value < 0:
        raise ValueError("Expected a nonnegative integer in " + field)
    return int(value)


def feasible_pairs(row):
    pairs = [(integer(row, "greedy_flowtime"), integer(row, "greedy_movements"))]
    for prefix in ("", "final_"):
        if row.get(prefix + "has_solution") == "1":
            pairs.append((integer(row, prefix + "flowtime"), integer(row, prefix + "movements")))
    return pairs


def percentage_gap(upper, lower):
    if not math.isfinite(lower) or lower < -1e-6 or lower > upper + 1e-6:
        raise ValueError(f"Lower bound {lower} is incompatible with upper bound {upper}")
    if upper == 0:
        return 0.0
    return 100 * (upper - min(upper, max(0, lower))) / upper


def improvement(greedy, best):
    if greedy == 0:
        if best != 0:
            raise ValueError("Cannot compute improvement relative to a zero greedy component")
        return 0.0
    return 100 * (greedy - best) / greedy


def component_metrics(row, best, lp):
    problem = parse_instance(row)
    coefficient = problem["flow_weight"]
    lx, ly = problem["Lx"], problem["Ly"]
    load_count = lx * ly - problem["escorts"]
    distances = [min(abs(x - ox) + abs(y - oy) for ox, oy in problem["outputs"])
                 for x, y in problem["targets"]]
    distance, maximum_distance = sum(distances), max(distances, default=0)
    if integer(row, "safe_movement_lower_bound") != distance:
        raise ValueError("Recorded distance lower bound differs from the initial instance")
    flow, moves = integer(row, "flowtime"), integer(row, "movements")
    if row["has_solution"] != "1" or row.get("error") or (flow, moves) < best:
        raise ValueError("Missing/error incumbent or invalid best-known reference")
    if integer(row, "scaled_objective") != coefficient * flow + moves:
        raise ValueError("Integer objective disagrees with R*FT+MV")
    if row.get("weighted_bound_consistent") != "1":
        raise ValueError("Inconsistent integer solver bound")
    bound = float(row["scaled_best_bound"]) if row["scaled_best_bound"] else 0.0
    if not math.isfinite(bound) or bound > coefficient * flow + moves + EPSILON:
        raise ValueError("Invalid initial-cutoff solver bound")
    # At a minimum-flow representative, T <= F-D+d_max and M <= K*T.
    # Thus L_Q <= (R+K)*F - K*(D-d_max). Both components are integers.
    flow_lb = max(distance, math.ceil((bound + load_count * (distance - maximum_distance)
                                     - EPSILON) / (coefficient + load_count)))
    movement_lb = max(distance, math.ceil(bound - coefficient * best[0] - EPSILON))
    flow_opt = row["phase1_flow_proven"] == "1"
    move_opt = row["weighted_proven"] == "1"
    if flow_opt:
        if flow != best[0]:
            raise ValueError("A best-known flow contradicts the initial FT proof")
        flow_lb = flow
    if move_opt:
        if (flow, moves) != best or not flow_opt:
            raise ValueError("A best-known pair contradicts the initial lexicographic proof")
        flow_lb, movement_lb = best
    total_cpu = float(row["weighted_cpu_time"])
    flow_cpu = float(row["first_flow_proof_cpu_time"]) if flow_opt else total_cpu
    if (not math.isfinite(total_cpu) or not math.isfinite(flow_cpu)
            or min(total_cpu, flow_cpu) < 0 or flow_cpu > total_cpu + 1e-6):
        raise ValueError("Invalid cutoff CPU times")
    if flow_opt and (integer(row, "first_flow_proof_flowtime") != flow
                     or float(row["first_flow_proof_runtime"]) > float(row["weighted_time_limit"]) + 1e-6):
        raise ValueError("A flow proof recorded after the cutoff cannot enter the comparison")
    validate_value(lp)
    if (lp["problem_sha256"] != problem["problem_sha256"]
            or reference_key(lp) != reference_key(problem)
            or int(lp["flow_weight"]) != coefficient
            or int(lp["physical_horizon"]) != problem["physical_horizon"]):
        raise ValueError("LP and integer instance, mode, coefficient, or horizon differ")
    reference = best[0] + best[1] / coefficient
    return dict(mode=problem["retrieval_mode"], layout=problem["layout"], outputs=len(problem["outputs"]),
                loads=problem["loads"], escorts=problem["escorts"], seed=problem["seed"], method=problem["method"],
                utilization=load_count / (lx * ly), naive_ft=distance,
                greedy_ft=integer(row, "greedy_flowtime"), greedy_mv=integer(row, "greedy_movements"),
                best_ft=best[0], best_mv=best[1], flow_weight=coefficient,
                ft_lower_bound=flow_lb, mv_lower_bound_at_best_ft=movement_lb,
                ft_opt=int(flow_opt), mv_opt=int(move_opt),
                ft_gap_pct=percentage_gap(best[0], flow_lb), mv_gap_pct=percentage_gap(best[1], movement_lb),
                lp_gap_pct=percentage_gap(reference, float(lp["lp_objective"])),
                ft_cpu_seconds=flow_cpu, total_cpu_seconds=total_cpu,
                ft_improvement_pct=improvement(integer(row, "greedy_flowtime"), best[0]),
                mv_improvement_pct=improvement(integer(row, "greedy_movements"), best[1]),
                cutoff_cpu_time_is_estimate=int(row["weighted_cpu_time_is_estimate"] == "1"))


def aggregate(records):
    groups = {}
    for row in records:
        key = tuple(row[k] for k in ("mode", "layout", "loads", "escorts", "method"))
        groups.setdefault(key, []).append(row)
    comparisons, solutions = [], []
    for key, rows in sorted(groups.items()):
        common = {field: rows[0][field] for field in ("mode", "layout", "outputs", "loads", "escorts", "utilization")}
        common["instances"] = len(rows)
        comparisons.append(dict(**common, method=key[-1],
            ft_opt_pct=100 * statistics.mean(r["ft_opt"] for r in rows),
            mv_opt_pct=100 * statistics.mean(r["mv_opt"] for r in rows),
            **{field: statistics.mean(r[field] for r in rows) for field in
               ("ft_gap_pct", "mv_gap_pct", "lp_gap_pct", "ft_cpu_seconds", "total_cpu_seconds")}))
        if key[-1] == "loadflow":
            solutions.append(dict(**common, **{field: statistics.mean(r[field] for r in rows) for field in
                ("naive_ft", "greedy_ft", "greedy_mv", "best_ft", "best_mv", "ft_improvement_pct", "mv_improvement_pct")}))
    return comparisons, solutions


def tabular_header(kind):
    if kind == 2:
        alignment = "lr|rrrrrrr|rrrrrrr"
        titles = [r"\multicolumn{2}{c|}{Config.} & \multicolumn{7}{c|}{Load-flow} & \multicolumn{7}{c}{Escort-flow} \\",
                  r"\midrule"]
        columns = [r"\makecell{PBS dim.\\(\# outs)}", r"\makecell{\#\\esc.}"]
        columns += [r"\makecell{" + title + "}" for _ in METHODS for title in
                    (r"FT\\opt.\\(\%)", r"MV\\opt.\\(\%)", r"FT\\gap\\(\%)", r"MV\\gap\\(\%)",
                     r"LP\\gap\\(\%)", r"FT\\CPU\\(s)", r"Total\\CPU\\(s)")]
    else:
        alignment = "lrr|rrrrrrr"
        titles = [r"\multicolumn{3}{c|}{Config.} & \multicolumn{7}{c}{Bounds and solutions} \\", r"\midrule"]
        columns = [r"\makecell{PBS dim.\\(\# outs)}", r"\makecell{\#\\esc.}", "Util."]
        columns += [r"\makecell{" + title + "}" for title in
                    (r"Naive\\FT LB", r"Greedy\\FT", r"Greedy\\MV", r"OPT\\FT", r"OPT\\MV",
                     r"FT\\improv.\\(\%)", r"MV\\improv.\\(\%)")]
    return [r"\scriptsize", r"\setlength{\tabcolsep}{2pt}", r"\renewcommand{\arraystretch}{1.08}",
            r"\begin{tabular*}{0.98\linewidth}{@{\extracolsep{\fill}}" + alignment + "@{}}",
            r"\toprule", *titles, " & ".join(columns) + r" \\", r"\midrule"]


def render_fragment(kind, rows):
    lines = tabular_header(kind)
    grouped = {}
    for row in rows:
        grouped.setdefault((row["layout"], row["escorts"]), {})[row.get("method", "shared")] = row
    previous = None
    layout_order = {f"{x}x{y}": index for index, (x, y, _, _) in enumerate(LAYOUTS)}
    for (layout, escorts), methods in sorted(grouped.items(), key=lambda x:
            (layout_order.get(x[0][0], len(LAYOUTS)), *map(int, x[0][0].split("x")), x[0][1])):
        if previous and previous != layout:
            lines.append(r"\midrule")
        previous = layout
        row = next(iter(methods.values()))
        cells = [layout.replace("x", r"$\times$") + f" ({row['outputs']})", str(escorts)]
        if kind == 2:
            if set(methods) != set(METHODS):
                raise ValueError("Incomplete method comparison")
            for method in METHODS:
                row = methods[method]
                cells += [f"{row['ft_opt_pct']:.2f}".rstrip("0").rstrip("."),
                          f"{row['mv_opt_pct']:.2f}".rstrip("0").rstrip("."),
                          "" if row["ft_opt_pct"] == 100 else f"{row['ft_gap_pct']:.2f}",
                          "" if row["mv_opt_pct"] == 100 else f"{row['mv_gap_pct']:.2f}", f"{row['lp_gap_pct']:.2f}"]
                for field in ("ft_cpu_seconds", "total_cpu_seconds"):
                    value = f"{row[field]:.2f}"
                    if round(row[field], 2) == min(round(r[field], 2) for r in methods.values()):
                        value = r"\textbf{" + value + "}"
                    cells.append(value)
        else:
            cells.append(f"{row['utilization']:.3f}")
            cells += [f"{row[field]:.2f}" for field in
                      ("naive_ft", "greedy_ft", "greedy_mv", "best_ft", "best_mv", "ft_improvement_pct", "mv_improvement_pct")]
        lines.append(" & ".join(cells) + r" \\")
    return "\n".join([*lines, r"\bottomrule", r"\end{tabular*}"]) + "\n"


def render_tables(directory, comparisons, solutions, seed_count):
    fragments = {}
    for kind, rows in ((2, comparisons), (3, solutions)):
        for mode in ("leave", "continue"):
            for loads in sorted({r["loads"] for r in rows if r["mode"] == mode}):
                chosen = [r for r in rows if r["mode"] == mode and r["loads"] == loads
                          and (kind != 2 or loads != 1 or 3 <= r["escorts"] <= 6)]
                name = f"table{kind}_{mode}_{'single' if loads == 1 else str(loads) + 'targets'}.tex"
                fragments[kind, mode, loads] = render_fragment(kind, chosen)
                (directory / name).write_text(f"% {seed_count} instances per configuration; generated by ReproducePaper.py\n" + fragments[kind, mode, loads])
        # Keep the current (a) single-target and (b) four-target leave panels
        # together in one portrait float with one shared caption.
        panels = []
        for loads, title in ((1, "Single target load"), (4, "Four target loads")):
            panels += [r"\subfloat[" + title + "]{%\n" + fragments[kind, "leave", loads] + "}"]
        caption = (f"Warm-started load-flow and escort-flow comparison over {seed_count} matched instances per row. "
                   "Optimality rates, integer bounds, and times use the initial cutoff. Component and LP gaps use "
                   "the common best-known solution, including extensions. Bold marks the lower time."
                   if kind == 2 else
                   f"Solution components and improvements over the greedy heuristic, averaged over {seed_count} "
                   "matched instances per row. OPT is the common best-known feasible pair over both formulations, including extensions.")
        (directory / f"table{kind}_leave.tex").write_text("\n".join([
            r"\begin{table}[p]", r"\centering", panels[0], r"\par\medskip", panels[1],
            r"\caption{" + caption + "}", r"\end{table}"]) + "\n")


def build_tables(root, manifest):
    root = Path(root)
    indexed = {}
    lp_values = {}
    inputs = []
    expected = {(c["mode"], f"{c['lx']}x{c['ly']}", c["loads"], c["escorts"], seed)
                for c in manifest["configurations"] for seed in manifest["seed_values"]}
    configurations = {(c["mode"], f"{c['lx']}x{c['ly']}", c["loads"], c["escorts"]): c
                      for c in manifest["configurations"]}
    source_hash = fingerprint({name: hashlib.sha256((root / "source" / name).read_bytes()).hexdigest() for name in MODEL_FILES})
    for mode in ("leave", "continue"):
        for method in METHODS:
            path = root / mode / f"{method}.csv"
            inputs.append(dict(file=str(path.relative_to(root)), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            for row in read_csv(path):
                key = instance_key(row)
                if key not in expected or row["formulation"] != method or key[0] != mode or (key, method) in indexed:
                    raise ValueError("Unexpected or duplicate integer instance")
                problem = parse_instance(row)
                config = configurations[key[:4]]
                if (problem["outputs"] != sorted(zip(config["outputs"][::2], config["outputs"][1::2]))
                        or row["warmstart"] != "1"
                        or any(float(row[name]) != manifest["settings"][name]
                               for name in ("weighted_time_limit", "extension_time_limit", "threads"))):
                    raise ValueError("Integer run does not match the saved experiment settings")
                indexed[key, method] = row
        path = root / mode / "lp_results.csv"
        inputs.append(dict(file=str(path.relative_to(root)), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        for row in read_csv(path):
            validate_value(row)
            layout, escorts, loads, seed, method, retrieval = reference_key(row)
            key = (retrieval, layout, loads, escorts, seed)
            if (key not in expected or (key, method) in lp_values or retrieval != mode
                    or row["source_sha256"] != source_hash):
                raise ValueError("Unexpected, duplicate, or changed-source LP result")
            lp_values[key, method] = row
    expected_keys = {(key, method) for key in expected for method in METHODS}
    if set(indexed) != expected_keys or set(lp_values) != expected_keys:
        raise ValueError("Missing integer or OPTIMAL LP results; cannot publish an incomplete table")
    shared_fields = ("flow_weight", "safe_flow_horizon", "weighted_physical_horizon", "threads", "warmstart",
                     "greedy_flowtime", "greedy_movements", "weighted_time_limit", "extension_time_limit")
    best, states = {}, {}
    for key in expected:
        lf, ef = (indexed[key, method] for method in METHODS)
        if coordinates(lf) != coordinates(ef) or any(lf[field] != ef[field] for field in shared_fields):
            raise ValueError("The formulations have different initial states or comparison settings")
        states[key] = coordinates(lf)
        best[key] = min(feasible_pairs(lf) + feasible_pairs(ef))
    references_used = 0
    for item in manifest["reference_inputs"]:
        for row in read_csv(root / item["file"]):
            key = instance_key(row)
            if key not in expected:
                continue
            if coordinates(row) != states[key]:
                raise ValueError("A historical reference has different coordinates for a selected seed")
            if (row.get("error") or row.get("movement_mode") != "BM"
                    or row.get("protocol") not in {"safe_integer_flow_timing_v4", "safe_integer_flow_timing_v5"}):
                raise ValueError("Failed or incompatible historical reference")
            best[key] = min(best[key], *feasible_pairs(row))
            references_used += 1
    for (key, method), row in indexed.items():
        if (row.get("final_flow_proven") == "1"
                and integer(row, "final_flowtime") != best[key][0]):
            raise ValueError("The common best-known flow contradicts an extension's FT proof")
        if (row.get("final_lexicographic_proven") == "1"
                and (integer(row, "final_flowtime"), integer(row, "final_movements")) != best[key]):
            raise ValueError("The common best-known pair contradicts an extension's lexicographic proof")
    records = [component_metrics(indexed[key, method], best[key], lp_values[key, method])
               for key in sorted(expected) for method in METHODS]
    comparisons, solutions = aggregate(records)
    directory = root / "tables"
    directory.mkdir(exist_ok=True)
    write_csv(directory / "instance_metrics.csv", records)
    write_csv(directory / "method_comparison.csv", comparisons)
    write_csv(directory / "bounds_and_solutions.csv", solutions)
    render_tables(directory, comparisons, solutions, len(manifest["seed_values"]))
    audit = dict(passed=True, matched_instances=len(expected), integer_records=len(indexed), optimal_lp_records=len(lp_values),
                 seed_values=manifest["seed_values"], reference_records_used=references_used, inputs=inputs,
                 initial_cutoff_seconds=manifest["settings"]["weighted_time_limit"],
                 gap_population="All selected instances, including zero gaps; average instance percentages",
                 ft_lower_bound="max(D, ceil((L_Q + K*(D-d_max) - 0.001)/(R+K))); cutoff proof overrides",
                 mv_lower_bound="max(D, ceil(L_Q - R*F_BK - 0.001)); conditional on best-known FT",
                 lp_gap="100*(F_BK+M_BK/R-LP)/(F_BK+M_BK/R); only OPTIMAL LP values",
                 improvement="100*(greedy_component-best_known_component)/greedy_component; same best pair for both components",
                 zero_reference_policy="Zero reference contributes zero percent only when its bound/component is also zero",
                 best_known="Common lexicographically smallest feasible (FT,MV) over greedy, both models, initial/final records, and supplied references",
                 timing="Initial-cutoff elapsed wall-clock time including model construction; unproved FT censored at cutoff",
                 cutoff_time_estimates=sum(r["cutoff_cpu_time_is_estimate"] for r in records),
                 best_known_not_certified=sum(not any(
                     indexed[key, method].get("final_lexicographic_proven") == "1" for method in METHODS) for key in expected))
    (directory / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit
