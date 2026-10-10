#!/usr/bin/env python3
"""Calculate both Table 2(b) LP relaxations from the copied four-target runs.

No arguments: read Experiment Oct2026/ beside this script; save and automatically
resume results_table2b_lp/. Rerun after copying additional integer results.
Recorded R, coordinates, retrieval mode, and each model's horizon are preserved.
See README.md and LP_REPRODUCIBILITY.md for Mac Studio commands and output files.
"""

import argparse
import csv
from datetime import datetime
import fcntl
import hashlib
import io
import json
import math
from pathlib import Path
import shlex
import shutil
import statistics
import subprocess
import sys

import RunStaticLP as lp
from RunStaticCampaign import LAYOUTS, automatic_threads

HERE = Path(__file__).resolve().parent
METHODS = ('loadflow', 'escortflow')
OUTPUTS = {f'{x}x{y}': sorted(zip(outs[::2], outs[1::2])) for x, y, _, outs in LAYOUTS}
SOURCE_FILES = ('RunTable2bLP.py', 'RunStaticLP.py', 'RunStaticCampaign.py',
                'static_integrated_lp.py', *lp.MODEL_FILES)
ENVIRONMENT_CHECK = '''import json,platform,shlex,sys
version=platform.python_version()
if sys.version_info < (3,10):
 raise SystemExit("Python " + version + " at " + sys.executable + " is too old. Use Python 3.10 or newer; installing another Python does not change this interpreter.")
try:
 import gurobipy as gp
except ImportError as exc:
 raise SystemExit("Cannot import gurobipy in " + sys.executable + ": " + str(exc) + "\\nInstall it in this interpreter with: " + shlex.quote(sys.executable) + " -m pip install gurobipy==13.0.3")
try:
 with gp.Env(empty=True) as env:
  env.setParam("OutputFlag",0); env.start()
  with gp.Model(env=env) as model:
   variables=model.addVars(2001,lb=0,obj=1)
   model.addConstr(variables.sum() >= 1)
   model.optimize()
   if model.Status != gp.GRB.OPTIMAL:
    raise SystemExit("Gurobi environment test did not solve to optimality.")
except gp.GurobiError as exc:
 if exc.errno == 10010:
  raise SystemExit("Gurobi's size-limited license is active. The paper LPs require your academic or commercial license. Configure that license for this interpreter, then rerun --check-environment.")
 raise SystemExit("Gurobi license/environment error " + str(exc.errno) + ": " + str(exc))
print(json.dumps(dict(host=platform.node(),platform=platform.platform(),python=sys.version,executable=sys.executable,gurobi=".".join(map(str,gp.gurobi.version())),license_check="passed_2001_variables")))
'''


def check_environment():
    """Check the actual interpreter and a license beyond the pip trial limit."""
    result = subprocess.run([sys.executable, '-c', ENVIRONMENT_CHECK], capture_output=True, text=True)
    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip() or 'No diagnostic output'
        raise ValueError(f'Environment check failed.\nPython executable: {sys.executable}\n'
                         f'Python version: {sys.version.split()[0]}\n{detail}\n'
                         'For a permanent Mac setup, run: bash RunTable2bLP.sh --setup python3.13')
    return json.loads(result.stdout.strip().splitlines()[-1])


def write_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def csv_text(rows, fields):
    text = io.StringIO(newline='')
    writer = csv.DictWriter(text, fieldnames=fields, lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return text.getvalue()


def save_csv(path, rows, fields):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(csv_text(rows, fields))
    temporary.replace(path)


def integer(row, field):
    value = float(row[field])
    if not math.isfinite(value) or value < 0 or int(value) != value:
        raise ValueError(f'Expected a nonnegative integer in {field}')
    return int(value)


def feasible_pairs(row):
    pairs = [(integer(row, 'greedy_flowtime'), integer(row, 'greedy_movements'))]
    for prefix in ('', 'final_'):
        if row.get(prefix + 'has_solution') == '1':
            pairs.append((integer(row, prefix + 'flowtime'), integer(row, prefix + 'movements')))
    return pairs


def read_inputs(args):
    """Read once into memory, so a solve always uses a consistent frozen snapshot."""
    records, problems, inputs, fields = {}, {}, [], []
    for layout in args.layouts:
        for method in METHODS:
            path = args.input_dir / f'table2b_{method}_{layout}.csv'
            if not path.exists():
                continue
            raw = path.read_bytes()
            reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig'), newline=''))
            if not reader.fieldnames:
                raise ValueError(f'Empty CSV: {path}')
            fields.extend(f for f in reader.fieldnames if f not in fields)
            selected = 0
            for line, row in enumerate(reader, 2):
                try:
                    if None in row or any(v is None for v in row.values()):
                        raise ValueError('Malformed or truncated CSV row')
                    problem = lp.parse_instance(row)
                    if (problem['layout'] != layout or problem['method'] != method
                            or problem['loads'] != 4 or problem['retrieval_mode'] != 'leave'
                            or problem['outputs'] != OUTPUTS[layout]):
                        raise ValueError('Expected the Table 2(b) four-target leave-mode configuration')
                    if problem['seed'] not in args.seed_values or problem['escorts'] not in args.escort_values:
                        continue
                    if row['protocol'] not in ('safe_integer_flow_timing_v4', 'safe_integer_flow_timing_v5'):
                        raise ValueError('Expected a recorded sufficient-integer-weight v4/v5 run')
                    if (row['warmstart'] != '1' or row['has_solution'] != '1'
                            or row['weighted_global_scope'] != '1' or row.get('error')):
                        raise ValueError('Integer run has an error, missing warm start, or insufficient global scope')
                    distance = sum(min(abs(x-a) + abs(y-b) for a, b in problem['outputs'])
                                   for x, y in problem['targets'])
                    bound = integer(row, 'safe_movement_bound')
                    expected_r = bound + 1 - (distance if row['protocol'].endswith('_v5') else 0)
                    safe_horizon = integer(row, 'safe_flow_horizon')
                    distances = [min(abs(x-a) + abs(y-b) for a, b in problem['outputs'])
                                 for x, y in problem['targets']]
                    if (bound != (problem['Lx'] * problem['Ly'] - problem['escorts']) * safe_horizon
                            or safe_horizon != integer(row, 'greedy_flowtime') - distance + max(distances)
                            or problem['flow_weight'] != expected_r or problem['physical_horizon'] < safe_horizon):
                        raise ValueError('Coefficient or physical horizon is inconsistent with the recorded protocol')
                    if any(min(pair) < distance for pair in feasible_pairs(row)):
                        raise ValueError('Feasible components contradict the Manhattan-distance lower bound')
                    key = lp.reference_key(problem)
                    if key in records:
                        raise ValueError(f'Duplicate instance label: {key}')
                    records[key], problems[key] = row, problem
                    selected += 1
                except (ValueError, KeyError, TypeError, SyntaxError) as exc:
                    raise ValueError(f'{path}:{line}: {exc}') from exc
            inputs.append(dict(file=str(path), sha256=hashlib.sha256(raw).hexdigest(), selected_rows=selected))
    if not records:
        raise ValueError(f'No selected Table 2(b) records in {args.input_dir}')
    # A one-period physical-horizon difference in archived v4 files is valid.
    # Preserve it, while requiring the same underlying instance and objective R.
    paired = {}
    for key, problem in problems.items():
        common = paired.setdefault(key[:4], problem)
        if any(problem[f] != common[f] for f in ('outputs', 'targets', 'escort_cells', 'flow_weight')):
            raise ValueError(f'The formulations have different coordinates or R at {key[:4]}')
    return records, problems, inputs, fields


def export_gaps(root, records, problems, args):
    """Use one best feasible pair per instance, including both extensions."""
    best = {}
    for key, row in records.items():
        candidate = min(feasible_pairs(row))
        best[key[:4]] = min(best.get(key[:4], candidate), candidate)
    result_path = root / 'lp_results.csv'
    saved = lp.read_reference(result_path) if result_path.exists() else {}
    source_hash = lp.fingerprint({name: lp.sha256(root / 'source' / name) for name in lp.MODEL_FILES})
    gaps = []
    for key, result in sorted(saved.items()):
        problem = problems.get(key)
        if (problem is None or result['problem_sha256'] != problem['problem_sha256']
                or result['source_sha256'] != source_hash
                or int(result['flow_weight']) != problem['flow_weight']
                or int(result['horizon']) != problem['horizon']
                or int(result['physical_horizon']) != problem['physical_horizon']):
            raise ValueError(f'Saved LP metadata differs from its recorded integer instance: {key}')
        flow, movements = best[key[:4]]
        reference = flow + movements / problem['flow_weight']
        relaxation = float(result['lp_objective'])
        if relaxation > reference + max(1e-6, 1e-8 * abs(reference)):
            raise ValueError(f'LP value exceeds the common best-known feasible objective: {key}')
        gap = max(0.0, 100 * (reference - max(0.0, relaxation)) / reference) if reference else 0.0
        gaps.append(dict(layout=key[0], escorts=key[1], loads=4, seed=key[3], method=key[4],
            flow_weight=problem['flow_weight'], horizon=problem['horizon'], physical_horizon=problem['physical_horizon'],
            best_ft=flow, best_mv=movements, best_objective=reference, lp_objective=relaxation,
            lp_flow=float(result['lp_flow']), lp_movements=float(result['lp_movements']), lp_gap_pct=gap))
    gap_fields = ('layout', 'escorts', 'loads', 'seed', 'method', 'flow_weight', 'horizon', 'physical_horizon',
                  'best_ft', 'best_mv', 'best_objective', 'lp_objective', 'lp_flow', 'lp_movements', 'lp_gap_pct')
    save_csv(root / 'lp_gaps.csv', gaps, gap_fields)
    summary = []
    for layout in args.layouts:
        for escorts in sorted(args.escort_values):
            for method in METHODS:
                mip_keys = [k for k in records if (k[0], k[1], k[4]) == (layout, escorts, method)]
                rows = [r for r in gaps if (r['layout'], r['escorts'], r['method']) == (layout, escorts, method)]
                complete = {r['seed'] for r in rows} == args.seed_values
                summary.append(dict(layout=layout, escorts=escorts, method=method,
                    expected_instances=len(args.seed_values), integer_instances_available=len(mip_keys),
                    lp_instances_complete=len(rows), complete_group=int(complete),
                    lp_gap_pct=statistics.mean(r['lp_gap_pct'] for r in rows) if complete else ''))
    save_csv(root / 'table2b_lp_summary.csv', summary, list(summary[0]))
    lines = [f'% Table 2(b) LP gaps; {len(args.seed_values)} selected seeds per complete group.',
             '% Empty cells mean a missing/incomplete integer or LP group.',
             r'\begin{tabular}{lr|rr}', r'\toprule',
             r'Dim. & \# esc. & Load-flow LP gap (\%) & Escort-flow LP gap (\%) \\', r'\midrule']
    for layout in args.layouts:
        for escorts in sorted(args.escort_values):
            values = [next(r['lp_gap_pct'] for r in summary if (r['layout'], r['escorts'], r['method'])
                           == (layout, escorts, method)) for method in METHODS]
            lines.append(' & '.join([layout.replace('x', r'$\times$'), str(escorts)] +
                                   [f'{v:.2f}' if v != '' else '' for v in values]) + r' \\')
    lines.extend([r'\bottomrule', r'\end{tabular}'])
    (root / 'table2b_lp_columns.tex').write_text('\n'.join(lines) + '\n')
    report = dict(integer_records_available=len(records), optimal_lp_records=len(gaps),
                  expected_lp_records=2 * len(args.layouts) * len(args.escort_values) * len(args.seed_values),
                  complete_groups=sum(r['complete_group'] for r in summary),
                  summary_groups=len(summary), best_known_reference='best lexicographic pair over both available initial/final integer runs and greedy',
                  averaging='per-instance LP percentages, then arithmetic mean over every selected seed; incomplete groups blank',
                  gap_formula='100*(FT_BK+MV_BK/R-Z_LP)/(FT_BK+MV_BK/R)',
                  updated=datetime.now().astimezone().isoformat())
    write_json(root / 'coverage.json', report)
    print(f"Saved {len(gaps)}/{len(records)} available LPs; {report['complete_groups']}/{len(summary)} complete groups.", flush=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input-dir', type=Path, default=HERE / 'Experiment Oct2026')
    parser.add_argument('--output-dir', type=Path, default=HERE / 'results_table2b_lp')
    parser.add_argument('--seeds', default='1-100', help='inclusive range or list, e.g. 1-10 or 1,7')
    parser.add_argument('--layouts', nargs='+', choices=list(OUTPUTS), default=list(OUTPUTS))
    parser.add_argument('--escorts', default='8,12,16', help='subset of 8,12,16, for a pilot')
    parser.add_argument('--workers', type=int, default=1, help='parallel LPs; one limits memory use')
    parser.add_argument('--threads', type=int, default=automatic_threads(), help='solver threads per LP; Mac performance cores by default')
    parser.add_argument('--time-limit', type=lp.positive_number, default=300)
    parser.add_argument('--retry-time-limit', type=lp.positive_number, default=600, help='barrier retry after TIME_LIMIT')
    parser.add_argument('--dry-run', action='store_true', help='validate and print coverage without writing or solving')
    parser.add_argument('--summarize-only', action='store_true', help='rebuild gaps from saved LPs without solving')
    parser.add_argument('--check-environment', action='store_true', help='check Python, gurobipy, and the full license using a tiny test; no campaign files or LP experiments')
    args = parser.parse_args(argv)
    try:
        args.seed_values, args.escort_values = lp.parse_range(args.seeds), lp.parse_range(args.escorts)
        if not args.seed_values or not args.escort_values or not args.escort_values <= {8, 12, 16}:
            raise ValueError('Select at least one seed and escorts from 8,12,16')
        if min(args.workers, args.threads) <= 0 or len(set(args.layouts)) != len(args.layouts):
            raise ValueError('Threads/workers must be positive and layouts distinct')
        args.input_dir, args.output_dir = args.input_dir.resolve(), args.output_dir.resolve()
        if args.output_dir == args.input_dir:
            raise ValueError('The LP output directory must be separate from the integer input directory')
    except ValueError as exc:
        parser.error(str(exc))
    return args


def forward_command(args, source):
    command = [sys.executable, '-u', str(source / 'RunTable2bLP.py'), '--input-dir', str(args.input_dir),
               '--output-dir', str(args.output_dir), '--seeds', args.seeds, '--escorts', args.escorts,
               '--layouts', *args.layouts, '--workers', str(args.workers), '--threads', str(args.threads),
               '--time-limit', str(args.time_limit), '--retry-time-limit', str(args.retry_time_limit)]
    if args.summarize_only:
        command.append('--summarize-only')
    return command


def main(argv=None):
    args = parse_args(argv)
    print(f'Python: {sys.executable} (version {sys.version.split()[0]})', flush=True)
    if args.check_environment:
        environment = check_environment()
        print(f"Gurobi {environment['gurobi']}: environment and license check passed.", flush=True)
        return 0
    root, source = args.output_dir, args.output_dir / 'source'
    selection = dict(seeds=sorted(args.seed_values), layouts=args.layouts, escorts=sorted(args.escort_values))
    manifest_path = root / 'manifest.json'
    if root.exists():
        if not manifest_path.exists():
            raise ValueError(f'Output folder exists without this runner\'s manifest: {root}; choose a fresh --output-dir')
        manifest = json.loads(manifest_path.read_text())
        if manifest['selection'] != selection:
            raise ValueError('Selected seeds/layouts/escorts differ from this saved run; use a fresh --output-dir')
        for name, expected in manifest['source_sha256'].items():
            if lp.sha256(source / name) != expected:
                raise ValueError('Frozen source changed: ' + name)
        if HERE != source and not args.dry_run:
            return subprocess.call(forward_command(args, source))
    elif args.summarize_only:
        raise ValueError('--summarize-only requires an existing LP output folder')
    records, problems, inputs, fields = read_inputs(args)
    print('Table 2(b): four targets, leave retrieval, both formulations.', flush=True)
    print('Inputs: ' + str(args.input_dir) + '\nLP results: ' + str(root), flush=True)
    for layout in args.layouts:
        counts = {m: sum(k[0] == layout and k[4] == m for k in records) for m in METHODS}
        print(f"  {layout}: LF {counts['loadflow']}, EF {counts['escortflow']} available integer records", flush=True)
    if args.dry_run:
        print(f'{len(records)} LP instances available; workers={args.workers}, threads={args.threads}. No files written or solves started.')
        return 0
    if not args.summarize_only:
        environment = check_environment()
    if not root.exists():
        root.mkdir(parents=True)
        source.mkdir()
        for name in SOURCE_FILES:
            shutil.copy2(HERE / name, source / name)
        manifest = dict(format_version=1, selection=selection,
                        source_sha256={name: lp.sha256(source / name) for name in SOURCE_FILES},
                        snapshots=[], sessions=[])
        write_json(manifest_path, manifest)
    with (root / '.runner.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('This LP output folder already has an active runner') from exc
        merged = csv_text([records[k] for k in sorted(records)], fields)
        digest = hashlib.sha256(merged.encode()).hexdigest()
        snapshot = root / 'inputs' / f'integer_{digest}.csv'
        snapshot.parent.mkdir(exist_ok=True)
        if snapshot.exists():
            if snapshot.read_text() != merged:
                raise ValueError('Saved input snapshot changed: ' + str(snapshot))
        else:
            snapshot.write_text(merged)
        if not any(item['sha256'] == digest for item in manifest['snapshots']):
            manifest['snapshots'].append(dict(file=str(snapshot.relative_to(root)), sha256=digest,
                                             inputs=inputs, records=len(records)))
        started = datetime.now().astimezone().isoformat()
        session = dict(started=started, summarize_only=args.summarize_only)
        if not args.summarize_only:
            session.update(environment=environment, workers=args.workers, threads=args.threads,
                           time_limit=args.time_limit, retry_time_limit=args.retry_time_limit)
        manifest['sessions'].append(session)
        write_json(manifest_path, manifest)
        exit_code = 0
        if not args.summarize_only:
            command = [sys.executable, '-u', str(source / 'RunStaticLP.py'), '--input', str(snapshot),
                       '--source-dir', str(source), '--workers', str(args.workers), '--threads', str(args.threads),
                       '--time-limit', str(args.time_limit), '--retry-time-limit', str(args.retry_time_limit),
                       '-f', str(root / 'lp_results.csv'), '--resume', '--extend']
            print('LP command: ' + shlex.join(command), flush=True)
            try:
                exit_code = subprocess.call(command)
            except KeyboardInterrupt:
                # subprocess.call waits for its interrupted child before returning.
                print('\nInterrupted; completed LPs remain saved. Rerun this command to continue.', flush=True)
                exit_code = 130
        export_gaps(root, records, problems, args)
        session.update(finished=datetime.now().astimezone().isoformat(), exit_code=exit_code)
        write_json(manifest_path, manifest)
        print('Gap summary: ' + str(root / 'table2b_lp_summary.csv'), flush=True)
        return exit_code


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError) as exc:
        print(f'Table 2(b) LP error: {exc}', file=sys.stderr)
        raise SystemExit(1)
