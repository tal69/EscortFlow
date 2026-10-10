#!/usr/bin/env python3
"""Generate all eight archived four-target LP batches independently of the MIPs.

No arguments: four paper layouts, both formulations, 8/12/16 escorts, seeds
1-100, leave mode, archived v4 R and horizons. Repeating the command resumes.
Copied integer CSVs are optional checks, never a filter on which LPs are run.
Uses the same instance/greedy/weight generator as the main runners' --lp mode.
"""

import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import RunStaticLP as lp
from RunStaticCampaign import LAYOUTS, automatic_threads
from RunTable2bLP import check_environment, save_csv, write_json
from static_generated_lp import SNAPSHOT_FILES, generated_row

HERE = Path(__file__).resolve().parent
METHODS = ('escortflow', 'loadflow')
SOURCE_FILES = tuple(dict.fromkeys(('RunFourTargetLP.py', 'RunTable2bLP.py',
                                   'RunStaticCampaign.py', *SNAPSHOT_FILES)))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--seeds', '-r', default='1-100', help='Inclusive range or comma-separated seeds (default: 1-100)')
    parser.add_argument('--escorts', '-e', default='8,12,16', help='Selected paper escort counts (default: 8,12,16)')
    parser.add_argument('--layouts', nargs='+', choices=[f'{x}x{y}' for x, y, _, _ in LAYOUTS],
                        default=[f'{x}x{y}' for x, y, _, _ in LAYOUTS])
    parser.add_argument('--methods', nargs='+', choices=METHODS, default=list(METHODS))
    parser.add_argument('--input-dir', type=Path, default=HERE/'Experiment Oct2026',
                        help='Optional copied integer results, used only to verify generated instances')
    parser.add_argument('--output-dir', type=Path, default=HERE/'results_four_target_lp',
                        help='Separate LP folder, automatically resumed (default: results_four_target_lp)')
    parser.add_argument('--workers', type=int, default=1, help='Parallel LPs within each part (default: 1)')
    parser.add_argument('--threads', type=int, default=automatic_threads(), help='Threads per LP (default: performance cores)')
    parser.add_argument('--time-limit', type=float, default=300, help='Initial LP solver limit in seconds')
    parser.add_argument('--retry-time-limit', type=float, default=600, help='Barrier retry limit after TIME_LIMIT')
    parser.add_argument('--dry-run', action='store_true', help='Show all eight parts and copied-result coverage without solving')
    parser.add_argument('--check-environment', action='store_true', help='Check this Python, NumPy, Gurobi and license, then exit')
    args = parser.parse_args(argv)
    try:
        args.seed_values = sorted(lp.parse_range(args.seeds))
        args.escort_values = sorted(lp.parse_range(args.escorts))
        if not args.seed_values or not args.escort_values or not set(args.escort_values) <= {8, 12, 16}:
            raise ValueError('Select seeds and escort counts from 8,12,16')
        if len(set(args.layouts)) != len(args.layouts) or len(set(args.methods)) != len(args.methods):
            raise ValueError('Layouts and methods cannot contain duplicates')
        if min(args.workers, args.threads) <= 0 or any(not math.isfinite(v) or v <= 0
                for v in (args.time_limit, args.retry_time_limit)):
            raise ValueError('Workers, threads and solver limits must be positive')
        args.input_dir, args.output_dir = args.input_dir.resolve(), args.output_dir.resolve()
        if (args.input_dir == args.output_dir or args.input_dir in args.output_dir.parents
                or args.output_dir in args.input_dir.parents):
            raise ValueError('Use a separate output directory outside the integer result directory')
    except ValueError as exc:
        parser.error(str(exc))
    return args


def batches(args):
    for lx, ly, _, flat in LAYOUTS:
        layout = f'{lx}x{ly}'
        if layout in args.layouts:
            for method in METHODS:
                if method in args.methods:
                    yield dict(name=f'{method}_{layout}', layout=layout, method=method,
                               lx=lx, ly=ly, outputs=sorted(zip(flat[::2], flat[1::2])))


def selection(args):
    return dict(protocol='v4', loads=4, retrieval_mode='leave', movement_mode='BM',
                seeds=args.seed_values, escorts=args.escort_values,
                parts=[b['name'] for b in batches(args)])


def environment():
    result = check_environment()
    numpy = subprocess.run([sys.executable, '-c', 'import numpy; print(numpy.__version__)'],
                           capture_output=True, text=True)
    if numpy.returncode:
        raise ValueError(f'NumPy is needed by the main instance/greedy generator.\n'
                         f'Install it in {sys.executable} with: python -m pip install numpy\n{numpy.stderr.strip()}')
    result['numpy'] = numpy.stdout.strip()
    return result


def read_references(args):
    """Check every available selected row, allowing any of the eight to be absent."""
    references, files = {}, []
    for batch in batches(args):
        path = args.input_dir/f"table2b_{batch['name']}.csv"
        if not path.is_file():
            continue
        raw = path.read_bytes()
        reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig'), newline=''))
        if not reader.fieldnames:
            raise ValueError(f'Empty integer CSV: {path}')
        selected = 0
        for line, row in enumerate(reader, 2):
            try:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError('Malformed or truncated integer CSV row; finish copying the file and rerun')
                if int(row['seed']) not in args.seed_values or int(row['# Escorts']) not in args.escort_values:
                    continue
                problem = lp.parse_instance(row)
                if (problem['layout'] != batch['layout'] or problem['method'] != batch['method']
                        or problem['loads'] != 4 or problem['retrieval_mode'] != 'leave'
                        or problem['outputs'] != batch['outputs']):
                    raise ValueError('Expected the four-target leave-mode paper configuration')
                if row.get('protocol') != 'safe_integer_flow_timing_v4':
                    raise ValueError('This launcher matches the currently running v4 campaign; the copied row uses another protocol')
                key = lp.reference_key(problem)
                if key in references:
                    raise ValueError(f'Duplicate integer instance: {key}')
                references[key] = problem
                selected += 1
            except (ValueError, KeyError, TypeError, SyntaxError) as exc:
                raise ValueError(f'{path}:{line}: {exc}') from exc
        files.append(dict(file=str(path), sha256=hashlib.sha256(raw).hexdigest(), selected_rows=selected))
    return references, files


def compare_references(rows, references):
    matched = 0
    for row in rows:
        problem = lp.parse_instance(row)
        recorded = references.get(lp.reference_key(problem))
        if recorded is not None:
            if problem != recorded:
                differences = [k for k in problem if k != 'problem_sha256' and problem[k] != recorded[k]]
                raise ValueError(f"Generated v4 instance differs from its copied result: {lp.reference_key(problem)}; "
                                 f"changed fields: {', '.join(differences)}. No LPs launched.")
            matched += 1
    return matched


def freeze(args):
    root = args.output_dir
    path = root/'campaign.json'
    if root.exists():
        if not path.is_file():
            raise ValueError('Output folder exists without its campaign manifest; choose a new --output-dir')
        saved = json.loads(path.read_text())
        if saved['selection'] != selection(args):
            raise ValueError('Seeds, escorts or parts changed; use a separate --output-dir for a different selection')
        for name, expected in saved['source_sha256'].items():
            if lp.sha256(root/'source'/name) != expected:
                raise ValueError('Frozen source changed: ' + name)
        for name, expected in saved['inputs_sha256'].items():
            if lp.sha256(root/'instances'/name) != expected:
                raise ValueError('Generated input snapshot changed: ' + name)
        return saved
    root.parent.mkdir(parents=True, exist_ok=True)
    pending = root.with_name(root.name+'.preparing')
    pending.mkdir()
    try:
        source = pending/'source'
        source.mkdir()
        hashes = {name: lp.sha256(HERE/name) for name in SOURCE_FILES}
        for name in SOURCE_FILES:
            shutil.copy2(HERE/name, source/name)
        if any(lp.sha256(source/name) != expected or lp.sha256(HERE/name) != expected
               for name, expected in hashes.items()):
            raise ValueError('Source changed while freezing; rerun after code updates finish')
        saved = dict(selection=selection(args), objective='FT + MV/R',
                     created_utc=datetime.now(timezone.utc).isoformat(),
                     source_sha256=hashes, inputs_sha256={})
        write_json(pending/'campaign.json', saved)
        pending.rename(root)
    finally:
        if pending.exists():
            shutil.rmtree(pending)
    return saved


def prepare(args, saved, references):
    root = args.output_dir
    (root/'instances').mkdir(exist_ok=True)
    (root/'parts').mkdir(exist_ok=True)
    prepared, checks = [], []
    total = len(args.seed_values)*len(args.escort_values)
    for batch in batches(args):
        filename = batch['name']+'.csv'
        path = root/'instances'/filename
        if filename in saved['inputs_sha256']:
            with path.open(newline='') as stream:
                rows = list(csv.DictReader(stream))
        else:
            if (root/'parts'/filename).exists() or Path(str(root/'parts'/filename)+'.manifest.json').exists():
                raise ValueError('LP results exist without a recognized instance snapshot: ' + filename)
            rows = []
            print(f"Preparing {batch['name']}: {total} instances from seeds (v4).", flush=True)
            for count in args.escort_values:
                for seed in args.seed_values:
                    rows.append(generated_row(batch['method'], batch['lx'], batch['ly'], batch['outputs'],
                                              count, 4, seed, 'leave', 'v4'))
                    if len(rows) % 25 == 0:
                        print(f"  {len(rows)}/{total} generated", flush=True)
            compare_references(rows, references)
            save_csv(path, rows, list(rows[0]))
            saved['inputs_sha256'][filename] = lp.sha256(path)
            write_json(root/'campaign.json', saved)
        matched = compare_references(rows, references)
        print(f"{batch['name']}: {len(rows)} generated; {matched} copied integer rows verified.", flush=True)
        checks.append(dict(part=batch['name'], generated=len(rows), integer_rows_verified=matched))
        prepared.append((batch, path, rows))
    return prepared, checks


def aggregate(args, prepared):
    root = args.output_dir
    source_hash = lp.fingerprint({name: lp.sha256(root/'source'/name) for name in lp.MODEL_FILES})
    results, coverage = [], []
    for batch, path, inputs in prepared:
        expected = {lp.reference_key(p): p for p in map(lp.parse_instance, inputs)}
        output = root/'parts'/path.name
        rows = lp.read_reference(output) if output.exists() else {}
        for key, row in rows.items():
            problem = expected.get(key)
            if (problem is None or row['problem_sha256'] != problem['problem_sha256']
                    or row['source_sha256'] != source_hash
                    or any(str(row[k]) != str(problem[k]) for k in
                           ('flow_weight', 'horizon', 'physical_horizon', 'retrieval_mode'))):
                raise ValueError(f'Saved LP metadata differs from the generated instance: {key}')
            results.append(row)
        coverage.append(dict(part=batch['name'], expected=len(expected), optimal=len(rows),
                             missing=len(expected)-len(rows)))
    save_csv(root/'lp_results.csv', results, lp.FIELDS)
    write_json(root/'coverage.json', dict(parts=coverage, expected=sum(c['expected'] for c in coverage),
                                         optimal=len(results), complete=all(c['missing'] == 0 for c in coverage)))
    return sum(c['missing'] for c in coverage)


def run(args, argv):
    print(f'Python: {sys.executable} ({sys.version.split()[0]})', flush=True)
    if args.check_environment:
        print(json.dumps(environment(), indent=2))
        return 0
    references, files = read_references(args)
    if args.dry_run:
        count = len(args.seed_values)*len(args.escort_values)
        for b in batches(args):
            available = sum(k[0] == b['layout'] and k[4] == b['method'] for k in references)
            print(f"{b['name']}: {count} LPs; {available} copied integer rows available for checking")
        print(f'Total: {len(list(batches(args)))*count} LPs. Protocol v4; four targets; leave mode.')
        print(f'Output: {args.output_dir}')
        return 0
    saved = freeze(args)
    source = args.output_dir/'source'
    if HERE != source:
        # Resume always uses the original generator and model code, even after git pull.
        return subprocess.call([sys.executable, '-u', str(source/'RunFourTargetLP.py'), *argv,
                                '--input-dir', str(args.input_dir), '--output-dir', str(args.output_dir)])
    info = environment()
    with (args.output_dir/'.runner.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('Another four-target LP launcher is using this output folder') from exc
        prepared, checks = prepare(args, saved, references)
        write_json(args.output_dir/'reference_checks.json', dict(integer_files=files, parts=checks))
        write_json(args.output_dir/'environment.json', info)
        status = 0
        try:
            for index, (batch, inputs, _) in enumerate(prepared, 1):
                output = args.output_dir/'parts'/inputs.name
                command = [sys.executable, '-u', str(source/'RunStaticLP.py'), '--input', str(inputs),
                           '-f', str(output), '--source-dir', str(source), '--workers', str(args.workers),
                           '--threads', str(args.threads), '--time-limit', str(args.time_limit),
                           '--retry-time-limit', str(args.retry_time_limit), '--resume']
                print(f"Part {index}/{len(prepared)}: {batch['name']}", flush=True)
                code = subprocess.call(command)
                status = max(status, int(code != 0))
                aggregate(args, prepared)
        finally:
            missing = aggregate(args, prepared)
        print(f'Saved: {args.output_dir/"lp_results.csv"}; {missing} LPs still missing.', flush=True)
        if missing:
            print('Repeat the same command to retry missing LPs and resume saved OPTIMAL values.', flush=True)
        return max(status, int(missing != 0))


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parse_args(argv)
    try:
        return run(args, argv)
    except (ValueError, KeyError, OSError) as exc:
        print(f'Four-target LP error: {exc}', file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print('Stopped. Repeat the same command to resume saved LPs.', file=sys.stderr)
        return 130


if __name__ == '__main__':
    raise SystemExit(main())
