"""LP bounds stored alongside their exact integer experiment instance."""

import csv
import math
from pathlib import Path
import time

import RunStaticLP as replay

LP_FIELDS = (
    'lp_requested', 'lp_relaxation_lower_bound', 'lp_flow', 'lp_movements',
    'lp_status', 'lp_elapsed_seconds', 'lp_threads', 'lp_time_limit',
    'lp_retry_time_limit', 'lp_problem_sha256', 'lp_source_sha256',
    'lp_solver_version', 'lp_algorithm', 'lp_error',
)
VALUE_FIELDS = ('lp_relaxation_lower_bound', 'lp_flow', 'lp_movements')


def source_fingerprint(source_dir):
    return replay.fingerprint({name: replay.sha256(Path(source_dir)/name) for name in replay.MODEL_FILES})


def result_fields(result):
    return dict(lp_requested=1, lp_relaxation_lower_bound=result['lp_objective'],
                lp_flow=result['lp_flow'], lp_movements=result['lp_movements'],
                lp_status=result['status'], lp_elapsed_seconds=result['elapsed_seconds'],
                lp_problem_sha256=result['problem_sha256'], lp_source_sha256=result['source_sha256'],
                lp_solver_version=result['solver_version'], lp_algorithm=result['algorithm'], lp_error='')


def embedded_value(row, expected_source=None):
    """Accept only an optimal, consistent LP bound for this exact integer row."""
    if str(row.get('lp_requested', '')) not in {'', '0', '1'}:
        raise ValueError('Invalid lp_requested flag')
    status = row.get('lp_status', '')
    if status != 'OPTIMAL':
        if any(row.get(key) not in ('', None) for key in VALUE_FIELDS):
            raise ValueError('Nonoptimal LP cannot supply a relaxation lower bound')
        if str(row.get('lp_requested', '')) == '1' and status in {'', 'NOT_RUN'}:
            raise ValueError('Requested LP has no status')
        return None
    if str(row.get('lp_requested', '')) != '1' or row.get('lp_error'):
        raise ValueError('Optimal LP has inconsistent request/error fields')
    problem = replay.parse_instance(row)
    if row.get('lp_problem_sha256') != problem['problem_sha256']:
        raise ValueError('Embedded LP has different coordinates, R, or horizon')
    if (not row.get('lp_source_sha256') or
            (expected_source is not None and row['lp_source_sha256'] != expected_source)):
        raise ValueError('Embedded LP model sources differ')
    result = {key: problem[key] for key in ('layout', 'escorts', 'loads', 'seed', 'method', 'flow_weight',
                'horizon', 'physical_horizon', 'retrieval_mode', 'problem_sha256')}
    result.update(lp_objective=row['lp_relaxation_lower_bound'], lp_flow=row['lp_flow'],
                  lp_movements=row['lp_movements'], status=status, elapsed_seconds=row['lp_elapsed_seconds'],
                  source_sha256=row['lp_source_sha256'], solver_version=row['lp_solver_version'],
                  algorithm=row['lp_algorithm'])
    replay.validate_value(result)
    elapsed = float(result['elapsed_seconds'])
    if not math.isfinite(elapsed) or elapsed < 0 or not result['solver_version'] or not result['algorithm']:
        raise ValueError('Invalid embedded LP timing or solver metadata')
    for prefix in ('greedy_', '', 'final_'):
        if prefix != 'greedy_' and str(row.get(prefix+'has_solution')) != '1':
            continue
        if row.get(prefix+'flowtime') in ('', None) or row.get(prefix+'movements') in ('', None):
            continue
        upper = float(row[prefix+'flowtime']) + float(row[prefix+'movements'])/problem['flow_weight']
        if float(result['lp_objective']) > upper + max(1e-6, abs(upper)*1e-8):
            raise ValueError('LP lower bound exceeds a known feasible integer objective')
    return result


def attach_lp(row, args, source_dir):
    """Run a separate continuous solve, without altering integer measurements."""
    start = time.perf_counter()
    row.update(lp_requested=1, lp_status='PENDING', lp_threads=args.lp_threads,
               lp_time_limit=args.lp_time_limit, lp_retry_time_limit=args.lp_retry_time_limit)
    try:
        problem = replay.parse_instance(row)
        source = source_fingerprint(source_dir)
        row.update(lp_problem_sha256=problem['problem_sha256'], lp_source_sha256=source)
        result = replay.solve_instance(problem, dict(threads=args.lp_threads, time_limit=args.lp_time_limit,
                      retry_time_limit=args.lp_retry_time_limit, source_sha256=source))
        row.update(result_fields(result))
        embedded_value(row, source)
    except Exception as exc:
        row.update({key: '' for key in VALUE_FIELDS})
        row.update(lp_status=getattr(exc, 'status', 'ERROR'), lp_elapsed_seconds=time.perf_counter()-start,
                   lp_error=f'{type(exc).__name__}: {exc}')
    print(f"seed={row['seed']} escorts={row['# Escorts']}: LP {row['lp_status']}, "
          f"lower bound={row['lp_relaxation_lower_bound']}, elapsed={row['lp_elapsed_seconds']:.3f}s", flush=True)


def read_embedded(paths, problems, source):
    values = {}
    for path in paths:
        with Path(path).open(newline='') as handle:
            for row in csv.DictReader(handle):
                if row.get('lp_status') != 'OPTIMAL':
                    continue
                problem = replay.parse_instance(row)
                key = problem['problem_sha256']
                if key not in problems:
                    continue
                if key in values:
                    raise ValueError('Duplicate embedded LP result')
                values[key] = embedded_value(row, source)
    return values


def fill_inputs(paths, completed, source, settings):
    """Complete new integrated CSVs after retry, preserving every integer field."""
    updated = 0
    for path in paths:
        path = Path(path)
        with path.open(newline='') as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            rows = list(reader)
        # Historical CSVs are read-only, even when this option is enabled.
        if not set(LP_FIELDS) <= set(fields or []):
            continue
        changed = False
        for row in rows:
            if str(row['lp_requested']) != '1':
                continue
            current = embedded_value(row, source)
            if current is not None:
                continue
            problem = replay.parse_instance(row)
            result = completed.get(problem['problem_sha256'])
            if result is None:
                continue
            if result['source_sha256'] != source:
                raise ValueError('Cannot fill an LP bound from different model sources')
            row.update(result_fields(result), lp_threads=settings['threads'], lp_time_limit=settings['time_limit'],
                       lp_retry_time_limit=settings['retry_time_limit'])
            embedded_value(row, source)
            changed = True
            updated += 1
        if changed:
            temporary = Path(str(path)+'.lp.tmp')
            with temporary.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            temporary.replace(path)
    return updated
