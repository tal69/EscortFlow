"""Read-only SSH inspection of the enrolled continue_v6 campaign."""

import argparse
import json
from pathlib import Path
import subprocess

REMOTE = r'''python3 - <<'REMOTE'
import csv, io, json, os, time
from pathlib import Path
root = Path.home() / 'EscortFlow/results_continue_v6'
env = json.loads((root / 'environment.json').read_text())
now = time.time()
processes = []
for path in Path('/proc').iterdir():
    if not path.name.isdigit() or int(path.name) == os.getpid():
        continue
    try:
        argv = [arg.decode(errors='replace') for arg in (path / 'cmdline').read_bytes().split(bytes([0])) if arg]
        cwd = (path / 'cwd').resolve()
        driver = (any(Path(arg).name in {'RunContinue.sh', 'RunContinueSixTargets.sh', 'RunContinueSixTargets.py'} for arg in argv)
                  and any(str((cwd / arg).resolve()) == str(root) for arg in argv if 'results_continue_v6' in arg))
        worker = any(arg == str(root / 'source' / name) for name in ('SolveStatic.py', 'RunStaticLP.py') for arg in argv)
        if not (driver or worker):
            continue
        raw = (path / 'stat').read_text()
        fields = raw[raw.rfind(')') + 2:].split()
        processes.append(dict(pid=int(path.name), ppid=int(fields[1]), state=fields[0],
                              kind='driver' if driver else 'worker', argv=argv, cwd=str(cwd)))
    except (OSError, PermissionError):
        continue
coverage, errors = [], []
seeds = set(env['seed_values'])
for targets in env['targets']:
    for layout in env['layouts']:
        for escorts in env['escorts']:
            for method in env['method_order']:
                name = f'{method}_{layout}_t{targets}_e{escorts}.csv'
                path = root / 'parts' / name
                rows, metadata = [], None
                if path.exists():
                    try:
                        stat = path.stat()
                        metadata = dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
                        rows = list(csv.DictReader(io.StringIO(path.read_text())))
                        if any(None in row or any(value is None for value in row.values()) for row in rows):
                            raise ValueError('Malformed or truncated CSV record')
                    except Exception as exc:
                        errors.append(dict(file=name, type='csv_parse', error=str(exc)))
                row_seeds = []
                good, optimal_lp = 0, 0
                for row in rows:
                    try:
                        row_seeds.append(int(row['seed']))
                    except (ValueError, KeyError, TypeError):
                        errors.append(dict(file=name, type='seed_parse', seed=row.get('seed')))
                    if row.get('error'):
                        errors.append(dict(file=name, type='integer', seed=row.get('seed'), error=row['error']))
                    elif row.get('has_solution') == '1' and row.get('final_has_solution') == '1':
                        good += 1
                    if row.get('lp_error'):
                        errors.append(dict(file=name, type='lp', seed=row.get('seed'), error=row['lp_error']))
                    if row.get('lp_status') == 'OPTIMAL' and not row.get('lp_error'):
                        optimal_lp += 1
                    if row.get('protocol') != env['protocol'] or row.get('retrieval_mode') != 'continue' or row.get('movement_mode') != 'BM':
                        errors.append(dict(file=name, type='protocol', seed=row.get('seed')))
                if len(row_seeds) != len(set(row_seeds)) or not set(row_seeds).issubset(seeds):
                    errors.append(dict(file=name, type='seed_coverage'))
                coverage.append(dict(file=name, targets=targets, layout=layout, escorts=escorts, method=method,
                                     records=len(rows), successful_records=good, optimal_lp_records=optimal_lp,
                                     integer_complete=good == len(seeds) and set(row_seeds) == seeds and len(row_seeds) == len(set(row_seeds)),
                                     metadata=metadata))
logs = sorted((root / 'logs').glob('*.log'), key=lambda path: path.stat().st_mtime_ns)
log_errors = []
for path in logs:
    lines = path.read_text(errors='replace').splitlines()[-60:]
    hits = [line for line in lines if any(term in line for term in
            ('Traceback (most recent call last)', 'MemoryError', 'Segmentation fault', 'Killed', 'GurobiError'))]
    if hits:
        log_errors.append(dict(file=path.name, lines=hits))
latest = None
if logs:
    path = logs[-1]
    stat = path.stat()
    latest = dict(file=path.name, size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                  age_seconds=now-stat.st_mtime, tail=path.read_text(errors='replace').splitlines()[-15:])
print(json.dumps(dict(checked_at_epoch=now, host=os.uname().nodename, processes=processes,
                     expected_files=len(coverage), expected_integer_records=env['integer_runs'],
                     present_files=sum(row['metadata'] is not None for row in coverage),
                     integer_records=sum(row['records'] for row in coverage),
                     successful_integer_records=sum(row['successful_records'] for row in coverage),
                     optimal_lp_records=sum(row['optimal_lp_records'] for row in coverage),
                     complete_integer_files=sum(row['integer_complete'] for row in coverage),
                     coverage=coverage, errors=errors, log_errors=log_errors, latest_log=latest), indent=2))
REMOTE'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--campaign', choices=('continue_v6', 'continue_v6_six_targets'), default='continue_v6')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    helper = Path.home() / '.codex/skills/lab-linux-access/scripts/lab_access.py'
    remote = REMOTE.replace('EscortFlow/results_continue_v6', 'EscortFlow/results_' + args.campaign)
    result = subprocess.run(['python3', str(helper), 'exec', 'bukhan1', '--command', remote],
                            capture_output=True, text=True, timeout=70)
    (args.output_dir / 'remote_status.stderr').write_text(result.stderr)
    (args.output_dir / 'remote_status.json').write_text(result.stdout)
    if result.returncode:
        print(result.stderr)
        print(result.stdout)
        raise SystemExit(result.returncode)
    status = json.loads(result.stdout)
    print(json.dumps({key: value for key, value in status.items() if key != 'coverage'}, indent=2))


if __name__ == '__main__':
    main()
