"""Direct paper LPs from the regular runners' dimensions, ranges, and seeds.

Generate coordinates and the greedy reference without solving an integer model.
Freeze these inputs and model sources, then use the existing continuous replay
engine. The v4 option reproduces the archived coefficient/horizon conventions.
"""

import csv
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

from PBSCom import GeneretaeRandomInstance, str2range
import RunStaticLP as lp

HERE = Path(__file__).resolve().parent
GENERATOR_FILES = ('static_generated_lp.py', 'PBSCom.py', 'OneStepHeuristic_v2.py',
                   'static_weighted_certification.py')
SNAPSHOT_FILES = tuple(dict.fromkeys((*lp.MODEL_FILES, *GENERATOR_FILES, 'RunStaticLP.py', 'static_integrated_lp.py')))


def generated_row(method, lx, ly, outputs, escorts, loads, seed, mode, protocol,
                  flow_weight=None, horizon=None):
    from OneStepHeuristic_v2 import SolveGreedy
    from static_weighted_certification import safe_weight_parameters
    targets, empty = GeneretaeRandomInstance(seed, sorted(itertools.product(range(lx), range(ly))), escorts, loads)
    trace = SolveGreedy(lx, ly, set(outputs), set(targets), set(empty), verbal=False,
                       max_steps=max(1, 4 * loads * (lx + ly - 2) + 1), retrieval_mode=mode)
    makespan, flow, movements = trace[:3]
    parameters = safe_weight_parameters(targets, outputs, lx * ly, escorts, flow)
    coefficient = (parameters['movement_bound'] + 1 if protocol == 'v4'
                   else parameters['flow_weight']) if flow_weight is None else flow_weight
    if coefficient < parameters['flow_weight']:
        raise ValueError('R is below the sufficient integer coefficient; use --legacy-lp for a historical weighted objective')
    safe_horizon = parameters['flow_horizon']
    physical = (safe_horizon + 1 if protocol == 'v4' and method == 'escortflow'
                else max(safe_horizon, makespan + 1))
    if horizon is None:
        horizon = physical - (method == 'escortflow')
    else:
        physical = horizon + (method == 'escortflow')
        if physical < safe_horizon:
            raise ValueError('Explicit horizon does not cover the sufficient physical horizon')
    return {'Lx x Ly': f'{lx}x{ly}', 'formulation': method, 'seed': seed,
            '#Loads': loads, '# Escorts': escorts, 'retrieval_mode': mode, 'movement_mode': 'BM',
            'IOs': repr(outputs), 'Target Loads': repr(targets), 'Escorts': repr(empty),
            'flow_weight': coefficient, 'movement_weight': 1/coefficient,
            'weighted_horizon': horizon, 'weighted_physical_horizon': physical,
            'protocol': 'generated_lp_' + protocol, 'safe_flow_horizon': safe_horizon,
            'greedy_makespan': makespan, 'greedy_flowtime': flow, 'greedy_movements': movements}


def run_standard(args, method):
    """Accept either regular EF/LF arguments or RunSafeWeightedStatic arguments."""
    raw_outputs = getattr(args, 'output_cells', getattr(args, 'outputs', None))
    if raw_outputs and isinstance(raw_outputs[0], tuple):
        outputs = sorted(raw_outputs)
    else:
        if not raw_outputs or len(raw_outputs) % 2:
            raise ValueError('Output coordinates must be x y pairs')
        outputs = sorted(zip(raw_outputs[::2], raw_outputs[1::2]))
    lx, ly = args.Lx, args.Ly
    loads = getattr(args, 'load_num', getattr(args, 'loads', 1))
    mode = getattr(args, 'retrieval_mode', 'leave')
    seeds = sorted(str2range(getattr(args, 'reps_range', getattr(args, 'seeds', '1-100'))))
    escorts = sorted(str2range(getattr(args, 'escorts_range', getattr(args, 'escorts', '3-8'))))
    protocol = getattr(args, 'lp_protocol', 'v5')
    if getattr(args, 'flow_weight', None) is not None and args.flow_weight <= 0:
        raise ValueError('R must be a positive integer')
    if getattr(args, 'horizon', None) is not None and args.horizon < 0:
        raise ValueError('The horizon must be nonnegative')
    output = Path(getattr(args, 'csv', getattr(args, 'output', 'lp.csv'))).resolve()
    resume = getattr(args, 'resume', False)
    if method not in {'escortflow', 'loadflow'} or protocol not in {'v4', 'v5'}:
        raise ValueError('Expected an EF/LF formulation and v4/v5 LP protocol')
    if (min(lx, ly, loads) <= 0 or not seeds or min(seeds) < 0 or not escorts or min(escorts) < 1
            or max(escorts) + loads > lx * ly):
        raise ValueError('Invalid grid, target/escort count, or seed range')
    if len(seeds) != len(set(seeds)) or len(escorts) != len(set(escorts)):
        raise ValueError('Seed and escort ranges cannot contain duplicates')
    if len(outputs) != len(set(outputs)) or any(not (0 <= x < lx and 0 <= y < ly) for x, y in outputs):
        raise ValueError('Outputs must be distinct cells inside the grid')
    if mode not in {'leave', 'continue'} or getattr(args, 'lm', False) or getattr(args, 'opl', False):
        raise ValueError('Direct paper LPs require Gurobi, BM, and leave/continue retrieval; use --legacy-lp for other legacy modes')
    for flag in ('warmstart', 'cutoff', 'lexicographic', 'greedy', 'naive', 'dp_file'):
        if getattr(args, flag, None):
            raise ValueError(f'--lp cannot combine the direct paper relaxation with {flag}')
    for flag in ('work_limit', 'phase1_time_limit'):
        if getattr(args, flag, None) is not None:
            raise ValueError(f'--lp cannot combine the direct paper relaxation with {flag}')
    if getattr(args, 'lazy', None) is not None or getattr(args, 'bnc', None) is not None:
        raise ValueError('Direct LPs cannot use lazy constraints or branch-and-cut')
    if (getattr(args, 'alpha', 0), getattr(args, 'beta', 1), getattr(args, 'gamma', .01)) != (0, 1, .01):
        raise ValueError('Direct --lp calculates R for F+M/R; use --legacy-lp to choose historical objective weights')
    threads = getattr(args, 'num_threads', getattr(args, 'threads', 1))
    if threads < 0:
        raise ValueError('Threads cannot be negative')
    # RunStaticLP requires a positive count; zero means automatic in the old CLI.
    if threads == 0:
        from RunStaticCampaign import automatic_threads
        threads = automatic_threads()
    limit = getattr(args, 'time_limit', getattr(args, 'weighted_time_limit', 300))
    if limit is None:
        limit = 300
    workers = getattr(args, 'lp_workers', 1)
    retry = getattr(args, 'lp_retry_time_limit', 600)
    if workers <= 0 or not math.isfinite(limit) or limit <= 0 or not math.isfinite(retry) or retry <= 0:
        raise ValueError('LP worker counts and time limits must be positive')
    if output.exists() and not resume:
        raise ValueError('LP output already exists; use a new -f filename or --resume')
    settings = dict(formulation=method, lx=lx, ly=ly, outputs=outputs, loads=loads, seeds=seeds,
                    escorts=escorts, mode=mode, protocol=protocol,
                    flow_weight=getattr(args, 'flow_weight', None), horizon=getattr(args, 'horizon', None))
    generator_hashes = {name: lp.sha256(HERE/name) for name in GENERATOR_FILES}
    metadata = json.loads(json.dumps(dict(settings=settings, generator_sha256=generator_hashes)))
    inputs = Path(str(output) + '.instances.csv')
    generation = Path(str(output) + '.generation.json')
    source = Path(str(output) + '.source')
    if generation.exists():
        saved = json.loads(generation.read_text())
        if saved['settings'] != metadata['settings']:
            raise ValueError('Generated LP settings changed; use a new output filename')
        if not inputs.is_file() or not source.is_dir():
            raise ValueError('Generated input/source snapshot is missing')
        if lp.sha256(inputs) != saved['inputs_sha256']:
            raise ValueError('Generated input snapshot changed')
        for name, expected in saved['source_sha256'].items():
            if lp.sha256(source/name) != expected:
                raise ValueError('Frozen LP source changed: ' + name)
    else:
        if inputs.exists() or source.exists():
            raise ValueError('An unrecognized input/source snapshot already exists; use a new output filename')
        rows = []
        for count in escorts:
            for seed in seeds:
                rows.append(generated_row(method, lx, ly, outputs, count, loads, seed, mode, protocol,
                                          settings['flow_weight'], settings['horizon']))
                if len(rows) % 25 == 0:
                    print(f'Prepared {len(rows)}/{len(escorts)*len(seeds)} direct LP instances.', flush=True)
        if any(lp.sha256(HERE/name) != expected for name, expected in generator_hashes.items()):
            raise ValueError('Generation code changed while preparing instances; restart with stable code')
        output.parent.mkdir(parents=True, exist_ok=True)
        source.mkdir()
        for name in SNAPSHOT_FILES:
            shutil.copy2(HERE/name, source/name)
        if any(lp.sha256(source/name) != expected for name, expected in generator_hashes.items()):
            shutil.rmtree(source)
            raise ValueError('Generation code changed while freezing sources; restart with stable code')
        with inputs.open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
        metadata.update(inputs_sha256=lp.sha256(inputs), source_sha256={name: lp.sha256(source/name) for name in SNAPSHOT_FILES})
        generation.write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'Direct {method} LPs: {len(seeds)*len(escorts)} instances generated from parameters; protocol {protocol}.', flush=True)
    print(f'Instance snapshot: {inputs}', flush=True)
    arguments = ['--input', str(inputs), '-f', str(output), '--source-dir', str(source),
                 '--workers', str(workers), '--threads', str(threads),
                 '--time-limit', str(limit), '--retry-time-limit', str(retry)]
    if resume:
        arguments.append('--resume')
    # The historical EF/LF entry points parse arguments at module scope. Start
    # the guarded replay entry point separately so macOS spawn workers cannot
    # import those entry points and repeat their generation/output operations.
    return subprocess.call([sys.executable, '-u', str(source/'RunStaticLP.py'), *arguments])
