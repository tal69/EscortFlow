"""Build multi-target continue tables from the saved v6 CSV snapshot. No compilation."""

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import sys
from zoneinfo import ZoneInfo

HERE = Path(__file__).resolve().parent
CODE = HERE.parent.parent
sys.path.insert(0, str(CODE))
import PaperTables as tables
from static_integrated_lp import embedded_value, source_fingerprint

MAIN = Path('/Users/talraviv/Library/CloudStorage/Dropbox/.metadata/Apps/Overleaf/'
            'Escort flow formulation for the PBS retrieval problem/static_escort_flow_IISE.tex')
RESULTS = CODE / 'results_continue_v6'
SIX_RESULTS = CODE / 'results_continue_v6_six_targets'
LAYOUTS = ('13x7', '10x10', '16x10', '27x10')
ESCORTS = (8, 12, 16)
ESCORTS_BY_TARGET = {2: ESCORTS, 4: ESCORTS, 6: (12, 16, 20)}
TARGET_COUNTS = (2, 4, 6)
OUTPUT_COUNTS = dict(zip(LAYOUTS, (1, 1, 2, 3)))
SHARED = ('flow_weight', 'safe_flow_horizon', 'weighted_physical_horizon', 'threads',
          'warmstart', 'greedy_flowtime', 'greedy_movements', 'greedy_makespan',
          'weighted_time_limit')


def campaign_root(loads):
    return SIX_RESULTS if loads == 6 else RESULTS


def extension_budget(loads, method):
    return 0 if loads == 6 and method == 'loadflow' else 300


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def panel(source, label):
    position = source.index(r'\label{' + label + '}')
    table_end = source.index(r'\end{table}', position)
    start = source.find(r'\scriptsize', position, table_end)
    if start < 0:
        # Standalone-table labels follow the tabular in the caption, unlike
        # subfloat labels, which precede their tabular.
        table_start = source.rfind(r'\begin{table}', 0, position)
        start = source.rfind(r'\scriptsize', table_start, position)
    assert start >= 0, ('Table formatting not found', label)
    stop = source.index(r'\end{tabular*}', start) + len(r'\end{tabular*}')
    return source[start:stop]


def table_rows(fragment):
    """Read rows with or without the newly added target-count column."""
    result = {}
    for line in fragment.splitlines():
        cells = line.removesuffix(r' \\').split(' & ')
        offset = int(cells[0] in {'2', '4', '6'})
        if len(cells) <= offset or not cells[offset].startswith(('13$', '10$', '16$', '27$')):
            continue
        loads = int(cells[0]) if offset else 4
        layout = cells[offset].split(' (')[0].replace(r'$\times$', 'x')
        key = (loads, layout, int(cells[offset + 1]))
        assert key not in result, ('Duplicate displayed configuration', key)
        result[key] = cells[offset + 2:]
    return result


def render_fragment(kind, template, displayed_targets, eligible, comparisons, solutions):
    """Keep the existing style and sort by targets, layout order, and escorts."""
    lines = template.splitlines()
    first_row = next(i for i, line in enumerate(lines) if line.startswith('13$'))
    header = '\n'.join(lines[:first_row])
    if kind == 2:
        header = header.replace('lr|rrrrrrr|rrrrrrr', 'rlr|rrrrrrr|rrrrrrr')
        header = header.replace(r'\multicolumn{2}{c|}{Config.}', r'\multicolumn{3}{c|}{Config.}')
    else:
        header = header.replace('lrr|rrrrrrr', 'rlrr|rrrrrrr')
        header = header.replace(r'\multicolumn{3}{c|}{Configuration}', r'\multicolumn{4}{c|}{Configuration}')
    for title in ('Dim.', 'PBS dim.'):
        header = header.replace(r'\makecell{' + title,
                                r'\makecell{\# of\\targets} & \makecell{' + title, 1)
    assert header.count(r'\makecell{\# of\\targets}') == 1
    output = [header]
    previous = None
    for loads in sorted(displayed_targets):
        for layout in LAYOUTS:
            if previous is not None:
                output.append(r'\midrule')
            previous = (loads, layout)
            for esc in ESCORTS_BY_TARGET[loads]:
                key = (loads, layout, esc)
                cells = [str(loads), layout.replace('x', r'$\times$') + f' ({OUTPUT_COUNTS[layout]})', str(esc)]
                if key not in eligible:
                    cells += [''] * (14 if kind == 2 else 8)
                elif kind == 2:
                    stats = {m: next(r for r in comparisons if (r['loads'], r['layout'], r['escorts'], r['method']) == (*key, m))
                             for m in tables.METHODS}
                    for method in tables.METHODS:
                        row = stats[method]
                        assert row['instances'] == 100
                        cells += [str(row['ft_opt_count']), str(row['mv_opt_count']),
                                  '' if row['ft_opt_pct'] == 100 else f"{row['ft_gap_pct']:.2f}",
                                  '' if row['mv_opt_pct'] == 100 else f"{row['mv_gap_pct']:.2f}", f"{row['lp_gap_pct']:.2f}"]
                        for field in ('ft_cpu_seconds', 'total_cpu_seconds'):
                            value = f'{row[field]:.2f}'
                            if round(row[field], 2) == min(round(r[field], 2) for r in stats.values()):
                                value = r'\textbf{' + value + '}'
                            cells.append(value)
                else:
                    row = next(r for r in solutions if (r['loads'], r['layout'], r['escorts']) == key)
                    cells += [f"{100 * row['utilization']:.1f}"]
                    for field in ('naive_ft', 'greedy_ft', 'greedy_mv', 'best_ft', 'best_mv', 'ft_improvement_pct', 'mv_improvement_pct'):
                        value = f'{row[field]:.2f}'
                        if (field == 'best_ft' and row['flow_unproved']) or (field == 'best_mv' and row['movements_unproved']):
                            value += r'\textsuperscript{*}'
                        cells.append(value)
                assert len(cells) == (17 if kind == 2 else 11)
                output.append(' & '.join(cells) + r' \\')
    return '\n'.join([*output, r'\bottomrule', r'\end{tabular*}'])


def main():
    global HERE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, help='New evidence directory for this refresh')
    args = parser.parse_args()
    if args.output_dir:
        HERE = args.output_dir.resolve()
        HERE.mkdir(parents=True, exist_ok=False)
    today = datetime.now(ZoneInfo('Asia/Jerusalem'))
    snapshot_date = today.strftime('%Y-%m-%d')
    display_date = f'{today:%B} {today.day}, {today.year}'
    snapshot = HERE / 'inputs'
    snapshot.mkdir(exist_ok=True)
    env = json.loads((RESULTS / 'environment.json').read_text())
    assert env['protocol'] == 'safe_integer_flow_timing_v6'
    assert env['retrieval_mode'] == 'continue'
    assert env['seed_values'] == list(range(1, 101))
    shutil.copy2(RESULTS / 'environment.json', snapshot / 'environment.json')
    source_hash = source_fingerprint(RESULTS / 'source')
    for name, expected in env['source_sha256'].items():
        assert sha(RESULTS / 'source' / name) == expected, ('Frozen source changed', name)
    campaign_environments = {str(RESULTS): dict(source_sha256=env['source_sha256'],
                                               model_source_sha256=source_hash)}
    if (SIX_RESULTS / 'environment.json').exists():
        six_env = json.loads((SIX_RESULTS / 'environment.json').read_text())
        assert six_env['protocol'] == env['protocol'] and six_env['retrieval_mode'] == 'continue'
        assert six_env['seed_values'] == list(range(1, 101))
        assert six_env['targets'] == [6] and six_env['escorts'] == [12, 16, 20]
        assert six_env['extension_time_limit_by_method'] == {'loadflow': 0, 'escortflow': 300}
        assert six_env['initial_time_limit'] == 300
        for name, expected in six_env['source_sha256'].items():
            assert sha(SIX_RESULTS / 'source' / name) == expected, ('Revised frozen source changed', name)
        six_hash = source_fingerprint(SIX_RESULTS / 'source')
        assert six_hash == source_hash, 'Revised campaign must retain the audited numerical model sources'
        shutil.copy2(SIX_RESULTS / 'environment.json', snapshot / 'environment_six_targets.json')
        campaign_environments[str(SIX_RESULTS)] = dict(source_sha256=six_env['source_sha256'],
                                                      model_source_sha256=six_hash)
    paired, inputs, eligible, coverage = {}, [], [], []
    for loads, layout, esc in ((n, l, e) for n in TARGET_COUNTS for l in LAYOUTS for e in ESCORTS_BY_TARGET[n]):
        complete = True
        for method in tables.METHODS:
            path = campaign_root(loads) / 'parts' / f'{method}_{layout}_t{loads}_e{esc}.csv'
            rows = tables.read_csv(path) if path.exists() else []
            if rows:
                assert str(campaign_root(loads)) in campaign_environments, 'Results lack a verified campaign environment'
            unique_seeds = {int(r['seed']) for r in rows}
            assert len(rows) == len(unique_seeds), ('Duplicate seeds', path)
            valid = (len(rows) == 100 and unique_seeds == set(range(1, 101))
                     and all(not r['error'] and r['lp_status'] == 'OPTIMAL' and not r['lp_error'] for r in rows))
            coverage.append(dict(loads=loads, layout=layout, escorts=esc, method=method, integer_records=len(rows),
                                 optimal_lp_records=sum(r['lp_status'] == 'OPTIMAL' and not r['lp_error'] for r in rows),
                                 errors=sum(bool(r['error'] or r['lp_error']) for r in rows), complete=valid))
            complete = complete and valid
        if complete:
            eligible.append((loads, layout, esc))
    assert eligible, 'No complete matched configurations'
    for loads, layout, esc in eligible:
            for method in tables.METHODS:
                path = campaign_root(loads) / 'parts' / f'{method}_{layout}_t{loads}_e{esc}.csv'
                saved = snapshot / path.name
                shutil.copy2(path, saved)
                rows = tables.read_csv(saved)
                assert len(rows) == 100 and {int(r['seed']) for r in rows} == set(range(1, 101))
                inputs.append(dict(file=str(saved.relative_to(HERE)), rows=100, sha256=sha(saved)))
                for row in rows:
                    assert row['protocol'] == env['protocol']
                    assert row['retrieval_mode'] == 'continue' and row['movement_mode'] == 'BM'
                    assert row['formulation'] == method and row['Lx x Ly'] == layout
                    assert int(row['#Loads']) == loads and int(row['# Escorts']) == esc
                    assert not row['error'] and row['phase1_snapshot_missing'] == '0'
                    assert all(row[f] == '1' for f in ('warmstart', 'has_solution', 'final_has_solution',
                                                     'weighted_global_scope', 'weighted_bound_consistent',
                                                     'final_bound_consistent'))
                    assert int(row['threads']) == 16
                    assert float(row['weighted_time_limit']) == 300
                    assert float(row['extension_time_limit']) == extension_budget(loads, method)
                    key = tables.instance_key(row)
                    assert method not in paired.setdefault(key, {})
                    paired[key][method] = row
    records, certificates = [], []
    for key, methods in sorted(paired.items()):
        assert set(methods) == set(tables.METHODS)
        lf, ef = (methods[m] for m in tables.METHODS)
        assert tables.coordinates(lf) == tables.coordinates(ef)
        assert all(lf[f] == ef[f] for f in SHARED)
        lx, ly = map(int, lf['Lx x Ly'].split('x'))
        distance = int(lf['safe_movement_lower_bound'])
        bound = (lx * ly - int(lf['# Escorts'])) * int(lf['safe_flow_horizon'])
        assert int(lf['flow_weight']) == bound - distance + 1
        candidates = [pair for row in methods.values() for pair in tables.feasible_pairs(row)]
        best = min(candidates)
        flow_proven = moves_proven = False
        for method, row in methods.items():
            for prefix in ('', 'final_'):
                pair = (tables.integer(row, prefix + 'flowtime'), tables.integer(row, prefix + 'movements'))
                assert tables.integer(row, prefix + 'scaled_objective') == int(row['flow_weight']) * pair[0] + pair[1]
            final = (int(row['final_flowtime']), int(row['final_movements']))
            if row['final_flow_proven'] == '1':
                assert final[0] == best[0]
                flow_proven = True
            if row['final_lexicographic_proven'] == '1':
                assert final == best
                moves_proven = True
            lp = embedded_value(row, source_hash)
            assert lp is not None and row['lp_threads'] == '1'
            records.append(tables.component_metrics(row, best, lp))
        certificates.append(dict(loads=key[2], layout=key[1], escorts=key[3], seed=key[4],
                                 flow_proven=int(flow_proven), movements_proven=int(moves_proven)))
    comparison, solutions = tables.aggregate(records)
    for row in solutions:
        matches = [r for r in certificates if (r['loads'], r['layout'], r['escorts']) == (row['loads'], row['layout'], row['escorts'])]
        assert len(matches) == 100
        row['flow_unproved'] = sum(not r['flow_proven'] for r in matches)
        row['movements_unproved'] = sum(not r['movements_proven'] for r in matches)
    assert len(records) == 200 * len(eligible) and len(paired) == 100 * len(eligible)
    assert len(comparison) == 2 * len(eligible) and len(solutions) == len(eligible)
    for name, rows in (('instance_metrics.csv', records), ('method_comparison.csv', comparison),
                       ('bounds_and_solutions.csv', solutions), ('certificates.csv', certificates)):
        tables.write_csv(HERE / name, rows)

    original = MAIN.read_text()
    assert r'\newcommand{\blind}{0}' in original
    displayed_targets = {key[0] for key in eligible}
    if '% BEGIN CONTINUE FOUR-TARGET TABLES' in original:
        for label in ('tab:continue-four-target-comparison', 'tab:continue-four-target-bounds'):
            for config, cells in table_rows(panel(original, label)).items():
                displayed_targets.add(config[0])
                if any(cell.strip() for cell in cells):
                    assert config in eligible, ('Previously populated row lost complete valid inputs', config)
    if not (HERE / 'before_main.tex').exists():
        (HERE / 'before_main.tex').write_text(original)
    fragments = {}
    for kind, label in ((2, 'tab: 4 target loads'), (3, 'tab:bounds-four-target')):
        fragments[kind] = render_fragment(kind, panel(original, label), displayed_targets, eligible, comparison, solutions)
        (HERE / f'table{kind}_continue_targets.tex').write_text(fragments[kind] + '\n')
    row_count = 12 * len(displayed_targets)
    introduction = (r'''Tables~\ref{tab:continue-four-target-comparison} and~\ref{tab:continue-four-target-bounds} report multi-target \emph{continue} retrieval under the corrected v6 specification. The SNAPSHOT_DATE snapshot contains 100 matched instances in each populated row. Both formulations use identical initial states, complete greedy warm starts, coefficients, physical horizons, and solver settings. The measures and averaging conventions follow Tables~\ref{tab:single-and-multi-target loads}(b) and~\ref{tab:bounds-and-improvement}(b): optimality rates, component lower bounds, and times use the initial 300-second cutoff, while common best-known solution components and certificates include the conditional extensions.''' .replace('SNAPSHOT_DATE', display_date)
                    + (' Incomplete configuration rows are left blank pending completion.' if len(eligible) < row_count else ''))
    introduction = introduction.replace('and solver settings.', 'and initial solver settings.')
    if 6 in displayed_targets:
        introduction += (' Six-target rows use 12, 16, and 20 escorts. Load-flow has no time extension in this campaign,'
                         ' while escort-flow retains a conditional 300-second extension; both initial comparison budgets remain 300 seconds.')
    captions = {
        2: r'''Multi-target \emph{continue}-mode comparison from the October 11, 2026 v6 snapshot. Each populated row uses 100 matched instances. Optimality rates, FT/MV incumbent components and lower bounds, and CPU times use each method's own initial 300-second cutoff. All gaps use $100(1-\sum\mathrm{LB}/\sum\mathrm{UB})$ from unrounded bounds; FT/MV upper bounds are that method's cutoff incumbent components, while LP upper bounds use the common best-known pair including extensions. MV bounds are conditional on that method's cutoff FT. FT/MV gap cells are blank when the corresponding optimality count is 100. Bold marks the lower mean time for each timing measure. The $16\times10$ and $27\times10$ result rows are pending. \label{tab:continue-four-target-comparison}''',
        3: r'''Solution components and improvements over the greedy heuristic for multi-target \emph{continue} retrieval from the October 11, 2026 v6 snapshot. Each populated row averages 100 matched instances. OPT values and improvements use the common best-known feasible pair across both formulations, including extensions. FT/MV improvements are calculated from the unrounded component means. An asterisk marks an OPT mean containing an unproved value for that component. The $16\times10$ and $27\times10$ result rows are pending. \label{tab:continue-four-target-bounds}''',
    }
    for kind in captions:
        captions[kind] = captions[kind].replace('October 11, 2026', display_date)
        captions[kind] = captions[kind].replace(r'The $16\times10$ and $27\times10$ result rows are pending.',
                                               'Blank result rows are pending.' if len(eligible) < row_count else '')
    floats = []
    for kind in (2, 3):
        floats.append('\n'.join([r'\begin{table}[p]', r'    \centering',
                                 r'    \renewcommand{\baselinestretch}{1}',
                                 r'    \captionsetup{font=scriptsize,skip=4pt}', fragments[kind],
                                 r'    \caption{' + captions[kind] + '}', r'\end{table}']))
    block = '% BEGIN CONTINUE FOUR-TARGET TABLES\n' + introduction + '\n\n' + '\n\n'.join(floats) + '\n% END CONTINUE FOUR-TARGET TABLES\n\n'
    start_marker = '% BEGIN CONTINUE FOUR-TARGET TABLES'
    end_marker = '% END CONTINUE FOUR-TARGET TABLES\n\n'
    if start_marker in original:
        start = original.index(start_marker)
        stop = original.index(end_marker, start) + len(end_marker)
        candidate = original[:start] + block + original[stop:]
    else:
        anchor = r'\section{Conclusions} \label{sec: Conclusions}'
        assert original.count(anchor) == 1
        candidate = original.replace(anchor, block + anchor)
    # All existing manuscript content is retained verbatim, including Tables 2-3.
    for label in ('tab:1 target load:feas', 'tab: 4 target loads', 'tab:bounds-single-target', 'tab:bounds-four-target'):
        assert panel(candidate, label) == panel(original, label)
    (HERE / 'candidate_main.tex').write_text(candidate)
    audit = dict(passed=True, snapshot_date=snapshot_date, protocol=env['protocol'],
                 source=str(RESULTS), model_source_sha256=source_hash, input_files=inputs,
                 matched_instances=len(paired), integer_records=len(records), optimal_lp_records=len(records),
                 complete_configurations=eligible, coverage=coverage,
                 escorts_by_target=ESCORTS_BY_TARGET, campaign_environments=campaign_environments,
                 configuration_campaigns={n: str(campaign_root(n)) for n in TARGET_COUNTS},
                 extension_time_limits={n: {m: extension_budget(n, m) for m in tables.METHODS} for n in TARGET_COUNTS},
                 displayed_target_counts=sorted(displayed_targets), displayed_row_count=row_count,
                 seed_values=list(range(1, 101)), existing_panels_preserved=True,
                 flow_unproved=sum(r['flow_unproved'] for r in solutions),
                 movements_unproved=sum(r['movements_unproved'] for r in solutions),
                 cutoff_time_estimates=sum(r['cutoff_cpu_time_is_estimate'] for r in records),
                 timing='Initial-cutoff elapsed wall time including model construction; unproved FT uses cutoff elapsed time',
                 percentages="FT/MV/LP gaps use 100*(1-sum(lower)/sum(upper)) from unrounded bounds; FT/MV use each method's own initial-cutoff incumbent components and bounds, with MV conditional on its cutoff FT; LP uses normalized F+M/R units and the common best-known initial/final pair; FT/MV improvements use the ratio of unrounded component means",
                 best_known='Lexicographically smallest feasible pair from the common greedy, initial, and final records of both methods',
                 compilation_performed=False, original_main_sha256=sha(MAIN))
    (HERE / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({k: v for k, v in audit.items() if k not in ('input_files', 'seed_values')}, indent=2))
    for kind in (2, 3):
        print(f'Generated table {kind} continue rows:')
        print('\n'.join(line for line in fragments[kind].splitlines() if line.startswith(('2 & ', '4 & ', '6 & '))))


if __name__ == '__main__':
    main()
