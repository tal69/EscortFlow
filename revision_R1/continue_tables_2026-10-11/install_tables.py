"""Validate and install the prepared tables without invoking a LaTeX compiler."""

import argparse
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path)
args = parser.parse_args()
if args.output_dir:
    HERE = args.output_dir.resolve()
audit = json.loads((HERE / 'audit.json').read_text())
main = Path('/Users/talraviv/Library/CloudStorage/Dropbox/.metadata/Apps/Overleaf/'
            'Escort flow formulation for the PBS retrieval problem/static_escort_flow_IISE.tex')
original = main.read_text()
candidate = (HERE / 'candidate_main.tex').read_text()
assert hashlib.sha256(main.read_bytes()).hexdigest() == audit['original_main_sha256'], 'Manuscript changed during preparation'
start = candidate.index('% BEGIN CONTINUE FOUR-TARGET TABLES')
end_marker = '% END CONTINUE FOUR-TARGET TABLES\n\n'
end = candidate.index(end_marker, start) + len(end_marker)
if '% BEGIN CONTINUE FOUR-TARGET TABLES' in original:
    old_start = original.index('% BEGIN CONTINUE FOUR-TARGET TABLES')
    old_end = original.index(end_marker, old_start) + len(end_marker)
    original_outside = original[:old_start] + original[old_end:]
else:
    original_outside = original
assert candidate[:start] + candidate[end:] == original_outside, 'Unrelated manuscript edit'
block = candidate[start:end]
slash = chr(92)
assert block.count(slash + 'begin{table}') == block.count(slash + 'end{table}') == 2
assert block.count(slash + 'begin{tabular*}') == block.count(slash + 'end{tabular*}') == 2
assert slash + 'newcommand{' + slash + 'blind}{0}' in candidate
labels = re.findall(re.escape(slash + 'label{') + r'([^}]+)}', candidate)
for label in ('tab:continue-four-target-comparison', 'tab:continue-four-target-bounds'):
    assert labels.count(label) == 1
    assert candidate.count(slash + 'ref{' + label + '}') >= 1
rows_by_table = []
for kind, columns in ((2, 17), (3, 11)):
    table = (HERE / f'table{kind}_continue_targets.tex').read_text()
    rows = [line for line in table.splitlines() if line.startswith(('2 & ', '4 & ', '6 & '))]
    assert len(rows) == audit['displayed_row_count']
    for row in rows:
        assert row.endswith(' ' + slash * 2)
        cells = row[:-3].split(' & ')
        assert len(cells) == columns, (kind, row, len(cells))
        layout = cells[1].split(' (')[0].replace('$' + slash + 'times$', 'x')
        if [int(cells[0]), layout, int(cells[2])] not in audit['complete_configurations']:
            assert all(not cell.strip() for cell in cells[3:]), 'Placeholder result was populated'
    complete_count = len(audit['complete_configurations'])
    expected_order = [(n, l, e) for n in audit['displayed_target_counts']
                      for l in ('13x7', '10x10', '16x10', '27x10')
                      for e in audit.get('escorts_by_target', {}).get(str(n), (8, 12, 16))]
    actual_order = [(int(c[0]), c[1].split(' (')[0].replace('$' + slash + 'times$', 'x'), int(c[2]))
                    for c in (row[:-3].split(' & ') for row in rows)]
    assert actual_order == expected_order, 'Incorrect target/layout/escort order'
    rows_by_table.append(dict(table_kind=kind, rows=len(rows), populated_rows=complete_count,
                              blank_rows=len(rows)-complete_count, columns=columns))
refs = re.findall(re.escape(slash + 'ref{') + r'([^}]+)}', block)
assert all(ref in labels for ref in refs)
uncommented = '\n'.join(line[:next((i for i, char in enumerate(line)
                                 if char == '%' and (i == 0 or line[i-1] != slash)), len(line))]
                        for line in block.splitlines())
depth = 0
for char in uncommented.replace(slash + '{', '').replace(slash + '}', ''):
    if char == '{':
        depth += 1
    elif char == '}':
        depth -= 1
        assert depth >= 0
assert depth == 0
stack = []
for action, env in re.findall(re.escape(slash) + r'(begin|end){([^}]+)}', uncommented):
    if action == 'begin':
        stack.append(env)
    else:
        assert stack.pop() == env
assert not stack
for item in audit['input_files']:
    assert hashlib.sha256((HERE / item['file']).read_bytes()).hexdigest() == item['sha256']
if candidate != original:
    main.write_text(candidate)
assert main.read_text() == candidate
validation = dict(passed=True, canonical_main=str(main), existing_content_preserved=True,
                  new_tables=rows_by_table, author_mode=True, new_labels_unique=True,
                  references_resolved=True, braces_and_environments_balanced=True,
                  input_snapshot_hashes_verified=True, compilation_performed=False,
                  final_main_sha256=hashlib.sha256(main.read_bytes()).hexdigest())
(HERE / 'source_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
print(json.dumps(validation, indent=2))
