#!/usr/bin/env bash
# Remember a working interpreter on this machine; no activation or PATH dependence.
set -euo pipefail
task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
task_config_dir="${XDG_CONFIG_HOME:-$HOME/.config}/escortflow"
task_config_file="$task_config_dir/table2b-python"

save_python() {
    local task_candidate="$1"
    "$task_candidate" -u "$task_dir/RunTable2bLP.py" --check-environment
    mkdir -p "$task_config_dir"
    "$task_candidate" -c 'import sys; print(sys.executable)' > "$task_config_file.tmp"
    mv "$task_config_file.tmp" "$task_config_file"
    printf 'Saved Python for future terminals and tmux sessions: '
    cat "$task_config_file"
    printf 'Start or resume with: bash "%s/RunTable2bLP.sh"\n' "$task_dir"
}

case "${1:-}" in
    --setup)
        if (( $# > 2 )); then
            printf 'Usage: bash RunTable2bLP.sh --setup [python3.13-or-absolute-path]\n' >&2
            exit 2
        fi
        task_base="$(command -v "${2:-python3.13}")" || {
            printf 'Cannot find %s. Pass the absolute path to your Python 3.13 installation.\n' "${2:-python3.13}" >&2
            exit 1
        }
        "$task_base" -c 'import sys; print("Setup Python:", sys.executable, sys.version.split()[0]); sys.exit(0 if (3,10) <= sys.version_info[:2] <= (3,14) else "Setup requires Python 3.10-3.14.")'
        task_environment="$HOME/.venvs/escortflow-table2b"
        if [[ ! -x "$task_environment/bin/python" ]]; then
            "$task_base" -m venv "$task_environment"
        fi
        "$task_environment/bin/python" -m pip install gurobipy==13.0.3
        save_python "$task_environment/bin/python"
        exit 0
        ;;
    --set-python)
        if (( $# != 2 )); then
            printf 'Usage: bash RunTable2bLP.sh --set-python /path/to/licensed/python\n' >&2
            exit 2
        fi
        task_candidate="$(command -v "$2")" || {
            printf 'Cannot find Python: %s\n' "$2" >&2
            exit 1
        }
        save_python "$task_candidate"
        exit 0
        ;;
esac

if [[ ! -f "$task_config_file" ]]; then
    printf 'Set up once: bash "%s/RunTable2bLP.sh" --setup python3.13\n' "$task_dir" >&2
    printf 'Or reuse a licensed environment: bash "%s/RunTable2bLP.sh" --set-python /path/to/python\n' "$task_dir" >&2
    exit 1
fi
IFS= read -r task_python < "$task_config_file"
if [[ ! -x "$task_python" ]]; then
    printf 'Saved Python is missing: %s. Run --setup or --set-python again.\n' "$task_python" >&2
    exit 1
fi
exec "$task_python" -u "$task_dir/RunTable2bLP.py" "$@"
