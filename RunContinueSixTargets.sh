#!/usr/bin/env bash
# Revised six-target campaign. The Python launcher freezes the audited v6 code.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec "${PYTHON:-python3}" -u "$script_dir/RunContinueSixTargets.py" "$@"
