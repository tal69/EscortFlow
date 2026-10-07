"""Run the static load-flow MILP with two-phase lexicographic objectives."""

from pathlib import Path
import runpy
import sys


def main():
    original_argv = sys.argv[:]
    try:
        # Prepend defaults so any explicit CSV argument supplied by the user wins.
        sys.argv = [
            original_argv[0],
            "--lexicographic",
            "--csv", "res_load_flow_lex.csv",
            *original_argv[1:],
        ]
        runpy.run_path(str(Path(__file__).with_name("LoadFlowStatic.py")), run_name="__main__")
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    main()
