"""Run static escort-flow experiments with two-phase lexicographic objectives."""

from pathlib import Path
import runpy
import sys


if __name__ == "__main__":
    # argparse keeps the last occurrence, so a caller's CSV argument wins.
    sys.argv[1:1] = ["--lexicographic", "--csv", "res_escort_flow_lex.csv"]
    runpy.run_path(str(Path(__file__).with_name("EscortFlowStatic.py")), run_name="__main__")
