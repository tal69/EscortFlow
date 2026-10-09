#!/usr/bin/env python3
"""Build a portable Table 3(a) LP reproduction bundle from the frozen study."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "revision_R1" / "table3_bounds_2026-10-09"
MODEL_FILES = (
    "escort_flow_static_gurobi.py", "load_flow_static_gurobi.py",
    "static_lexicographic.py", "static_weighted_certification.py",
    "static_safe_weighted_search.py",
)

PACKAGE_README = """# Table 3(a): continuous LP relaxation reproduction

This package contains the exact 4,800 model inputs, archived v4 model sources,
and reference LP objectives used in Table 3(a). It covers four layouts, escorts
3 through 8, 100 seeds, and both formulations. FT and MV are fractional LP
components. The objective is FT + MV/R using each recorded R and horizon.

Use Python 3.10 or later, Gurobi 13.0.3, and a valid Gurobi license for these
model sizes. Install the Python API if needed:

```bash
python3 -m pip install -r requirements.txt
```

From this directory, reproduce the entire batch:

```bash
python3 RunStaticLP.py --input 'inputs/table2a_*.csv' --source-dir source \\
  --workers 4 --threads 1 --check-against expected_lp_results.csv \\
  -f lp_reproduced.csv
```

For an eight-case check, add `--seeds 1 --escorts 3` and use a different output.
To resume an interrupted batch, repeat its identical command with `--resume`.
Existing results are reused only after source, input, and optimality checks.
Every objective must match the reference within absolute tolerance 1e-6 and
relative tolerance 1e-8. Alternative optimal FT/MV decompositions are permitted.

The first solve has a 300-second limit. TIME_LIMIT cases are retried using
barrier without crossover and a 600-second limit. Only OPTIMAL results are
saved; an incomplete batch exits with a nonzero status. Each run writes a
manifest, per-worker solver logs, and failure diagnostics when needed.

`MANIFEST.json` records the checksums of all packaged files. Verify them with:

```bash
python3 -c 'import hashlib,json,pathlib; m=json.loads(pathlib.Path("MANIFEST.json").read_text()); assert all(hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()==h for p,h in m["files"].items()); print("Checksums verified")'
```

The archived weight and horizon must be retained. Recomputing the current v5
coefficient or using the same numeric T for both methods changes the LP being
reproduced. EF physical horizon is T+1; LF physical horizon is T.
See `LP_REPRODUCIBILITY.md` for metric definitions and standard-runner options.
The standard runners are available in the full code supplement; this standalone
bundle uses the frozen model sources directly.
"""


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "reproducibility" / "table3_lp")
    args = parser.parse_args(argv)
    destination = args.output_dir.resolve()
    zip_path = destination.parent / (destination.name + ".zip")
    if destination.exists() or zip_path.exists():
        parser.error("Destination or ZIP already exists; choose a fresh destination")
    inputs = sorted((ARCHIVE / "inputs").glob("table2a_*.csv"))
    if len(inputs) != 8:
        parser.error("Expected eight frozen Table 3 input files")
    total = 0
    for path in inputs:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != 600:
            parser.error(f"Expected 600 frozen rows in {path.name}")
        total += len(rows)
    with (ARCHIVE / "lp_results.csv").open(newline="") as handle:
        reference = list(csv.DictReader(handle))
    if len(reference) != total or any(row["status"] != "OPTIMAL" for row in reference):
        parser.error("The frozen reference must contain 4,800 OPTIMAL LP values")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="lp_package_", dir=destination.parent) as temporary:
        stage = Path(temporary) / destination.name
        (stage / "inputs").mkdir(parents=True)
        (stage / "source").mkdir()
        for path in inputs:
            shutil.copyfile(path, stage / "inputs" / path.name)
        for name in MODEL_FILES:
            shutil.copyfile(ARCHIVE / "source" / name, stage / "source" / name)
        for name in ("RunStaticLP.py", "LP_REPRODUCIBILITY.md"):
            shutil.copyfile(HERE / name, stage / name)
        shutil.copyfile(ARCHIVE / "lp_results.csv", stage / "expected_lp_results.csv")
        (stage / "requirements.txt").write_text("gurobipy==13.0.3\n")
        (stage / "README.md").write_text(PACKAGE_README)
        manifest = dict(objective="FT + MV/R", coefficient_protocol="archived v4; retain recorded R",
                        instances=2400, LP_solves=total, solver="Gurobi 13.0.3",
                        files={str(p.relative_to(stage)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(stage.rglob("*")) if p.is_file()})
        (stage / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
        staged_zip = Path(temporary) / zip_path.name
        with zipfile.ZipFile(staged_zip, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(stage.rglob("*")):
                if path.is_file():
                    archive.write(path, str(Path(destination.name) / path.relative_to(stage)))
        stage.rename(destination)
        staged_zip.rename(zip_path)
    print(f"Created {zip_path} with {total} frozen LP inputs and model sources.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
