#!/usr/bin/env python3
"""
run_grid.py -- simple algorithm to set up (and optionally launch) 
an YxY grid of 3D-PDR models, varying G0 and the cosmic-ray ionization rate (CRIR).

Structure of params.dat (confirmed by reading the file, not guessed):
  - It is read by Fortran list-directed I/O (`read(12,*) ...`) in
    src/input_parameters.F90, so column alignment does not matter for
    parsing -- but we preserve it anyway for human readability.
  - Line 7  (1-indexed, as in a text editor) = Output prefix
  - Line 11                                  = G0 (Draine field units)
  - Line 12                                  = Cosmic-ray ionization rate (s^-1)
  - indir (line 4) and outdir (line 6) are "ics" and "sims" respectively,
    and are shared/relative to wherever the 3DPDR executable is run from.
    Multiple models can therefore coexist by using distinct output
    prefixes -- there is no need for separate per-model working
    directories; the code already supports many params files pointing at
    one shared sims/ output directory (see src/read_command_line.F90,
    which accepts `-p=<paramfile>` to pick which params file to read).

This script only ever touches its own generated copies -- the original
params.dat (the template) is never modified.

Usage:
    python3 run_grid.py                # setup-only (default): write the 9
                                        # params files and print what would
                                        # be run -- does NOT launch 3DPDR
    python3 run_grid.py --setup-only   # same as above, explicit
    python3 run_grid.py --execute      # actually run all models,
                                        # sequentially, in the foreground
                                        # (slow -- see README/RAM notes)
"""

import argparse
import subprocess
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
BASE_PARAMS = BASE_DIR / "params.dat"
GRID_DIR = BASE_DIR / "grid_runs"
EXECUTABLE = BASE_DIR / "3DPDR"

# 1-indexed line numbers (as counted in a text editor) that we edit.
LINE_OUTPUT_PREFIX = 7
LINE_G0 = 11
LINE_CRIR = 12

G0_VALUES = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
CRIR_VALUES = [1e-17, 2e-17, 5e-17, 1e-16, 2e-16, 5e-16, 1e-15, 2e-15, 5e-15, 1e-14]


def make_prefix(g0, crir):
    """Filesystem-safe, sortable, descriptive prefix, e.g. 'G0_100_CR_1e-15'."""
    crir_str = f"{crir:.0e}"  # '1e-17', '1e-16', '1e-15'
    return f"G0_{g0}_CR_{crir_str}"


def replace_value(line, new_value):
    """
    Replace the value portion of a params.dat line, keeping the existing
    '!' comment (and its column position) exactly as in the original file.
    """
    stripped = line.rstrip("\n")
    bang_idx = stripped.index("!")
    comment = stripped[bang_idx:]
    return new_value.ljust(bang_idx) + comment + "\n"


def build_params_lines(g0, crir, prefix):
    base_lines = BASE_PARAMS.read_text().splitlines(keepends=True)

    new_lines = list(base_lines)  # shallow copy; only 3 entries get replaced
    new_lines[LINE_OUTPUT_PREFIX - 1] = replace_value(
        base_lines[LINE_OUTPUT_PREFIX - 1], prefix
    )
    new_lines[LINE_G0 - 1] = replace_value(
        base_lines[LINE_G0 - 1], str(g0)
    )
    new_lines[LINE_CRIR - 1] = replace_value(
        base_lines[LINE_CRIR - 1], f"{crir:.3E}"
    )
    return new_lines


def setup_grid():
    """Write all params files into grid_runs/. Returns list of (prefix, path)."""
    GRID_DIR.mkdir(exist_ok=True)
    models = []

    for g0 in G0_VALUES:
        for crir in CRIR_VALUES:
            prefix = make_prefix(g0, crir)
            lines = build_params_lines(g0, crir, prefix)

            out_path = GRID_DIR / f"params_{prefix}.dat"
            out_path.write_text("".join(lines))

            print(f"Wrote {out_path.relative_to(BASE_DIR)}")
            print(f"  line {LINE_OUTPUT_PREFIX:>2} : {lines[LINE_OUTPUT_PREFIX - 1].rstrip()}")
            print(f"  line {LINE_G0:>2} : {lines[LINE_G0 - 1].rstrip()}")
            print(f"  line {LINE_CRIR:>2} : {lines[LINE_CRIR - 1].rstrip()}")
            print()

            models.append((prefix, out_path))

    return models


def run_model(prefix, params_path):
    """Actually launch 3DPDR for one model, in the foreground, from BASE_DIR."""
    rel_params = params_path.relative_to(BASE_DIR)
    cmd = [str(EXECUTABLE), f"-p={rel_params}"]
    print(f"[RUN] {prefix}: {' '.join(cmd)}  (cwd={BASE_DIR})")
    subprocess.run(cmd, cwd=BASE_DIR, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually launch all 3DPDR runs sequentially (slow). "
             "Without this flag, the script only writes the params files.",
    )
    parser.add_argument(
        "--setup-only",
        action="store_true",
        help="Explicit no-op flag (this is the default behaviour) -- "
             "only write the params.dat files, do not run anything.",
    )
    args = parser.parse_args()

    if not BASE_PARAMS.exists():
        sys.exit(f"Base params file not found: {BASE_PARAMS}")

    print(f"Base template : {BASE_PARAMS}")
    print(f"Grid dir      : {GRID_DIR}")
    print(f"G0 values     : {G0_VALUES}")
    print(f"CRIR values   : {CRIR_VALUES}")
    print()

    models = setup_grid()

    print(f"Generated {len(models)} params files in {GRID_DIR}")

    if args.execute:
        print("\n--execute given: launching all models sequentially...\n")
        for prefix, params_path in models:
            run_model(prefix, params_path)
    else:
        print("\nSetup-only (default). Nothing was run.")
        print("To actually launch all models sequentially, re-run with --execute.")
        print("Example single-model command (run manually from this directory):")
        example_prefix, example_path = models[0]
        print(f"  ./3DPDR -p={example_path.relative_to(BASE_DIR)}")


if __name__ == "__main__":
    main()
