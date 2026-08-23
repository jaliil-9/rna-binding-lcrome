"""
Run the MATLAB/Octave LCRFinder implementation from a standard terminal.

"""

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def matlab_quote(value):
    """Quote a path/string for a MATLAB or Octave expression."""
    return str(value).replace("'", "''").replace("\\", "/")


def find_engine(requested):
    """Return the requested executable, or discover Octave/MATLAB on PATH."""
    if requested:
        executable = shutil.which(requested) or requested
        return executable

    for executable in ("octave-cli", "octave", "matlab"):
        if shutil.which(executable):
            return executable

    raise FileNotFoundError(
        "Neither Octave nor MATLAB was found on PATH. "
        "Install GNU Octave or pass --engine with the executable path."
    )


def build_command(engine, lcrfinder_dir, fasta_path, output_dir, resolution):
    """Build a non-interactive Octave or MATLAB command."""
    expression = (
        f"addpath('{matlab_quote(lcrfinder_dir)}'); "
        f"run_lcrfinder_batch('{matlab_quote(fasta_path)}', "
        f"'{matlab_quote(output_dir)}', '{resolution}');"
    )

    if Path(engine).name.lower().startswith("matlab"):
        return [engine, "-batch", expression]

    return [engine, "--quiet", "--eval", expression]


def main():
    parser = argparse.ArgumentParser(
        description="Run LCRFinder and produce LCRFinder.csv and LCRFinder.fasta."
    )
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--lcrfinder-dir",
        type=Path,
        required=True,
        help="Directory containing run_lcrfinder_batch.m and source/.",
    )
    parser.add_argument(
        "--resolution",
        choices=["Low", "Medium", "High"],
        default="Medium",
        help="LCRFinder recurrence resolution (default: Medium).",
    )
    parser.add_argument(
        "--engine",
        help="Octave or MATLAB executable; auto-detected if omitted.",
    )
    args = parser.parse_args()

    fasta_path = args.input_fasta.resolve()
    output_dir = args.output_dir.resolve()
    lcrfinder_dir = args.lcrfinder_dir.resolve()
    driver = lcrfinder_dir / "run_lcrfinder_batch.m"
    source_dir = lcrfinder_dir / "source"

    if not fasta_path.is_file():
        sys.exit(f"Input FASTA not found: {fasta_path}")
    if not driver.is_file() or not source_dir.is_dir():
        sys.exit(
            "Invalid --lcrfinder-dir. It must contain run_lcrfinder_batch.m and source/."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    engine = find_engine(args.engine)
    command = build_command(
        engine, lcrfinder_dir, fasta_path, output_dir, args.resolution
    )

    manifest = {
        "tool": "LCRFinder",
        "resolution": args.resolution,
        "input_fasta": str(fasta_path),
        "lcrfinder_dir": str(lcrfinder_dir),
        "engine": engine,
        "command": command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    (output_dir / "lcrfinder_run_manifest.json").write_text(
        json.dumps(manifest, indent=2)
    )

    print("Running LCRFinder...")
    print("Command:", " ".join(command))
    completed = subprocess.run(command)
    if completed.returncode != 0:
        sys.exit(f"LCRFinder failed with exit code {completed.returncode}.")

    expected = [output_dir / "LCRFinder.csv", output_dir / "LCRFinder.fasta"]
    missing = [str(path) for path in expected if not path.is_file()]
    if missing:
        sys.exit("LCRFinder ended but expected output is missing: " + ", ".join(missing))

    print("Finished successfully.")
    print(f"CSV:   {expected[0]}")
    print(f"FASTA: {expected[1]}")


if __name__ == "__main__":
    main()
