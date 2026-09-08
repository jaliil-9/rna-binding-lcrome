"""
Run AlcoR mapper on a protein FASTA and write AlcoR.csv and AlcoR.fasta.

"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_MODELS = "6:10:0:0:10:0.9/3:1:0.9"


def read_fasta(path):
    records = []
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence found before a FASTA header in {path}")
                sequence.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(sequence)))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def find_alcor(requested):
    candidate = Path(requested).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    discovered = shutil.which(requested)
    if discovered:
        return Path(discovered).resolve()
    raise FileNotFoundError(
        f"AlcoR executable was not found: {requested}. "
        "Pass --alcor with the path to the compiled AlcoR executable."
    )


def masked_intervals(sequence):
    start = None
    for index, residue in enumerate(sequence, start=1):
        if residue.islower() and start is None:
            start = index
        elif not residue.islower() and start is not None:
            yield start, index - 1
            start = None
    if start is not None:
        yield start, len(sequence)


def protein_id(header):
    """UniProt accession from a FASTA header.

    '>sp|A0A087X1C5|...' or '>tr|...|...' -> 'A0A087X1C5'
    '>Q5T200|ensembl_protein_id=...'     -> 'Q5T200'
    """
    text = header[1:].strip()
    if "|" in text:
        parts = text.split("|")
        if parts[0] in {"sp", "tr"} and len(parts) > 1:
            return parts[1].strip()
        return parts[0].strip()
    return text.split()[0]


def description(fragment):
    counts = Counter(fragment.upper())
    residue, _ = max(counts.items(), key=lambda item: (item[1], item[0]))
    return f"{residue} rich region"


def main():
    parser = argparse.ArgumentParser(
        description="Run AlcoR mapper and convert its lowercase mask into CSV and FASTA LCR calls."
    )
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--alcor",
        default="./build/AlcoR",
        help="AlcoR executable path or executable on PATH (default: ./build/AlcoR).",
    )
    parser.add_argument(
        "--models",
        default=DEFAULT_MODELS,
        help="AlcoR protein models passed to -m (default: %(default)s).",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=5,
        help="AlcoR moving-average window passed to -w (default: %(default)s).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        help="Optional AlcoR segmentation threshold passed to -t.",
    )
    parser.add_argument(
        "--ignore",
        type=int,
        help="Optional minimum AlcoR region length passed to -i.",
    )
    args = parser.parse_args()

    input_fasta = args.input_fasta.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not input_fasta.is_file():
        sys.exit(f"Input FASTA not found: {input_fasta}")
    if args.window < 1:
        sys.exit("--window must be at least 1.")
    if args.ignore is not None and args.ignore < 1:
        sys.exit("--ignore must be at least 1.")

    alcor = find_alcor(args.alcor)
    output_dir.mkdir(parents=True, exist_ok=True)
    masked_fasta = output_dir / "AlcoR_masked.fasta"
    csv_path = output_dir / "AlcoR.csv"
    lcr_fasta_path = output_dir / "AlcoR.fasta"
    manifest_path = output_dir / "alcor_run_manifest.json"

    command = [
        str(alcor), "mapper", "-v", "-w", str(args.window), "-m", args.models,
        "-k", "-o", str(masked_fasta), str(input_fasta),
    ]
    if args.threshold is not None:
        command[2:2] = ["-t", str(args.threshold)]
    if args.ignore is not None:
        command[2:2] = ["-i", str(args.ignore)]

    manifest = {
        "tool": "AlcoR",
        "subcommand": "mapper",
        "input_fasta": str(input_fasta),
        "alcor_executable": str(alcor),
        "models": args.models,
        "window": args.window,
        "threshold": args.threshold,
        "ignore": args.ignore,
        "command": command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("Running AlcoR mapper...")
    print("Command:", " ".join(command))
    completed = subprocess.run(command)
    if completed.returncode != 0:
        sys.exit(f"AlcoR failed with exit code {completed.returncode}.")
    if not masked_fasta.is_file():
        sys.exit(f"AlcoR finished but did not create the masked FASTA: {masked_fasta}")

    original_records = read_fasta(input_fasta)
    masked_records = read_fasta(masked_fasta)
    if len(original_records) != len(masked_records):
        sys.exit(
            "AlcoR masked FASTA has a different number of records than the input; "
            "cannot safely assign positions to the original headers."
        )

    rows = []
    with lcr_fasta_path.open("w", encoding="utf-8", newline="\n") as fasta_out:
        for (header, original), (_, masked) in zip(original_records, masked_records):
            if len(original) != len(masked):
                sys.exit(f"Masked sequence length differs from input for {protein_id(header)}.")
            for start, end in masked_intervals(masked):
                fragment = original[start - 1:end].upper()
                length = end - start + 1
                rows.append({
                    "protein_id": protein_id(header),
                    "header": header,
                    "method": "AlcoR",
                    "start": start,
                    "end": end,
                    "length": length,
                    "sequence": fragment,
                    "description": description(fragment),
                })
                fasta_out.write(
                    f"{header} method=AlcoR start={start} end={end} length={length}\n"
                )
                fasta_out.write(fragment + "\n")

    fieldnames = ["protein_id", "header", "method", "start", "end", "length", "sequence", "description"]
    with csv_path.open("w", encoding="utf-8", newline="") as csv_out:
        writer = csv.DictWriter(csv_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["lcr_count"] = len(rows)
    manifest["masked_fasta"] = str(masked_fasta)
    manifest["csv"] = str(csv_path)
    manifest["lcr_fasta"] = str(lcr_fasta_path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("Finished successfully.")
    print(f"LCR calls: {len(rows)}")
    print(f"CSV:         {csv_path}")
    print(f"LCR FASTA:   {lcr_fasta_path}")
    print(f"Masked FASTA:{masked_fasta}")


if __name__ == "__main__":
    main()
