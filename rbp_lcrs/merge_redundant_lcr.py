"""
merge_redundant_lcr.py

Merge redundant LCR calls within the same protein and same detection method,
then populate the sequence column by slicing from reference sequences.

Rules applied (per protein_id x method group, never across methods/proteins):
  T0  exact duplicate coordinates          -> keep one copy
  T1  containment (B inside A)             -> keep the outer interval
  T2  overlap >= OVERLAP_THRESHOLD of the
      shorter interval                     -> merge into union span

Sequences: every merged row's sequence is sliced from the reference
(uniprot_sequence) using 1-based inclusive coordinates, seq[start-1:end].
This is uniform across singletons and all merge tiers.

Input : INPUT_XLSX (sheet "all_results", raw calls)
        REF_XLSX   (sheet "Combined", reference sequences)
Output: OUTPUT_XLSX with
  "all_results"  -> the MERGED table with sequences filled
  "summary"      -> per-method merge statistics
  one sheet per method -> merged calls for that method only, with
                    provenance columns
                    (n_members, member_intervals, merge_tier, raw_row_ids)

raw_row_ids are 0-based row indices into the INPUT file's all_results sheet.
"""

import sys
import pandas as pd

INPUT_XLSX = "rbp_lcrs/lcr_methods_combined.xlsx"
REF_XLSX = "datasets/combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx"
OUTPUT_XLSX = "rbp_lcrs/lcr_methods_combined_merged.xlsx"

RAW_SHEET = "Combined"
REF_SHEET = "Combined"
REF_SEQ_COL = "uniprot_sequence"
REF_ACC_COL = "uniprot_accession"
OVERLAP_THRESHOLD = 0.7



def load_reference_sequences(path):
    """Read the reference sheet and return {accession: sequence}."""
    ref = pd.read_excel(path, sheet_name=REF_SHEET)

    acc_col = REF_ACC_COL
    if acc_col is None:
        lower = {c.lower(): c for c in ref.columns}
        acc_col = lower.get(REF_ACC_COL.lower())
        if acc_col is None:
            raise KeyError(
                f"No accession column found in {path}; columns: {list(ref.columns)}"
            )

    seqs = (
        ref[[acc_col, REF_SEQ_COL]]
        .dropna()
        .assign(**{REF_SEQ_COL: lambda d: d[REF_SEQ_COL].str.upper().str.strip()})
        .drop_duplicates(subset=acc_col)
        .set_index(acc_col)[REF_SEQ_COL]
        .to_dict()
    )
    print(f"Reference: {len(seqs)} sequences (key column: '{acc_col}')")
    return seqs


def slice_sequence(ref_seqs, protein_id, start, end):
    """Slice seq[start-1:end] (1-based inclusive) for the given protein."""
    acc = str(protein_id).split("|")[0]
    seq = ref_seqs.get(acc)
    if seq is None:
        return None, "missing_accession"
    if end > len(seq):
        return None, "end_beyond_sequence"
    return seq[start - 1:end], "ok"


def overlap_fraction_of_shorter(a_start, a_end, b_start, b_end):
    """Overlap length as a fraction of the shorter of the two intervals."""
    ov = min(a_end, b_end) - max(a_start, b_start) + 1
    if ov <= 0:
        return 0.0
    shorter = min(a_end - a_start + 1, b_end - b_start + 1)
    return ov / shorter


def merge_group(grp):
    """Merge one (protein_id, method) group, sorted by (start, end).

    Returns a list of dicts: merged coordinates + member info + tiers.
    """
    records = [
        {"start": int(r.start), "end": int(r.end), "raw_id": int(r.raw_id)}
        for r in grp.itertuples()
    ]

    merged = []
    for rec in records:
        # T0: exact duplicate of the current merged interval
        if merged and rec["start"] == merged[-1]["start"] and rec["end"] == merged[-1]["end"]:
            merged[-1]["members"].append(rec)
            merged[-1]["tiers"].add("T0")
            continue

        if merged:
            cur = merged[-1]
            contains = rec["start"] >= cur["start"] and rec["end"] <= cur["end"]
            ov = overlap_fraction_of_shorter(cur["start"], cur["end"], rec["start"], rec["end"])

            if contains:
                # T1: keep the outer interval, record the contained one
                cur["members"].append(rec)
                cur["tiers"].add("T1")
            elif ov >= OVERLAP_THRESHOLD:
                # T2: extend to the union span
                cur["start"] = min(cur["start"], rec["start"])
                cur["end"] = max(cur["end"], rec["end"])
                cur["members"].append(rec)
                cur["tiers"].add("T2")
            else:
                merged.append({"start": rec["start"], "end": rec["end"],
                               "members": [rec], "tiers": set()})
        else:
            merged.append({"start": rec["start"], "end": rec["end"],
                           "members": [rec], "tiers": set()})

    out = []
    for m in merged:
        out.append({
            "start": m["start"],
            "end": m["end"],
            "n_members": len(m["members"]),
            "member_intervals": ";".join(
                f"{x['start']}-{x['end']}" for x in m["members"]
            ),
            "merge_tier": "".join(sorted(m["tiers"])) if m["tiers"] else "none",
            "raw_row_ids": ";".join(str(x["raw_id"]) for x in m["members"]),
        })
    return out


def main():
    raw = pd.read_excel(INPUT_XLSX, sheet_name=RAW_SHEET)
    raw = raw.reset_index().rename(columns={"index": "raw_id"})
    ref_seqs = load_reference_sequences(REF_XLSX)

    # sanity: 1-based inclusive coordinates
    bad = raw[raw["end"] < raw["start"]]
    if not bad.empty:
        raise ValueError(f"{len(bad)} rows have end < start; check coordinate convention")

    pieces = []
    seq_status = {"ok": 0, "missing_accession": 0, "end_beyond_sequence": 0}
    qc_mismatch = 0  # singleton rows whose sliced sequence differs from input

    for (source_method, pid), grp in raw.groupby(["source_method", "protein_id"], sort=False):
        grp = grp.sort_values(["start", "end"])
        for m in merge_group(grp):
            # representative row: the raw member that exactly spans the merged
            # interval (the surviving outer call); fall back to longest member
            member_ids = set(map(int, m["raw_row_ids"].split(";")))
            members = [r for r in grp.itertuples() if int(r.raw_id) in member_ids] # type: ignore
            exact = [r for r in members
                     if int(r.start) == m["start"] and int(r.end) == m["end"] # type: ignore
                     and m["n_members"] > 1]
            rep = exact[0] if exact else max(members, key=lambda r: r.end - r.start) # type: ignore

            seq, status = slice_sequence(ref_seqs, pid, m["start"], m["end"])
            seq_status[status] += 1

            # QC: untouched rows should reproduce their input sequence exactly
            if status == "ok" and m["merge_tier"] == "none" \
                    and str(rep.sequence) != seq:
                qc_mismatch += 1

            pieces.append({
                "source_method": rep.source_method,
                "protein_id": pid,
                "header": rep.header,
                "source_method": rep.source_method,
                "start": m["start"],
                "end": m["end"],
                "length": m["end"] - m["start"] + 1,
                "sequence": seq,
                "description": rep.description,
                "n_members": m["n_members"],
                "member_intervals": m["member_intervals"],
                "merge_tier": m["merge_tier"],
                "raw_row_ids": m["raw_row_ids"],
            })

    merged_df = pd.DataFrame(pieces)

    # per-method summary statistics
    summary = (
        merged_df.groupby("source_method")
        .agg(merged_calls=("start", "size"),
             multi_member=("n_members", lambda s: int((s > 1).sum())),
             raw_calls_absorbed=("n_members", lambda s: int(s.sum())))
    )
    summary["pct_redundant"] = (
        100 * (summary["raw_calls_absorbed"] - summary["merged_calls"])
        / summary["raw_calls_absorbed"]
    ).round(1)
    summary = summary.reset_index()

    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as xw:
        merged_df.to_excel(xw, sheet_name=RAW_SHEET, index=False)
        summary.to_excel(xw, sheet_name="summary", index=False)
        for method, sub in merged_df.groupby("source_method", sort=False):
            sheet = str(method)[:31]  # Excel sheet-name limit
            sub.to_excel(xw, sheet_name=sheet, index=False)

    print(summary.to_string(index=False))
    print(f"\nSequence slicing: {seq_status}")
    if qc_mismatch:
        print(f"WARNING: {qc_mismatch} singleton rows differ from input sequence "
              f"(check coordinate convention or isoform mismatch)")
    print(f"Wrote {OUTPUT_XLSX}")


if __name__ == "__main__":
    sys.exit(main())
