"""
Create a simple LCR-caller comparison workbook and four figures.

Usage:
  python lcr_tool_comparison.py --lcr lcr_methods_combined.xlsx \
      --proteins combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx \
      --out lcr_tool_comparison

Expected LCR columns (case-insensitive aliases are accepted): protein ID, start, end.
Optional: tool/method and LCR sequence. Each sheet is treated as one tool when no
Tool/Method column is present. Expected protein columns: protein ID and sequence.
"""
import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from openpyxl import load_workbook
from openpyxl.drawing.image import Image as XLImage

AA = list("ACDEFGHIKLMNPQRSTVWY")
ALIASES = {
    "protein_id": [
        "ensemblproteinid",
        "ensembl_protein_id",
        "protein_id",
        "protein id",
        "protein",
        "accession",
        "uniprot",
        "uniprot_id",
        "uniprotaccession",
        "entry",
        "id",
        "sequence_id",
        "seq_id"
    ],
    "start": ["start", "from", "begin", "lcr_start", "start_position", "start position"],
    "end": ["end", "to", "stop", "lcr_end", "end_position", "end position"],
    "sequence": [
        "uniprotsequence",
        "uniprot_sequence",
        "sequence",
        "protein_sequence",
        "protein sequence",
        "seq",
        "aa_sequence",
        "amino_acid_sequence"
    ],
    "lcr_sequence": ["lcr_sequence", "lcr sequence", "region_sequence", "region sequence", "subsequence", "subseq"],
    "tool": ["tool", "method", "caller", "lcr_method", "lcr method"]
}

def canonical_id(x):
    x = str(x).strip()

    return x.split("|")[0].split()[0]

def find_col(df, key, required=True):
    norm = {str(c).strip().lower(): c for c in df.columns}
    for alias in ALIASES[key]:
        if alias in norm: return norm[alias]
    if required: raise ValueError(f"Cannot find {key!r} column. Available: {list(df.columns)}")
    return None

def clean_sequence(x):
    return re.sub(r"[^A-Za-z]", "", str(x).upper())

def entropy(seq):
    seq = clean_sequence(seq)
    if not seq: return np.nan
    p = np.array([seq.count(a) / len(seq) for a in AA]); p = p[p > 0]
    return float(-(p * np.log2(p)).sum())

def read_proteins(path):
    def normalize_header(x):
        return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())

    book = pd.ExcelFile(path)

    for sheet in book.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)

        normalized_columns = {
            normalize_header(column): column
            for column in df.columns
        }

        accession_col = normalized_columns.get("uniprotaccession")
        sequence_col = normalized_columns.get("uniprotsequence")

        if accession_col is None or sequence_col is None:
            continue

        print(f"Using protein sheet: {sheet}")
        print(f"Using ID column: {accession_col}")
        print(f"Using sequence column: {sequence_col}")

        out = df[[accession_col, sequence_col]].copy()
        out.columns = ["protein_id", "protein_sequence"]

        out["protein_id"] = (
            out["protein_id"]
            .astype(str)
            .str.strip()
            .str.upper()
        )

        out["protein_sequence"] = out["protein_sequence"].map(clean_sequence)

        out = (
            out[out["protein_sequence"].str.len() > 0]
            .drop_duplicates("protein_id")
            .reset_index(drop=True)
        )

        return out

    available = {
        sheet: pd.read_excel(path, sheet_name=sheet, nrows=0).columns.tolist()
        for sheet in book.sheet_names
    }

    raise ValueError(
        "Could not find UniProt accession and sequence columns. "
        f"Available columns by sheet: {available}"
    )

def read_calls(path):
    book = pd.ExcelFile(path)

    sheet_lookup = {
        re.sub(r"[^a-z0-9]+", "", str(sheet).lower()): sheet
        for sheet in book.sheet_names
    }

    if "allresults" not in sheet_lookup:
        raise ValueError(
            f"Could not find the 'all_results' sheet. "
            f"Available sheets: {book.sheet_names}"
        )

    sheet = sheet_lookup["allresults"]
    df = pd.read_excel(path, sheet_name=sheet)

    pid = find_col(df, "protein_id")
    start = find_col(df, "start")
    end = find_col(df, "end")
    toolcol = find_col(df, "tool", required=False)
    lcrseq = find_col(df, "lcr_sequence", required=False)

    selected_columns = [pid, start, end]

    if toolcol:
        selected_columns.append(toolcol)

    if lcrseq:
        selected_columns.append(lcrseq)

    out = df[selected_columns].copy()

    out.columns = (
        ["protein_id", "start", "end"]
        + (["tool"] if toolcol else [])
        + (["input_lcr_sequence"] if lcrseq else [])
    )

    if not toolcol:
        out["tool"] = sheet

    out["protein_id"] = (
        out["protein_id"]
        .map(canonical_id)
        .astype(str)
        .str.strip()
        .str.upper()
    )

    out["tool"] = (
        out["tool"]
        .astype(str)
        .str.strip()
    )

    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"] = pd.to_numeric(out["end"], errors="coerce")

    out = out.dropna(subset=["protein_id", "start", "end", "tool"]).copy()
    out[["start", "end"]] = out[["start", "end"]].astype(int)

    raw_n = len(out)

    # Remove repeated copies of exactly the same call.
    # Overlapping calls with different coordinates remain separate.
    out = out.drop_duplicates(
        subset=["tool", "protein_id", "start", "end"]
    ).copy()

    duplicate_n = raw_n - len(out)

    print(f"Using LCR sheet: {sheet}")
    print(f"Raw LCR rows: {raw_n}")
    print(f"Exact duplicate rows removed: {duplicate_n}")
    print("\nUnique LCR calls per tool:")
    print(out.groupby("tool").size().sort_values(ascending=False).to_string())

    return out

def build_table(calls, proteins):
    d = calls.merge(proteins, on="protein_id", how="left", validate="many_to_one")
    if d["protein_sequence"].notna().sum() == 0:
        raise ValueError(
            "No LCR calls matched protein sequences. "
            "Check that LCR protein IDs and protein-workbook IDs use the same identifier type."
        )
    missing = d.protein_sequence.isna().sum()
    if missing: print(f"Warning: {missing} calls lack a matched whole-protein sequence and are omitted.")
    d = d.dropna(subset=["protein_sequence"]).copy()
    d["protein_length"] = d.protein_sequence.str.len()
    d = d[(d.start >= 1) & (d.end >= d.start) & (d.end <= d.protein_length)].copy()
    d["lcr_length"] = d.end - d.start + 1
    d["lcr_sequence"] = [s[a-1:b] for s,a,b in zip(d.protein_sequence, d.start, d.end)]
    d["shannon_entropy_bits"] = d.lcr_sequence.map(entropy)
    d["location"] = pd.cut((d.start + d.end) / 2 / d.protein_length,
                           bins=[0, 1/3, 2/3, 1], labels=["N-terminal", "Internal", "C-terminal"], include_lowest=True)
    for a in AA: d[f"freq_{a}"] = d.lcr_sequence.map(lambda s, aa=a: s.count(aa) / len(s))
    keep = ["tool", "protein_id", "start", "end", "protein_length", "lcr_length", "lcr_sequence", "shannon_entropy_bits", "location"] + [f"freq_{a}" for a in AA]
    return d[keep].sort_values(["tool", "protein_id", "start", "end"])

def save_figures(d, proteins, outdir):
    sns.set_theme(style="whitegrid", context="talk")
    order = sorted(d.tool.unique())
    palette = dict(zip(order, sns.color_palette("tab10", len(order))))
    figs = []
    n_total_proteins = proteins["protein_id"].nunique()

    summary = (
        d.groupby("tool")
        .agg(
            LCRs=("tool", "size"),
            proteins_with_LCR=("protein_id", "nunique"),
            mean_lcr_length=("lcr_length", "mean"),
            mean_shannon_entropy_bits=("shannon_entropy_bits", "mean")
        )
        .reindex(order)
    )

    summary["mean_lcrs_per_protein"] = (
        summary["LCRs"] / n_total_proteins
    )

    summary["mean_lcrs_per_positive_protein"] = (
        summary["LCRs"] / summary["proteins_with_LCR"]
    )

    summary["mean_lcr_length"] = summary["mean_lcr_length"].round(2)
    summary["mean_shannon_entropy_bits"] = summary["mean_shannon_entropy_bits"].round(3)
    summary["mean_lcrs_per_protein"] = summary["mean_lcrs_per_protein"].round(3)
    summary["mean_lcrs_per_positive_protein"] = (
        summary["mean_lcrs_per_positive_protein"].round(3)
    )
    fig, ax = plt.subplots(figsize=(10, 5))

    x = np.arange(len(summary))
    width = 0.38

    bars_lcrs = ax.bar(
        x - width / 2,
        summary["LCRs"],
        width,
        label="LCRs",
        color="#4C78A8"
    )

    bars_proteins = ax.bar(
        x + width / 2,
        summary["proteins_with_LCR"],
        width,
        label="Proteins with ≥1 LCR",
        color="#F58518"
    )

    ax.set_xlabel("Tool")
    ax.set_ylabel("Count")
    ax.set_xticks(x)
    ax.set_xticklabels(summary.index)
    ax.legend()

    ax.bar_label(
        bars_lcrs,
        labels=[f"{value:,.0f}" for value in summary["LCRs"]],
        padding=3,
        fontsize=9
    )

    ax.bar_label(
        bars_proteins,
        labels=[f"{value:,.0f}" for value in summary["proteins_with_LCR"]],
        padding=3,
        fontsize=9
    )

    ax.margins(y=0.12)
    fig.tight_layout()

    figs.append((fig, "figure_1_detection_abundance.png"))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5)); sns.boxplot(data=d, x="tool", y="lcr_length", order=order, palette=palette, showfliers=False, ax=axes[0])
    axes[0].set(yscale="log", xlabel="Tool", ylabel="LCR length (aa, log scale)")
    sns.boxplot(data=d, x="tool", y="shannon_entropy_bits", order=order, palette=palette, showfliers=False, ax=axes[1])
    axes[1].set(xlabel="Tool", ylabel="Shannon entropy (bits)")
    for ax in axes: ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    figs.append((fig, "figure_2_length_entropy.png"))
    loc = pd.crosstab(d.tool, d.location, normalize="index").reindex(order).reindex(columns=["N-terminal","Internal","C-terminal"], fill_value=0)*100
    summary = summary.join(
        loc.rename(columns={
            "N-terminal": "pct_N_terminal",
            "Internal": "pct_internal",
            "C-terminal": "pct_C_terminal"
        })
    )

    summary[["pct_N_terminal", "pct_internal", "pct_C_terminal"]] = (
        summary[["pct_N_terminal", "pct_internal", "pct_C_terminal"]].round(2)
    )
    fig, ax = plt.subplots(figsize=(9,5)); loc.plot.bar(stacked=True, ax=ax, color=["#4C78A8", "#BAB0AC", "#F58518"])
    ax.set(xlabel="Tool", ylabel="LCRs (%)"); ax.legend(title="Location", bbox_to_anchor=(1.02,1), loc="upper left"); plt.xticks(rotation=0); fig.tight_layout(); figs.append((fig, "figure_3_location.png"))
    bgseq = "".join(proteins.protein_sequence.tolist()); bg = np.array([bgseq.count(a)/len(bgseq) for a in AA])
    enrich = pd.DataFrame({tool: np.log2((d[d.tool==tool][[f"freq_{a}" for a in AA]].mean().to_numpy()+1e-6)/(bg+1e-6)) for tool in order}, index=AA)
    mean_aa_frequency = (
        d.groupby("tool")[[f"freq_{aa}" for aa in AA]]
        .mean()
        .reindex(order)
    )

    mean_aa_frequency.columns = [
        f"mean_freq_{column.replace('freq_', '')}"
        for column in mean_aa_frequency.columns
    ]

    summary = summary.join(mean_aa_frequency)

    aa_enrichment = enrich.T.copy()
    aa_enrichment.columns = [
        f"log2_enrichment_{aa}"
        for aa in aa_enrichment.columns
    ]

    summary = summary.join(aa_enrichment)

    for aa in AA:
        summary[f"mean_freq_{aa}"] = summary[f"mean_freq_{aa}"].round(4)
        summary[f"log2_enrichment_{aa}"] = summary[f"log2_enrichment_{aa}"].round(3)
    fig, ax = plt.subplots(figsize=(9,8)); sns.heatmap(enrich, cmap="vlag", center=0, linewidths=.4, cbar_kws={"label":"log2 enrichment vs full dataset"}, ax=ax)
    ax.set(xlabel="Tool", ylabel="Amino acid"); fig.tight_layout(); figs.append((fig, "figure_4_amino_acid_enrichment.png"))
    paths=[]
    for fig, name in figs:
        p=outdir/name; fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig); paths.append(p)
    return summary.reset_index(), paths

def main():
    parser=argparse.ArgumentParser(); 
    parser.add_argument("--lcr", default="rbp_lcrs/lcr_methods_combined_merged.xlsx") 
    parser.add_argument("--proteins", default="datasets/combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx")
    parser.add_argument("--out", default="rbp_lcrs/lcr_tool_comparison/after_merge")

    args=parser.parse_args(); outdir=Path(args.out) 
    outdir.mkdir(parents=True, exist_ok=True)

    proteins=read_proteins(args.proteins)
    d=build_table(read_calls(args.lcr), proteins) 
    summary, figs=save_figures(d, proteins, outdir)

    xlsx=outdir/"lcr_tool_comparison.xlsx"

    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        d.to_excel(w, sheet_name="Combined_results", index=False)
        summary.to_excel(w, sheet_name="Summary", index=False)

        for tool, sub in d.groupby("tool", sort=True): 
            sub.to_excel(w, sheet_name=re.sub(r"[\\/*?:\[\]]", "_", str(tool))[:31], index=False) # type: ignore

    wb=load_workbook(xlsx)
    ws=wb.create_sheet("Figures")

    for i,p in enumerate(figs):
        img=XLImage(str(p)) 
        img.width=720 
        img.height=400
        ws.add_image(img, f"A{1+i*23}")
    wb.save(xlsx)

    print(f"Wrote {xlsx} and {len(figs)} PNG figures to {outdir}")

if __name__ == "__main__": main()
