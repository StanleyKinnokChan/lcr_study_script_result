#!/usr/bin/env python3
"""
Amino acid identity of terminal vs internal LCRs.

fLPS 2.0 (SINGLE-type rows) names the dominant amino acid in col 8 as {X}.
This script re-parses raw fLPS output to extract residue identity, classifies
LCRs as terminal (bin 1 or 20) vs internal, and computes amino acid frequency
distributions across phyla.

Outputs:
  results/aa_composition.tsv  — per-phylum × residue × location counts
  figures/fig5_aa_composition.pdf
"""

import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from config import (
    PROJECT_DIR, RESULTS_DIR, FIGURES_DIR, FLPS_DIR,
    N_BINS, TERMINAL_BINS, PURITY_THRESHOLD, MIN_LCR_LENGTH,
    PHYLUM_ORDER,
    AA_ORDER, AA_COLOURS,
)


def parse_flps_with_residue(filepath: Path) -> list[dict]:
    records = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "SINGLE":
                continue
            try:
                protein_len = int(cols[1])
                start       = int(cols[4])
                end         = int(cols[5])
                aa_count    = int(cols[6])
                residue_raw = cols[8]                     # e.g. "{P}"
            except (ValueError, IndexError):
                continue

            lcr_len = end - start + 1
            if lcr_len < MIN_LCR_LENGTH:
                continue
            if aa_count / lcr_len < PURITY_THRESHOLD:
                continue

            residue = residue_raw.strip("{}").upper()
            if len(residue) != 1 or residue not in AA_ORDER:
                continue                                  # skip non-standard

            midpoint = (start + end) / 2
            norm_mid = midpoint / protein_len
            bin_num  = min(int(norm_mid * N_BINS) + 1, N_BINS)

            records.append({
                "residue":     residue,
                "bin":         bin_num,
                "is_terminal": bin_num in TERMINAL_BINS,
            })
    return records


def main():
    manifest_path = RESULTS_DIR / "species_manifest.tsv"
    if not manifest_path.exists():
        print("ERROR: species_manifest.tsv missing — run 01_download_proteomes.py first.")
        return

    manifest = pd.read_csv(manifest_path, sep="\t")
    all_rows = []

    for _, row in manifest.iterrows():
        sp_key  = row["species_key"]
        phylum  = row["phylum"]
        flps_f  = FLPS_DIR / f"{sp_key}.flps.txt"
        if not flps_f.exists():
            continue
        recs = parse_flps_with_residue(flps_f)
        for r in recs:
            r["phylum"] = phylum
        all_rows.extend(recs)
        print(f"  [{phylum}] {sp_key}: {len(recs)} LCRs parsed")

    if not all_rows:
        print("No records found.")
        return

    df = pd.DataFrame(all_rows)

    # ── Frequency table: phylum × location × residue ─────────────────────────
    df["location"] = df["is_terminal"].map({True: "terminal", False: "internal"})
    counts = (
        df.groupby(["phylum", "location", "residue"])
        .size()
        .reset_index(name="count")
    )
    # Normalise to fraction within each phylum × location
    totals = counts.groupby(["phylum", "location"])["count"].transform("sum")
    counts["fraction"] = counts["count"] / totals

    out_tsv = RESULTS_DIR / "aa_composition.tsv"
    counts.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nAA composition table: {out_tsv}  ({len(counts)} rows)")

    # ── Global summary: which residues are enriched at termini? ──────────────
    global_counts = (
        df.groupby(["location", "residue"])
        .size()
        .reset_index(name="count")
    )
    global_totals = global_counts.groupby("location")["count"].transform("sum")
    global_counts["fraction"] = global_counts["count"] / global_totals
    pivot_global = global_counts.pivot(index="residue", columns="location", values="fraction").fillna(0)
    pivot_global["enrichment_ratio"] = (
        pivot_global.get("terminal", 0) / pivot_global.get("internal", 0.001)
    )
    pivot_global = pivot_global.sort_values("enrichment_ratio", ascending=False)
    print("\nGlobal terminal vs internal amino acid enrichment (terminal/internal ratio):")
    print(pivot_global[["terminal", "internal", "enrichment_ratio"]].round(3).to_string())

    # ── Figure 5: stacked bar per phylum, terminal vs internal ───────────────
    phyla_present = [p for p in PHYLUM_ORDER if p in counts["phylum"].values]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

    for ax, loc in zip(axes, ["terminal", "internal"]):
        pivot = (
            counts[counts["location"] == loc]
            .pivot_table(index="phylum", columns="residue", values="fraction", aggfunc="mean")
            .reindex(phyla_present)
            .fillna(0)
        )
        # Reorder columns by AA_ORDER, keep only those present
        cols = [aa for aa in AA_ORDER if aa in pivot.columns]
        pivot = pivot[cols]

        bottom = np.zeros(len(pivot))
        for aa in cols:
            vals = pivot[aa].values
            ax.bar(
                range(len(pivot)), vals, bottom=bottom,
                color=AA_COLOURS.get(aa, "#cccccc"),
                label=aa, width=0.8
            )
            bottom += vals

        ax.set_xticks(range(len(pivot)))
        ax.set_xticklabels(phyla_present, rotation=30, ha="right", fontsize=9)
        ax.set_title(f"{'Terminal (bins 1 & 20)' if loc == 'terminal' else 'Internal (bins 2–19)'}\nLCR amino acid composition", fontsize=11)
        ax.set_ylabel("Fraction of LCRs" if loc == "terminal" else "")
        ax.set_ylim(0, 1)

    # Shared legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=AA_COLOURS.get(aa, "#ccc"), label=aa)
               for aa in AA_ORDER]
    fig.legend(handles=handles, title="Residue", loc="center right",
               bbox_to_anchor=(1.0, 0.5), fontsize=8, ncol=2)
    fig.suptitle("Amino acid composition of terminal vs internal LCRs across phyla\n"
                 "(purity ≥ 70%, SINGLE-residue LCRs only)", fontsize=12, y=1.01)
    plt.tight_layout()

    out = FIGURES_DIR / "fig5_aa_composition.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nFigure 5 saved: {out}")
    plt.close()

    # ── Figure 5b: enrichment ratio bar (global) ─────────────────────────────
    top_aas = pivot_global.head(10)
    fig, ax = plt.subplots(figsize=(7, 4))
    colours = [AA_COLOURS.get(aa, "#ccc") for aa in top_aas.index]
    ax.bar(top_aas.index, top_aas["enrichment_ratio"], color=colours,
           edgecolor="black", linewidth=0.7)
    ax.axhline(1.0, color="red", linestyle=":", linewidth=1.2,
               label="No enrichment (ratio = 1)")
    ax.set_ylabel("Terminal / internal frequency ratio", fontsize=11)
    ax.set_title("Amino acids enriched at protein termini (all species pooled)", fontsize=12)
    ax.legend(fontsize=9)
    plt.tight_layout()

    out2 = FIGURES_DIR / "fig5b_aa_enrichment_ratio.pdf"
    fig.savefig(out2, dpi=300, bbox_inches="tight")
    fig.savefig(str(out2).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 5b saved: {out2}")
    plt.close()


if __name__ == "__main__":
    main()
