#!/usr/bin/env python3
"""
Per-domain amino acid composition of terminal vs internal LCRs.

Groups species into four broad domains and compares which amino acids dominate
terminal LCRs in each — a direct mechanistic test of whether the same residues
are enriched at termini across all domains or whether prokaryotes show a
distinct compositional signature.

Domains:
  Bacteria     — phylum == "Bacteria"
  Archaea      — phylum == "Archaea"
  Viridiplantae— phylum == "Viridiplantae"
  Metazoa      — all metazoan phyla
  Other_Euk    — remaining eukaryotes (SAR, Excavata, Fungi, etc.)

Expected result under universal mechanism: same residues (C, E, K) enriched
at termini in all domains.  If prokaryote terminal LCRs are dominated by
different residues (e.g. P, G, A), that supports distinct mechanisms.

Outputs:
  results/domain_aa_composition.tsv  — Supp Table S7
  figures/suppfig_domain_aa.pdf
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from config import (
    RESULTS_DIR, FIGURES_DIR, FLPS_DIR,
    N_BINS, TERMINAL_BINS, PURITY_THRESHOLD, MIN_LCR_LENGTH,
    SUPERGROUP_OF, AA_ORDER, AA_COLOURS,
)

METAZOAN_PHYLA = {
    p for p, sg in SUPERGROUP_OF.items() if sg == "Metazoa"
}

DOMAIN_MAP: dict[str, str] = {
    "Bacteria": "Bacteria",
    "Archaea":  "Archaea",
    "Viridiplantae": "Viridiplantae",
}


def classify_domain(phylum: str) -> str:
    if phylum in DOMAIN_MAP:
        return DOMAIN_MAP[phylum]
    if phylum in METAZOAN_PHYLA:
        return "Metazoa"
    return "Other_Euk"


def parse_flps_residue(filepath: Path) -> list[dict]:
    records = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "SINGLE":
                continue
            try:
                protein_len = int(cols[1])
                start       = int(cols[4])
                end         = int(cols[5])
                aa_count    = int(cols[6])
                residue_raw = cols[8]
            except (ValueError, IndexError):
                continue

            lcr_len = end - start + 1
            if lcr_len < MIN_LCR_LENGTH:
                continue
            if aa_count / lcr_len < PURITY_THRESHOLD:
                continue

            residue = residue_raw.strip("{}").upper()
            if len(residue) != 1 or residue not in AA_ORDER:
                continue

            midpoint = (start + end) / 2
            bin_num  = min(int((midpoint / protein_len) * N_BINS) + 1, N_BINS)

            records.append({
                "residue":     residue,
                "bin":         bin_num,
                "is_terminal": bin_num in TERMINAL_BINS,
                "location":    "terminal" if bin_num in TERMINAL_BINS else "internal",
            })
    return records


def enrichment_ratio(df: pd.DataFrame) -> pd.DataFrame:
    """Per-residue terminal/internal frequency ratio."""
    counts = df.groupby(["location", "residue"]).size().reset_index(name="count")
    totals = counts.groupby("location")["count"].transform("sum")
    counts["freq"] = counts["count"] / totals
    pivot = counts.pivot(index="residue", columns="location", values="freq").fillna(0)
    pivot["ratio"] = pivot.get("terminal", 0) / (pivot.get("internal", 0) + 1e-9)
    return pivot.sort_values("ratio", ascending=False)


def main():
    manifest_path = RESULTS_DIR / "species_manifest.tsv"
    if not manifest_path.exists():
        print(f"ERROR: {manifest_path} missing — run 01_download_proteomes.py first.")
        return

    manifest = pd.read_csv(manifest_path, sep="\t")
    print(f"Loaded manifest: {len(manifest)} species\n")

    # ── Collect per-residue records by domain ─────────────────────────────────
    all_rows = []
    for _, row in manifest.iterrows():
        phylum = row["phylum"]
        domain = classify_domain(phylum)
        flps_f = FLPS_DIR / f"{row['species_key']}.flps.txt"
        if not flps_f.exists():
            continue
        recs = parse_flps_residue(flps_f)
        for r in recs:
            r["domain"] = domain
        all_rows.extend(recs)
        print(f"  [{domain}] {row['species_key']}: {len(recs)} LCRs")

    if not all_rows:
        print("No records found. Check fLPS output files.")
        return

    df = pd.DataFrame(all_rows)
    domains = ["Bacteria", "Archaea", "Viridiplantae", "Metazoa", "Other_Euk"]

    # ── Build Supp Table S7: per-domain × residue × location ─────────────────
    counts = (
        df.groupby(["domain", "location", "residue"])
        .size()
        .reset_index(name="count")
    )
    totals = counts.groupby(["domain", "location"])["count"].transform("sum")
    counts["freq"] = (counts["count"] / totals).round(4)

    # Add terminal / internal ratio per domain × residue
    ratio_rows = []
    for domain in domains:
        sub = df[df["domain"] == domain]
        if sub.empty:
            continue
        rat = enrichment_ratio(sub)
        for residue, row_ in rat.iterrows():
            ratio_rows.append({
                "domain":   domain,
                "residue":  residue,
                "freq_terminal": round(row_.get("terminal", 0), 4),
                "freq_internal": round(row_.get("internal", 0), 4),
                "ratio_terminal_over_internal": round(row_["ratio"], 3),
            })

    ratio_df = pd.DataFrame(ratio_rows)
    out_tsv = RESULTS_DIR / "domain_aa_composition.tsv"
    ratio_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nDomain AA composition table (Supp Table S7): {out_tsv}")

    # ── Console summary: top-5 enriched residues per domain ──────────────────
    print("\nTop 5 residues enriched at termini per domain (terminal/internal ratio):")
    for domain in domains:
        sub = ratio_df[ratio_df["domain"] == domain].nlargest(5, "ratio_terminal_over_internal")
        if sub.empty:
            continue
        residues = ", ".join(
            f"{r['residue']}({r['ratio_terminal_over_internal']:.2f})"
            for _, r in sub.iterrows()
        )
        print(f"  {domain:<18}: {residues}")

    # ── Figure: heatmap of ratio per domain × residue ─────────────────────────
    pivot = ratio_df.pivot(index="domain", columns="residue",
                           values="ratio_terminal_over_internal").fillna(1.0)
    # Reorder rows and columns
    pivot = pivot.reindex(
        [d for d in domains if d in pivot.index]
    )[[aa for aa in AA_ORDER if aa in pivot.columns]]

    fig, ax = plt.subplots(figsize=(14, 4))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                   vmin=0.5, vmax=2.0)
    plt.colorbar(im, ax=ax, label="terminal / internal frequency ratio")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.axvline(x=-0.5, color="white", linewidth=0)
    ax.set_title(
        "Amino acid enrichment at termini by domain\n"
        "(ratio > 1 = enriched at termini; ratio < 1 = depleted)", fontsize=12
    )
    # Mark ratio=1 contour
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            text = f"{val:.2f}"
            colour = "white" if (val > 1.6 or val < 0.7) else "black"
            ax.text(j, i, text, ha="center", va="center",
                    fontsize=7, color=colour)

    plt.tight_layout()
    out_fig = FIGURES_DIR / "suppfig_domain_aa.pdf"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    fig.savefig(str(out_fig).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nDomain AA composition figure: {out_fig}")
    plt.close()

    # ── Side-by-side barplot: top enriched residues per domain ───────────────
    n_domains = len([d for d in domains if d in ratio_df["domain"].values])
    fig, axes = plt.subplots(1, n_domains, figsize=(4 * n_domains, 5), sharey=False)
    if n_domains == 1:
        axes = [axes]

    for ax, domain in zip(axes, [d for d in domains if d in ratio_df["domain"].values]):
        sub = ratio_df[ratio_df["domain"] == domain].set_index("residue")
        sub = sub.reindex([aa for aa in AA_ORDER if aa in sub.index])
        colours = [AA_COLOURS.get(aa, "#ccc") for aa in sub.index]
        ax.barh(
            range(len(sub)), sub["ratio_terminal_over_internal"],
            color=colours, edgecolor="black", linewidth=0.5
        )
        ax.axvline(1.0, color="red", linestyle=":", linewidth=1.2)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub.index)
        ax.set_xlabel("Term / internal ratio")
        ax.set_title(domain, fontsize=11)

    plt.suptitle("Terminal LCR amino acid enrichment per domain\n"
                 "(red line = no enrichment; ratio > 1 = enriched at termini)",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    out_fig2 = FIGURES_DIR / "suppfig_domain_aa_barplot.pdf"
    fig.savefig(out_fig2, dpi=300, bbox_inches="tight")
    fig.savefig(str(out_fig2).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Domain AA barplot: {out_fig2}")
    plt.close()


if __name__ == "__main__":
    main()
