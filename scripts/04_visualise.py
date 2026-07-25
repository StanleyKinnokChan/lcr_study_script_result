#!/usr/bin/env python3
"""
Produce two publication-ready figures:

Figure 1 — Bin heatmap
  Rows = species grouped by phylum, columns = bins 1–20.
  Cell colour = fraction of LCRs in that bin.
  Highlights terminal enrichment pattern visually.

Figure 2 — Terminal enrichment bar chart
  % terminal LCRs per phylum (invertebrates) with Teekas Tetrapoda range shaded.
  The headline comparison figure for the paper.
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
from pathlib import Path

from config import (
    PROJECT_DIR, RESULTS_DIR, FIGURES_DIR,
    N_BINS,
    TETRAPODA_TERMINAL_LOW, TETRAPODA_TERMINAL_HIGH,
    PHYLUM_ORDER, PHYLUM_COLOURS,
)


def fig1_bin_heatmap(pos_df: pd.DataFrame):
    """Heatmap: species × bin, colour = fraction of LCRs in that bin."""
    # Pivot: fraction per bin per species
    species_order = (
        pos_df.groupby(["phylum", "display_name"])
        .size()
        .reset_index()
        .sort_values(["phylum", "display_name"])
        ["display_name"]
        .tolist()
    )

    pivot = (
        pos_df.groupby(["display_name", "bin"])
        .size()
        .unstack(fill_value=0)
    )
    # Ensure all 20 bins present
    for b in range(1, N_BINS + 1):
        if b not in pivot.columns:
            pivot[b] = 0
    pivot = pivot[sorted(pivot.columns)]
    # Normalise to fraction per species
    pivot = pivot.div(pivot.sum(axis=1), axis=0)
    pivot = pivot.reindex([s for s in species_order if s in pivot.index])

    fig, ax = plt.subplots(figsize=(12, max(4, len(pivot) * 0.55)))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="YlOrRd",
        vmin=0,
        vmax=0.15,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Fraction of LCRs in bin"},
    )
    ax.set_xlabel("Normalized protein position (bin, 1=N-term, 20=C-term)", fontsize=11)
    ax.set_ylabel("")
    ax.set_title("LCR positional distribution across invertebrate proteomes\n"
                 f"(purity ≥ 70%, bins 1 & {N_BINS} are terminal)", fontsize=12)
    ax.axvline(1, color="blue", linewidth=1.5, linestyle="--", alpha=0.6)
    ax.axvline(N_BINS, color="blue", linewidth=1.5, linestyle="--", alpha=0.6)

    plt.tight_layout()
    out = FIGURES_DIR / "fig1_bin_heatmap.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 1 saved: {out}")
    plt.close()


def fig2_terminal_barchart(enr_df: pd.DataFrame, phylum_df: pd.DataFrame):
    """Bar chart: % terminal LCRs per phylum with Tetrapoda range shaded."""
    phylum_df = phylum_df.copy()
    # Drop any phyla absent from PHYLUM_ORDER (avoids NaN category → Matplotlib TypeError)
    unknown = set(phylum_df["phylum"].astype(str)) - set(PHYLUM_ORDER)
    if unknown:
        print(f"  [WARN] fig2: unknown phyla dropped from bar chart: {sorted(unknown)}")
    phylum_df = phylum_df[phylum_df["phylum"].isin(PHYLUM_ORDER)].copy()
    ordered_cats = [p for p in PHYLUM_ORDER if p in phylum_df["phylum"].values]
    phylum_df["phylum"] = pd.Categorical(
        phylum_df["phylum"].astype(str),
        categories=ordered_cats,
        ordered=True,
    )
    phylum_df = phylum_df.sort_values("phylum").dropna(subset=["phylum", "pct_terminal"])

    fig, ax = plt.subplots(figsize=(8, 5))

    # Tetrapoda reference band
    ax.axhspan(
        TETRAPODA_TERMINAL_LOW,
        TETRAPODA_TERMINAL_HIGH,
        color="grey",
        alpha=0.18,
        label=f"Tetrapoda range (Teekas 2024): {TETRAPODA_TERMINAL_LOW}–{TETRAPODA_TERMINAL_HIGH}%",
    )

    bars = ax.bar(
        phylum_df["phylum"],
        phylum_df["pct_terminal"],
        color=[PHYLUM_COLOURS.get(p, "#aaaaaa") for p in phylum_df["phylum"]],
        edgecolor="black",
        linewidth=0.7,
        width=0.6,
    )

    # Species-level dots overlaid
    phylum_list = list(phylum_df["phylum"])
    for _, row in enr_df.iterrows():
        if pd.isna(row["pct_terminal"]):
            continue
        if row["phylum"] not in phylum_list:
            continue
        ax.scatter(
            phylum_list.index(row["phylum"]),
            row["pct_terminal"],
            color="black",
            s=22,
            zorder=5,
            alpha=0.7,
        )

    # Null expectation line (2/20 = 10%)
    ax.axhline(10.0, color="red", linestyle=":", linewidth=1.2,
               label="Null expectation (2/20 bins = 10%)")

    ax.set_ylabel("% LCRs in terminal bins (bins 1 & 20)", fontsize=11)
    ax.set_xlabel("Phylum", fontsize=11)
    ax.set_title("Terminal LCR enrichment in invertebrates vs. Tetrapoda\n"
                 "(black dots = individual species)", fontsize=12)
    ax.set_ylim(0, max(phylum_df["pct_terminal"].max() * 1.2, TETRAPODA_TERMINAL_HIGH + 5))
    ax.legend(fontsize=9, loc="upper right")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()

    out = FIGURES_DIR / "fig2_terminal_enrichment.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 2 saved: {out}")
    plt.close()


def fig3_bin_profile(pos_df: pd.DataFrame):
    """
    U-shaped bin profile per phylum — the visual evidence that terminal enrichment
    is consistent in shape, not just magnitude, across all phyla.
    Each line = fraction of LCRs in each bin, averaged across species in that phylum.
    """
    phyla = [p for p in PHYLUM_ORDER if p in pos_df["phylum"].values]
    fig, ax = plt.subplots(figsize=(10, 5))

    for phylum in phyla:
        sub = pos_df[pos_df["phylum"] == phylum]
        bin_counts = sub.groupby("bin").size()
        bin_frac = bin_counts / bin_counts.sum()
        # Ensure all 20 bins represented
        all_bins = pd.Series(0.0, index=range(1, N_BINS + 1))
        all_bins.update(bin_frac)
        ax.plot(
            all_bins.index,
            all_bins.values * 100,
            marker="o",
            markersize=3,
            linewidth=1.5,
            color=PHYLUM_COLOURS.get(phylum, "#888888"),
            label=phylum,
        )

    # Null expectation
    ax.axhline(100 / N_BINS, color="black", linestyle=":", linewidth=1,
               label="Null (uniform = 5%)")
    ax.axvline(1.5, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(N_BINS - 0.5, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)

    ax.set_xlabel("Normalized protein position (bin 1 = N-term, 20 = C-term)", fontsize=11)
    ax.set_ylabel("% of LCRs in bin", fontsize=11)
    ax.set_title("U-shaped LCR positional profile is conserved across all invertebrate phyla\n"
                 "(dashed lines mark terminal bins)", fontsize=12)
    ax.set_xticks(range(1, N_BINS + 1))
    ax.legend(fontsize=9, loc="upper center", ncol=4)
    plt.tight_layout()

    out = FIGURES_DIR / "fig3_bin_profile.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 3 saved: {out}")
    plt.close()


def main():
    pos_path    = RESULTS_DIR / "lcr_positions.tsv"
    enr_path    = RESULTS_DIR / "enrichment.tsv"
    phylum_path = RESULTS_DIR / "phylum_summary.tsv"

    for p in [pos_path, enr_path, phylum_path]:
        if not p.exists():
            print(f"ERROR: missing {p} — run 03_analyse.py first.")
            return

    pos_df    = pd.read_csv(pos_path,    sep="\t")
    enr_df    = pd.read_csv(enr_path,    sep="\t")
    phylum_df = pd.read_csv(phylum_path, sep="\t")

    print(f"Loaded {len(pos_df)} LCRs across {enr_df['species_key'].nunique()} species")
    fig1_bin_heatmap(pos_df)
    fig2_terminal_barchart(enr_df, phylum_df)
    fig3_bin_profile(pos_df)
    print("\nAll figures saved to:", FIGURES_DIR)


if __name__ == "__main__":
    main()
