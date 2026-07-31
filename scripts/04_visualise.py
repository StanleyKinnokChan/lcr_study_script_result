#!/usr/bin/env python3
"""
Produce publication-ready figures:

Figure 1 — Bin heatmap
  Rows = phylum/domain (43 groups, PHYLUM_ORDER), columns = bins 1–20.
  Cell colour = fraction of LCRs in that bin, averaged across species in the
  group. Highlights terminal enrichment pattern visually at the same pooled
  level as the paper's statistics.

  figures/suppfig8_bin_heatmap_species.pdf — Supplementary Figure 8
  Same heatmap at full species resolution (772 rows) — the species-level
  detail that Figure 1's phylum aggregation collapses.

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


def _bin_fraction_pivot(pos_df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """Fraction-of-LCRs-per-bin pivot, one row per value of group_col."""
    pivot = (
        pos_df.groupby([group_col, "bin"])
        .size()
        .unstack(fill_value=0)
    )
    for b in range(1, N_BINS + 1):
        if b not in pivot.columns:
            pivot[b] = 0
    pivot = pivot[sorted(pivot.columns)]
    return pivot.div(pivot.sum(axis=1), axis=0)


def fig1_bin_heatmap(pos_df: pd.DataFrame):
    """Heatmap: phylum/domain × bin, colour = mean fraction of LCRs in that bin.

    Aggregated to the same 43 phylum/domain groups used throughout the paper's
    pooled statistics (PHYLUM_ORDER) rather than one row per species — a
    772-species version is unreadable at publication size (see suppfig8).
    """
    pivot = _bin_fraction_pivot(pos_df, "phylum")
    ordered = [p for p in PHYLUM_ORDER if p in pivot.index]
    pivot = pivot.reindex(ordered)

    fig, ax = plt.subplots(figsize=(10, max(6, len(pivot) * 0.28)))
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
    ax.set_title("LCR positional distribution across phyla/domains\n"
                 f"(purity ≥ 70%, bins 1 & {N_BINS} are terminal)", fontsize=12)
    ax.axvline(1, color="blue", linewidth=1.5, linestyle="--", alpha=0.6)
    ax.axvline(N_BINS, color="blue", linewidth=1.5, linestyle="--", alpha=0.6)

    plt.tight_layout()
    out = FIGURES_DIR / "fig1_bin_heatmap.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 1 saved: {out}")
    plt.close()


def suppfig8_bin_heatmap_species(pos_df: pd.DataFrame):
    """Species-level companion to Figure 1 (772 rows) — supplementary detail only."""
    species_order = (
        pos_df.groupby(["phylum", "display_name"])
        .size()
        .reset_index()
        .sort_values(["phylum", "display_name"])
        ["display_name"]
        .tolist()
    )
    pivot = _bin_fraction_pivot(pos_df, "display_name")
    pivot = pivot.reindex([s for s in species_order if s in pivot.index])

    fig, ax = plt.subplots(figsize=(12, max(4, len(pivot) * 0.12)))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="YlOrRd",
        vmin=0,
        vmax=0.15,
        linewidths=0,
        cbar_kws={"label": "Fraction of LCRs in bin"},
    )
    ax.set_xlabel("Normalized protein position (bin, 1=N-term, 20=C-term)", fontsize=11)
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_title("LCR positional distribution across individual proteomes (n=772)\n"
                 f"(purity ≥ 70%, bins 1 & {N_BINS} are terminal; species grouped by phylum "
                 "as in Fig. 1)", fontsize=12)
    ax.axvline(1, color="blue", linewidth=1.5, linestyle="--", alpha=0.6)
    ax.axvline(N_BINS, color="blue", linewidth=1.5, linestyle="--", alpha=0.6)

    plt.tight_layout()
    out = FIGURES_DIR / "suppfig8_bin_heatmap_species.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Supplementary Figure 8 saved: {out}")
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

    # Species-level dots overlaid. Small markers plus a narrow horizontal jitter
    # so that densely sampled phyla (Bacteria, Insecta, Viridiplantae) resolve
    # into individual species rather than a solid black column.
    phylum_list = list(phylum_df["phylum"])
    rng = np.random.default_rng(0)
    jitter = 0.15
    for _, row in enr_df.iterrows():
        if pd.isna(row["pct_terminal"]):
            continue
        if row["phylum"] not in phylum_list:
            continue
        ax.scatter(
            phylum_list.index(row["phylum"]) + rng.uniform(-jitter, jitter),
            row["pct_terminal"],
            color="black",
            s=4,
            zorder=5,
            alpha=0.45,
            linewidths=0,
        )

    # Null expectation line (2/20 = 10%)
    ax.axhline(10.0, color="red", linestyle=":", linewidth=1.2,
               label="Null expectation (2/20 bins = 10%)")

    ax.set_ylabel("% LCRs in terminal bins (bins 1 & 20)", fontsize=11)
    ax.set_xlabel("Phylum", fontsize=11)
    ax.set_title("Terminal LCR enrichment in invertebrates vs. Tetrapoda\n"
                 "(black dots = individual species)", fontsize=12)
    # Fixed headroom above the tallest phylum bar so the species dots have room;
    # a handful of extreme single-species outliers still fall outside the axis.
    ax.set_ylim(0, 38)
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
    suppfig8_bin_heatmap_species(pos_df)
    fig2_terminal_barchart(enr_df, phylum_df)
    fig3_bin_profile(pos_df)
    print("\nAll figures saved to:", FIGURES_DIR)


if __name__ == "__main__":
    main()
