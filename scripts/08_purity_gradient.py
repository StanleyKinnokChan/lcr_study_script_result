#!/usr/bin/env python3
"""
Purity gradient analysis.

Tests whether terminal LCRs have higher compositional purity than internal ones.
Higher purity at termini would indicate stronger selection for amino acid homogeneity
at protein ends — a distinct functional signal beyond mere positional enrichment.

Also tests whether the purity difference varies across phyla.

Outputs:
  results/purity_gradient.tsv  — per-phylum Mann-Whitney U results
  figures/fig7_purity_gradient.pdf
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from pathlib import Path

from config import (
    PROJECT_DIR, RESULTS_DIR, FIGURES_DIR,
    PHYLUM_ORDER, PHYLUM_COLOURS,
)

def main():
    pos_path = RESULTS_DIR / "lcr_positions.tsv"
    if not pos_path.exists():
        print(f"ERROR: {pos_path} missing — run 03_analyse.py first.")
        return

    pos_df = pd.read_csv(pos_path, sep="\t")

    # ── Global test ───────────────────────────────────────────────────────────
    term   = pos_df[pos_df["is_terminal"]]["purity"]
    intern = pos_df[~pos_df["is_terminal"]]["purity"]
    stat, p = mannwhitneyu(term, intern, alternative="greater")
    print(f"Global Mann-Whitney U (terminal > internal purity):")
    print(f"  Terminal  mean purity = {term.mean():.4f} ± {term.std():.4f}")
    print(f"  Internal  mean purity = {intern.mean():.4f} ± {intern.std():.4f}")
    print(f"  U={stat:.0f}, p={p:.2e}  {'***' if p < 0.05 else 'ns'}\n")

    # ── Per-phylum ────────────────────────────────────────────────────────────
    rows = []
    print(f"{'Phylum':<22}  term_mean  int_mean  delta   U_p")
    for phylum in PHYLUM_ORDER:
        sub = pos_df[pos_df["phylum"] == phylum]
        if sub.empty:
            continue
        t = sub[sub["is_terminal"]]["purity"]
        i = sub[~sub["is_terminal"]]["purity"]
        if len(t) < 5 or len(i) < 5:
            continue
        stat_ph, p_ph = mannwhitneyu(t, i, alternative="greater")
        delta = t.mean() - i.mean()
        rows.append({
            "phylum":       phylum,
            "n_terminal":   len(t),
            "n_internal":   len(i),
            "mean_purity_terminal": round(t.mean(), 4),
            "mean_purity_internal": round(i.mean(), 4),
            "delta_purity": round(delta, 4),
            "mwu_stat":     round(stat_ph, 0),
            "pvalue":       round(p_ph, 6),
            "sig":          p_ph < 0.05,
        })
        print(f"  {phylum:<20}  {t.mean():.4f}     {i.mean():.4f}   "
              f"{delta:+.4f}  {'***' if p_ph < 0.05 else 'ns'}")

    out_df = pd.DataFrame(rows)
    out_tsv = RESULTS_DIR / "purity_gradient.tsv"
    out_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nPurity gradient table: {out_tsv}")

    phyla_present = [p for p in PHYLUM_ORDER if p in out_df["phylum"].values]

    # ── Figure 7a: Violin plot global terminal vs internal ────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    data_violin = [
        pos_df[pos_df["is_terminal"]]["purity"].values,
        pos_df[~pos_df["is_terminal"]]["purity"].values,
    ]
    parts = ax.violinplot(data_violin, positions=[1, 2], showmedians=True,
                          showextrema=False)
    parts["bodies"][0].set_facecolor("#d6604d")
    parts["bodies"][0].set_alpha(0.7)
    parts["bodies"][1].set_facecolor("#4393c3")
    parts["bodies"][1].set_alpha(0.7)
    parts["cmedians"].set_color("black")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Terminal\n(bins 1 & 20)", "Internal\n(bins 2–19)"])
    ax.set_ylabel("LCR purity (dominant aa / LCR length)", fontsize=11)
    ax.set_title(f"Terminal LCRs have higher purity\n(Mann-Whitney U, p={p:.2e})", fontsize=11)
    ax.set_ylim(0.65, 1.05)

    # ── Figure 7b: delta purity per phylum ───────────────────────────────────
    ax2 = axes[1]
    colours = [PHYLUM_COLOURS.get(p, "#aaa") for p in phyla_present]
    out_present = out_df.set_index("phylum").reindex(phyla_present)

    bars = ax2.bar(
        range(len(phyla_present)),
        out_present["delta_purity"],
        color=colours, edgecolor="black", linewidth=0.7
    )
    ax2.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax2.set_xticks(range(len(phyla_present)))
    ax2.set_xticklabels(phyla_present, rotation=60, ha="right", fontsize=9)
    ax2.set_ylabel("Δ purity (terminal − internal)", fontsize=11)
    ax2.set_title("Per-phylum terminal purity excess\n"
                  "(positive = terminal LCRs are purer)", fontsize=11)

    # Significance markers
    for i, phylum in enumerate(phyla_present):
        row_ = out_present.loc[phylum]
        if row_["sig"]:
            ypos = row_["delta_purity"] + 0.0008 if row_["delta_purity"] >= 0 else row_["delta_purity"] - 0.002
            ax2.text(i, ypos, "***", ha="center", va="bottom", fontsize=8)

    plt.suptitle("Terminal LCR purity gradient across metazoan phyla\n"
                 "(purity = fraction of LCR composed of dominant amino acid)", fontsize=12, y=1.01)
    plt.tight_layout()

    out = FIGURES_DIR / "fig7_purity_gradient.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 7 saved: {out}")
    plt.close()


if __name__ == "__main__":
    main()
