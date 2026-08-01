#!/usr/bin/env python3
"""
Purity gradient analysis.

Tests whether terminal LCRs have higher compositional purity than internal ones.
Higher purity at termini would indicate stronger selection for amino acid homogeneity
at protein ends — a distinct functional signal beyond mere positional enrichment.

Also tests whether the purity difference varies across phyla.

Outputs:
  results/supp_table_S5_purity_gradient.tsv  — Supp Table S5: per-phylum
                                                Mann-Whitney U results
  figures/suppfig6_purity_gradient.pdf  — Supplementary Figure 6 (the 9 phyla
                                           with a significant terminal/internal
                                           purity difference)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from pathlib import Path

from config import (
    PROJECT_DIR, RESULTS_DIR, FIGURES_DIR,
    PHYLUM_ORDER,
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
    out_tsv = RESULTS_DIR / "supp_table_S5_purity_gradient.tsv"
    out_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nPurity gradient table: {out_tsv}")

    # ── Supplementary Figure 6: terminal vs internal purity, 9 significant phyla ─
    NINE_SIG_PHYLA = [
        "Apicomplexa", "Rhizaria", "Euglenozoa", "Viridiplantae", "Rotifera",
        "Annelida", "Chelicerata", "Echinodermata", "Hemichordata",
    ]
    sig_df = out_df.set_index("phylum")
    order = [p for p in NINE_SIG_PHYLA if p in sig_df.index]

    sub9 = pos_df[pos_df["phylum"].isin(order)].copy()
    sub9["location"] = sub9["is_terminal"].map({True: "Terminal", False: "Internal"})

    fig, ax = plt.subplots(figsize=(12, 6.5))
    sns.violinplot(
        data=sub9, x="phylum", y="purity", hue="location", order=order,
        split=True, inner="quartile", ax=ax,
        palette={"Terminal": "#e41a1c", "Internal": "#377eb8"},
        cut=0, linewidth=0.8,
    )
    for i, phylum in enumerate(order):
        p_ph = sig_df.loc[phylum, "pvalue"]
        label = "p<0.001" if p_ph < 0.001 else f"p={p_ph:.2g}"
        ax.text(i, 1.03, label, ha="center", va="bottom", fontsize=8,
                transform=ax.get_xaxis_transform())
    ax.set_ylim(0.6, 1.14)
    ax.set_xlabel("")
    ax.set_ylabel("LCR purity (dominant aa / LCR length)", fontsize=11)
    ax.set_title("Terminal vs internal LCR purity, nine phyla with significant\n"
                 "Mann-Whitney U difference (bins 1+20 vs bins 2-19)", fontsize=12, y=1.08)
    ax.legend(title="", fontsize=9, loc="lower right")
    plt.tight_layout()

    out = FIGURES_DIR / "suppfig6_purity_gradient.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Supplementary Figure 6 saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
