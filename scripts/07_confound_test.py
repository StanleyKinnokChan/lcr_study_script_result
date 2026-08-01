#!/usr/bin/env python3
"""
Protein-length confound test.

Short proteins mechanically inflate terminal bin counts: an LCR at position 1–5 aa
of a 50-aa protein hits bin 1, whereas the same LCR in a 5000-aa protein is deep
internal. This script stratifies proteins by length quartile and shows that terminal
enrichment is significant within every quartile, ruling out length confounding.

Outputs:
  results/supp_table_S4_length_stratified.tsv  — Supp Table S4: Fisher's exact
                                                  results per length quartile × phylum
  figures/fig5_prokaryote_length_quartile.pdf  — Figure 5 (main text)
  figures/suppfig2_length_heatmap.pdf          — Supplementary Figure 2
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.proportion import proportion_confint
from pathlib import Path

from config import (
    PROJECT_DIR, RESULTS_DIR, FIGURES_DIR,
    N_BINS, NULL_TERMINAL,
    PHYLUM_ORDER,
    QUARTILE_LABELS,
)


def enrich_stats(n_terminal: int, total: int) -> dict:
    if total == 0:
        return {"pct_terminal": None, "odds_ratio": None, "pvalue": None, "sig": False}
    exp_t = total * NULL_TERMINAL
    exp_i = total - exp_t
    table = [[n_terminal, total - n_terminal], [exp_t, exp_i]]
    or_, p = fisher_exact(table, alternative="greater")
    return {
        "pct_terminal": round(n_terminal / total * 100, 2),
        "odds_ratio":   round(or_, 3),
        "pvalue":       round(p, 6),
        "sig":          p < 0.05,
    }


def main():
    pos_path = RESULTS_DIR / "lcr_positions.tsv"
    if not pos_path.exists():
        print(f"ERROR: {pos_path} missing — run 03_analyse.py first.")
        return

    pos_df = pd.read_csv(pos_path, sep="\t")

    # Assign length quartile labels globally (across all species)
    quartile_breaks = pos_df["protein_len"].quantile([0.25, 0.50, 0.75]).values
    q1, q2, q3 = quartile_breaks

    def _quartile(length):
        if length <= q1:
            return 1
        if length <= q2:
            return 2
        if length <= q3:
            return 3
        return 4

    pos_df["length_quartile"] = pos_df["protein_len"].apply(_quartile)

    print(f"Protein length quartile breaks: Q1={q1:.0f} aa, Q2={q2:.0f} aa, Q3={q3:.0f} aa\n")

    rows = []
    # ── Pooled across all species (main result) ───────────────────────────────
    print("Pooled (all species):")
    for q in [1, 2, 3, 4]:
        sub = pos_df[pos_df["length_quartile"] == q]
        n_term = sub["is_terminal"].sum()
        stats  = enrich_stats(n_term, len(sub))
        label  = f"Q{q}"
        rows.append({"group": "ALL", "quartile": label, "n_lcr": len(sub),
                     "n_terminal": n_term, **stats})
        print(f"  {label}: {n_term}/{len(sub)} terminal  "
              f"({stats['pct_terminal']}%)  OR={stats['odds_ratio']}  "
              f"p={stats['pvalue']}  {'***' if stats['sig'] else 'ns'}")

    # ── Per-phylum breakdown ──────────────────────────────────────────────────
    print("\nPer-phylum breakdown:")
    for phylum in PHYLUM_ORDER:
        sub_ph = pos_df[pos_df["phylum"] == phylum]
        if sub_ph.empty:
            continue
        print(f"  {phylum}:")
        for q in [1, 2, 3, 4]:
            sub = sub_ph[sub_ph["length_quartile"] == q]
            if len(sub) < 10:
                continue
            n_term = sub["is_terminal"].sum()
            stats  = enrich_stats(n_term, len(sub))
            rows.append({"group": phylum, "quartile": f"Q{q}", "n_lcr": len(sub),
                         "n_terminal": n_term, **stats})
            print(f"    Q{q}: {n_term}/{len(sub)} ({stats['pct_terminal']}%)  "
                  f"OR={stats['odds_ratio']}  p={stats['pvalue']}  "
                  f"{'***' if stats['sig'] else 'ns'}")

    out_df = pd.DataFrame(rows)
    out_tsv = RESULTS_DIR / "supp_table_S4_length_stratified.tsv"
    out_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nLength confound table: {out_tsv}")

    # ── Figure 5 (main text): prokaryote length-quartile enrichment ──────────
    prok = out_df[out_df["group"].isin(["Bacteria", "Archaea"])].copy()
    prok["quartile"] = pd.Categorical(prok["quartile"], ["Q1", "Q2", "Q3", "Q4"], ordered=True)
    prok = prok.sort_values(["group", "quartile"])
    lo, hi = proportion_confint(prok["n_terminal"], prok["n_lcr"], method="wilson")
    prok["ci_lo"] = prok["pct_terminal"] - lo * 100
    prok["ci_hi"] = hi * 100 - prok["pct_terminal"]

    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.arange(4)
    width = 0.35
    domain_colours = {"Bacteria": "#1f78b4", "Archaea": "#e31a1c"}
    for i, domain in enumerate(["Bacteria", "Archaea"]):
        d = prok[prok["group"] == domain]
        offset = (i - 0.5) * width
        ax.bar(
            x + offset, d["pct_terminal"], width, label=domain,
            color=domain_colours[domain], edgecolor="black", linewidth=0.7,
            yerr=[d["ci_lo"], d["ci_hi"]], capsize=4, error_kw={"linewidth": 1.1},
        )
    ax.axhline(10.0, color="grey", linestyle=":", linewidth=1.2, label="10% null")
    ax.set_xticks(x)
    ax.set_xticklabels(QUARTILE_LABELS)
    ax.set_ylabel("% LCRs in terminal bins (bins 1 & 20)", fontsize=11)
    ax.set_xlabel("Protein-length quartile (boundaries global, defined once per domain)", fontsize=10)
    ax.set_title("Protein-length stratified terminal LCR enrichment in prokaryotes", fontsize=12)
    ax.legend(fontsize=9)
    plt.tight_layout()

    out5 = FIGURES_DIR / "fig5_prokaryote_length_quartile.pdf"
    fig.savefig(out5, dpi=300, bbox_inches="tight")
    fig.savefig(str(out5).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 5 saved: {out5}")
    plt.close(fig)

    # ── Supplementary Figure 2: 4 quartiles x 43 groups heatmap ──────────────
    sub = out_df[out_df["group"] != "ALL"].copy()
    pivot = sub.pivot(index="group", columns="quartile", values="pct_terminal")
    pivot = pivot[["Q1", "Q2", "Q3", "Q4"]]
    order = [p for p in PHYLUM_ORDER if p in pivot.index]
    pivot = pivot.reindex(order)

    fig, ax = plt.subplots(figsize=(5, max(8, len(pivot) * 0.32)))
    sns.heatmap(
        pivot, ax=ax, cmap="YlOrRd", vmin=0, vmax=pivot.to_numpy(dtype=float).max(),
        linewidths=0.4, linecolor="white",
        cbar_kws={"label": "% LCRs in terminal bins"},
    )
    ax.set_xlabel("Protein-length quartile", fontsize=11)
    ax.set_ylabel("")
    ax.set_title("Length-stratified terminal LCR enrichment,\nall 43 phyla/groups", fontsize=12)
    plt.tight_layout()

    out2 = FIGURES_DIR / "suppfig2_length_heatmap.pdf"
    fig.savefig(out2, dpi=300, bbox_inches="tight")
    fig.savefig(str(out2).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Supplementary Figure 2 saved: {out2}")
    plt.close(fig)


if __name__ == "__main__":
    main()
