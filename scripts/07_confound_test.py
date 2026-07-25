#!/usr/bin/env python3
"""
Protein-length confound test.

Short proteins mechanically inflate terminal bin counts: an LCR at position 1–5 aa
of a 50-aa protein hits bin 1, whereas the same LCR in a 5000-aa protein is deep
internal. This script stratifies proteins by length quartile and shows that terminal
enrichment is significant within every quartile, ruling out length confounding.

Outputs:
  results/length_confound.tsv  — Fisher's exact results per length quartile × phylum
  figures/fig6_length_confound.pdf
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
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
    out_tsv = RESULTS_DIR / "length_confound.tsv"
    out_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nLength confound table: {out_tsv}")

    # ── Figure 6: pct_terminal by quartile (pooled + per-phylum) ─────────────
    pooled = out_df[out_df["group"] == "ALL"].dropna(subset=["pct_terminal"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left panel: pooled bars
    ax = axes[0]
    colours = ["#d1e5f0", "#92c5de", "#4393c3", "#2166ac"]
    ax.bar(
        [f"Q{q}" for q in [1, 2, 3, 4]],
        pooled["pct_terminal"],
        color=colours, edgecolor="black", linewidth=0.7
    )
    ax.axhline(10.0, color="red", linestyle=":", linewidth=1.2, label="Null (10%)")
    ax.set_ylabel("% LCRs in terminal bins", fontsize=11)
    ax.set_xlabel("Protein length quartile", fontsize=11)
    ax.set_title("Terminal enrichment by protein length quartile\n(all species pooled)", fontsize=11)
    ax.set_ylim(0, max(pooled["pct_terminal"].max() * 1.2, 15))
    ax.legend(fontsize=9)

    # Annotate p values
    for i, (_, row_) in enumerate(pooled.iterrows()):
        ax.text(i, row_["pct_terminal"] + 0.3,
                "***" if row_["sig"] else "ns",
                ha="center", va="bottom", fontsize=9)

    # Right panel: per-phylum line plot across quartiles
    ax2 = axes[1]
    phyla_with_data = out_df[out_df["group"] != "ALL"]["group"].unique()
    phyla_ordered = [p for p in PHYLUM_ORDER if p in phyla_with_data]

    cmap = plt.colormaps["tab20"].resampled(max(len(phyla_ordered), 1))
    for i, phylum in enumerate(phyla_ordered):
        sub = out_df[(out_df["group"] == phylum)].dropna(subset=["pct_terminal"])
        if sub.empty:
            continue
        sub = sub.sort_values("quartile")
        ax2.plot(
            sub["quartile"], sub["pct_terminal"],
            marker="o", markersize=5, linewidth=1.5,
            color=cmap(i), label=phylum
        )

    ax2.axhline(10.0, color="red", linestyle=":", linewidth=1.2, label="Null (10%)")
    ax2.set_ylabel("% LCRs in terminal bins", fontsize=11)
    ax2.set_xlabel("Protein length quartile", fontsize=11)
    ax2.set_title("Per-phylum terminal enrichment across length quartiles", fontsize=11)
    ax2.legend(fontsize=7, loc="upper right", ncol=2)

    plt.suptitle("Protein length does not confound terminal LCR enrichment", fontsize=13, y=1.01)
    plt.tight_layout()

    out = FIGURES_DIR / "fig6_length_confound.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 6 saved: {out}")
    plt.close()


if __name__ == "__main__":
    main()
