#!/usr/bin/env python3
"""
N-terminal vs C-terminal asymmetry analysis.

Teekas et al. (2024) reported terminal enrichment as a combined bin-1+bin-20 metric.
This script separates the two termini and tests each independently against a null of
5% (1/20 bins), then asks whether the two termini are enriched equally.

Outputs:
  results/asymmetry.tsv  — per-species N/C enrichment
  figures/fig4_asymmetry.pdf
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, wilcoxon
from pathlib import Path

from config import (
    PROJECT_DIR, RESULTS_DIR, FIGURES_DIR,
    N_BINS, NULL_SINGLE,
    PHYLUM_ORDER, PHYLUM_COLOURS,
)


def asymmetry_stats(df_sp: pd.DataFrame) -> dict:
    total = len(df_sp)
    if total == 0:
        return {}

    n_nterm = (df_sp["bin"] == 1).sum()
    n_cterm = (df_sp["bin"] == N_BINS).sum()

    def _fisher(n_term, n_total):
        expected_term = n_total * NULL_SINGLE
        expected_other = n_total - expected_term
        table = [[n_term, n_total - n_term], [expected_term, expected_other]]
        or_, p = fisher_exact(table, alternative="greater")
        return round(or_, 3), round(p, 6)

    nterm_or, nterm_p = _fisher(n_nterm, total)
    cterm_or, cterm_p = _fisher(n_cterm, total)

    return {
        "n_lcr": total,
        "n_nterm": int(n_nterm),
        "n_cterm": int(n_cterm),
        "pct_nterm": round(n_nterm / total * 100, 2),
        "pct_cterm": round(n_cterm / total * 100, 2),
        "nterm_or": nterm_or,
        "nterm_p": nterm_p,
        "nterm_sig": nterm_p < 0.05,
        "cterm_or": cterm_or,
        "cterm_p": cterm_p,
        "cterm_sig": cterm_p < 0.05,
        "asymmetry_ratio":     round(n_nterm / n_cterm, 3) if n_cterm > 0 else None,
        "log_odds_asymmetry":  round(
            float(np.log((n_nterm + 0.5) / (n_cterm + 0.5))), 4
        ),
        "low_conf_asymmetry":  n_cterm < 10,
    }


def main():
    pos_path = RESULTS_DIR / "lcr_positions.tsv"
    if not pos_path.exists():
        print(f"ERROR: {pos_path} missing — run 03_analyse.py first.")
        return

    pos_df = pd.read_csv(pos_path, sep="\t")
    rows = []

    for (sp_key, display, phylum), grp in pos_df.groupby(
        ["species_key", "display_name", "phylum"]
    ):
        stats = asymmetry_stats(grp)
        if not stats:
            continue
        rows.append({"species_key": sp_key, "display_name": display,
                     "phylum": phylum, **stats})

    asym_df = pd.DataFrame(rows)
    out_tsv = RESULTS_DIR / "asymmetry.tsv"
    asym_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Asymmetry table: {out_tsv}  ({len(asym_df)} species)\n")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("Species  N-term%  C-term%  ratio  N-sig  C-sig")
    for _, r in asym_df.iterrows():
        print(f"  {r.display_name:<40} {r.pct_nterm:>6.1f}%  "
              f"{r.pct_cterm:>6.1f}%  {r.asymmetry_ratio or 'NA':>6}  "
              f"{'***' if r.nterm_sig else 'ns ':>4}   {'***' if r.cterm_sig else 'ns'}")

    # Wilcoxon signed-rank: is N-term enrichment different from C-term?
    paired = asym_df.dropna(subset=["pct_nterm", "pct_cterm"])
    if len(paired) >= 5:
        stat, p = wilcoxon(paired["pct_nterm"], paired["pct_cterm"])
        print(f"\nWilcoxon signed-rank (N-term vs C-term): stat={stat:.2f}, p={p:.4f}")
        print("  → N- and C-termini are " +
              ("asymmetrically enriched" if p < 0.05 else "symmetrically enriched"))

    # Phylum-level averages
    phylum_asym = (
        asym_df.groupby("phylum")[["pct_nterm", "pct_cterm"]]
        .mean()
        .reindex([p for p in PHYLUM_ORDER if p in asym_df["phylum"].values])
    )
    print(f"\nPhylum averages:\n{phylum_asym.round(2).to_string()}")

    # ── Figure 4: N vs C enrichment per phylum ───────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    phyla = phylum_asym.index.tolist()
    x = np.arange(len(phyla))
    w = 0.35

    bars_n = ax.bar(
        x - w / 2, phylum_asym["pct_nterm"], w,
        label="N-terminal (bin 1)", color="#4393c3", edgecolor="black", linewidth=0.7
    )
    bars_c = ax.bar(
        x + w / 2, phylum_asym["pct_cterm"], w,
        label="C-terminal (bin 20)", color="#d6604d", edgecolor="black", linewidth=0.7
    )

    # Species-level dots. Small markers plus a narrow horizontal jitter so that
    # densely sampled phyla (Bacteria, Insecta, Viridiplantae) resolve into
    # individual species rather than a solid black column.
    rng = np.random.default_rng(0)
    jitter = w * 0.30
    for i, phylum in enumerate(phyla):
        sub = asym_df[asym_df["phylum"] == phylum]
        for offset, col in ((-w / 2, "pct_nterm"), (w / 2, "pct_cterm")):
            ax.scatter(
                i + offset + rng.uniform(-jitter, jitter, len(sub)), sub[col],
                color="black", s=4, zorder=5, alpha=0.45, linewidths=0
            )

    ax.axhline(5.0, color="red", linestyle=":", linewidth=1.2,
               label="Null (1/20 bins = 5%)")
    ax.set_xticks(x)
    ax.set_xticklabels(phyla, rotation=60, ha="right")
    ax.set_ylabel("% LCRs in terminal bin", fontsize=11)
    ax.set_title("N-terminal vs C-terminal LCR enrichment per phylum\n"
                 "(black dots = individual species)", fontsize=12)
    ax.legend(fontsize=9)
    plt.tight_layout()

    out = FIGURES_DIR / "fig4_asymmetry.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nFigure 4 saved: {out}")
    plt.close()


if __name__ == "__main__":
    main()
