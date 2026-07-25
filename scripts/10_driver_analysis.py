#!/usr/bin/env python3
"""
Multiple-LCR protein driver analysis.

Tests whether the terminal enrichment signal is driven by a specific class of
proteins — those carrying multiple LCRs (potential disordered scaffolds,
transcription factors, stress-response proteins).

Rationale: if multi-LCR proteins disproportionately place LCRs at termini, the
enrichment is class-specific; if singleton-LCR proteins show equal enrichment, the
mechanism is generic across the proteome.

Outputs:
  results/driver_analysis.tsv   — singleton vs multi-LCR enrichment per phylum
  figures/fig9_driver.pdf
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu
from pathlib import Path

from config import (
    PROJECT_DIR, RESULTS_DIR, FIGURES_DIR,
    N_BINS, NULL_TERMINAL,
    PHYLUM_ORDER, PHYLUM_COLOURS,
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

    # Count LCRs per protein per species (to classify singleton vs multi)
    lcr_counts = (
        pos_df.groupby(["species_key", "protein_id"])
        .size()
        .reset_index(name="n_lcr_in_protein")
    )
    pos_df = pos_df.merge(lcr_counts, on=["species_key", "protein_id"])
    pos_df["protein_class"] = pos_df["n_lcr_in_protein"].apply(
        lambda n: "singleton" if n == 1 else "multi"
    )

    # ── Global summary ────────────────────────────────────────────────────────
    print("Global LCR class breakdown:")
    for cls in ["singleton", "multi"]:
        sub = pos_df[pos_df["protein_class"] == cls]
        n_t = sub["is_terminal"].sum()
        stats = enrich_stats(n_t, len(sub))
        print(f"  {cls:12s}: {len(sub):>6} LCRs  terminal={stats['pct_terminal']}%  "
              f"OR={stats['odds_ratio']}  p={stats['pvalue']}  "
              f"{'***' if stats['sig'] else 'ns'}")

    # Compare pct_terminal between singleton and multi at species level
    sp_rows = []
    for (sp_key, display, phylum), grp in pos_df.groupby(
        ["species_key", "display_name", "phylum"]
    ):
        for cls in ["singleton", "multi"]:
            sub = grp[grp["protein_class"] == cls]
            if len(sub) < 5:
                continue
            n_t = sub["is_terminal"].sum()
            stats = enrich_stats(n_t, len(sub))
            sp_rows.append({
                "species_key": sp_key, "display_name": display,
                "phylum": phylum, "protein_class": cls,
                "n_lcr": len(sub), **stats,
            })

    sp_df = pd.DataFrame(sp_rows)

    # Mann-Whitney: is multi enrichment > singleton?
    sing_pct = sp_df[sp_df["protein_class"] == "singleton"]["pct_terminal"].dropna()
    mult_pct = sp_df[sp_df["protein_class"] == "multi"]["pct_terminal"].dropna()
    if len(sing_pct) >= 5 and len(mult_pct) >= 5:
        stat, p_mwu = mannwhitneyu(mult_pct, sing_pct, alternative="greater")
        print(f"\nMann-Whitney U (multi > singleton pct_terminal):")
        print(f"  multi mean={mult_pct.mean():.1f}%  singleton mean={sing_pct.mean():.1f}%")
        print(f"  U={stat:.0f}, p={p_mwu:.4f}  "
              f"{'→ Multi-LCR proteins drive enrichment ***' if p_mwu < 0.05 else '→ Both classes contribute equally'}")

    # ── Per-phylum breakdown ──────────────────────────────────────────────────
    rows = []
    print(f"\n{'Phylum':<22}  singleton%   multi%   driver")
    for phylum in PHYLUM_ORDER:
        sub_ph = pos_df[pos_df["phylum"] == phylum]
        if sub_ph.empty:
            continue
        row = {"phylum": phylum}
        for cls in ["singleton", "multi"]:
            sub = sub_ph[sub_ph["protein_class"] == cls]
            if len(sub) < 5:
                continue
            n_t = sub["is_terminal"].sum()
            stats = enrich_stats(n_t, len(sub))
            row[f"pct_terminal_{cls}"]  = stats["pct_terminal"]
            row[f"odds_ratio_{cls}"]    = stats["odds_ratio"]
            row[f"pvalue_{cls}"]        = stats["pvalue"]
            row[f"sig_{cls}"]           = stats["sig"]
        rows.append(row)
        s_pct = row.get("pct_terminal_singleton", float("nan"))
        m_pct = row.get("pct_terminal_multi", float("nan"))
        driver = ("multi" if not np.isnan(m_pct) and not np.isnan(s_pct) and m_pct > s_pct + 2
                  else "shared" if not np.isnan(m_pct) and not np.isnan(s_pct)
                  else "n/a")
        print(f"  {phylum:<22}  {s_pct:>6.1f}%    {m_pct:>6.1f}%   {driver}")

    out_df = pd.DataFrame(rows)
    out_tsv = RESULTS_DIR / "driver_analysis.tsv"
    out_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nDriver analysis table: {out_tsv}")

    # ── Figure 9: singleton vs multi per phylum ───────────────────────────────
    phyla_present = [p for p in PHYLUM_ORDER if p in out_df["phylum"].values]
    x = np.arange(len(phyla_present))
    w = 0.35

    out_plot = out_df.set_index("phylum").reindex(phyla_present)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(x - w / 2, out_plot["pct_terminal_singleton"], w,
           label="Singleton-LCR proteins", color="#4393c3",
           edgecolor="black", linewidth=0.7)
    ax.bar(x + w / 2, out_plot["pct_terminal_multi"], w,
           label="Multi-LCR proteins (≥2 LCRs)", color="#d6604d",
           edgecolor="black", linewidth=0.7)

    # Species-level dots
    for i, phylum in enumerate(phyla_present):
        for j, cls in enumerate(["singleton", "multi"]):
            sub = sp_df[(sp_df["phylum"] == phylum) &
                        (sp_df["protein_class"] == cls)].dropna(subset=["pct_terminal"])
            xpos = i + (j - 0.5) * w
            ax.scatter(
                np.full(len(sub), xpos), sub["pct_terminal"],
                color="black", s=18, zorder=5, alpha=0.7
            )

    ax.axhline(10.0, color="red", linestyle=":", linewidth=1.2,
               label="Null expectation (10%)")
    ax.set_xticks(x)
    ax.set_xticklabels(phyla_present, rotation=22, ha="right", fontsize=9)
    ax.set_ylabel("% LCRs in terminal bins", fontsize=11)
    ax.set_title("Terminal LCR enrichment: singleton vs multi-LCR proteins\n"
                 "(black dots = individual species)", fontsize=12)
    ax.legend(fontsize=9)
    ax.set_ylim(0, max(
        out_plot["pct_terminal_singleton"].dropna().max(),
        out_plot["pct_terminal_multi"].dropna().max()
    ) * 1.2 + 2)

    # Annotate MWU result
    if len(sing_pct) >= 5 and len(mult_pct) >= 5:
        ax.text(0.99, 0.97,
                f"multi vs singleton MWU p={p_mwu:.3f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    plt.tight_layout()
    out = FIGURES_DIR / "fig9_driver.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Figure 9 saved: {out}")
    plt.close()


if __name__ == "__main__":
    main()
