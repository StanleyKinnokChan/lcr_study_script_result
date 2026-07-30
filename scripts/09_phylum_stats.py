#!/usr/bin/env python3
"""
Between-phylum and within-phylum statistical analysis.

Tests:
  1. Kruskal-Wallis across all phyla — is phylum a significant predictor of
     terminal enrichment magnitude?
  2. Pairwise Mann-Whitney U between phyla (Holm-Bonferroni corrected).
  3. Within-phylum variance for phyla with ≥3 species — is enrichment
     consistent within a phylum?
  4. Nematode anomaly test — is Nematoda significantly higher than the rest?

Outputs:
  results/phylum_stats.tsv         — pairwise p-values matrix
  results/within_phylum_cv.tsv     — coefficient of variation per phylum
  figures/suppfig1_phylum_distribution.pdf  — Supplementary Figure 1
  figures/suppfig4_cv_vs_nspecies.pdf       — Supplementary Figure 4
"""

import itertools
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms
from scipy.stats import kruskal, mannwhitneyu
from pathlib import Path

from config import (
    PROJECT_DIR, RESULTS_DIR, FIGURES_DIR,
    PHYLUM_ORDER, PHYLUM_COLOURS,
)


def _place_labels_no_overlap(fig, ax, xs, ys, labels, fontsize=7, pad_px=2.0,
                              marker_size_pt2=55):
    """
    Annotate each (x, y) point with its label, picking the first candidate
    offset (in points, tried in increasing distance from the point) whose
    rendered bounding box doesn't overlap a already-placed label or any
    scatter marker. Falls back to the least-bad candidate if every option
    collides — dense clusters (e.g. several low-n phyla at similar CV) would
    otherwise stack labels on top of each other or under a neighbouring dot
    (whose higher z-order then paints over the label text) with a fixed
    offset.
    """
    fig.canvas.draw()  # establish a renderer so get_window_extent() works
    renderer = fig.canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer)  # plot area in display pixels

    candidates = [
        (6, 4, "left"), (6, -11, "left"), (-6, 4, "right"), (-6, -11, "right"),
        (6, 13, "left"), (-6, 13, "right"), (6, -20, "left"), (-6, -20, "right"),
        (16, 4, "left"), (-16, 4, "right"), (16, -11, "left"), (-16, -11, "right"),
        (0, 20, "center"), (0, -22, "center"),
        (26, 4, "left"), (-26, 4, "right"), (0, 30, "center"), (0, -32, "center"),
    ]

    # Every scatter marker is a permanent obstacle so a label never lands on
    # top of a (higher z-order) dot and gets visually clipped by it.
    marker_radius_px = ((marker_size_pt2 / 3.14159) ** 0.5) * fig.dpi / 72.0 + pad_px
    placed_boxes = []  # display-coordinate (x0, y0, x1, y1) of obstacles (markers + placed labels)
    for x, y in zip(xs, ys):
        cx, cy = ax.transData.transform((x, y))
        placed_boxes.append(matplotlib.transforms.Bbox(
            [[cx - marker_radius_px, cy - marker_radius_px],
             [cx + marker_radius_px, cy + marker_radius_px]]))

    def _overlaps(b1, b2):
        return not (b1.x1 + pad_px < b2.x0 or b2.x1 + pad_px < b1.x0 or
                    b1.y1 + pad_px < b2.y0 or b2.y1 + pad_px < b1.y0)

    def _inside_axes(bbox):
        return (bbox.x0 >= axes_box.x0 and bbox.x1 <= axes_box.x1 and
                bbox.y0 >= axes_box.y0 and bbox.y1 <= axes_box.y1)

    for x, y, label in zip(xs, ys, labels):
        best_dx, best_dy, best_ha = candidates[0][0], candidates[0][1], candidates[0][2]
        best_score = None  # (out_of_bounds: bool, n_overlap: int) — lower is better
        for dx, dy, ha in candidates:
            txt = ax.annotate(
                label, xy=(x, y), xytext=(dx, dy), textcoords="offset points",
                ha=ha, va="center", fontsize=fontsize,
            )
            bbox = txt.get_window_extent(renderer)
            out_of_bounds = not _inside_axes(bbox)
            n_overlap = sum(_overlaps(bbox, b) for b in placed_boxes)
            score = (out_of_bounds, n_overlap)
            if score == (False, 0):
                placed_boxes.append(bbox)
                best_score = score
                break
            txt.remove()
            if best_score is None or score < best_score:
                best_score = score
                best_dx, best_dy, best_ha = dx, dy, ha
        else:
            # Nothing both in-bounds and collision-free: keep the best-scoring candidate.
            txt = ax.annotate(
                label, xy=(x, y), xytext=(best_dx, best_dy),
                textcoords="offset points", ha=best_ha, va="center",
                fontsize=fontsize,
            )
            placed_boxes.append(txt.get_window_extent(renderer))


def holm_bonferroni(pvalues: list[float]) -> list[float]:
    """Holm-Bonferroni step-down correction with monotonicity enforcement."""
    n = len(pvalues)
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [0.0] * n
    running_max = 0.0
    for rank, (orig_idx, p) in enumerate(indexed):
        adj = min(p * (n - rank), 1.0)
        running_max = max(running_max, adj)   # step-down monotonicity
        adjusted[orig_idx] = running_max
    return adjusted


def main():
    enr_path = RESULTS_DIR / "enrichment.tsv"
    if not enr_path.exists():
        print(f"ERROR: {enr_path} missing — run 03_analyse.py first.")
        return

    enr_df = pd.read_csv(enr_path, sep="\t").dropna(subset=["pct_terminal"])
    phyla_present = [p for p in PHYLUM_ORDER if p in enr_df["phylum"].values]

    # ── 1. Kruskal-Wallis ────────────────────────────────────────────────────
    groups = [enr_df[enr_df["phylum"] == p]["pct_terminal"].values
              for p in phyla_present if len(enr_df[enr_df["phylum"] == p]) >= 2]
    if len(groups) >= 3:
        h, p_kw = kruskal(*groups)
        print(f"Kruskal-Wallis across phyla: H={h:.3f}, p={p_kw:.4f}")
        print(f"  → Phylum {'is' if p_kw < 0.05 else 'is NOT'} a significant "
              f"predictor of terminal enrichment magnitude\n")
    else:
        print("Not enough multi-species phyla for Kruskal-Wallis.\n")
        p_kw = None

    # ── 2. Pairwise Mann-Whitney U (Holm-corrected) ──────────────────────────
    # Only compare phyla with ≥2 species; single-species phyla cannot be tested.
    multi_phyla = [p for p in phyla_present
                   if len(enr_df[enr_df["phylum"] == p]) >= 2]
    single_phyla = [p for p in phyla_present if p not in multi_phyla]
    if single_phyla:
        print(f"Skipping single-species phyla from pairwise test: {', '.join(single_phyla)}")

    pairs = list(itertools.combinations(multi_phyla, 2))
    raw_ps = []
    pair_labels = []
    for p1, p2 in pairs:
        g1 = enr_df[enr_df["phylum"] == p1]["pct_terminal"].values
        g2 = enr_df[enr_df["phylum"] == p2]["pct_terminal"].values
        _, p = mannwhitneyu(g1, g2, alternative="two-sided")
        raw_ps.append(p)
        pair_labels.append((p1, p2))

    adj_ps = holm_bonferroni(raw_ps) if raw_ps else []
    sig_pairs = [(p1, p2, rp, ap)
                 for (p1, p2), rp, ap in zip(pair_labels, raw_ps, adj_ps)
                 if ap < 0.05]

    print(f"\nPairwise comparisons between {len(multi_phyla)} multi-species phyla "
          f"(Holm-corrected, {len(pairs)} pairs tested):")
    if sig_pairs:
        for p1, p2, rp, ap in sig_pairs:
            print(f"  {p1} vs {p2}:  p_raw={rp:.4f}  p_adj={ap:.4f}  ***")
    else:
        print("  No pairs significantly different after Holm correction.")

    # Save a flat table (not the all-1.0 matrix that mixed single+multi phyla)
    pair_rows = []
    for (p1, p2), rp, ap in zip(pair_labels, raw_ps, adj_ps):
        pair_rows.append({"phylum_1": p1, "phylum_2": p2,
                          "pvalue_raw": round(rp, 6),
                          "pvalue_holm": round(ap, 6),
                          "sig": ap < 0.05})
    pmat_path = RESULTS_DIR / "phylum_stats.tsv"
    pd.DataFrame(pair_rows).to_csv(pmat_path, sep="\t", index=False)
    print(f"\nPairwise p-value table (multi-species phyla only): {pmat_path}")

    # ── 3. Within-phylum CV ──────────────────────────────────────────────────
    cv_rows = []
    print("\nWithin-phylum consistency (CV of pct_terminal):")
    for phylum in phyla_present:
        sub = enr_df[enr_df["phylum"] == phylum]["pct_terminal"]
        if len(sub) < 2:
            continue
        mean_ = sub.mean()
        std_  = sub.std()
        cv    = std_ / mean_ * 100 if mean_ > 0 else None
        cv_rows.append({"phylum": phylum, "n_species": len(sub),
                        "mean_pct_terminal": round(mean_, 2),
                        "std_pct_terminal": round(std_, 2),
                        "cv_pct": round(cv, 1) if cv else None})
        print(f"  {phylum:<22}  n={len(sub)}  mean={mean_:.1f}%  "
              f"SD={std_:.1f}  CV={cv:.1f}%" if cv else
              f"  {phylum:<22}  n={len(sub)}  mean={mean_:.1f}%")

    cv_df = pd.DataFrame(cv_rows)
    cv_path = RESULTS_DIR / "within_phylum_cv.tsv"
    cv_df.to_csv(cv_path, sep="\t", index=False)
    print(f"\nWithin-phylum CV table: {cv_path}")

    # ── Supplementary Figure 4: within-phylum CV vs. sampling depth ─────────
    cv_plot = cv_df[cv_df["n_species"] >= 2].dropna(subset=["cv_pct"])
    fig4, ax4 = plt.subplots(figsize=(9, 7))
    colours4 = [PHYLUM_COLOURS.get(p, "#888") for p in cv_plot["phylum"]]
    ax4.scatter(cv_plot["n_species"], cv_plot["cv_pct"], c=colours4, s=55,
                edgecolors="black", linewidths=0.6, zorder=5)
    ax4.set_xscale("log")
    # Extra headroom beyond the data range so edge-point labels (e.g. the
    # leftmost single-digit-species phyla, or Insecta at the far right) have
    # somewhere to go instead of being pushed past the plot border.
    x_min, x_max = cv_plot["n_species"].min(), cv_plot["n_species"].max()
    ax4.set_xlim(x_min / 1.8, x_max * 2.5)
    y_min, y_max = cv_plot["cv_pct"].min(), cv_plot["cv_pct"].max()
    y_pad = (y_max - y_min) * 0.08
    ax4.set_ylim(y_min - y_pad, y_max + y_pad)
    _place_labels_no_overlap(fig4, ax4, cv_plot["n_species"].values,
                              cv_plot["cv_pct"].values, cv_plot["phylum"].values)
    ax4.set_xlabel("Number of species (log scale)", fontsize=11)
    ax4.set_ylabel("Within-phylum CV of % terminal (species-level)", fontsize=11)
    ax4.set_title("Within-phylum coefficient of variation vs. sampling depth\n"
                  "(phyla/groups with n ≥ 2 species)", fontsize=12)
    plt.tight_layout()
    out4 = FIGURES_DIR / "suppfig4_cv_vs_nspecies.pdf"
    fig4.savefig(out4, dpi=300, bbox_inches="tight")
    fig4.savefig(str(out4).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"Supplementary Figure 4 saved: {out4}")
    plt.close(fig4)

    # ── 3b. Holm-Bonferroni on phylum-level pooled Fisher p-values ──────────
    phylum_path = RESULTS_DIR / "phylum_summary.tsv"
    if phylum_path.exists():
        phylum_df = pd.read_csv(phylum_path, sep="\t")
        if "pooled_pvalue" in phylum_df.columns:
            phylum_ps = phylum_df["pooled_pvalue"].tolist()
            adj_phylum_ps = holm_bonferroni(phylum_ps)
            phylum_df["pooled_pvalue_holm"] = [round(p, 6) for p in adj_phylum_ps]
            phylum_df["pooled_sig_holm"] = phylum_df["pooled_pvalue_holm"] < 0.05
            phylum_df.to_csv(phylum_path, sep="\t", index=False)
            n_sig_raw  = (phylum_df["pooled_pvalue"] < 0.05).sum()
            n_sig_holm = phylum_df["pooled_sig_holm"].sum()
            print(f"\nPhylum-level pooled Fisher p-values (Holm-Bonferroni across "
                  f"{len(phylum_df)} phyla):")
            print(f"  Significant before correction: {n_sig_raw}/{len(phylum_df)}")
            print(f"  Significant after  correction: {n_sig_holm}/{len(phylum_df)}")
            provisional = phylum_df[
                (phylum_df["pooled_pvalue"] < 0.05) &
                (~phylum_df["pooled_sig_holm"])
            ]["phylum"].tolist()
            if provisional:
                print(f"  Reclassified as provisional (0.05 > p_holm): "
                      f"{', '.join(provisional)}")
            print(f"  Updated {phylum_path} with pooled_pvalue_holm column")

    # ── 4. Nematode anomaly test ──────────────────────────────────────────────
    nema = enr_df[enr_df["phylum"] == "Nematoda"]["pct_terminal"].values
    rest = enr_df[enr_df["phylum"] != "Nematoda"]["pct_terminal"].values
    if len(nema) >= 2:
        _, p_nema = mannwhitneyu(nema, rest, alternative="greater")
        print(f"\nNematoda vs all others (Mann-Whitney, one-sided greater):")
        print(f"  Nematoda mean = {nema.mean():.1f}%  vs  others mean = {rest.mean():.1f}%")
        print(f"  p = {p_nema:.4f}  {'→ Nematoda significantly higher ***' if p_nema < 0.05 else '→ ns'}")

    # ── Figure 8: box + strip plot per phylum ────────────────────────────────
    phyla_multi = [p for p in phyla_present
                   if len(enr_df[enr_df["phylum"] == p]) >= 2]

    fig, ax = plt.subplots(figsize=(12, 5))

    positions = list(range(len(phyla_present)))
    for i, phylum in enumerate(phyla_present):
        sub = enr_df[enr_df["phylum"] == phylum]["pct_terminal"].values
        colour = PHYLUM_COLOURS.get(phylum, "#aaa")

        if len(sub) >= 2:
            bp = ax.boxplot(
                sub, positions=[i], widths=0.5,
                patch_artist=True,
                boxprops=dict(facecolor=colour, alpha=0.6),
                medianprops=dict(color="black", linewidth=1.5),
                whiskerprops=dict(linewidth=0.8),
                capprops=dict(linewidth=0.8),
                flierprops=dict(marker="o", markersize=4, alpha=0.5),
                showfliers=False,
            )
        # Individual species dots
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(sub))
        ax.scatter(
            np.full(len(sub), i) + jitter, sub,
            color=colour, edgecolor="black", s=35, zorder=5, linewidth=0.6
        )

    ax.axhline(10.0, color="red", linestyle=":", linewidth=1.2,
               label="Null expectation (10%)")
    ax.axhspan(15, 25, color="grey", alpha=0.15,
               label="Teekas Tetrapoda range (15–25%)")

    ax.set_xticks(positions)
    ax.set_xticklabels(phyla_present, rotation=60, ha="right", fontsize=9)
    ax.set_ylabel("% LCRs in terminal bins (bins 1 & 20)", fontsize=11)
    ax.set_title("Terminal LCR enrichment distribution per phylum\n"
                 "(box = IQR, dots = individual species)", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_ylim(0, max(enr_df["pct_terminal"].max() * 1.15, 28))

    # Annotate KW result
    if p_kw is not None:
        ax.text(0.99, 0.97, f"Kruskal-Wallis p={p_kw:.3f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    plt.tight_layout()
    out = FIGURES_DIR / "suppfig1_phylum_distribution.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nSupplementary Figure 1 saved: {out}")
    plt.close()


if __name__ == "__main__":
    main()
