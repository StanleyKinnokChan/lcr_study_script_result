#!/usr/bin/env python3
"""
Holm-Bonferroni multiple-testing correction across all phylum-level Fisher's tests.

Reads phylum_summary.tsv (output of 03_analyse.py) and applies Holm-Bonferroni
correction to the pooled Fisher's exact p-values. Writes a corrected summary and
flags phyla that do not survive correction. This is the single owner of this
correction — do not duplicate it in another script (09_phylum_stats.py used to;
that copy was removed since nothing downstream read it).

Backs manuscript Table 1/Table 2 and the "41/43 phyla significant, only
Acanthocephala and Nematomorpha provisional" claim (manuscript_v10.md). The
provisional set will differ if the species set or filters change — re-run this
script and re-check the printed "Provisional" line against the manuscript
whenever upstream data changes; do not hardcode which phyla it names.

Outputs:
  results/table1_2_phylum_lineage_enrichment.tsv  — backs main-text Table 1
                                                      (Metazoa) and Table 2
                                                      (non-metazoan eukaryotes);
                                                      split by phylum membership
                                                      in config.SUPERGROUP_OF
                                                      ("Metazoa" vs other), no
                                                      metazoan flag column of its
                                                      own
  results/multiple_testing_report.txt   — human-readable summary
"""

import pandas as pd
import numpy as np
from pathlib import Path

from config import RESULTS_DIR, PHYLUM_ORDER


def holm_bonferroni(pvalues: list[float]) -> list[float]:
    """
    Step-down Holm-Bonferroni correction.
    Returns adjusted p-values in the same order as input.
    """
    n = len(pvalues)
    order = np.argsort(pvalues)           # ascending
    adjusted = np.empty(n)
    running_max = 0.0
    for rank, idx in enumerate(order):
        adj = pvalues[idx] * (n - rank)
        running_max = max(running_max, adj)
        adjusted[idx] = min(running_max, 1.0)
    return adjusted.tolist()


def sig_stars(p_corrected: float) -> str:
    if p_corrected < 0.001:
        return "***"
    if p_corrected < 0.01:
        return "**"
    if p_corrected < 0.05:
        return "*"
    return "ns"


def main():
    summary_path = RESULTS_DIR / "phylum_summary.tsv"
    if not summary_path.exists():
        print(f"ERROR: {summary_path} missing — run 03_analyse.py first.")
        return

    df = pd.read_csv(summary_path, sep="\t")
    df = df.dropna(subset=["pooled_pvalue"]).copy()

    n_tests = len(df)
    print(f"Applying Holm-Bonferroni correction across {n_tests} phylum-level tests.\n")

    raw_ps = df["pooled_pvalue"].tolist()
    adj_ps = holm_bonferroni(raw_ps)

    df["pooled_pvalue_holm"] = [round(p, 6) for p in adj_ps]
    df["pooled_sig_corrected"] = df["pooled_pvalue_holm"] < 0.05
    df["sig_stars"] = df["pooled_pvalue_holm"].apply(sig_stars)
    df["provisional"] = ~df["pooled_sig_corrected"]

    # Order by PHYLUM_ORDER where possible
    phylum_rank = {p: i for i, p in enumerate(PHYLUM_ORDER)}
    df["_rank"] = df["phylum"].map(lambda p: phylum_rank.get(p, 999))
    df = df.sort_values("_rank").drop(columns="_rank")

    out_path = RESULTS_DIR / "table1_2_phylum_lineage_enrichment.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Corrected summary: {out_path}\n")

    # ── Human-readable report ─────────────────────────────────────────────────
    report_lines = [
        "Holm-Bonferroni Multiple-Testing Correction Report",
        f"Tests performed: {n_tests} (one pooled Fisher's exact test per phylum/group)",
        f"Correction method: Holm-Bonferroni (step-down)",
        "",
        f"{'Phylum':<25} {'N_spp':>5}  {'%Term':>6}  {'p_raw':>10}  {'p_holm':>10}  {'Sig':>4}  {'Status'}",
        "─" * 85,
    ]
    n_sig = 0
    provisional_list = []
    for _, row in df.iterrows():
        status = "PROVISIONAL" if row["provisional"] else "significant"
        if not row["provisional"]:
            n_sig += 1
        else:
            provisional_list.append(row["phylum"])
        report_lines.append(
            f"  {row['phylum']:<23} {int(row['n_species']):>5}  "
            f"{row['pct_terminal']:>5.1f}%  "
            f"{row['pooled_pvalue']:>10.4e}  "
            f"{row['pooled_pvalue_holm']:>10.4e}  "
            f"{row['sig_stars']:>4}  {status}"
        )

    provisional_note = (
        f"  - {', '.join(provisional_list)} {'is' if len(provisional_list) == 1 else 'are'} not "
        f"significant after Holm-Bonferroni and must be described as 'provisional' "
        f"pending additional species."
        if provisional_list else
        "  - All phyla are significant after Holm-Bonferroni; none are provisional."
    )
    report_lines += [
        "─" * 85,
        f"",
        f"Significant after correction (p_holm < 0.05): {n_sig} / {n_tests}",
        f"Provisional (does not survive correction): {', '.join(provisional_list) if provisional_list else 'none'}",
        "",
        "Interpretation for manuscript:",
        "  - All text referring to significance should cite corrected p-values.",
        provisional_note,
    ]

    report_text = "\n".join(report_lines)
    print(report_text)

    report_path = RESULTS_DIR / "multiple_testing_report.txt"
    report_path.write_text(report_text + "\n")
    print(f"\nReport saved: {report_path}")

    # ── Console summary ──────────────────────────────────────────────────────
    print(f"\n→ {n_sig}/{n_tests} phyla significant after Holm-Bonferroni correction.")
    if provisional_list:
        print(f"→ Provisional (reclassified): {', '.join(provisional_list)}")


if __name__ == "__main__":
    main()
