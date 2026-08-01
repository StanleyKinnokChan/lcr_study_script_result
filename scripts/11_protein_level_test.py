#!/usr/bin/env python3
"""
Protein-level sensitivity analysis.

Backs the "Protein-level confirmation" main-text paragraph; not currently a
numbered supplementary table in the manuscript.

Addresses the concern that LCR-level Fisher's exact tests treat multiple LCRs
within the same protein as independent observations.

For each species: compute the fraction of proteins bearing ≥1 terminal LCR
vs. the fraction bearing ≥1 internal-only LCR. Test against the 10% null
using a one-sided binomial test.

Also computes: fraction of proteins with any terminal LCR vs all proteins
with LCRs — a more intuitive protein-centric enrichment measure.

Outputs:
  results/protein_level_enrichment.tsv
"""

import pandas as pd
from scipy.stats import binomtest
from pathlib import Path

from config import PROJECT_DIR, RESULTS_DIR, N_BINS, NULL_TERMINAL


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
        # Per-protein: does it have at least one terminal LCR?
        protein_has_terminal = grp.groupby("protein_id")["is_terminal"].any()
        n_proteins_total    = len(protein_has_terminal)
        n_proteins_terminal = protein_has_terminal.sum()
        pct_proteins_terminal = round(n_proteins_terminal / n_proteins_total * 100, 2) if n_proteins_total > 0 else None

        # Binomial test: observed proportion vs null 10%
        result = binomtest(int(n_proteins_terminal), n_proteins_total, NULL_TERMINAL,
                           alternative="greater")
        p = round(result.pvalue, 6)

        rows.append({
            "species_key":                sp_key,
            "display_name":               display,
            "phylum":                     phylum,
            "n_proteins_with_lcr":        n_proteins_total,
            "n_proteins_with_terminal_lcr": int(n_proteins_terminal),
            "pct_proteins_terminal":      pct_proteins_terminal,
            "binomial_p":                 p,
            "significant":                p < 0.05,
        })
        print(f"[{phylum}] {display}: {n_proteins_terminal}/{n_proteins_total} proteins "
              f"have terminal LCR ({pct_proteins_terminal}%)  p={p}  "
              f"{'***' if p < 0.05 else 'ns'}")

    out_df = pd.DataFrame(rows)
    out_path = RESULTS_DIR / "protein_level_enrichment.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nProtein-level enrichment table: {out_path}")

    n_sig = out_df["significant"].sum()
    print(f"Significant: {n_sig}/{len(out_df)} species")
    print(f"Overall: {out_df['n_proteins_with_terminal_lcr'].sum()} / "
          f"{out_df['n_proteins_with_lcr'].sum()} proteins have ≥1 terminal LCR")


if __name__ == "__main__":
    main()
