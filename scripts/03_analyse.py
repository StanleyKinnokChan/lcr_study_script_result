#!/usr/bin/env python3
"""
Parse fLPS 2.0 output, compute terminal LCR enrichment per species/phylum,
and run statistical tests replicating Teekas et al. (2024) Open Biology.

Key metric: fraction of LCRs falling in terminal bins (bins 1 or 20 of 20),
compared to random expectation (2/20 = 10%).

Outputs:
  results/lcr_positions.tsv     — per-LCR position table
  results/enrichment.tsv        — per-species enrichment + Fisher's exact p
  results/phylum_summary.tsv    — aggregated by phylum with pooled Fisher's test
"""

import re
import csv
from pathlib import Path
from collections import defaultdict
from scipy.stats import fisher_exact, wilcoxon
import pandas as pd

from config import (
    PROJECT_DIR, RESULTS_DIR, FLPS_DIR, MANIFEST,
    N_BINS, TERMINAL_BINS, NULL_TERMINAL, PURITY_THRESHOLD, MIN_LCR_LENGTH,
)

# ── fLPS output parser ────────────────────────────────────────────────────────

def parse_flps(filepath: Path) -> list[dict]:
    """
    Parse fLPS 2.0 output into a list of LCR records.

    fLPS2 actual output format (tab-separated, no header lines):
        protein_id  protein_len  TYPE  rank  start  end  aa_count  pvalue  {residue}  [flag  lambda]

    TYPE is SINGLE (single-residue bias), MULTIPLE (multi-residue), or WHOLE (entire sequence).
    We keep only SINGLE rows — these are classic single-amino-acid LCRs.
    Purity = aa_count / lcr_length (computed here; Teekas threshold ≥ 70%).
    """
    records = []

    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            if cols[2] != "SINGLE":
                continue

            try:
                protein_id  = cols[0]
                protein_len = int(cols[1])
                start       = int(cols[4])
                end         = int(cols[5])
                aa_count    = int(cols[6])
            except (ValueError, IndexError):
                continue

            lcr_len = end - start + 1
            if lcr_len < MIN_LCR_LENGTH:
                continue

            purity = aa_count / lcr_len
            if purity < PURITY_THRESHOLD:
                continue

            midpoint = (start + end) / 2
            norm_mid = midpoint / protein_len
            bin_num  = min(int(norm_mid * N_BINS) + 1, N_BINS)  # 1–20

            records.append({
                "protein_id":    protein_id,
                "protein_len":   protein_len,
                "lcr_start":     start,
                "lcr_end":       end,
                "lcr_len":       lcr_len,
                "purity":        round(purity, 4),
                "norm_midpoint": round(norm_mid, 4),
                "bin":           bin_num,
                "is_terminal":   bin_num in TERMINAL_BINS,
            })

    return records


# ── Enrichment stats ──────────────────────────────────────────────────────────

def enrichment_stats(records: list[dict]) -> dict:
    """
    Fisher's exact test: terminal LCRs vs internal LCRs vs null expectation.
    Null: terminal bins = 2/20 = 10% of all bins.
    """
    n_terminal = sum(1 for r in records if r["is_terminal"])
    n_internal = len(records) - n_terminal
    total = len(records)
    if total == 0:
        return {"n_lcr": 0, "n_terminal": 0, "pct_terminal": None,
                "odds_ratio": None, "pvalue": None}

    # Expected terminal count under null (2/20 of all LCRs)
    expected_terminal = total * (2 / N_BINS)
    expected_internal = total - expected_terminal

    # 2x2 contingency: [[observed_terminal, observed_internal],
    #                    [expected_terminal, expected_internal]]
    table = [
        [n_terminal,          n_internal],
        [expected_terminal,   expected_internal],
    ]
    odds_ratio, pvalue = fisher_exact(table, alternative="greater")

    return {
        "n_lcr":        total,
        "n_terminal":   n_terminal,
        "pct_terminal": round(n_terminal / total * 100, 2),
        "odds_ratio":   round(odds_ratio, 3),
        "pvalue":       round(pvalue, 6),
        "significant":  pvalue < 0.05,
    }


def pooled_fisher(n_terminal: int, n_total: int) -> tuple[float, float]:
    """
    Fisher's exact test on pooled LCRs from an entire group (e.g. all Bacteria).
    Same test as species-level but applied to aggregate counts.
    Returns (odds_ratio, pvalue).
    """
    if n_total == 0:
        return (None, None)
    exp_term = n_total * (2 / N_BINS)
    exp_int  = n_total - exp_term
    n_int    = n_total - n_terminal
    table = [[n_terminal, n_int], [exp_term, exp_int]]
    or_, p = fisher_exact(table, alternative="greater")
    return (round(or_, 3), round(p, 6))


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if not MANIFEST.exists():
        print(f"ERROR: manifest not found at {MANIFEST}")
        print("Run 01_download_proteomes.py first.")
        return

    manifest = pd.read_csv(MANIFEST, sep="\t")
    print(f"Loaded manifest: {len(manifest)} species\n")

    # Authoritative domain per phylum, taken from the manifest's `domain` column
    # (written by 01a/01b from NCBI Taxonomy). This covers protist phyla the legacy
    # DOMAIN_MAP below never enumerated. Falls back to that map for older manifests.
    if "domain" in manifest.columns:
        phylum_domain = (
            manifest.dropna(subset=["domain"])
            .groupby("phylum")["domain"]
            .agg(lambda s: s.mode().iat[0] if not s.mode().empty else s.iloc[0])
            .to_dict()
        )
    else:
        phylum_domain = {}

    all_records  = []
    enrichment_rows = []

    for _, row in manifest.iterrows():
        sp_key   = row["species_key"]
        phylum   = row["phylum"]
        display  = row["display_name"]
        flps_out = FLPS_DIR / f"{sp_key}.flps.txt"

        if not flps_out.exists():
            print(f"[SKIP] No fLPS output for {display} — run 02_run_flps.sh")
            continue

        print(f"[{phylum}] {display}")
        records = parse_flps(flps_out)
        print(f"  LCRs parsed (purity ≥{PURITY_THRESHOLD}): {len(records)}")

        for r in records:
            r["species_key"] = sp_key
            r["display_name"] = display
            r["phylum"] = phylum
        all_records.extend(records)

        stats = enrichment_stats(records)
        enrichment_rows.append({
            "species_key":  sp_key,
            "display_name": display,
            "phylum":       phylum,
            **stats,
        })
        sig = "***" if stats.get("significant") else "ns"
        low_conf = stats["n_lcr"] < 50
        conf_tag = "  [LOW CONFIDENCE: n_lcr<50, test underpowered]" if low_conf else ""
        print(f"  Terminal LCRs: {stats['n_terminal']}/{stats['n_lcr']} "
              f"({stats['pct_terminal']}%)  OR={stats['odds_ratio']}  "
              f"p={stats['pvalue']}  {sig}{conf_tag}")
        print()

    if not all_records:
        print("No records to analyse. Check fLPS output files.")
        return

    # ── Per-LCR position table
    pos_df = pd.DataFrame(all_records)
    pos_path = RESULTS_DIR / "lcr_positions.tsv"
    pos_df.to_csv(pos_path, sep="\t", index=False)
    print(f"Position table: {pos_path}  ({len(pos_df)} LCRs)")

    # ── Per-species enrichment table
    enr_df = pd.DataFrame(enrichment_rows)
    enr_path = RESULTS_DIR / "enrichment.tsv"
    enr_df.to_csv(enr_path, sep="\t", index=False)
    print(f"Enrichment table: {enr_path}")

    # ── Phylum-level summary with pooled Fisher's exact test ─────────────────
    # The pooled test applies the same 70% purity filter and same Fisher's test
    # as the species-level analysis, but pools all LCRs within a phylum/domain.
    # This is the statistically appropriate way to test prokaryote domains where
    # individual species are underpowered (n_lcr < 50).
    phylum_agg = (
        enr_df.groupby("phylum")
        .agg(
            n_species=("species_key", "count"),
            total_lcr=("n_lcr", "sum"),
            total_terminal=("n_terminal", "sum"),
        )
        .assign(pct_terminal=lambda d: (d.total_terminal / d.total_lcr * 100).round(2))
        .reset_index()
    )

    pooled_rows = []
    for _, prow in phylum_agg.iterrows():
        or_, p = pooled_fisher(int(prow.total_terminal), int(prow.total_lcr))
        pooled_rows.append({
            **prow.to_dict(),
            "pooled_OR":      or_,
            "pooled_pvalue":  p,
            "pooled_sig":     (p < 0.05) if p is not None else None,
            "single_species": prow.n_species == 1,
        })

    phylum_summary = pd.DataFrame(pooled_rows)
    phylum_path = RESULTS_DIR / "phylum_summary.tsv"
    phylum_summary.to_csv(phylum_path, sep="\t", index=False)

    print(f"\nPhylum summary (with pooled Fisher's exact test):")
    print(phylum_summary.to_string(index=False))
    print(f"\nPhylum table: {phylum_path}")

    # ── Key comparison vs Teekas (2024) Tetrapoda baseline ──────────────────────
    # Fallback map for manifests predating the `domain` column; phylum_domain
    # (from the manifest) takes precedence and covers any phylum NCBI labelled.
    DOMAIN_MAP = {
        "Bacteria": "Bacteria", "Archaea": "Archaea",
        "Fungi": "Non-metazoan Eukaryota", "Viridiplantae": "Non-metazoan Eukaryota",
        "Amoebozoa": "Non-metazoan Eukaryota", "Apicomplexa": "Non-metazoan Eukaryota",
        "Euglenozoa": "Non-metazoan Eukaryota", "Metamonada": "Non-metazoan Eukaryota",
        "Chlorophyta": "Non-metazoan Eukaryota", "Heterolobosea": "Non-metazoan Eukaryota",
    }
    print("\n── Comparison to Teekas (2024) Tetrapoda baseline ──")
    print("Tetrapoda clades reported ~15–25% of LCRs in terminal bins.")
    print(f"{'Phylum':<25} {'n_sp':>4}  {'%Term':>6}  {'pooled_p':>10}  {'vs Tetrapoda range':}")
    for _, r in phylum_summary.iterrows():
        domain = phylum_domain.get(r.phylum) or DOMAIN_MAP.get(r.phylum, "Metazoa")
        if r.pct_terminal is None:
            continue
        flag = "↑ above tetrapod range" if r.pct_terminal > 25 else \
               "↓ below tetrapod range" if r.pct_terminal < 15 else \
               "≈ within tetrapod range"
        p_str = f"{r.pooled_pvalue:.2e}" if r.pooled_pvalue is not None else "n/a"
        sig_str = " ***" if r.pooled_sig else ""
        print(f"  [{domain:<26}] {r.phylum:<22} n={int(r.n_species):<3}  "
              f"{r.pct_terminal:>5.1f}%  p={p_str}{sig_str}  {flag}")


if __name__ == "__main__":
    main()
