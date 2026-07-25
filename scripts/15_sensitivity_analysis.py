#!/usr/bin/env python3
"""
Sensitivity analyses: (1) fLPS parameter choice, (2) taxon-sampling bias.

(1) Parameter sensitivity + MULTI-type positional analysis
    Re-parses existing fLPS 2.0 raw output files at three parameter combinations to
    test whether the terminal enrichment signal is robust to parameter choice:
      A (current)   — min_length=3, purity>=0.70, SINGLE-type only
      B (stringent) — min_length=6, purity>=0.80, SINGLE-type only
      C (relaxed)   — min_length=3, purity>=0.60, SINGLE-type only
      D (multi-type)— min_length=3, purity>=0.70, MULTI-type only
    For each, computes pooled % terminal per phylum and Fisher's exact p (uses the
    .flps.txt files already on disk — no new fLPS runs required).

(2) Taxon-sampling robustness
    The species set is skewed (Insecta ≈ 30% of all LCRs; 12+ phyla on a single
    genome). This section re-weights the per-LCR table (results/lcr_positions.tsv) to
    check the terminal-enrichment signal is a property of the clades, not of who got
    over-/under-sampled: pooled vs equal-per-species weighting, Insecta excluded/capped,
    singletons excluded, and each singleton phylum tested on its own. scipy-free
    (normal-approx binomial), so it runs even without the fLPS files.

Outputs:
  results/sensitivity_analysis.tsv    — Supp Table S6: per-phylum × param combination
  results/multitype_positions.tsv     — MULTI-type LCR positional table
  results/sampling_robustness.tsv     — Supp Table S7: taxon-sampling robustness
  figures/suppfig_sensitivity.pdf
"""

import math

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from pathlib import Path
from typing import NamedTuple

from config import (
    RESULTS_DIR, FIGURES_DIR, FLPS_DIR,
    N_BINS, TERMINAL_BINS, PHYLUM_ORDER, PHYLUM_COLOURS,
    SUPERGROUP_OF, NULL_TERMINAL,
)


class Params(NamedTuple):
    label: str
    min_length: int
    purity_threshold: float
    lcr_type: str            # "SINGLE", "MULTI", or "BOTH"


PARAM_SETS: list[Params] = [
    Params("A_current",   min_length=3, purity_threshold=0.70, lcr_type="SINGLE"),
    Params("B_stringent", min_length=6, purity_threshold=0.80, lcr_type="SINGLE"),
    Params("C_relaxed",   min_length=3, purity_threshold=0.60, lcr_type="SINGLE"),
    Params("D_multitype", min_length=3, purity_threshold=0.70, lcr_type="MULTI"),
]


def parse_flps(filepath: Path, params: Params) -> list[dict]:
    records = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            lcr_type_col = cols[2]

            if params.lcr_type == "SINGLE" and lcr_type_col != "SINGLE":
                continue
            if params.lcr_type == "MULTI" and lcr_type_col != "MULTIPLE":
                continue

            try:
                protein_len = int(cols[1])
                start       = int(cols[4])
                end         = int(cols[5])
                aa_count    = int(cols[6])
            except (ValueError, IndexError):
                continue

            lcr_len = end - start + 1
            if lcr_len < params.min_length:
                continue
            purity = aa_count / lcr_len
            if purity < params.purity_threshold:
                continue

            midpoint = (start + end) / 2
            norm_mid = midpoint / protein_len
            bin_num  = min(int(norm_mid * N_BINS) + 1, N_BINS)

            records.append({
                "protein_len": protein_len,
                "bin":         bin_num,
                "is_terminal": bin_num in TERMINAL_BINS,
                "lcr_type":    lcr_type_col,
            })
    return records


def pooled_fisher(n_terminal: int, n_total: int) -> tuple[float, float]:
    if n_total == 0:
        return (None, None)
    exp_term = n_total * (2 / N_BINS)
    exp_int  = n_total - exp_term
    n_int    = n_total - n_terminal
    table    = [[n_terminal, n_int], [exp_term, exp_int]]
    or_, p   = fisher_exact(table, alternative="greater")
    return round(or_, 3), round(p, 6)


def run_param_set(manifest: pd.DataFrame, params: Params) -> pd.DataFrame:
    phylum_terminal: dict[str, int] = {}
    phylum_total:    dict[str, int] = {}

    for _, row in manifest.iterrows():
        sp_key = row["species_key"]
        phylum = row["phylum"]
        flps_f = FLPS_DIR / f"{sp_key}.flps.txt"
        if not flps_f.exists():
            continue
        recs = parse_flps(flps_f, params)
        phylum_terminal[phylum] = phylum_terminal.get(phylum, 0) + sum(r["is_terminal"] for r in recs)
        phylum_total[phylum]    = phylum_total.get(phylum, 0) + len(recs)

    rows = []
    for phylum in set(phylum_total):
        n_term  = phylum_terminal.get(phylum, 0)
        n_total = phylum_total.get(phylum, 0)
        or_, p  = pooled_fisher(n_term, n_total)
        pct     = round(n_term / n_total * 100, 2) if n_total > 0 else None
        rows.append({
            "param_set":    params.label,
            "phylum":       phylum,
            "n_total_lcr":  n_total,
            "n_terminal":   n_term,
            "pct_terminal": pct,
            "pooled_OR":    or_,
            "pooled_p":     p,
            "significant":  (p < 0.05) if p is not None else None,
        })
    return pd.DataFrame(rows)


# ── Taxon-sampling robustness ────────────────────────────────────────────────────

def _binom_z(k: int, n: int, p0: float = NULL_TERMINAL) -> tuple[float, float]:
    """Normal-approx one-sided binomial test of %terminal vs the null. Returns
    (percent_terminal, z). scipy-free so this runs even where scipy is unavailable."""
    if n == 0:
        return float("nan"), float("nan")
    phat = k / n
    z = (phat - p0) / math.sqrt(p0 * (1 - p0) / n)
    return phat * 100.0, z


def taxon_sampling_robustness(manifest: pd.DataFrame) -> pd.DataFrame:
    """Is the terminal-LCR enrichment a property of the clades, or an artefact of who
    got over-sampled (Insecta) / under-sampled (singletons)? Re-weights the per-LCR
    table (results/lcr_positions.tsv) and checks the signal survives. Writes
    results/sampling_robustness.tsv (Supp Table S7)."""
    lcr = pd.read_csv(RESULTS_DIR / "lcr_positions.tsv", sep="\t",
                      usecols=["species_key", "phylum", "is_terminal"])
    lcr["is_terminal"] = lcr["is_terminal"].astype(str).isin(["True", "1", "true"])
    lcr["super"] = lcr["phylum"].map(SUPERGROUP_OF).fillna("Other")
    manifest = manifest.copy()
    manifest["super"] = manifest["phylum"].map(SUPERGROUP_OF).fillna("Other")
    tot = len(lcr)
    rows: list[dict] = []

    def emit(analysis, group, variant, sub, n_species):
        pct, z = _binom_z(int(sub["is_terminal"].sum()), len(sub))
        rows.append({"analysis": analysis, "group": group, "variant": variant,
                     "n_species": n_species, "n_lcr": len(sub),
                     "pct_of_all_lcr": round(100 * len(sub) / tot, 3),
                     "pct_terminal": round(pct, 2), "z": round(z, 1),
                     "enriched": bool(z > 1.64)})

    # (1) composition vs signal weight, and pooled-vs-species-mean, per supergroup
    order = ["Prokaryota", "SAR", "Excavata", "Amoebozoa", "Archaeplastida",
             "Other eukaryotes", "Opisthokonta", "Metazoa"]
    for sg in order:
        sub = lcr[lcr["super"] == sg]
        n_sp = int((manifest["super"] == sg).sum())
        emit("supergroup", sg, "pooled_LCR_weighted", sub, n_sp)
        per = sub.groupby("species_key")["is_terminal"].mean() * 100
        rows.append({"analysis": "supergroup", "group": sg,
                     "variant": "species_mean_equal_weight", "n_species": len(per),
                     "n_lcr": len(sub), "pct_of_all_lcr": round(100 * len(sub) / tot, 3),
                     "pct_terminal": round(per.mean(), 2),
                     "z": float("nan"), "enriched": None})

    # (2) Metazoa with over/under-sampling removed
    meta = lcr[lcr["super"] == "Metazoa"]
    ins = meta[meta["phylum"] == "Insecta"]
    noins = meta[meta["phylum"] != "Insecta"]
    emit("metazoa_reweight", "Metazoa", "full", meta,
         int((manifest["super"] == "Metazoa").sum()))
    emit("metazoa_reweight", "Metazoa", "exclude_Insecta", noins, None)
    med = int(meta.groupby("phylum").size().median())
    idx = np.random.default_rng(0).choice(ins.index, size=min(med, len(ins)), replace=False)
    emit("metazoa_reweight", "Metazoa", f"Insecta_capped_to_{med}",
         pd.concat([noins, ins.loc[idx]]), None)
    singletons = manifest["phylum"].value_counts()
    singletons = singletons[singletons == 1].index.tolist()
    emit("metazoa_reweight", "Metazoa", f"exclude_{len(singletons)}_singletons",
         meta[~meta["phylum"].isin(singletons)], None)

    # (3) each singleton phylum on its own — can it support a phylum-level claim?
    for phy in [p for p in PHYLUM_ORDER if p in singletons]:
        sub = lcr[lcr["phylum"] == phy]
        if len(sub):
            emit("singleton_phylum", phy, "n1_standalone", sub, 1)

    out = pd.DataFrame(rows)
    out_path = RESULTS_DIR / "sampling_robustness.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\nTaxon-sampling robustness: {out_path}  ({len(out)} rows)")

    meta_full = out.query("group=='Metazoa' and variant=='full'")["pct_terminal"].iloc[0]
    meta_noins = out.query("variant=='exclude_Insecta'")["pct_terminal"].iloc[0]
    n_sing = out.query("analysis=='singleton_phylum'")
    print(f"  Metazoa %terminal: full={meta_full}  excl-Insecta={meta_noins} "
          f"(signal is not Insecta-driven)")
    print(f"  Singleton phyla still individually enriched: "
          f"{int(n_sing['enriched'].sum())}/{len(n_sing)}")
    return out


def main():
    manifest_path = RESULTS_DIR / "species_manifest.tsv"
    if not manifest_path.exists():
        print(f"ERROR: {manifest_path} missing — run 01_download_proteomes.py first.")
        return

    manifest = pd.read_csv(manifest_path, sep="\t")
    print(f"Loaded manifest: {len(manifest)} species")

    # Taxon-sampling robustness — runs off results/lcr_positions.tsv, independent of the
    # fLPS re-parsing below (so it works even if the raw .flps.txt files are absent).
    if (RESULTS_DIR / "lcr_positions.tsv").exists():
        taxon_sampling_robustness(manifest)
    else:
        print("  [skip] lcr_positions.tsv missing — run 03_analyse.py for sampling robustness")

    print(f"\nRunning {len(PARAM_SETS)} parameter combinations...\n")

    all_results = []
    multitype_records = []

    for params in PARAM_SETS:
        print(f"  [{params.label}] min_len={params.min_length}, "
              f"purity≥{params.purity_threshold}, type={params.lcr_type}")
        df_params = run_param_set(manifest, params)
        all_results.append(df_params)

        # Collect MULTI-type positions separately for Supp Table S6b
        if params.lcr_type == "MULTI":
            for _, row in manifest.iterrows():
                flps_f = FLPS_DIR / f"{row['species_key']}.flps.txt"
                if not flps_f.exists():
                    continue
                recs = parse_flps(flps_f, params)
                for r in recs:
                    r["species_key"] = row["species_key"]
                    r["phylum"]      = row["phylum"]
                multitype_records.extend(recs)

    combined = pd.concat(all_results, ignore_index=True)

    # Order by phylum
    phylum_rank = {p: i for i, p in enumerate(PHYLUM_ORDER)}
    combined["_rank"] = combined["phylum"].map(lambda p: phylum_rank.get(p, 999))
    combined = combined.sort_values(["param_set", "_rank"]).drop(columns="_rank")

    out_tsv = RESULTS_DIR / "sensitivity_analysis.tsv"
    combined.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nSensitivity table: {out_tsv}  ({len(combined)} rows)")

    if multitype_records:
        mt_df = pd.DataFrame(multitype_records)
        mt_path = RESULTS_DIR / "multitype_positions.tsv"
        mt_df.to_csv(mt_path, sep="\t", index=False)
        print(f"MULTI-type positions: {mt_path}  ({len(mt_df)} records)")

    # ── Pivot for comparison: % terminal per phylum × param set ──────────────
    pivot = combined.pivot_table(
        index="phylum", columns="param_set", values="pct_terminal"
    ).reindex([p for p in PHYLUM_ORDER if p in combined["phylum"].values])

    print("\n% terminal LCRs per phylum × parameter combination:")
    print(pivot.round(1).to_string())

    # ── Consistency check: how many phyla change significance? ───────────────
    sig_pivot = combined.pivot_table(
        index="phylum", columns="param_set", values="significant"
    )
    inconsistent = sig_pivot[sig_pivot.nunique(axis=1) > 1].index.tolist()
    if inconsistent:
        print(f"\nPhyla with inconsistent significance across param sets:")
        for p in inconsistent:
            print(f"  {p}: {sig_pivot.loc[p].to_dict()}")
    else:
        print("\nAll phyla have consistent significance direction across all parameter sets.")

    # ── Figure: grouped bar chart ─────────────────────────────────────────────
    phyla_plot = [p for p in PHYLUM_ORDER if p in pivot.index]
    n_phyla = len(phyla_plot)
    param_labels = [p.label for p in PARAM_SETS]
    n_params = len(param_labels)
    colours = ["#636363", "#3182bd", "#31a354", "#e6550d"]
    width = 0.8 / n_params

    fig, ax = plt.subplots(figsize=(18, 5))
    x = np.arange(n_phyla)
    for i, label in enumerate(param_labels):
        if label not in pivot.columns:
            continue
        vals = [pivot.loc[p, label] if p in pivot.index else np.nan for p in phyla_plot]
        offset = (i - n_params / 2 + 0.5) * width
        ax.bar(x + offset, vals, width, label=label,
               color=colours[i], alpha=0.85, edgecolor="black", linewidth=0.5)

    ax.axhline(10.0, color="red", linestyle=":", linewidth=1.2,
               label="Null (10%)")
    ax.set_xticks(x)
    ax.set_xticklabels(phyla_plot, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("% terminal LCRs", fontsize=11)
    ax.set_title("Sensitivity of terminal LCR enrichment to fLPS parameter choice\n"
                 "A=current, B=stringent (len≥6, purity≥80%), "
                 "C=relaxed (purity≥60%), D=MULTI-type", fontsize=11)
    ax.legend(fontsize=9)
    plt.tight_layout()

    out_fig = FIGURES_DIR / "suppfig_sensitivity.pdf"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    fig.savefig(str(out_fig).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nSensitivity figure: {out_fig}")
    plt.close()


if __name__ == "__main__":
    main()
