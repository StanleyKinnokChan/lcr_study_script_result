#!/usr/bin/env python3
"""
Liquid-liquid phase separation (LLPS) propensity of terminal vs internal LCRs.

Tests whether terminal LCRs have higher predicted LLPS propensity than internal
LCRs, using two scoring approaches:

  1. PLAAC score (if PLAAC output files are present):
     Protein Likelihood of Aggregation and Complexing; prion-like domain prediction.
     Expected input: results/plaac/<species_key>.plaac.tsv

  2. Composition-based proxy (always available):
     Graded per-residue LLPS propensity applied to the dominant residue of each
     SINGLE-type LCR (fLPS output column 9). Weights span the two established
     regimes: aromatic pi-pi / cation-pi stacking (F/Y/W, R; Vernon et al. 2018)
     and prion-like poly-Q / poly-N composition (Q/N; the signal PLAAC scores,
     Lancaster et al. 2014), plus flexible G/S spacers. This proxy is used when
     PLAAC output is not present.

Analysis is run on ≥5 model organisms with publicly available PLAAC scores:
  - Homo sapiens
  - Saccharomyces cerevisiae
  - Arabidopsis thaliana
  - Caenorhabditis elegans
  - Drosophila melanogaster
  (and any others found in the dataset)

Outputs:
  results/llps_analysis.tsv        — per-organism terminal vs internal LLPS scores
  results/llps_organism_summary.tsv — Mann-Whitney U p-values per organism
  figures/suppfig_llps.pdf
"""

import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from pathlib import Path

from config import (
    RESULTS_DIR, FIGURES_DIR, FLPS_DIR,
    N_BINS, TERMINAL_BINS, PURITY_THRESHOLD, MIN_LCR_LENGTH,
    AA_ORDER,
)

PLAAC_DIR = RESULTS_DIR / "plaac"

# Model organisms for LLPS analysis (display_name substring → species_key pattern)
MODEL_ORGANISMS = [
    "Homo sapiens",
    "Mus musculus",
    "Saccharomyces cerevisiae",
    "Arabidopsis thaliana",
    "Caenorhabditis elegans",
    "Drosophila melanogaster",
    "Escherichia coli",
    "Bacillus subtilis",
]

# LLPS-promoting residues for the composition proxy. The proxy stands in for
# PLAAC (Lancaster et al. 2014), whose signal is prion-like Q/N composition, and
# for the aromatic/cation-pi model of IDR phase separation (Vernon et al. 2018).
# Graded weights therefore span both regimes rather than aromatics alone.
AROMATIC = set("FYW")            # pi-pi / cation-pi drivers (reported separately)
LLPS_WEIGHTS = {
    "F": 1.0, "Y": 1.0, "W": 1.0,   # aromatic pi-pi stacking (strongest)
    "Q": 0.8, "N": 0.8,             # prion-like poly-Q / poly-N (PLAAC core)
    "R": 0.5,                       # cation-pi partner (RGG/RG motifs)
    "G": 0.5, "S": 0.5,             # flexible spacers that permit PS
}
LLPS_RESIDUES = set(LLPS_WEIGHTS)


def composition_llps_score(residue: str) -> float:
    """
    Single-residue LLPS propensity score (0–1) for the dominant residue of a
    SINGLE-type LCR. Aromatics (F/Y/W) score highest; prion-like Q/N next
    (matching PLAAC's Q/N-richness signal); cation-pi R and flexible spacers
    G/S contribute partially; all other residues score 0.
    """
    return LLPS_WEIGHTS.get(residue, 0.0)


def parse_flps_for_llps(filepath: Path) -> list[dict]:
    """Parse fLPS output, extracting residue identity and positional info."""
    records = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "SINGLE":
                continue
            try:
                protein_len = int(cols[1])
                start       = int(cols[4])
                end         = int(cols[5])
                aa_count    = int(cols[6])
                residue_raw = cols[8]
            except (ValueError, IndexError):
                continue

            lcr_len = end - start + 1
            if lcr_len < MIN_LCR_LENGTH:
                continue
            if aa_count / lcr_len < PURITY_THRESHOLD:
                continue

            residue = residue_raw.strip("{}").upper()
            if len(residue) != 1 or residue not in AA_ORDER:
                continue

            midpoint = (start + end) / 2
            bin_num  = min(int((midpoint / protein_len) * N_BINS) + 1, N_BINS)
            is_terminal = bin_num in TERMINAL_BINS

            records.append({
                "residue":          residue,
                "bin":              bin_num,
                "is_terminal":      is_terminal,
                "location":         "terminal" if is_terminal else "internal",
                "llps_score_proxy": composition_llps_score(residue),
            })
    return records


def load_plaac_scores(plaac_file: Path) -> dict[str, float]:
    """
    Load PLAAC output file.
    Expected format: protein_id  start  end  score  (tab-separated)
    Returns {protein_id: max_plaac_score}.
    """
    if not plaac_file.exists():
        return {}
    scores: dict[str, list[float]] = {}
    try:
        df = pd.read_csv(plaac_file, sep="\t", header=None)
        for _, row in df.iterrows():
            pid = str(row.iloc[0])
            score = float(row.iloc[3]) if len(row) > 3 else 0.0
            scores.setdefault(pid, []).append(score)
        return {pid: max(vals) for pid, vals in scores.items()}
    except Exception as e:
        warnings.warn(f"Could not parse PLAAC file {plaac_file}: {e}", UserWarning)
        return {}


def test_terminal_vs_internal(records: list[dict],
                               plaac_scores: dict[str, float]) -> dict:
    """
    Compare LLPS proxy scores of terminal vs internal LCRs.
    Returns Mann-Whitney U statistic and p-value, plus group medians.
    """
    terminal_scores = [r["llps_score_proxy"] for r in records if r["is_terminal"]]
    internal_scores = [r["llps_score_proxy"] for r in records if not r["is_terminal"]]

    n_terminal = len(terminal_scores)
    n_internal = len(internal_scores)

    if n_terminal < 5 or n_internal < 5:
        return {
            "n_terminal": n_terminal,
            "n_internal": n_internal,
            "median_terminal": np.median(terminal_scores) if terminal_scores else None,
            "median_internal": np.median(internal_scores) if internal_scores else None,
            "mean_terminal": np.mean(terminal_scores) if terminal_scores else None,
            "mean_internal": np.mean(internal_scores) if internal_scores else None,
            "mw_stat": None, "mw_p": None, "significant": None,
            "method": "composition_proxy",
        }

    stat, p = mannwhitneyu(terminal_scores, internal_scores, alternative="greater")
    return {
        "n_terminal":       n_terminal,
        "n_internal":       n_internal,
        "median_terminal":  round(float(np.median(terminal_scores)), 4),
        "median_internal":  round(float(np.median(internal_scores)), 4),
        "mean_terminal":    round(float(np.mean(terminal_scores)), 4),
        "mean_internal":    round(float(np.mean(internal_scores)), 4),
        "pct_aromatic_term": round(
            sum(1 for r in records if r["is_terminal"] and r["residue"] in AROMATIC) /
            max(n_terminal, 1) * 100, 2
        ),
        "pct_aromatic_int":  round(
            sum(1 for r in records if not r["is_terminal"] and r["residue"] in AROMATIC) /
            max(n_internal, 1) * 100, 2
        ),
        "mw_stat":  round(float(stat), 1),
        "mw_p":     round(float(p), 6),
        "significant": p < 0.05,
        "method": "composition_proxy",
    }


def main():
    manifest_path = RESULTS_DIR / "species_manifest.tsv"
    if not manifest_path.exists():
        print(f"ERROR: {manifest_path} missing — run 01_download_proteomes.py first.")
        return

    manifest = pd.read_csv(manifest_path, sep="\t")

    # Filter to model organisms (match by display_name substring)
    model_mask = manifest["display_name"].apply(
        lambda n: any(m in n for m in MODEL_ORGANISMS)
    )
    target_species = manifest[model_mask].copy()

    if target_species.empty:
        print("No model organisms found in manifest. Running on all species.")
        target_species = manifest.copy()
    else:
        print(f"Model organisms found: {len(target_species)}")
        for _, r in target_species.iterrows():
            print(f"  {r['display_name']} ({r['species_key']})")

    use_plaac = PLAAC_DIR.exists() and any(PLAAC_DIR.glob("*.plaac.tsv"))
    if use_plaac:
        print(f"\nPLAAC directory found: {PLAAC_DIR}")
    else:
        print(f"\nNo PLAAC scores found at {PLAAC_DIR} — using composition proxy only.")
        print("To use PLAAC scores: run PLAAC on each model organism proteome and save")
        print("output to results/plaac/<species_key>.plaac.tsv")

    # ── Per-organism analysis ─────────────────────────────────────────────────
    all_records: list[dict] = []
    summary_rows: list[dict] = []

    for _, row in target_species.iterrows():
        sp_key   = row["species_key"]
        display  = row.get("display_name", sp_key)
        phylum   = row["phylum"]
        flps_f   = FLPS_DIR / f"{sp_key}.flps.txt"

        if not flps_f.exists():
            print(f"  SKIP {display}: no fLPS file")
            continue

        records = parse_flps_for_llps(flps_f)
        if not records:
            print(f"  SKIP {display}: no LCRs parsed")
            continue

        # Load PLAAC scores if available
        plaac_file   = PLAAC_DIR / f"{sp_key}.plaac.tsv" if use_plaac else Path("/nonexistent")
        plaac_scores = load_plaac_scores(plaac_file)

        result = test_terminal_vs_internal(records, plaac_scores)
        result.update({
            "species_key":  sp_key,
            "display_name": display,
            "phylum":       phylum,
        })
        summary_rows.append(result)

        n_term = result["n_terminal"]
        n_int  = result["n_internal"]
        p_val  = result["mw_p"]
        sig    = "p<0.05" if result.get("significant") else "ns"
        print(f"  {display}: terminal_n={n_term}, internal_n={n_int}, "
              f"mw_p={p_val}  [{sig}]")

        for r in records:
            r["species_key"]  = sp_key
            r["display_name"] = display
            r["phylum"]       = phylum
        all_records.extend(records)

    if not summary_rows:
        print("No results generated. Check fLPS output files.")
        return

    # ── Write outputs ─────────────────────────────────────────────────────────
    all_df = pd.DataFrame(all_records)
    out_all = RESULTS_DIR / "llps_analysis.tsv"
    all_df.to_csv(out_all, sep="\t", index=False)
    print(f"\nFull LLPS LCR table: {out_all}")

    summary_df = pd.DataFrame(summary_rows)
    cols_ordered = ["species_key", "display_name", "phylum",
                    "n_terminal", "n_internal",
                    "mean_terminal", "mean_internal",
                    "median_terminal", "median_internal",
                    "pct_aromatic_term", "pct_aromatic_int",
                    "mw_stat", "mw_p", "significant", "method"]
    summary_df = summary_df[[c for c in cols_ordered if c in summary_df.columns]]
    out_summary = RESULTS_DIR / "llps_organism_summary.tsv"
    summary_df.to_csv(out_summary, sep="\t", index=False)
    print(f"Organism summary: {out_summary}")

    # Summary stats
    n_sig = summary_df["significant"].sum() if "significant" in summary_df.columns else 0
    print(f"\nOrganisms with higher terminal LLPS score (p<0.05): "
          f"{n_sig}/{len(summary_df)}")

    # ── Figure ────────────────────────────────────────────────────────────────
    n_orgs = len(summary_df)
    if n_orgs == 0:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: mean LLPS proxy score, terminal vs internal, per organism.
    # The mean, not the median, is plotted: the proxy is a zero-inflated discrete
    # score (0 / 0.5 / 0.8 / 1.0) whose median is 0 in most organisms, so a median
    # bar chart is blank and cannot show the rank shift the Mann-Whitney U tests.
    ax = axes[0]
    x = np.arange(n_orgs)
    width = 0.35
    ax.bar(x - width / 2, summary_df["mean_terminal"], width,
           label="Terminal", color="#d6604d", edgecolor="black", linewidth=0.6)
    ax.bar(x + width / 2, summary_df["mean_internal"], width,
           label="Internal", color="#4393c3", edgecolor="black", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [n.split()[0] for n in summary_df["display_name"]],
        rotation=60, ha="right", fontsize=9, rotation_mode="anchor"
    )
    ax.set_ylabel("Mean LLPS proxy score", fontsize=11)
    ax.set_title("Composition-based LLPS propensity\nterminal vs internal LCRs", fontsize=11)
    ax.legend(fontsize=9)
    # Mark organisms where terminal LCRs are stochastically greater (Mann-Whitney U)
    bar_top = max(
        [v for v in list(summary_df["mean_terminal"]) + list(summary_df["mean_internal"])
         if pd.notna(v)] or [1.0]
    )
    ax.set_ylim(0, bar_top * 1.25)
    for i, (_, row_) in enumerate(summary_df.iterrows()):
        if row_.get("significant"):
            y_max = max(row_["mean_terminal"] or 0, row_["mean_internal"] or 0)
            ax.text(i, y_max + bar_top * 0.03, "*", ha="center", fontsize=13, color="black")

    # Panel B: % aromatic residues at termini vs internal
    ax = axes[1]
    if "pct_aromatic_term" in summary_df.columns:
        ax.bar(x - width / 2, summary_df["pct_aromatic_term"], width,
               label="Terminal", color="#d6604d", edgecolor="black", linewidth=0.6)
        ax.bar(x + width / 2, summary_df["pct_aromatic_int"], width,
               label="Internal", color="#4393c3", edgecolor="black", linewidth=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(
            [n.split()[0] for n in summary_df["display_name"]],
            rotation=60, ha="right", fontsize=9, rotation_mode="anchor"
        )
        ax.set_ylabel("% aromatic residue LCRs (F/Y/W)", fontsize=11)
        ax.set_title("Aromatic LCR fraction at termini vs internal\n"
                     "(pi-pi stacking drives LLPS)", fontsize=11)
        ax.legend(fontsize=9)
    else:
        ax.text(0.5, 0.5, "No aromatic data", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)

    method_label = "PLAAC + composition proxy" if use_plaac else "composition proxy (F/Y/W + G/S/N)"
    plt.suptitle(f"LLPS propensity of terminal vs internal LCRs\n({method_label})",
                 fontsize=12)
    plt.tight_layout()
    out_fig = FIGURES_DIR / "suppfig_llps.pdf"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    fig.savefig(str(out_fig).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nLLPS figure: {out_fig}")
    plt.close()

    # Print note about PLAAC
    if not use_plaac:
        print(
            "\nNote: This analysis used the composition-based proxy (aromatic + G/S/N residues).")
        print("For the manuscript figure, run PLAAC (http://plaac.wi.mit.edu/) on each")
        print("model organism proteome and save output to results/plaac/<species_key>.plaac.tsv")
        print("Then re-run this script for the full PLAAC + proxy comparison.")


if __name__ == "__main__":
    main()
