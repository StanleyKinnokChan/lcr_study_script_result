#!/usr/bin/env python3
"""
Signal peptide stratification of bacterial terminal LCRs.

Tests whether N-terminal LCR enrichment in bacteria is explained by signal
peptides at the N-terminus of secreted proteins, or is a general property of
all bacterial proteins.

Data source: UniProt flat-file annotations embedded in proteome FASTA headers
or a pre-downloaded signal peptide annotation TSV.  The script supports two
input modes:
  MODE A (preferred): a TSV file with columns [protein_id, has_signal_peptide]
          produced by running `scripts/fetch_uniprot_annotations.py` (see below)
          or downloaded manually from UniProt.
  MODE B (fallback):  naive N-terminal methionine + small residue heuristic
          (predicts likely signal-peptide proteins by +2 residue MAP rule).

Signal peptide TSV format (MODE A):
  protein_id  has_signal_peptide   subcellular_location
  Q9EXJ4      True                 Secreted
  P00274      False                Cytoplasm

If no annotation file is found, the script runs in MODE B and prints a warning.

Outputs:
  results/signal_peptide_stratification.tsv  — Supp Table S8
  figures/suppfig_signal_peptide.pdf

To generate the annotation TSV from UniProt:
  For each bacterial species in UniProt reference proteomes, download the
  .txt flat file and grep for "FT   SIGNAL" lines.  A helper script is
  provided at: scripts/fetch_uniprot_annotations.py
"""

import sys
import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu
from pathlib import Path

from config import (
    RESULTS_DIR, FIGURES_DIR, FLPS_DIR,
    N_BINS, TERMINAL_BINS, PURITY_THRESHOLD, MIN_LCR_LENGTH,
)

ANNOTATION_FILE = RESULTS_DIR / "uniprot_signal_peptides.tsv"
SIGNAL_PEPTIDE_MARKER = "SIGNAL"      # UniProt flat file feature key
PROTEOMES_DIR = RESULTS_DIR.parent / "data" / "proteomes"


def accession(protein_id: str) -> str:
    """
    Reduce a fLPS protein identifier to a UniProt accession for annotation
    lookup. UniProt FASTA headers are `sp|ACC|NAME` / `tr|ACC|NAME`; anything
    else (Ensembl/NCBI-sourced ids) is returned unchanged and simply won't
    match a UniProt annotation.
    """
    parts = protein_id.split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr"):
        return parts[1]
    return protein_id


def parse_flps_bacterial(filepath: Path) -> list[dict]:
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
            if aa_count / lcr_len < PURITY_THRESHOLD:
                continue

            midpoint = (start + end) / 2
            bin_num  = min(int((midpoint / protein_len) * N_BINS) + 1, N_BINS)
            is_nterm = bin_num == 1
            is_cterm = bin_num == N_BINS

            records.append({
                "protein_id":  protein_id,
                "protein_len": protein_len,
                "bin":         bin_num,
                "is_terminal": bin_num in TERMINAL_BINS,
                "is_nterm":    is_nterm,
                "is_cterm":    is_cterm,
            })
    return records


def fisher_nterm(n_nterm: int, n_total: int) -> tuple[float, float]:
    if n_total == 0:
        return (None, None)
    exp_term = n_total * (1 / N_BINS)
    exp_other = n_total - exp_term
    n_other = n_total - n_nterm
    table = [[n_nterm, n_other], [exp_term, exp_other]]
    or_, p = fisher_exact(table, alternative="greater")
    return round(or_, 3), round(p, 6)


def stratify_and_test(records: list[dict],
                      annot: dict[str, bool]) -> dict:
    """
    Split records into signal-peptide vs non-signal-peptide proteins.
    Test N-terminal LCR enrichment in each class.
    """
    with_sp = [r for r in records if annot.get(accession(r["protein_id"]), False)]
    without_sp = [r for r in records
                  if accession(r["protein_id"]) in annot
                  and not annot[accession(r["protein_id"])]]
    unannotated = [r for r in records if accession(r["protein_id"]) not in annot]

    def stats(subset: list[dict], label: str) -> dict:
        total = len(subset)
        n_nterm = sum(r["is_nterm"] for r in subset)
        or_, p = fisher_nterm(n_nterm, total)
        return {
            "class":         label,
            "n_lcr":         total,
            "n_nterm":       n_nterm,
            "pct_nterm":     round(n_nterm / total * 100, 2) if total > 0 else None,
            "nterm_OR":      or_,
            "nterm_p":       p,
            "nterm_sig":     (p < 0.05) if p is not None else None,
        }

    return {
        "with_signal_peptide":    stats(with_sp,    "with_signal_peptide"),
        "without_signal_peptide": stats(without_sp, "without_signal_peptide"),
        "unannotated":            stats(unannotated, "unannotated"),
    }


def mode_b_heuristic(fasta_path: Path) -> dict[str, bool]:
    """
    Naive heuristic for MODE B: predict signal peptides from FASTA sequence.
    A protein is flagged as likely-secreted if residue 2 is large/charged
    (preventing MAP methionine excision) AND the N-terminal region is hydrophobic.
    This is a rough approximation; MODE A (UniProt annotations) is preferred.
    """
    warnings.warn(
        "MODE B: Using sequence heuristic for signal peptide prediction. "
        "Results are approximate. Download UniProt annotations for MODE A.",
        UserWarning, stacklevel=2
    )
    predictions: dict[str, bool] = {}
    large_residues = set("RKDEFHYWLIMV")   # block MAP if at position 2
    current_id = None
    seq_chars: list[str] = []

    def _classify(pid: str, seq: list[str]) -> None:
        if not pid or len(seq) < 15:
            return
        # Key by accession so MODE B matches the same lookup key as MODE A.
        key = accession(pid)
        if len(seq) >= 2 and seq[1] in large_residues:
            # Check hydrophobic stretch in first 15 aa
            hydro = sum(1 for aa in seq[:15] if aa in "ILMFVWA")
            predictions[key] = hydro >= 6

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    _classify(current_id, seq_chars)
                current_id = line[1:].split()[0]
                seq_chars = []
            else:
                seq_chars.extend(line.upper())
    if current_id:
        _classify(current_id, seq_chars)

    return predictions


def main():
    manifest_path = RESULTS_DIR / "species_manifest.tsv"
    if not manifest_path.exists():
        print(f"ERROR: {manifest_path} missing — run 01_download_proteomes.py first.")
        return

    manifest = pd.read_csv(manifest_path, sep="\t")
    bacterial = manifest[manifest["phylum"] == "Bacteria"]

    if bacterial.empty:
        print("No bacterial species found in manifest.")
        return

    print(f"Bacterial species: {len(bacterial)}")

    # ── Load or generate signal peptide annotations ───────────────────────────
    use_mode_b = False
    if ANNOTATION_FILE.exists():
        print(f"MODE A: Loading signal peptide annotations from {ANNOTATION_FILE}")
        sp_annot_df = pd.read_csv(ANNOTATION_FILE, sep="\t")
        # Parse booleans robustly: read_csv may deliver the column as strings
        # ("True"/"False"), and str→bool via astype(bool) is truthy for both.
        truthy = {"true", "1", "yes", "t"}
        has_sp = (sp_annot_df["has_signal_peptide"].astype(str)
                  .str.strip().str.lower().isin(truthy))
        global_annot: dict[str, bool] = dict(
            zip(sp_annot_df["protein_id"].astype(str), has_sp)
        )
        print(f"  Annotations loaded: {len(global_annot):,} proteins "
              f"({sum(global_annot.values()):,} with signal peptide)")
    else:
        print(
            f"WARNING: {ANNOTATION_FILE} not found.\n"
            "  Falling back to MODE B (sequence heuristic).\n"
            "  To use UniProt annotations (MODE A), run:\n"
            "    python3 scripts/fetch_uniprot_annotations.py\n"
            "  and re-run this script."
        )
        use_mode_b = True
        global_annot = {}

    # ── Per-species stratification ────────────────────────────────────────────
    result_rows = []

    for _, row in bacterial.iterrows():
        sp_key = row["species_key"]
        display = row.get("display_name", sp_key)
        flps_f = FLPS_DIR / f"{sp_key}.flps.txt"

        if not flps_f.exists():
            continue

        records = parse_flps_bacterial(flps_f)
        if not records:
            continue

        if use_mode_b:
            # Find FASTA file (proteomes live under data/proteomes/)
            fasta_candidates = list(PROTEOMES_DIR.glob(f"{sp_key}*.fa*"))
            if fasta_candidates:
                sp_annot = mode_b_heuristic(fasta_candidates[0])
            else:
                sp_annot = {}
        else:
            sp_annot = global_annot

        strat = stratify_and_test(records, sp_annot)
        for class_label, stats in strat.items():
            result_rows.append({
                "species_key": sp_key,
                "display_name": display,
                **stats,
            })
        print(f"  {display}: total={len(records)} LCRs  "
              f"with_SP={strat['with_signal_peptide']['n_lcr']}  "
              f"without_SP={strat['without_signal_peptide']['n_lcr']}")

    if not result_rows:
        print("No results generated. Check fLPS bacterial output files.")
        return

    results_df = pd.DataFrame(result_rows)
    out_tsv = RESULTS_DIR / "signal_peptide_stratification.tsv"
    results_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nSignal peptide stratification table (Supp Table S8): {out_tsv}")

    # ── Pooled comparison ─────────────────────────────────────────────────────
    print("\nPooled N-terminal enrichment by signal-peptide class:")
    for class_label in ["with_signal_peptide", "without_signal_peptide"]:
        sub = results_df[results_df["class"] == class_label]
        total = sub["n_lcr"].sum()
        n_nterm = sub["n_nterm"].sum()
        pct = round(n_nterm / total * 100, 2) if total > 0 else 0
        or_, p = fisher_nterm(int(n_nterm), int(total))
        sig = "***" if (p is not None and p < 0.001) else \
              "**" if (p is not None and p < 0.01) else \
              "*" if (p is not None and p < 0.05) else "ns"
        print(f"  {class_label:<30}: {n_nterm}/{total} "
              f"({pct}% N-terminal)  OR={or_}  p={p}  {sig}")

    # ── Mann-Whitney: pct_nterm distributions by class ────────────────────────
    with_sp_pcts = results_df[
        (results_df["class"] == "with_signal_peptide") &
        (results_df["pct_nterm"].notna())
    ]["pct_nterm"].values
    without_sp_pcts = results_df[
        (results_df["class"] == "without_signal_peptide") &
        (results_df["pct_nterm"].notna())
    ]["pct_nterm"].values

    if len(with_sp_pcts) >= 3 and len(without_sp_pcts) >= 3:
        stat, p_mw = mannwhitneyu(with_sp_pcts, without_sp_pcts, alternative="two-sided")
        print(f"\nMann-Whitney U (with_SP vs without_SP % N-terminal): "
              f"U={stat:.1f}, p={p_mw:.4f}")
        if p_mw < 0.05:
            direction = ("higher" if np.median(with_sp_pcts) > np.median(without_sp_pcts)
                         else "lower")
            print(f"  → Signal-peptide proteins have significantly {direction} "
                  f"N-terminal LCR fraction")
        else:
            print(f"  → No significant difference between classes (supports "
                  f"mechanism beyond signal peptides alone)")

    # ── Figure ────────────────────────────────────────────────────────────────
    classes = ["with_signal_peptide", "without_signal_peptide"]
    class_labels = ["With signal\npeptide", "Without signal\npeptide"]
    colours = ["#3182bd", "#de2d26"]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Panel A: pooled pct_nterm per class
    ax = axes[0]
    pooled_pcts = []
    for cl in classes:
        sub = results_df[results_df["class"] == cl]
        total = sub["n_lcr"].sum()
        n_nterm = sub["n_nterm"].sum()
        pooled_pcts.append(n_nterm / total * 100 if total > 0 else 0)

    ax.bar(class_labels, pooled_pcts, color=colours,
           edgecolor="black", linewidth=0.7)
    ax.axhline(5.0, color="red", linestyle=":", linewidth=1.2,
               label="Null (5% = 1/20 bins)")
    ax.set_ylabel("% LCRs in N-terminal bin (bin 1)", fontsize=11)
    ax.set_title("Pooled N-terminal LCR enrichment\nby protein class", fontsize=11)
    ax.legend(fontsize=9)

    # Panel B: per-species violin
    ax = axes[1]
    data_by_class = [
        results_df[(results_df["class"] == cl) & results_df["pct_nterm"].notna()]["pct_nterm"].values
        for cl in classes
    ]
    parts = ax.violinplot(data_by_class, positions=[1, 2],
                          showmedians=True, showextrema=True)
    for pc, colour in zip(parts["bodies"], colours):
        pc.set_facecolor(colour)
        pc.set_alpha(0.7)

    for i, (data, cl) in enumerate(zip(data_by_class, class_labels)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(data))
        ax.scatter(np.full(len(data), i + 1) + jitter, data,
                   color="black", s=20, alpha=0.6, zorder=5)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(class_labels)
    ax.axhline(5.0, color="red", linestyle=":", linewidth=1.2)
    ax.set_ylabel("Species-level % N-terminal LCRs")
    ax.set_title("Species distribution of\nN-terminal LCR enrichment", fontsize=11)
    if len(with_sp_pcts) >= 3 and len(without_sp_pcts) >= 3:
        p_label = f"p={p_mw:.3f}"
        ax.text(0.5, 0.95, p_label, transform=ax.transAxes,
                ha="center", va="top", fontsize=10,
                bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    mode_str = "MODE A (UniProt annotations)" if not use_mode_b else "MODE B (sequence heuristic)"
    plt.suptitle(f"Signal peptide stratification of bacterial N-terminal LCRs\n({mode_str})",
                 fontsize=12)
    plt.tight_layout()
    out_fig = FIGURES_DIR / "suppfig_signal_peptide.pdf"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    fig.savefig(str(out_fig).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nFigure: {out_fig}")
    plt.close()


if __name__ == "__main__":
    main()
