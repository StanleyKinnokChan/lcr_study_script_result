#!/usr/bin/env python3
"""
Main pipeline orchestrator for the invertebrate terminal LCR study.

Runs all scripts in logical order, skipping steps whose outputs already exist.
Each phase is independent — you can run from any phase after earlier outputs exist.

Usage:
  python3 main.py                      # run all phases
  python3 main.py --phase 1            # phase 1 only (download)
  python3 main.py --from-phase 3       # phases 3–5
  python3 main.py --skip-download      # skip 01a/01b/01c (use existing proteomes)
  python3 main.py --dry-run            # print what would run, don't execute
  python3 main.py --force              # re-run all steps even if outputs exist
  python3 main.py --list               # show status of every step then exit

Phases:
  1 — Data acquisition  (01a, 01b, 01c, 02_run_flps)
  2 — Core analysis     (03_analyse, 04_visualise)
  3 — Extended analysis (05_asymmetry … 11_protein_level_test)
  4 — GO enrichment     (12_go_enrichment)
  5 — Robustness & mechanism (14_multiple_testing … 20_llps_analysis, with
      fetch_uniprot_annotations before 17 and build_viridiplantae_backbone_tree
      before 19 as data-prep prerequisites)

Note: 02_run_flps.sh and the fLPS binary are Unix-only — run this on Mac/Linux.
"""

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

PROJECT_DIR = Path(__file__).parent
RESULTS_DIR = PROJECT_DIR / "results"
FIGURES_DIR = PROJECT_DIR / "figures"

# ── Pipeline definition ───────────────────────────────────────────────────────
# Each step: (command, label, sentinel_output, is_shell_cmd)
# sentinel_output: relative path whose existence means "already done". None = always run.

PHASES = {
    1: {
        "name": "Data acquisition",
        "steps": [
            (
                "python3 scripts/01a_download_metazoa.py",
                "Download metazoan proteomes (Ensembl Metazoa r63)",
                "results/supp_table_S1_species_list.tsv",
                False,
            ),
            (
                "python3 scripts/01b_download_outgroups.py",
                "Download outgroup proteomes (UniProt: bacteria/archaea/fungi/plants)",
                None,   # appends to manifest — always check internally
                False,
            ),
            (
                "python3 scripts/01c_download_viridiplantae_expansion.py",
                "Download Viridiplantae tier expansion (ferns, gymnosperms, charophytes)",
                None,   # always checks internally per-species
                False,
            ),
            (
                "bash scripts/02_run_flps.sh",
                "Run fLPS 2.0 on all proteomes (skips existing per-species)",
                None,   # always check — skips per-species internally
                True,
            ),
        ],
    },
    2: {
        "name": "Core analysis",
        "steps": [
            (
                "python3 scripts/03_analyse.py",
                "Parse fLPS output → LCR positions, enrichment stats, phylum summary",
                "results/enrichment.tsv",
                False,
            ),
            (
                "python3 scripts/04_visualise.py",
                "Core figures: bin heatmap, terminal bar chart, U-profile",
                "figures/fig1_bin_heatmap.pdf",
                False,
            ),
        ],
    },
    3: {
        "name": "Extended analyses",
        "steps": [
            (
                "python3 scripts/05_asymmetry.py",
                "N-terminal vs C-terminal enrichment asymmetry",
                "results/asymmetry.tsv",
                False,
            ),
            (
                "python3 scripts/06_aa_composition.py",
                "Amino acid identity of terminal vs internal LCRs",
                "results/aa_composition.tsv",
                False,
            ),
            (
                "python3 scripts/07_confound_test.py",
                "Protein-length confound test (quartile stratification)",
                "results/supp_table_S4_length_stratified.tsv",
                False,
            ),
            (
                "python3 scripts/08_purity_gradient.py",
                "LCR purity gradient: terminal vs internal",
                "results/supp_table_S5_purity_gradient.tsv",
                False,
            ),
            (
                "python3 scripts/09_phylum_stats.py",
                "Between-phylum Kruskal-Wallis + within-phylum variance",
                "results/phylum_stats.tsv",
                False,
            ),
            (
                "python3 scripts/10_driver_analysis.py",
                "Singleton vs multi-LCR protein driver analysis",
                "results/supp_table_S2_driver_analysis.tsv",
                False,
            ),
            (
                "python3 scripts/11_protein_level_test.py",
                "Protein-level binomial test (sensitivity analysis, main-text only)",
                "results/protein_level_enrichment.tsv",
                False,
            ),
        ],
    },
    4: {
        "name": "GO enrichment",
        "steps": [
            (
                "python3 scripts/12_go_enrichment.py",
                "GO term enrichment in terminal-LCR proteins (Ensembl BioMart)",
                "results/go_enrichment.tsv",
                False,
            ),
        ],
    },
    5: {
        "name": "Extended robustness & mechanism analyses",
        "steps": [
            (
                "python3 scripts/14_multiple_testing.py",
                "Holm-Bonferroni correction across 43 phylum-level Fisher tests",
                "results/table1_2_phylum_lineage_enrichment.tsv",
                False,
            ),
            (
                "python3 scripts/15_sensitivity_analysis.py",
                "fLPS parameter sensitivity: 4 param sets (Supp Table S6)",
                "results/supp_table_S6_flps_sensitivity.tsv",
                False,
            ),
            (
                "python3 scripts/16_domain_aa_composition.py",
                "Per-domain AA composition of terminal vs internal LCRs (Supp Table S7)",
                "results/supp_table_S7_domain_aa_composition.tsv",
                False,
            ),
            (
                # Prerequisite for step 17: fetches UniProtKB SIGNAL annotations.
                "python3 scripts/fetch_uniprot_annotations.py",
                "Fetch UniProt signal-peptide annotations (prereq for step 17)",
                "results/uniprot_signal_peptides.tsv",
                False,
            ),
            (
                "python3 scripts/17_signal_peptide_stratification.py",
                "Signal peptide stratification of bacterial terminal LCRs (Supp Table S8)",
                "results/supp_table_S8_signal_peptide_stratification.tsv",
                False,
            ),
            (
                "python3 scripts/18_mixed_effects_model.py",
                "Bootstrap CI + cluster-robust GEE for prokaryote terminal enrichment",
                "results/bootstrap_ci_domain.tsv",
                False,
            ),
            (
                # Prerequisite for step 19: builds the conservative dated backbone
                # tree; without it, 19 degrades to OLS-on-tier-rank.
                "python3 scripts/build_viridiplantae_backbone_tree.py",
                "Build Viridiplantae backbone phylogeny (prereq for step 19 PGLS)",
                "results/viridiplantae_timetree.nwk",
                False,
            ),
            (
                "python3 scripts/19_pgls_viridiplantae.py",
                "PGLS regression: N/C asymmetry vs evolutionary tier in land plants",
                "results/pgls_viridiplantae.tsv",
                False,
            ),
            (
                "python3 scripts/20_llps_analysis.py",
                "LLPS propensity of terminal vs internal LCRs (model organisms)",
                "results/supp_table_S9_llps_organism_summary.tsv",
                False,
            ),
        ],
    },
}

# ── Terminal colours ──────────────────────────────────────────────────────────

RESET  = "\033[0m"
BOLD   = "\033[1m"
GREEN  = "\033[32m"
YELLOW = "\033[33m"
RED    = "\033[31m"
CYAN   = "\033[36m"
GREY   = "\033[90m"


def banner(text: str, colour: str = CYAN):
    w = 72
    print(f"\n{colour}{BOLD}{'─' * w}{RESET}")
    print(f"{colour}{BOLD}  {text}{RESET}")
    print(f"{colour}{BOLD}{'─' * w}{RESET}")


def log(symbol: str, msg: str, colour: str = ""):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"{GREY}[{ts}]{RESET} {colour}{symbol}{RESET}  {msg}", flush=True)


def sentinel_done(sentinel: str | None, force: bool) -> bool:
    if force or sentinel is None:
        return False
    return (PROJECT_DIR / sentinel).exists()


def run_step(cmd: str, label: str, sentinel: str | None,
             is_shell: bool, dry_run: bool, force: bool) -> bool:
    """Execute one pipeline step. Returns True on success or skip."""
    if sentinel_done(sentinel, force):
        log("✓", f"{label}", GREY)
        log(" ", f"  → skipping ({sentinel} already exists)", GREY)
        return True

    log("▶", label, YELLOW)

    if dry_run:
        log(" ", f"  {GREY}[dry-run] {cmd}{RESET}")
        return True

    t0 = time.time()
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}

    if is_shell:
        proc = subprocess.run(cmd, shell=True, cwd=PROJECT_DIR, env=env)
    else:
        proc = subprocess.run(cmd.split(), cwd=PROJECT_DIR, env=env)

    elapsed = time.time() - t0

    if proc.returncode == 0:
        log("✓", f"{label}  {GREY}({elapsed:.0f}s){RESET}", GREEN)
        return True
    else:
        log("✗", f"{label}  [exit code {proc.returncode}]", RED)
        return False


def status_table():
    """Print a status table of all steps and whether their outputs exist."""
    banner("Pipeline status", CYAN)
    for ph_num, phase in PHASES.items():
        print(f"\n{BOLD}Phase {ph_num}: {phase['name']}{RESET}")
        for cmd, label, sentinel, _ in phase["steps"]:
            done = sentinel_done(sentinel, force=False)
            if done:
                sym = f"{GREEN}✓{RESET}"
                note = f"{GREY}done ({sentinel}){RESET}"
            elif sentinel is None:
                sym = f"{YELLOW}~{RESET}"
                note = f"{GREY}always re-checks internally{RESET}"
            else:
                sym = f"{YELLOW}○{RESET}"
                note = f"{GREY}pending{RESET}"
            print(f"  {sym}  {label:<65}  {note}")


# ── Argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Orchestrate the invertebrate terminal LCR study pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    g = p.add_mutually_exclusive_group()
    g.add_argument("--phase", type=int, choices=[1, 2, 3, 4, 5],
                   help="Run only this phase.")
    g.add_argument("--from-phase", type=int, choices=[1, 2, 3, 4, 5], dest="from_phase",
                   help="Run from this phase to phase 5.")
    p.add_argument("--skip-download", action="store_true",
                   help="Skip 01a/01b/01c downloads (use existing proteomes).")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without running them.")
    p.add_argument("--force", action="store_true",
                   help="Re-run all steps even if outputs exist.")
    p.add_argument("--list", action="store_true",
                   help="Show step status and exit.")
    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if args.list:
        status_table()
        return

    if args.phase:
        phases_to_run = [args.phase]
    elif args.from_phase:
        phases_to_run = list(range(args.from_phase, 6))
    else:
        phases_to_run = [1, 2, 3, 4, 5]

    banner("Invertebrate Terminal LCR Study — Pipeline Orchestrator", CYAN)
    print(f"  Phases to run : {phases_to_run}  (1–5)")
    print(f"  Skip download : {args.skip_download}")
    print(f"  Dry run       : {args.dry_run}")
    print(f"  Force rerun   : {args.force}")
    print(f"  Working dir   : {PROJECT_DIR}")

    step_results: list[tuple[str, bool]] = []
    abort = False

    for ph_num in phases_to_run:
        phase = PHASES[ph_num]
        banner(f"Phase {ph_num} — {phase['name']}", YELLOW)

        for cmd, label, sentinel, is_shell in phase["steps"]:
            # Skip download scripts when requested
            if args.skip_download and ("01a" in cmd or "01b" in cmd or "01c" in cmd):
                log("–", f"Skipping (--skip-download): {label}", GREY)
                step_results.append((label, True))
                continue

            ok = run_step(cmd, label, sentinel, is_shell, args.dry_run, args.force)
            step_results.append((label, ok))

            if not ok:
                log("!", f"Step failed — stopping pipeline.", RED)
                print(f"\n  Fix the error above, then resume:")
                print(f"    python3 main.py --from-phase {ph_num}")
                print(f"  Or force-rerun from this phase:")
                print(f"    python3 main.py --from-phase {ph_num} --force\n")
                abort = True
                break

        if abort:
            break

    # ── Summary ──────────────────────────────────────────────────────────────
    banner("Summary", CYAN)
    n_ok   = sum(1 for _, ok in step_results if ok)
    n_fail = sum(1 for _, ok in step_results if not ok)

    for label, ok in step_results:
        sym = f"{GREEN}✓{RESET}" if ok else f"{RED}✗{RESET}"
        print(f"  {sym}  {label}")

    print()
    if n_fail == 0:
        print(f"{GREEN}{BOLD}  All {n_ok} steps completed successfully.{RESET}\n")
    else:
        print(f"{RED}{BOLD}  {n_fail} step(s) failed. See above for details.{RESET}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
