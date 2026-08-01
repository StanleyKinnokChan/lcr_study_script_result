#!/usr/bin/env python3
"""
Shared configuration for all analysis scripts.

Single source of truth for:
  - Project paths
  - Analysis constants (bins, thresholds, null expectations)
  - Canonical phylum order and colour scheme (43 groups, derived from
    supp_table_S1_species_list.tsv / phylum_summary.tsv for the 772-species dataset)
  - Supergroup membership
  - Amino acid ordering and colours
"""

from pathlib import Path

# ── Project paths ─────────────────────────────────────────────────────────────
PROJECT_DIR = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIGURES_DIR = PROJECT_DIR / "figures"
FLPS_DIR    = RESULTS_DIR / "flps"
MANIFEST    = RESULTS_DIR / "supp_table_S1_species_list.tsv"
FIGURES_DIR.mkdir(exist_ok=True)

# ── Analysis constants ────────────────────────────────────────────────────────
N_BINS            = 20
TERMINAL_BINS     = {1, N_BINS}
NULL_TERMINAL     = 2 / N_BINS   # 10%  — combined two-terminus null
NULL_SINGLE       = 1 / N_BINS   # 5%   — single-terminus null
MIN_LCR_LENGTH    = 3
# ≥70% purity filter applied to ALL species including Bacteria and Archaea.
# Prokaryotes yield fewer LCRs under this filter — that is a biological result,
# not a methodological failure.  Do NOT lower it for prokaryotes; doing so
# would make their results incomparable to Teekas et al. and to the eukaryote
# data in this study.
PURITY_THRESHOLD  = 0.70

# Teekas et al. (2024) Tetrapoda reference range
TETRAPODA_TERMINAL_LOW  = 15.0
TETRAPODA_TERMINAL_HIGH = 25.0

# Length-quartile labels used by confound-test and related plots
QUARTILE_LABELS = ["Q1\n(shortest)", "Q2", "Q3", "Q4\n(longest)"]

# ── Canonical phylum order (43 groups, phylogenetically arranged) ─────────────
# Contains every phylum/group present in phylum_summary.tsv for the 772-species
# dataset.  Scripts filter to `phyla_present` at runtime so extra entries here
# are harmless — but a phylum present in the data and MISSING from this list is
# silently dropped from any script that uses PHYLUM_ORDER as a filter (not just
# an order), e.g. `[p for p in PHYLUM_ORDER if p in df["phylum"].values]`. This
# happened to Heterolobosea (fixed 2026-08-01) and dropped it from Figures 1-3
# and Supp Tables S2/S4/S5 despite passing Holm-Bonferroni in the main analysis.
PHYLUM_ORDER = [
    # ── Prokaryotes ───────────────────────────────────────────────────────────
    "Bacteria",
    "Archaea",
    # ── SAR: Stramenopiles ────────────────────────────────────────────────────
    "Oomycota",
    "Bacillariophyta",
    # ── SAR: Alveolata ────────────────────────────────────────────────────────
    "Apicomplexa",
    "Ciliophora",
    "Perkinsozoa",
    # ── SAR: Rhizaria ─────────────────────────────────────────────────────────
    "Rhizaria",
    # ── Excavata ──────────────────────────────────────────────────────────────
    "Euglenozoa",
    "Metamonada",
    "Heterolobosea",
    # ── Amoebozoa ─────────────────────────────────────────────────────────────
    "Amoebozoa",
    # ── Archaeplastida ────────────────────────────────────────────────────────
    "Viridiplantae",
    "Chlorophyta",
    "Rhodophyta",
    # ── Other/unplaced eukaryotes ─────────────────────────────────────────────
    "Haptophyta",
    "Cryptophyta",
    "Protist",
    # ── Opisthokonta (non-metazoan) ───────────────────────────────────────────
    "Fungi",
    # ── Metazoa: basal lineages ───────────────────────────────────────────────
    "Porifera",
    "Ctenophora",
    "Placozoa",
    "Xenacoelomorpha",
    "Cnidaria",
    # ── Metazoa: Lophotrochozoa ───────────────────────────────────────────────
    "Platyhelminthes",
    "Rotifera",
    "Nemertea",
    "Annelida",
    "Brachiopoda",
    "Mollusca",
    "Acanthocephala",
    # ── Metazoa: Ecdysozoa ────────────────────────────────────────────────────
    "Nematomorpha",
    "Nematoda",
    "Priapulida",
    "Tardigrada",
    "Collembola",
    "Crustacea",
    "Myriapoda",
    "Chelicerata",
    "Insecta",
    # ── Metazoa: Deuterostomia ────────────────────────────────────────────────
    "Echinodermata",
    "Hemichordata",
    "Chordata",
]

# ── Colour scheme (one colour per phylum/group) ───────────────────────────────
PHYLUM_COLOURS = {
    # Prokaryotes — greys
    "Bacteria":        "#636363",
    "Archaea":         "#969696",
    # SAR Stramenopiles — teals
    "Oomycota":        "#80cdc1",
    "Bacillariophyta": "#35978f",
    # SAR Alveolata — greens
    "Apicomplexa":     "#a6d96a",
    "Ciliophora":      "#66c2a5",
    "Perkinsozoa":     "#abdda4",
    # SAR Rhizaria — blue
    "Rhizaria":        "#3288bd",
    # Excavata — dark greens
    "Euglenozoa":      "#66bd63",
    "Metamonada":      "#1a9850",
    "Heterolobosea":   "#78c679",
    # Amoebozoa — yellow-green
    "Amoebozoa":       "#d9ef8b",
    # Archaeplastida — greens and red
    "Viridiplantae":   "#4dac26",
    "Chlorophyta":     "#238443",
    "Rhodophyta":      "#d73027",
    # Other eukaryotes — oranges / neutral
    "Haptophyta":      "#fdae61",
    "Cryptophyta":     "#f46d43",
    "Protist":         "#bababa",
    # Opisthokonta (non-metazoan) — yellow
    "Fungi":           "#fee08b",
    # Basal Metazoa — light blues
    "Porifera":        "#c6dbef",
    "Ctenophora":      "#9ecae1",
    "Placozoa":        "#a6cee3",
    "Xenacoelomorpha": "#6baed6",
    "Cnidaria":        "#7fc97f",
    # Lophotrochozoa — purples and pinks
    "Platyhelminthes": "#cab2d6",
    "Rotifera":        "#e7d4e8",
    "Nemertea":        "#d4b9da",
    "Annelida":        "#fb9a99",
    "Brachiopoda":     "#fdbf6f",
    "Mollusca":        "#fdc086",
    "Acanthocephala":  "#b2abd2",
    # Ecdysozoa — blues and browns
    "Nematomorpha":    "#c7e9b4",
    "Nematoda":        "#beaed4",
    "Priapulida":      "#b2df8a",
    "Tardigrada":      "#41b6c4",
    "Collembola":      "#dfc27d",
    "Crustacea":       "#e6ab02",
    "Myriapoda":       "#bf812d",
    "Chelicerata":     "#1f78b4",
    "Insecta":         "#386cb0",
    # Deuterostomia — reds and greens
    "Echinodermata":   "#f0027f",
    "Hemichordata":    "#33a02c",
    "Chordata":        "#bf5b17",
}

# ── Supergroup membership ─────────────────────────────────────────────────────
SUPERGROUP_OF: dict[str, str] = {
    "Bacteria":        "Prokaryota",
    "Archaea":         "Prokaryota",
    "Oomycota":        "SAR",
    "Bacillariophyta": "SAR",
    "Apicomplexa":     "SAR",
    "Ciliophora":      "SAR",
    "Perkinsozoa":     "SAR",
    "Rhizaria":        "SAR",
    "Euglenozoa":      "Excavata",
    "Metamonada":      "Excavata",
    "Heterolobosea":   "Excavata",
    "Amoebozoa":       "Amoebozoa",
    "Viridiplantae":   "Archaeplastida",
    "Chlorophyta":     "Archaeplastida",
    "Rhodophyta":      "Archaeplastida",
    "Haptophyta":      "Other eukaryotes",
    "Cryptophyta":     "Other eukaryotes",
    "Protist":         "Other eukaryotes",
    "Fungi":           "Opisthokonta",
    **{p: "Metazoa" for p in [
        "Porifera", "Ctenophora", "Placozoa", "Xenacoelomorpha", "Cnidaria",
        "Platyhelminthes", "Rotifera", "Nemertea", "Annelida", "Brachiopoda",
        "Mollusca", "Acanthocephala", "Nematomorpha", "Nematoda", "Priapulida",
        "Tardigrada", "Collembola", "Crustacea", "Myriapoda", "Chelicerata",
        "Insecta", "Echinodermata", "Hemichordata", "Chordata",
    ]},
}

# ── Amino acid ordering and colours (for composition plots) ───────────────────
AA_ORDER = list("QNHKRDESTAGPVILMFYWC")

AA_COLOURS: dict[str, str] = {
    "Q": "#4dac26", "N": "#72c147",                         # polar uncharged
    "H": "#0571b0", "K": "#4393c3", "R": "#92c5de",         # positive
    "D": "#ca0020", "E": "#f4a582",                         # negative
    "S": "#f1a340", "T": "#fdae61",                         # polar hydroxyl
    "A": "#d9d9d9", "G": "#bdbdbd", "P": "#969696",         # small/special
    "V": "#e7d4e8", "I": "#c2a5cf", "L": "#9970ab",         # aliphatic
    "M": "#762a83",                                          # sulfur
    "F": "#1b7837", "Y": "#5aae61", "W": "#a6dba0",         # aromatic
    "C": "#fee08b",                                          # cysteine
}
