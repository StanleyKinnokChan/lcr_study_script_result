#!/usr/bin/env python3
"""
PGLS regression: N/C asymmetry ratio vs evolutionary tier in Viridiplantae.

Tests whether N-terminal LCR dominance increases progressively from algae →
charophytes → bryophytes → lycophytes → ferns → gymnosperms →
basal angiosperms → eudicots → grasses.

Taxonomic tier scheme (10 tiers, ordered by divergence time from the root):
  0  — Chlorophyta      (Chlamydomonas, Ostreococcus)
  1  — Charophytes      (Chara, Klebsormidium, Mesotaenium)
  2  — Bryophytes       (Physcomitrium, Marchantia, Sphagnum, Anthoceros)
  3  — Lycophytes       (Selaginella, Isoetes)
  4  — Ferns            (Azolla, Salvinia, Ceratopteris)          ← expansion
  5  — Gymnosperms      (Ginkgo, Picea, Pinus)                    ← expansion
  6  — Basal angiosperms (ANA grade; Amborella, Nymphaea, Magnolia)
  7  — Eudicots         (Arabidopsis, Beta, Brassica, Glycine, Solanum …)
  8  — Non-grass monocots (Asparagus, Ananas, Dioscorea, Musa)
  9  — Grasses          (Oryza, Triticum, Hordeum, Brachypodium …)

Data sources:
  EnsemblPlants r63 covers tiers 0, 1 (Chara only), 2 (partial), 3 (Selaginella
  only), 6 (partial), 7, 8, 9.  Tiers 4 and 5 are entirely absent; tiers 1-3
  are sparse (n=1-2 species each).

  Run scripts/01c_download_viridiplantae_expansion.py FIRST to download:
    - Charophytes: Klebsormidium nitens, Mesotaenium endlicherianum
    - Bryophytes:  Sphagnum fallax, Anthoceros agrestis
    - Lycophytes:  Isoetes taiwanensis
    - Ferns:       Azolla filiculoides, Salvinia cucullata, Ceratopteris richardii
    - Gymnosperms: Ginkgo biloba, Picea abies, Pinus taeda (large genomes)
    - Basal angiosperms: Magnolia biondii
  Then re-run 02_run_flps.sh and main.py --from-phase 2.

De-duplication:
  Ensembl Plants has many cultivars (50+ Hordeum, 20+ Triticum).  By default
  (--one-per-genus) one representative per genus is chosen (highest n_lcr) to
  avoid cultivar pseudoreplication in the OLS regression.

Method:
  Primary: PGLS under Brownian motion (requires Biopython + TimeTree Newick tree
           at results/viridiplantae_timetree.nwk).
  Fallback: Spearman rank correlation + OLS on tier rank when tree is absent.

Outputs:
  results/pgls_viridiplantae.tsv    — per-species data with tier assignment
  results/pgls_tier_summary.tsv     — mean ± SD asymmetry ratio per tier
  results/pgls_regression_stats.tsv — regression coefficients and p-value
  figures/suppfig_pgls.pdf
"""

import argparse
import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from pathlib import Path

from config import RESULTS_DIR, FIGURES_DIR

TREE_FILE = RESULTS_DIR / "viridiplantae_timetree.nwk"
ASYM_FILE = RESULTS_DIR / "asymmetry.tsv"

# Phyla included in this analysis.
# 01b_download_outgroups.py labels Chlorophyta separately from "Viridiplantae",
# so we filter on both.  Rhodophyta (Cyanidioschyzon, Galdieria, Chondrus) are
# excluded — they receive their own phylum label in the manifest.
VIRIDIPLANTAE_PHYLA = {"Viridiplantae", "Chlorophyta"}

# ── Tier assignment ───────────────────────────────────────────────────────────
# Keys are genera (first word of display_name, case-insensitive).
# Tier values match the scheme in the module docstring.
# Note: Rhodophyta genera (Cyanidioschyzon, Galdieria, Chondrus) are NOT
# Viridiplantae; they are present in EnsemblPlants but receive a different
# phylum label and will not appear in this analysis.

GENUS_TIER: dict[str, int] = {
    # Tier 0 — Chlorophyta (green algae; phylum="Chlorophyta" in manifest)
    "chlamydomonas":     0,
    "ostreococcus":      0,
    "volvox":            0,   # 01c expansion
    "micromonas":        0,   # 01c expansion
    "chlorella":         0,   # 01c expansion
    "coccomyxa":         0,   # 01c expansion

    # Tier 1 — Charophytes (streptophyte algae, sister to land plants)
    "chara":             1,
    "klebsormidium":     1,   # 01c expansion
    "closterium":        1,   # 01c expansion (Zygnematophyceae)
    "mesotaenium":       1,   # 01c (portal-only; kept for manual adds)

    # Tier 2 — Bryophytes
    "physcomitrium":     2,
    "marchantia":        2,
    "sphagnum":          2,
    "ceratodon":         2,   # 01c expansion
    "anthoceros":        2,   # 01c (hornwort; portal-only, kept for manual adds)
    "takakia":           2,   # 01c (portal-only, kept for manual adds)

    # Tier 3 — Lycophytes
    "selaginella":       3,
    "isoetes":           3,   # 01c (portal-only, kept for manual adds)
    "huperzia":          3,
    "diphasiastrum":     3,   # 01c expansion

    # Tier 4 — Ferns (Polypodiopsida)
    # Added via 01c_download_viridiplantae_expansion.py
    "azolla":            4,   # portal-only, kept for manual adds
    "salvinia":          4,   # portal-only, kept for manual adds
    "ceratopteris":      4,   # 01c expansion
    "adiantum":          4,   # 01c expansion

    # Tier 5 — Gymnosperms
    # Added via 01c_download_viridiplantae_expansion.py
    # (absent from Ensembl Plants r63; large genomes)
    "ginkgo":            5,
    "picea":             5,
    "pinus":             5,
    "abies":             5,
    "larix":             5,
    "cycas":             5,
    "taxus":             5,
    "cryptomeria":       5,   # 01c expansion

    # Tier 6 — Basal angiosperms (ANA grade)
    "amborella":         6,
    "nymphaea":          6,
    "magnolia":          6,   # from 01c expansion
    "victoria":          6,
    "liriodendron":      6,
    "cinnamomum":        6,   # 01c expansion (magnoliid)
    "aristolochia":      6,   # 01c expansion (magnoliid)

    # Tier 7 — Eudicots
    "actinidia":         7,
    "arabidopsis":       7,
    "arabis":            7,
    "arachis":           7,
    "beta":              7,
    "brassica":          7,
    "cajanus":           7,
    "camelina":          7,
    "cannabis":          7,
    "capsicum":          7,
    "chenopodium":       7,
    "citrullus":         7,
    "citrus":            7,
    "coffea":            7,
    "corchorus":         7,
    "corylus":           7,
    "corymbia":          7,
    "cucumis":           7,
    "cynara":            7,
    "daucus":            7,
    "eucalyptus":        7,
    "eutrema":           7,
    "ficus":             7,
    "fraxinus":          7,
    "glycine":           7,
    "gossypium":         7,
    "helianthus":        7,
    "ipomoea":           7,
    "juglans":           7,
    "kalanchoe":         7,
    "lablab":            7,
    "lactuca":           7,
    "lathyrus":          7,
    "lupinus":           7,
    "malus":             7,
    "manihot":           7,
    "medicago":          7,
    "nicotiana":         7,
    "olea":              7,
    "papaver":           7,
    "phaseolus":         7,
    "pistacia":          7,
    "pisum":             7,
    "populus":           7,
    "prunus":            7,
    "quercus":           7,
    "rosa":              7,
    "sesamum":           7,
    "solanum":           7,
    "sphenostylis":      7,
    "theobroma":         7,
    "trifolium":         7,
    "vicia":             7,
    "vigna":             7,
    "vitis":             7,

    # Tier 8 — Non-grass monocots
    "ananas":            8,
    "asparagus":         8,
    "dioscorea":         8,
    "musa":              8,
    "phoenix":           8,
    "vanilla":           8,

    # Tier 9 — Grasses (Poaceae)
    "aegilops":          9,
    "avena":             9,
    "brachypodium":      9,
    "digitaria":         9,
    "echinochloa":       9,
    "eragrostis":        9,
    "hordeum":           9,
    "leersia":           9,
    "lolium":            9,
    "oryza":             9,
    "panicum":           9,
    "saccharum":         9,
    "secale":            9,
    "setaria":           9,
    "sorghum":           9,
    "triticum":          9,
    "zea":               9,
}

TIER_LABELS = {
    0: "Chlorophyta",
    1: "Charophytes",
    2: "Bryophytes",
    3: "Lycophytes",
    4: "Ferns",
    5: "Gymnosperms",
    6: "Basal\nangiosperms",
    7: "Eudicots",
    8: "Non-grass\nmonocots",
    9: "Grasses",
}

TIER_COLOURS = {
    0: "#cccccc",
    1: "#a1d99b",
    2: "#74c476",
    3: "#41ab5d",
    4: "#238b45",
    5: "#005a32",
    6: "#9ecae1",
    7: "#4292c6",
    8: "#f16913",
    9: "#d94801",
}


def assign_tier(display_name: str) -> int | None:
    genus = display_name.split()[0].lower()
    return GENUS_TIER.get(genus)


def select_representatives(df: pd.DataFrame) -> pd.DataFrame:
    """Select one species per genus to avoid cultivar pseudoreplication."""
    df = df.copy()
    df["genus"] = df["display_name"].str.split().str[0].str.lower()
    # Pick the species with most LCRs per genus (most data = most reliable asymmetry ratio)
    reps = (
        df.sort_values("n_lcr", ascending=False)
        .drop_duplicates(subset="genus", keep="first")
    )
    return reps.drop(columns="genus")


def pgls_brownian(trait: np.ndarray, predictor: np.ndarray,
                  vcv: np.ndarray) -> dict:
    """PGLS under Brownian motion using the phylogenetic variance-covariance matrix."""
    n = len(trait)
    try:
        from scipy.linalg import inv
        V = np.asarray(vcv, dtype=float)
        V_inv = inv(V)
        X = np.column_stack([np.ones(n), predictor])
        XtVi = X.T @ V_inv
        beta = np.linalg.solve(XtVi @ X, XtVi @ trait)
        residuals = trait - X @ beta
        sigma2 = (residuals @ V_inv @ residuals) / (n - 2)
        cov_beta = sigma2 * np.linalg.inv(XtVi @ X)
        se_slope = np.sqrt(cov_beta[1, 1])
        t_stat = beta[1] / se_slope
        from scipy.stats import t as t_dist
        p_value = 2 * t_dist.sf(abs(t_stat), df=n - 2)
        ss_res = residuals @ V_inv @ residuals
        y_mean = np.sum(V_inv @ trait) / np.sum(V_inv)
        ss_tot = (trait - y_mean) @ V_inv @ (trait - y_mean)
        r2 = 1 - ss_res / ss_tot
        return {
            "method": "PGLS_Brownian",
            "intercept": round(float(beta[0]), 4),
            "slope": round(float(beta[1]), 4),
            "se_slope": round(float(se_slope), 4),
            "t_stat": round(float(t_stat), 3),
            "p_value": round(float(p_value), 6),
            "r2_pgls": round(float(r2), 4),
            "n": n,
        }
    except Exception as e:
        warnings.warn(f"PGLS failed: {e}", RuntimeWarning)
        return {}


def load_tree_vcv(tree_file: Path, plants_df: pd.DataFrame):
    """
    Build the Brownian-motion variance-covariance matrix (VCV) for the analysis
    species from a Newick tree (e.g. a TimeTree "Load a List of Species" export
    of results/pgls_species_list.txt).

    Tree leaves are matched to analysis species by genus — the analysis uses one
    representative per genus, so a leaf named 'Avena_sterilis' maps to the Avena
    representative regardless of the exact species/cultivar TimeTree returned.

    For a dated (ultrametric) tree the correct Brownian VCV is
        V[i,j] = root-to-MRCA(i,j) path length   (shared evolutionary history)
        V[i,i] = root-to-tip path length,
    which is the proper covariance structure for PGLS — not the raw pairwise
    distance matrix.

    Returns (vcv, positions) where positions are row indices into
    plants_df.reset_index(drop=True) for the matched species (VCV rows/cols are
    in that order), or (None, None) if the tree is absent/unreadable or matches
    fewer than 5 genera.
    """
    if not tree_file.exists():
        return None, None
    try:
        from Bio import Phylo
        tree = Phylo.read(str(tree_file), "newick")
    except ImportError:
        warnings.warn("Biopython not installed — tree-based PGLS unavailable.", UserWarning)
        return None, None
    except Exception as e:
        warnings.warn(f"Could not read tree {tree_file}: {e}", UserWarning)
        return None, None

    root = tree.root
    leaf_by_genus: dict[str, object] = {}
    for leaf in tree.get_terminals():
        genus = (leaf.name or "").replace(" ", "_").split("_")[0].lower()
        if genus and genus not in leaf_by_genus:
            leaf_by_genus[genus] = leaf

    df = plants_df.reset_index(drop=True)
    genera = df["display_name"].str.split().str[0].str.lower().tolist()
    matched = [(i, leaf_by_genus[g]) for i, g in enumerate(genera) if g in leaf_by_genus]
    if len(matched) < 5:
        warnings.warn(
            f"Tree matched only {len(matched)} analysis genera — too few for PGLS; "
            "check that tree leaf names are genus-resolvable.", UserWarning)
        return None, None

    positions = [i for i, _ in matched]
    leaves    = [lf for _, lf in matched]
    n = len(leaves)
    depth = {id(lf): tree.distance(root, lf) for lf in leaves}
    vcv = np.zeros((n, n))
    for a in range(n):
        vcv[a, a] = depth[id(leaves[a])]
        for b in range(a + 1, n):
            mrca = tree.common_ancestor([leaves[a], leaves[b]])
            shared = tree.distance(root, mrca)
            vcv[a, b] = vcv[b, a] = shared
    print(f"  Tree matched {n}/{len(df)} analysis species by genus.")
    return vcv, positions


def main(one_per_genus: bool = True):
    if not ASYM_FILE.exists():
        print(f"ERROR: {ASYM_FILE} missing — run 05_asymmetry.py first.")
        return

    asym_df = pd.read_csv(ASYM_FILE, sep="\t")
    plants = asym_df[asym_df["phylum"].isin(VIRIDIPLANTAE_PHYLA)].copy()

    if plants.empty:
        print("No Viridiplantae found in asymmetry.tsv.")
        return

    # Assign tiers
    plants["tier"] = plants["display_name"].apply(assign_tier)
    unrecognised = plants[plants["tier"].isna()]["display_name"].tolist()
    if unrecognised:
        print(f"WARNING: {len(unrecognised)} species without tier assignment "
              f"(will be excluded):")
        for name in unrecognised[:20]:
            print(f"  {name}")
        if len(unrecognised) > 20:
            print(f"  ... and {len(unrecognised) - 20} more")

    plants = plants.dropna(subset=["tier", "asymmetry_ratio"])
    plants["tier"] = plants["tier"].astype(int)

    # Remove species with unstable asymmetry ratio (too few C-terminal LCRs)
    MIN_CTERM = 10
    plants = plants[plants["n_cterm"] >= MIN_CTERM].copy()
    print(f"\nViridiplantae species with ≥{MIN_CTERM} C-terminal LCRs and tier assignment: "
          f"{len(plants)}")

    # Tier composition before de-duplication
    print("\nTier composition (all species):")
    for tier in sorted(plants["tier"].unique()):
        sp_in_tier = plants[plants["tier"] == tier]
        print(f"  Tier {tier} — {TIER_LABELS.get(tier, '?'):<22}: "
              f"{len(sp_in_tier):3d} species   "
              f"(genera: {', '.join(sorted(sp_in_tier['display_name'].str.split().str[0].unique()[:6]))})")

    # De-duplicate: one representative per genus
    if one_per_genus:
        plants_ols = select_representatives(plants)
        print(f"\nAfter selecting one representative per genus: {len(plants_ols)} species")
        print("Tier composition (de-duplicated):")
        for tier in sorted(plants_ols["tier"].unique()):
            sp_in_tier = plants_ols[plants_ols["tier"] == tier]
            genera = sorted(sp_in_tier["display_name"].str.split().str[0].unique())
            print(f"  Tier {tier} — {TIER_LABELS.get(tier, '?'):<22}: "
                  f"{len(sp_in_tier):3d} species   "
                  f"({', '.join(genera)})")
    else:
        plants_ols = plants.copy()
        print(f"\nUsing all {len(plants_ols)} species (no genus de-duplication).")

    if len(plants_ols) < 5:
        print(f"\nToo few species for regression (n={len(plants_ols)}). "
              f"Check that 05_asymmetry.py ran on Viridiplantae.")
        return

    predictor = plants_ols["tier"].values
    trait = plants_ols["asymmetry_ratio"].values

    # ── Attempt PGLS ──────────────────────────────────────────────────────────
    vcv, matched_pos = load_tree_vcv(TREE_FILE, plants_ols)
    use_pgls = vcv is not None

    if use_pgls:
        plants_pgls = plants_ols.reset_index(drop=True).iloc[matched_pos]
        predictor_pgls = plants_pgls["tier"].values.astype(float)
        trait_pgls = plants_pgls["asymmetry_ratio"].values
        print(f"\nUsing phylogenetic tree for PGLS ({len(matched_pos)} species on tree).")
        stats = pgls_brownian(trait_pgls, predictor_pgls, vcv)
        print(f"  PGLS slope={stats.get('slope')}, p={stats.get('p_value')}, "
              f"R²={stats.get('r2_pgls')}")
    else:
        # ── OLS + Spearman fallback ───────────────────────────────────────────
        if not TREE_FILE.exists():
            print(f"\nTree file not found at {TREE_FILE}.")
        print("Falling back to OLS on taxonomic tier.")
        r, p_r = pearsonr(predictor, trait)
        rho, p_rho = spearmanr(predictor, trait)
        coeffs = np.polyfit(predictor, trait, 1)
        stats = {
            "method": "OLS_taxonomic_tier_fallback",
            "slope": round(float(coeffs[0]), 4),
            "intercept": round(float(coeffs[1]), 4),
            "pearson_r": round(float(r), 4),
            "pearson_p": round(float(p_r), 6),
            "spearman_rho": round(float(rho), 4),
            "spearman_p": round(float(p_rho), 6),
            "n": len(plants_ols),
            "one_per_genus": one_per_genus,
            "note": (
                "Approximate — provide results/viridiplantae_timetree.nwk for full PGLS. "
                "Upload results/pgls_species_list.txt to timetree.org "
                "('Load a List of Species'), download the Newick, and save it as "
                "results/viridiplantae_timetree.nwk"
            ),
        }
        print(f"  Pearson  r  = {r:.3f}, p = {p_r:.5f}")
        print(f"  Spearman ρ  = {rho:.3f}, p = {p_rho:.5f}")
        print(f"\n  Note: gymnosperms are absent from Ensembl Plants r63.")
        print("  To add full PGLS: upload results/pgls_species_list.txt to")
        print("  timetree.org ('Load a List of Species'), then save the downloaded")
        print("  Newick as results/viridiplantae_timetree.nwk and re-run this script.\n")

    # ── Tier-level summary ────────────────────────────────────────────────────
    tier_summary_rows = []
    for tier in sorted(plants_ols["tier"].unique()):
        sub = plants_ols[plants_ols["tier"] == tier]["asymmetry_ratio"]
        tier_summary_rows.append({
            "tier":          tier,
            # TIER_LABELS carries newlines for two-line plot labels; flatten them
            # so the TSV stays one row per tier.
            "tier_label":    TIER_LABELS.get(tier, f"Tier{tier}").replace("\n", " "),
            "n_species":     len(sub),
            "mean_ratio":    round(float(sub.mean()), 4),
            "sd_ratio":      round(float(sub.std()), 4),
            "median_ratio":  round(float(sub.median()), 4),
        })
    tier_df = pd.DataFrame(tier_summary_rows)
    tier_df.to_csv(RESULTS_DIR / "pgls_tier_summary.tsv", sep="\t", index=False)
    print(f"Tier summary: {RESULTS_DIR / 'pgls_tier_summary.tsv'}")
    print(tier_df.to_string(index=False))

    # ── Write full per-species table ──────────────────────────────────────────
    plants_out = plants_ols[["species_key", "display_name", "n_lcr",
                              "n_nterm", "n_cterm", "pct_nterm", "pct_cterm",
                              "asymmetry_ratio", "tier"]].copy()
    plants_out["tier_label"] = (plants_out["tier"].map(TIER_LABELS)
                                .str.replace("\n", " ", regex=False))
    plants_out.to_csv(RESULTS_DIR / "pgls_viridiplantae.tsv", sep="\t", index=False)

    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(RESULTS_DIR / "pgls_regression_stats.tsv", sep="\t", index=False)
    print(f"\nPer-species data: {RESULTS_DIR / 'pgls_viridiplantae.tsv'}")
    print(f"Regression stats: {RESULTS_DIR / 'pgls_regression_stats.tsv'}")

    # ── Figure ────────────────────────────────────────────────────────────────
    tiers_present = sorted(plants_ols["tier"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: individual species scatter + OLS line
    ax = axes[0]
    for tier in tiers_present:
        sub = plants_ols[plants_ols["tier"] == tier]
        jitter = np.random.default_rng(tier).uniform(-0.18, 0.18, len(sub))
        ax.scatter(
            np.full(len(sub), tier) + jitter,
            sub["asymmetry_ratio"],
            c=TIER_COLOURS.get(tier, "#aaa"),
            s=60, alpha=0.8, edgecolors="black", linewidths=0.4, zorder=5,
        )
    # Regression line
    x_line = np.linspace(min(tiers_present), max(tiers_present), 200)
    if not use_pgls:
        ax.plot(x_line, np.polyval(coeffs, x_line),
                color="black", linewidth=1.5, linestyle="--",
                label=f"OLS (Pearson r={stats['pearson_r']}, p={stats['pearson_p']:.4f})")
    else:
        b0, b1 = stats["intercept"], stats["slope"]
        ax.plot(x_line, b0 + b1 * x_line,
                color="black", linewidth=1.5, linestyle="--",
                label=f"PGLS slope={b1}, p={stats.get('p_value', 'NA')}")
    ax.axhline(1.0, color="red", linestyle=":", linewidth=1.2, label="Symmetry (ratio=1)")
    ax.set_xticks(tiers_present)
    ax.set_xticklabels(
        [TIER_LABELS.get(t, str(t)).replace("\n", " ") for t in tiers_present],
        rotation=30, ha="right", fontsize=8,
    )
    ax.set_ylabel("N/C asymmetry ratio (pct_nterm / pct_cterm)", fontsize=10)
    ax.set_title("Individual species (one per genus)", fontsize=10)
    ax.legend(fontsize=8)

    # Panel B: mean ± SD per tier (bar chart)
    ax = axes[1]
    x = np.arange(len(tier_df))
    bar_colours = [TIER_COLOURS.get(t, "#aaa") for t in tier_df["tier"]]
    ax.bar(x, tier_df["mean_ratio"], yerr=tier_df["sd_ratio"],
           color=bar_colours, edgecolor="black", linewidth=0.6,
           capsize=4, error_kw={"linewidth": 1.2, "ecolor": "black"})
    for i, row in tier_df.iterrows():
        ax.text(i, 0.02, f"n={row['n_species']}", ha="center", va="bottom",
                fontsize=7, color="white", fontweight="bold")
    ax.axhline(1.0, color="red", linestyle=":", linewidth=1.2, label="Symmetry (ratio=1)")
    ax.set_xticks(x)
    ax.set_xticklabels(
        [TIER_LABELS.get(t, str(t)).replace("\n", " ") for t in tier_df["tier"]],
        rotation=30, ha="right", fontsize=8,
    )
    ax.set_ylabel("Mean N/C asymmetry ratio", fontsize=10)
    ax.set_title("Tier means ± SD", fontsize=10)
    ax.legend(fontsize=8)

    method_note = ("PGLS" if use_pgls
                   else "OLS on tier rank (PGLS requires viridiplantae_timetree.nwk)")
    n_tiers = len(tier_df)
    has_ferns = 4 in tier_df["tier"].values
    has_gymno = 5 in tier_df["tier"].values
    missing = []
    if not has_ferns:
        missing.append("ferns")
    if not has_gymno:
        missing.append("gymnosperms")
    missing_note = (f"run 01c to add: {', '.join(missing)}" if missing
                    else "all 10 tiers present")
    plt.suptitle(
        f"N-terminal LCR asymmetry across Viridiplantae evolutionary tiers "
        f"({n_tiers}/10 tiers present)\n"
        f"({method_note}; {missing_note})",
        fontsize=11,
    )
    plt.tight_layout()

    out_fig = FIGURES_DIR / "suppfig_pgls.pdf"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    fig.savefig(str(out_fig).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nPGLS figure: {out_fig}")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all-species", action="store_true",
        help="Use all species without de-duplicating to one per genus. "
             "Not recommended for PGLS due to cultivar pseudoreplication.",
    )
    args = parser.parse_args()
    main(one_per_genus=not args.all_species)
