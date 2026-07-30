#!/usr/bin/env python3
"""
GO term enrichment analysis for terminal-LCR-bearing proteins.

For each target species, downloads protein → GO term mappings from Ensembl
Metazoa BioMart, then tests which GO terms (Biological Process, Molecular
Function, Cellular Component) are enriched among proteins with ≥1 terminal
LCR vs all proteins that carry any LCR.

TARGET_SPECIES covers every Ensembl Metazoa species in our pipeline that plausibly
has BioMart GO annotation. BioMart dataset names follow the convention:
  {first_letter_genus}{species_epithet}_eg_gene
For species with assembly suffixes in the species key (e.g. aedes_aegypti_lvpagwg),
only the first two name parts are used for the dataset name (aaegypti_eg_gene).
Species not present in BioMart are skipped gracefully via the diagnostic checks.

GO data is cached locally so re-runs are fast; delete data/go_cache/ to force refresh.

Outputs:
  results/go_enrichment.tsv      — per-species enriched GO terms (FDR ≤ 0.05)
  results/go_consensus.tsv       — GO terms enriched in ≥3 species
  figures/go_enrichment.pdf  (supports the Discussion text on GO enrichment,
                               not a numbered figure)
"""

import time
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import fisher_exact

try:
    import requests
except ImportError:
    raise SystemExit("pip install requests")

from config import PROJECT_DIR, RESULTS_DIR, FIGURES_DIR

CACHE_DIR = PROJECT_DIR / "data" / "go_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

BIOMART_URL = "https://metazoa.ensembl.org/biomart/martservice"

# All Ensembl Metazoa species in our pipeline with plausible BioMart datasets.
# Dataset name convention: {first_letter_genus}{species_name}_eg_gene
# Species absent from BioMart are skipped at runtime via diagnostic check.
TARGET_SPECIES = {
    # ── Insecta ──────────────────────────────────────────────────────────────
    "drosophila_melanogaster":              ("dmelanogaster_eg_gene",    "Insecta"),
    "anopheles_gambiae":                    ("agambiae_eg_gene",         "Insecta"),
    "aedes_aegypti_lvpagwg":              ("aaegypti_eg_gene",          "Insecta"),
    "culex_quinquefasciatus_jhu":           ("cquinquefasciatus_eg_gene","Insecta"),
    "musca_domestica_gca030504385v2rs":   ("mdomestica_eg_gene",        "Insecta"),
    "glossina_morsitans_morsitans":         ("gmorsitans_eg_gene",       "Insecta"),
    "bombyx_mori":                          ("bmori_eg_gene",            "Insecta"),
    "danaus_plexippus":                     ("dplexippus_eg_gene",       "Insecta"),
    "heliconius_melpomene":                 ("hmelpomene_eg_gene",       "Insecta"),
    "manduca_sexta_gca014839925v2":         ("msexta_eg_gene",           "Insecta"),
    "tribolium_castaneum":                  ("tcastaneum_eg_gene",       "Insecta"),
    "onthophagus_taurus_gca002045055v1":    ("otaurus_eg_gene",          "Insecta"),
    "apis_mellifera":                       ("amellifera_eg_gene",       "Insecta"),
    "nasonia_vitripennis":                  ("nvitripennis_eg_gene",     "Insecta"),
    "bombus_terrestris_gca000214255v1":     ("bterrestris_eg_gene",      "Insecta"),
    "linepithema_humile_gca000217595v1":    ("lhumile_eg_gene",          "Insecta"),
    "harpegnathos_saltator_gca000217535v2": ("hsaltator_eg_gene",        "Insecta"),
    "pediculus_humanus":                    ("phumanus_eg_gene",         "Insecta"),
    "acyrthosiphon_pisum":                  ("apisum_eg_gene",           "Insecta"),
    "nilaparvata_lugens_gca014356525v1rs":  ("nlugens_eg_gene",         "Insecta"),
    "rhodnius_prolixus":                    ("rprolixus_eg_gene",        "Insecta"),
    "zootermopsis_nevadensis":              ("znevadensis_eg_gene",      "Insecta"),
    "locusta_migratoria_gca023397235v2":    ("lmigratoria_eg_gene",      "Insecta"),
    "frankliniella_occidentalis_gca000696945v1rs": ("foccidentalis_eg_gene", "Insecta"),
    "blattella_germanica_gca009764705v2":   ("bgermanica_eg_gene",       "Insecta"),

    # ── Chelicerata ───────────────────────────────────────────────────────────
    "ixodes_scapularis_ise6":               ("iscapularis_eg_gene",      "Chelicerata"),
    "rhipicephalus_microplus_gca013339725v2": ("rmicroplus_eg_gene",     "Chelicerata"),
    "tetranychus_urticae":                  ("turticae_eg_gene",         "Chelicerata"),
    "varroa_destructor_gca002443255v1":     ("vdestructor_eg_gene",      "Chelicerata"),
    "limulus_polyphemus_gca000517525v1":    ("lpolyphemus_eg_gene",      "Chelicerata"),
    "parasteatoda_tepidariorum_gca043381705v1rs": ("ptepidariorum_eg_gene", "Chelicerata"),
    "stegodyphus_dumicola_gca010614865v2":  ("sdumicola_eg_gene",        "Chelicerata"),
    "latrodectus_hesperus_gca000350385v1":  ("lhesperus_eg_gene",        "Chelicerata"),
    "centruroides_sculpturatus_gca000671375v2": ("csculpturatus_eg_gene","Chelicerata"),
    "mesobuthus_martensii_gca000695325v1":  ("mmartensii_eg_gene",       "Chelicerata"),

    # ── Crustacea ─────────────────────────────────────────────────────────────
    "daphnia_pulex_gca021134715v1rs":       ("dpulex_eg_gene",           "Crustacea"),
    "artemia_franciscana_gca030216415v1":   ("afranciscana_eg_gene",     "Crustacea"),
    "penaeus_vannamei_gca042767895v1rs":    ("pvannamei_eg_gene",        "Crustacea"),
    "lepeophtheirus_salmonis_gca000181255v2": ("lsalmonis_eg_gene",      "Crustacea"),

    # ── Myriapoda ─────────────────────────────────────────────────────────────
    "strigamia_maritima":                   ("smaritima_eg_gene",        "Myriapoda"),

    # ── Nematoda ──────────────────────────────────────────────────────────────
    "caenorhabditis_elegans":               ("celegans_eg_gene",         "Nematoda"),
    "caenorhabditis_briggsae":              ("cbriggsae_eg_gene",        "Nematoda"),
    "caenorhabditis_brenneri":              ("cbrenneri_eg_gene",        "Nematoda"),
    "caenorhabditis_remanei":               ("cremanei_eg_gene",         "Nematoda"),
    "pristionchus_pacificus":               ("ppacificus_eg_gene",       "Nematoda"),
    "haemonchus_contortus":                 ("hcontortus_eg_gene",       "Nematoda"),
    "strongyloides_ratti":                  ("sratti_eg_gene",           "Nematoda"),
    "brugia_malayi":                        ("bmalayi_eg_gene",          "Nematoda"),
    "ascaris_suum":                         ("asuum_eg_gene",            "Nematoda"),
    "trichinella_spiralis":                 ("tspiralis_eg_gene",        "Nematoda"),
    "meloidogyne_hapla":                    ("mhapla_eg_gene",           "Nematoda"),
    "globodera_pallida_gca000724045v1":     ("gpallida_eg_gene",         "Nematoda"),

    # ── Platyhelminthes ───────────────────────────────────────────────────────
    "schistosoma_mansoni":                  ("smansoni_eg_gene",         "Platyhelminthes"),
    "schistosoma_haematobium_gca000699445v3": ("shaematobium_eg_gene",   "Platyhelminthes"),
    "echinococcus_granulosus":              ("egranulosus_eg_gene",      "Platyhelminthes"),
    "echinococcus_multilocularis":          ("emultilocularis_eg_gene",  "Platyhelminthes"),
    "clonorchis_sinensis_gca000236345v2":   ("csinensis_eg_gene",        "Platyhelminthes"),
    "opisthorchis_viverrini_gca000715545v1": ("oviverrini_eg_gene",     "Platyhelminthes"),
    "schmidtea_mediterranea_gca045838265v1cm": ("smediterranea_eg_gene", "Platyhelminthes"),
    "macrostomum_lignano_gca002269645v3":   ("mlignano_eg_gene",         "Platyhelminthes"),
    "taenia_solium":                        ("tsolium_eg_gene",          "Platyhelminthes"),
    "hymenolepis_microstoma_gca000469805v3": ("hmicrostoma_eg_gene",    "Platyhelminthes"),

    # ── Annelida ──────────────────────────────────────────────────────────────
    "capitella_teleta":                     ("cteleta_eg_gene",          "Annelida"),
    "platynereis_dumerilii_gca026936325v1cm": ("pdumerilii_eg_gene",    "Annelida"),
    "helobdella_robusta":                   ("hrobusta_eg_gene",         "Annelida"),
    "hirudo_medicinalis_gca900006815v1":    ("hmedicinalis_eg_gene",     "Annelida"),

    # ── Mollusca ──────────────────────────────────────────────────────────────
    "lottia_gigantea":                      ("lgigantea_eg_gene",        "Mollusca"),
    "aplysia_californica_gca000002075v2":   ("acalifornica_eg_gene",     "Mollusca"),
    "biomphalaria_glabrata_gca000457545v3": ("bglabrata_eg_gene",        "Mollusca"),
    "patella_vulgata_gca900773765v1":       ("pvulgata_eg_gene",         "Mollusca"),
    "magallana_gigas_gca963853765v1rs":     ("mgigas_eg_gene",           "Mollusca"),
    "crassostrea_virginica_gca002022765v4": ("cvirginica_eg_gene",       "Mollusca"),
    "mizuhopecten_yessoensis_gca002113885v2": ("myessoensis_eg_gene",   "Mollusca"),
    "pecten_maximus_gca902652985v1":        ("pmaximus_eg_gene",         "Mollusca"),
    "octopus_bimaculoides_gca001194135v2rs": ("obimaculoides_eg_gene",  "Mollusca"),
    "sepia_officinalis_gca027596085v2":     ("sofficinalis_eg_gene",     "Mollusca"),

    # ── Echinodermata ─────────────────────────────────────────────────────────
    "strongylocentrotus_purpuratus":        ("spurpuratus_eg_gene",      "Echinodermata"),
    "lytechinus_variegatus_gca018143015v1": ("lvariegatus_eg_gene",      "Echinodermata"),
    "lytechinus_pictus_gca037042905v1rs":   ("lpictus_eg_gene",          "Echinodermata"),
    "patiria_miniata_gca015706575v1":       ("pminiata_eg_gene",         "Echinodermata"),
    "asterias_rubens_gca902459465v3":       ("arubens_eg_gene",          "Echinodermata"),
    "acanthaster_planci_gca001949145v1":    ("aplanci_eg_gene",          "Echinodermata"),
    "apostichopus_japonicus_gca001835935v2": ("ajaponicus_eg_gene",      "Echinodermata"),
    "ophioderma_brevispinum_gca016038115v1": ("obrevispinum_eg_gene",    "Echinodermata"),

    # ── Cnidaria ──────────────────────────────────────────────────────────────
    "nematostella_vectensis":               ("nvectensis_eg_gene",       "Cnidaria"),
    "hydra_vulgaris_gca038396675v1rs":      ("hvulgaris_eg_gene",        "Cnidaria"),
    "acropora_millepora_gca013753865v1":    ("amillepora_eg_gene",       "Cnidaria"),
    "acropora_digitifera_gca000222465v2":   ("adigitfera_eg_gene",       "Cnidaria"),
    "orbicella_faveolata_gca002042975v1":   ("ofaveolata_eg_gene",       "Cnidaria"),
    "pocillopora_damicornis_gca003704095v1": ("pdamicornis_eg_gene",     "Cnidaria"),
    "stylophora_pistillata_gca002571385v1": ("spistillata_eg_gene",      "Cnidaria"),

    # ── Placozoa ──────────────────────────────────────────────────────────────
    "trichoplax_adhaerens":                 ("tadhaerens_eg_gene",       "Placozoa"),

    # ── Porifera ──────────────────────────────────────────────────────────────
    "amphimedon_queenslandica":             ("aqueenslandica_eg_gene",   "Porifera"),

    # ── Ctenophora ────────────────────────────────────────────────────────────
    "mnemiopsis_leidyi":                    ("mleidyi_eg_gene",          "Ctenophora"),

    # ── Hemichordata ──────────────────────────────────────────────────────────
    "saccoglossus_kowalevskii_gca000003605v1": ("skowalevskii_eg_gene", "Hemichordata"),

    # ── Cephalochordata ───────────────────────────────────────────────────────
    "branchiostoma_lanceolatum":            ("blanceolatum_eg_gene",     "Cephalochordata"),
    "branchiostoma_floridae_gca000003815v2": ("bfloridae_eg_gene",      "Cephalochordata"),

    # ── Brachiopoda ───────────────────────────────────────────────────────────
    "lingula_anatina_gca001039355v2":       ("lanatina_eg_gene",         "Brachiopoda"),
}

# GO namespace short labels
NS_LABELS = {
    "biological_process":  "BP",
    "molecular_function":  "MF",
    "cellular_component":  "CC",
}

FDR_THRESHOLD   = 0.05
MIN_TERM_COUNT  = 5    # minimum proteins annotated to a GO term to test it


def biomart_query(dataset: str) -> pd.DataFrame | None:
    """
    Download protein_stable_id → go_id + go_name + namespace from BioMart.
    Returns a DataFrame or None on failure. Caches to disk.
    """
    cache_file = CACHE_DIR / f"{dataset}_go.tsv"
    if cache_file.exists():
        df = pd.read_csv(cache_file, sep="\t")
        print(f"  Loaded from cache: {cache_file.name}  ({len(df)} rows)")
        return df

    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="metazoa_mart" formatter="TSV" header="1" uniqueRows="1" count="">
  <Dataset name="{dataset}" interface="default">
    <Attribute name="ensembl_peptide_id" />
    <Attribute name="go_id" />
    <Attribute name="name_1006" />
    <Attribute name="namespace_1003" />
  </Dataset>
</Query>"""

    try:
        r = requests.get(
            BIOMART_URL,
            params={"query": xml},
            timeout=300,
            headers={"User-Agent": "lcr-invertebrate-study/1.0"},
        )
        r.raise_for_status()
    except Exception as e:
        print(f"  [WARN] BioMart request failed for {dataset}: {e}")
        return None

    text = r.text.strip()
    if not text or text.startswith("ERROR"):
        print(f"  [WARN] BioMart returned error for {dataset}: {text[:200]}")
        return None

    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        print(f"  [WARN] No data returned for {dataset}")
        return None

    rows = []
    header = lines[0].split("\t")
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) == len(header):
            rows.append(parts)

    df = pd.DataFrame(rows, columns=header)
    df.columns = ["protein_id", "go_id", "go_name", "namespace"]
    df = df[df["go_id"].str.startswith("GO:", na=False)]
    df = df.dropna(subset=["protein_id", "go_id"])
    df.to_csv(cache_file, sep="\t", index=False)
    print(f"  Downloaded {len(df)} protein→GO rows, cached to {cache_file.name}")
    time.sleep(1)  # be polite to BioMart
    return df


def bh_correction(pvalues: pd.Series) -> pd.Series:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvalues)
    ranked = pvalues.rank(method="first")
    fdr = pvalues * n / ranked
    # Enforce monotonicity from right
    fdr_values = fdr.values.copy()
    for i in range(n - 2, -1, -1):
        fdr_values[i] = min(fdr_values[i], fdr_values[i + 1])
    return pd.Series(np.minimum(fdr_values, 1.0), index=pvalues.index)


def enrich_go(go_df: pd.DataFrame, terminal_proteins: set, all_lcr_proteins: set) -> pd.DataFrame:
    """
    Fisher's exact test per GO term:
      foreground = terminal-LCR proteins with this GO term
      background = all-LCR proteins with this GO term
    """
    go_df = go_df[go_df["protein_id"].isin(all_lcr_proteins)].copy()
    go_df["is_terminal"] = go_df["protein_id"].isin(terminal_proteins)

    n_terminal_total = len(terminal_proteins & set(go_df["protein_id"]))
    n_internal_total = len(all_lcr_proteins & set(go_df["protein_id"])) - n_terminal_total

    rows = []
    for (go_id, go_name, namespace), grp in go_df.groupby(["go_id", "go_name", "namespace"]):
        n_term_with = grp["is_terminal"].sum()
        n_all_with  = len(grp)
        if n_all_with < MIN_TERM_COUNT:
            continue
        n_term_without = n_terminal_total - n_term_with
        n_int_with     = n_all_with - n_term_with
        n_int_without  = n_internal_total - n_int_with

        table = [[n_term_with, n_term_without],
                 [n_int_with,  n_int_without]]
        try:
            or_, p = fisher_exact(table, alternative="greater")
        except Exception:
            continue

        rows.append({
            "go_id":         go_id,
            "go_name":       go_name,
            "namespace":     namespace,
            "ns_short":      NS_LABELS.get(namespace, namespace),
            "n_terminal":    int(n_term_with),
            "n_all_lcr":     int(n_all_with),
            "pct_terminal":  round(n_term_with / n_all_with * 100, 1) if n_all_with else None,
            "odds_ratio":    round(or_, 3),
            "pvalue":        round(p, 6),
        })

    if not rows:
        return pd.DataFrame()

    result = pd.DataFrame(rows).sort_values("pvalue")
    result["fdr"] = bh_correction(result["pvalue"]).round(4).values
    return result[result["fdr"] <= FDR_THRESHOLD].sort_values("fdr")


def main():
    pos_path = RESULTS_DIR / "lcr_positions.tsv"
    if not pos_path.exists():
        print(f"ERROR: {pos_path} missing — run 03_analyse.py first.")
        return

    pos_df = pd.read_csv(pos_path, sep="\t")
    # Only attempt species that are actually in the positions table
    pipeline_species = set(pos_df["species_key"].unique())
    species_to_run = {k: v for k, v in TARGET_SPECIES.items() if k in pipeline_species}
    skipped = set(TARGET_SPECIES) - pipeline_species
    if skipped:
        print(f"Skipping {len(skipped)} species not in lcr_positions.tsv "
              f"(not yet downloaded/analysed):")
        for s in sorted(skipped):
            print(f"  {s}")

    print(f"\nRunning GO enrichment for {len(species_to_run)} species\n")

    all_results  = []
    go_term_hits = {}   # go_id → list of species where it's enriched

    for sp_key, (dataset, phylum) in species_to_run.items():
        print(f"\n[{phylum}] {sp_key}")

        sp_pos = pos_df[pos_df["species_key"] == sp_key]

        # Protein sets
        terminal_proteins = set(sp_pos[sp_pos["is_terminal"]]["protein_id"])
        all_lcr_proteins  = set(sp_pos["protein_id"])
        print(f"  Proteins with ≥1 LCR: {len(all_lcr_proteins)}  "
              f"(terminal: {len(terminal_proteins)})")

        # Fetch GO annotations
        go_df = biomart_query(dataset)
        if go_df is None or go_df.empty:
            print(f"  Skipping GO enrichment — no annotations available.")
            continue

        # Diagnostic: check protein ID overlap between BioMart and LCR positions
        biomart_ids = set(go_df["protein_id"].dropna())
        overlap = biomart_ids & all_lcr_proteins
        print(f"  BioMart protein IDs: {len(biomart_ids)} unique  |  "
              f"LCR proteins: {len(all_lcr_proteins)}  |  "
              f"Overlap: {len(overlap)}")
        if len(overlap) == 0:
            print(f"  [WARN] Zero protein ID overlap — BioMart IDs don't match FASTA "
                  f"headers. Example BioMart ID: {next(iter(biomart_ids), 'none')!r}  "
                  f"Example LCR ID: {next(iter(all_lcr_proteins), 'none')!r}")
            print(f"  Skipping GO enrichment for {sp_key}.")
            continue

        # Run enrichment
        enr = enrich_go(go_df, terminal_proteins, all_lcr_proteins)
        if enr.empty:
            print(f"  No GO terms enriched at FDR ≤ {FDR_THRESHOLD}.")
            continue

        enr["species_key"] = sp_key
        enr["phylum"]      = phylum
        all_results.append(enr)

        print(f"  Enriched GO terms (FDR ≤ {FDR_THRESHOLD}): {len(enr)}")
        print(f"  Top 5:")
        for _, row in enr.head(5).iterrows():
            print(f"    [{row.ns_short}] {row.go_id} {row.go_name[:50]:<50}  "
                  f"OR={row.odds_ratio}  FDR={row.fdr}")

        for go_id in enr["go_id"]:
            go_term_hits.setdefault(go_id, []).append(sp_key)

    if not all_results:
        print("\nNo enrichment results to write. Check BioMart connectivity.")
        return

    full_df = pd.concat(all_results, ignore_index=True)
    out_tsv = RESULTS_DIR / "go_enrichment.tsv"
    full_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nGO enrichment table: {out_tsv}  ({len(full_df)} rows)")

    # ── Cross-species consensus GO terms ─────────────────────────────────────
    consensus_rows = []
    for go_id, species_list in go_term_hits.items():
        if len(species_list) < 3:
            continue
        rep = full_df[full_df["go_id"] == go_id].iloc[0]
        consensus_rows.append({
            "go_id":        go_id,
            "go_name":      rep["go_name"],
            "ns_short":     rep["ns_short"],
            "n_species":    len(species_list),
            "species":      ", ".join(species_list),
            "mean_or":      round(full_df[full_df["go_id"] == go_id]["odds_ratio"].mean(), 3),
            "min_fdr":      round(full_df[full_df["go_id"] == go_id]["fdr"].min(), 4),
        })

    if consensus_rows:
        cons_df = pd.DataFrame(consensus_rows).sort_values("n_species", ascending=False)
        cons_path = RESULTS_DIR / "go_consensus.tsv"
        cons_df.to_csv(cons_path, sep="\t", index=False)
        print(f"\nCross-species consensus GO terms (≥3 species):")
        print(cons_df.head(20).to_string(index=False))
        print(f"\nConsensus table: {cons_path}")
    else:
        print("\nNo GO terms enriched in ≥3 species.")
        cons_df = pd.DataFrame()

    # ── Figure 10: Dot plot of top GO terms per species ──────────────────────
    top_per_species = (
        full_df.sort_values("fdr")
        .groupby("species_key")
        .head(10)
    )

    if top_per_species.empty:
        print("No data for figure.")
        return

    pivot = (
        top_per_species
        .assign(neg_log_fdr=lambda d: -np.log10(d["fdr"].clip(lower=1e-10)))
        .pivot_table(index="go_name", columns="species_key",
                     values="neg_log_fdr", aggfunc="max")
        .fillna(0)
    )

    fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 1.5),
                                    max(6, len(pivot) * 0.35)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", vmin=0)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(
        [c.replace("_", " ").title() for c in pivot.columns],
        rotation=35, ha="right", fontsize=9
    )
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    plt.colorbar(im, ax=ax, label="-log₁₀(FDR)", shrink=0.6)
    ax.set_title("GO terms enriched in terminal-LCR proteins\n"
                 "(top 10 per species, FDR ≤ 0.05)", fontsize=12)
    plt.tight_layout()

    out = FIGURES_DIR / "go_enrichment.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nGO-enrichment figure saved: {out}")
    plt.close()


if __name__ == "__main__":
    main()
