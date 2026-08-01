# Terminal Low-Complexity Regions — a cross-domain architectural study

Tests whether the terminal enrichment of low-complexity regions (LCRs) first described in
yeast (Coletta et al. 2010) and Tetrapoda (Teekas et al. 2024) is a **general
architectural property of proteins** across all three domains of life — and dissects its
mechanism through amino-acid composition and N/C polarity.

The manuscript, working notes, and archived analyses are maintained separately and are not
part of this repository — this repo tracks only what reproduces the results below: the
pipeline scripts, the result tables, and the figures.

## Research question

Are LCRs — protein segments dominated by one or a few amino acid types — enriched at
protein **termini** (N- and C-terminal ends), and is that enrichment a domain-independent
feature of protein architecture rather than a lineage-specific trait? The project asks
whether the pattern (1) holds across all Metazoa incl. the most basal phyla, (2) spans all
major eukaryotic supergroups, and (3) — the primary novel claim — extends to **prokaryotes
(Bacteria, Archaea)**, and what the amino-acid identity and N/C polarity of terminal LCRs
reveal about mechanism.

## Dataset

**772 proteomes · 43 phyla/groups · all 3 domains of life · 1,004,572 LCRs.**

| Domain / grouping | Species | Groups |
|---|---|---|
| Bacteria | 92 | 1 |
| Archaea | 27 | 1 |
| Metazoa | 367 | 24 phyla |
| Other eukaryotes (all supergroups: SAR, Archaeplastida, Excavata, Amoebozoa, Opisthokonta-Fungi, Haptophyta, Cryptophyta) | 286 | 17 |
| **Total** | **772** | **43** |

Full per-species list: `results/supp_table_S1_species_list.tsv` (manuscript Supp Table S1).

## Pipeline (orchestrated by `main.py`)

One-time setup (not orchestrated): `bash scripts/00_setup.sh` — fLPS2 binary + Python deps.

```
Phase 1 — Data acquisition
  01a_download_metazoa.py                 Ensembl Metazoa (metazoan invertebrates)
  01b_download_outgroups.py               UniProt reference proteomes (Bacteria/Archaea/Fungi/protists)
  01c_download_viridiplantae_expansion.py Ensembl Plants (Viridiplantae expansion)
  02_run_flps.sh                          Run fLPS 2.0 on all *.longest.fa (skips existing)

Phase 2 — Core analysis
  03_analyse.py                           Parse fLPS output → per-LCR table, enrichment stats
  04_visualise.py                         Core figures: bin heatmap, terminal bar chart, U-profile

Phase 3 — Extended analyses
  05_asymmetry.py                         N- vs C-terminal split (log-odds asymmetry)
  06_aa_composition.py                    Amino acid identity, terminal vs internal
  07_confound_test.py                     Protein-length quartile stratification
  08_purity_gradient.py                   Terminal vs internal purity (Mann-Whitney U)
  09_phylum_stats.py                      Between-phylum KW test + within-phylum CV
  10_driver_analysis.py                   Singleton- vs multi-LCR protein enrichment
  11_protein_level_test.py                Protein-level binomial test (sensitivity)

Phase 4 — GO enrichment
  12_go_enrichment.py                     GO term enrichment of terminal-LCR proteins (BioMart)

Phase 5 — Robustness & mechanism analyses
  14_multiple_testing.py                  Holm-Bonferroni across 43 phylum tests
  15_sensitivity_analysis.py              fLPS parameter sensitivity (4 param sets) + MULTI-type
  16_domain_aa_composition.py             Per-domain terminal AA composition fingerprint
  fetch_uniprot_annotations.py            PREREQ for 17 — fetch UniProtKB SIGNAL annotations
  17_signal_peptide_stratification.py     Bacterial N-term LCRs vs signal-peptide presence
  18_mixed_effects_model.py               Bootstrap CI + cluster-robust GEE (prokaryote pooling)
  build_viridiplantae_backbone_tree.py    PREREQ for 19 — conservative dated backbone tree
  19_pgls_viridiplantae.py                PGLS of the plant N/C polarity gradient
  20_llps_analysis.py                     LLPS propensity (PLAAC proxy), terminal vs internal
```

Not orchestrated: `build_timetree_from_api.py` (optional timetree.org tree, future within-tier
PGLS). A step 13 (deep-time dated phylogeny) existed in earlier pipeline versions and was
dropped with the mechanism-first pivot; it is not part of this repository.

### Run order (on Mac / Linux)

```bash
bash scripts/00_setup.sh          # one-time: fLPS binary + deps
python3 main.py                   # full pipeline (phases 1–5), skips done steps
```

Common options:
```bash
python3 main.py --list            # status of every step
python3 main.py --dry-run         # preview without executing
python3 main.py --phase 5         # only phase 5
python3 main.py --from-phase 3    # phase 3 onward
python3 main.py --skip-download   # skip 01a/01b/01c, use existing proteomes
python3 main.py --force           # re-run even if outputs exist
```

`main.py` skips any step whose sentinel output already exists; download/fLPS steps also skip
already-processed species internally. **After adding species:** `python3 main.py --skip-download`
reruns fLPS → analysis on the updated set.

## Key outputs

Every table in `results/` — what it is, which manuscript table/figure it backs (if any), and
which script generates it — is catalogued in **[`results/TABLE_MANIFEST.md`](results/TABLE_MANIFEST.md)**.
Keep that file, not this README, up to date when outputs change.

Figures follow the naming convention documented in `CLAUDE.md` (`figN_*.pdf` / `suppfigN_*.pdf`,
N = manuscript figure number): `figures/fig1–fig10*`, `figures/suppfig_*`.

## Method (mirrors Teekas et al. 2024)

- **LCR detection:** fLPS 2.0 (`-m 3`); SINGLE-type rows only; purity ≥70% (dominant residue).
- **Terminal definition:** 20 equal-length positional bins; bins 1 & 20 = terminal; null = 10%.
- **Primary statistic:** one-sided Fisher's exact vs the 10% null (species + pooled group level).
- **Length-confound control:** length quartile breaks computed globally across the pooled
  dataset (not per domain/group); enrichment then tested within each quartile, per group
  (rising Q1→Q3 = evidence against a length artefact).
- **Robustness:** parameter sensitivity (4 sets), cluster-robust GEE + species bootstrap
  (prokaryotes), Holm-Bonferroni across all 43 group tests (41/43 survive correction; the two
  that don't — Acanthocephala, Nematomorpha — are both single-species groups).
- **Mechanism:** per-domain terminal AA composition; N/C polarity (log-odds asymmetry); bacterial
  signal-peptide stratification; PGLS of the plant polarity gradient; LLPS proxy.

## What is novel here

- **First systematic, proteome-wide positional (terminal) LCR analysis in prokaryotes** — Bacteria
  27.6% and Archaea 25.7% terminal, the highest pooled fractions among the multi-species groups
  (single-species groups such as Heterolobosea, 31.0% at n=1, are reported as data points, not
  robust comparisons), with a decisive length-confound control. (Prior work noted C-terminal LCRs
  only in the ribosomal-protein subset; Ntountoumi et al. 2019.)
- **A mechanistic fingerprint:** leucine is the most terminally-biased residue in every domain and
  methionine second in eukaryotes (N-terminal signal/anchor + initiator-Met processing); N/C
  polarity is lineage-specific (only land plants robustly N-dominant; bacteria balanced), which is
  itself incompatible with a single universal directional mechanism.

## Citation to anchor

Teekas L, Sharma S, Vijay N (2024). Terminal regions of a protein are a hotspot for low complexity
regions and selection. *Open Biology* 14(6):230439. https://doi.org/10.1098/rsob.230439
