# results/ table manifest

Inventory of `results/` and what each file backs in `manuscript/manuscript_v10.md`.

**Update this file whenever a filename changes or a script starts/stops backing a manuscript
table.** If a row's "Manuscript reference" is wrong, the file's *name* is the first thing to
fix — see the naming rule below.

## Naming rule

Any file that backs a numbered item in the manuscript is named after that item, so the name
alone tells you what it is:
- `supp_table_SN_description.tsv` → Supplementary Table N (`manuscript_v10.md`'s own numbering)
- `tableN_description.tsv` → main-text Table N
- Figures follow the same rule already (`figN_*.pdf`, `suppfigN_*.pdf`) — see `CLAUDE.md`

Everything else (intermediates, caches, exploratory outputs not cited in the manuscript) keeps
a plain descriptive name with **no number**, so it can never be mistaken for a citable table.
Each table has exactly **one** owning script — never add a second script that recomputes the
same statistic into a different file, and never add a copier/export script that regenerates a
manuscript-facing table from another script's output. The fix for a wrong or missing number
always belongs in the script that owns that analysis.

## Manuscript-cited tables

| File | Manuscript reference | Purpose | Generating script |
|---|---|---|---|
| `supp_table_S1_species_list.tsv` | Supp Table S1 | Full 772-species list: name, phylum, source, proteome ID, protein/LCR counts, pct_terminal | `01a_download_metazoa.py` + `01b_download_outgroups.py` + `01c_download_viridiplantae_expansion.py` (cumulative) |
| `supp_table_S2_driver_analysis.tsv` | Supp Table S2 | Singleton- vs multi-LCR terminal enrichment per phylum | `10_driver_analysis.py` |
| `supp_table_S3_within_phylum_cv.tsv` | Supp Table S3 | Within-phylum coefficient of variation of pct_terminal | `09_phylum_stats.py` |
| `supp_table_S4_length_stratified.tsv` | Supp Table S4 | Fisher's exact per protein-length quartile × phylum | `07_confound_test.py` |
| `supp_table_S5_purity_gradient.tsv` | Supp Table S5 | Terminal vs internal LCR purity, Mann-Whitney U per phylum | `08_purity_gradient.py` |
| `supp_table_S6_flps_sensitivity.tsv` | Supp Table S6 | fLPS parameter sensitivity (3 param sets) per phylum | `15_sensitivity_analysis.py` |
| `supp_table_S6_multitype_positions.tsv` | Supp Table S6 (MULTI-type component) | Raw MULTI-type LCR positional records | `15_sensitivity_analysis.py` |
| `supp_table_S7_domain_aa_composition.tsv` | Supp Table S7 | Per-domain (Bacteria/Archaea/Viridiplantae/Metazoa/Other_Euk) terminal AA enrichment ratios | `16_domain_aa_composition.py` |
| `supp_table_S8_signal_peptide_stratification.tsv` | Supp Table S8 | Bacterial N-terminal LCRs by UniProt SIGNAL annotation | `17_signal_peptide_stratification.py` |
| `supp_table_S9_llps_organism_summary.tsv` | Supp Table S9 | LLPS proxy score, terminal vs internal, per model organism | `20_llps_analysis.py` |
| `table1_2_phylum_lineage_enrichment.tsv` | Table 1 (Metazoa rows) + Table 2 (non-metazoan rows) | Pooled % terminal + Holm-Bonferroni corrected significance per phylum/lineage — **the single canonical Holm-Bonferroni correction**, do not duplicate elsewhere | `14_multiple_testing.py` (reads `phylum_summary.tsv` from `03_analyse.py`) |
| `multiple_testing_report.txt` | Backs the "41/43 significant, only Acanthocephala + Nematomorpha provisional" claim throughout the manuscript | Human-readable Holm-Bonferroni report | `14_multiple_testing.py` |

**Table 1 and Table 2 are hand-copied into the manuscript markdown as inline tables, not
generated files.** After any pipeline re-run, diff `table1_2_phylum_lineage_enrichment.tsv`
against the manuscript's Table 1/2 rows before treating a "no regression" claim as verified —
nothing enforces this sync automatically, so it is the single easiest place for the manuscript
to silently drift from the data.

## Cited only in prose (no numbered table/figure of their own)

| File | Manuscript reference | Purpose | Generating script |
|---|---|---|---|
| `protein_level_enrichment.tsv` | "Protein-level confirmation" paragraph (main text) | Binomial test: fraction of proteins with ≥1 terminal LCR per species | `11_protein_level_test.py` |
| `go_enrichment.tsv` | "Driver analysis and GO enrichment" paragraph (main text) | GO term enrichment in terminal-LCR proteins, 5 BioMart-compatible species | `12_go_enrichment.py` |
| `bootstrap_ci_domain.tsv`, `mixed_model_prokaryote.tsv` | Prokaryote clustering-robustness check (main text) | GEE / bootstrap CI on domain-level % terminal, species as cluster | `18_mixed_effects_model.py` |
| `pgls_regression_stats.tsv`, `pgls_tier_summary.tsv`, `pgls_viridiplantae.tsv` | Supp Figure 7 + PGLS prose (r=0.63, PGLS slope=0.253) | PGLS regression of N/C asymmetry vs. evolutionary tier, Viridiplantae | `19_pgls_viridiplantae.py` |
| `aa_composition.tsv` | Supp Figure 5 | Terminal vs internal AA composition, pooled across Metazoa | `06_aa_composition.py` |

## Intermediate / not cited in the manuscript

Legitimate current outputs — not orphaned, just not manuscript-facing. Keep plain names.

| File | Purpose | Generating script |
|---|---|---|
| `lcr_positions.tsv` | Per-LCR raw table (bin, purity, location) — primary intermediate, read by 9+ downstream scripts | `03_analyse.py` |
| `enrichment.tsv` | Per-species Fisher's exact results | `03_analyse.py` |
| `phylum_summary.tsv` | Per-phylum pooled Fisher's test, pre-Holm-correction — input to `14_multiple_testing.py` only | `03_analyse.py` |
| `phylum_stats.tsv` | Pairwise Mann-Whitney U between phyla (595 pairs, Holm-corrected) — not cited in `manuscript_v10.md` | `09_phylum_stats.py` |
| `asymmetry.tsv` | Per-species N/C ratio + log-odds asymmetry | `05_asymmetry.py` |
| `llps_analysis.tsv` | Raw per-LCR LLPS proxy scores (S9's summary file is the cited one, this is its intermediate) | `20_llps_analysis.py` |
| `sampling_robustness.tsv` | Taxon-sampling robustness (Insecta-excluded, equal-weighting checks) — not currently a numbered supp table | `15_sensitivity_analysis.py` |

## Inputs / caches (not results, do not delete)

| File | Purpose | Generating script |
|---|---|---|
| `download_status_metazoa.tsv`, `download_status_outgroups.tsv` | Per-species download success/failure log | `01a_download_metazoa.py`, `01b_download_outgroups.py` |
| `taxonomy_cache.tsv` | Shared NCBI taxon-ID cache | `01a_download_metazoa.py` + `01b_download_outgroups.py` (both append) |
| `timetree_pairwise_cache.json` | TimeTree API pairwise divergence-date cache | `build_timetree_from_api.py` |
| `uniprot_signal_peptides.tsv` | UniProtKB SIGNAL annotation cache, prereq for Supp Table S8 | `fetch_uniprot_annotations.py` (read by `17_signal_peptide_stratification.py`) |
| `viridiplantae_timetree.nwk` | Time-calibrated backbone phylogeny, prereq for PGLS | `build_viridiplantae_backbone_tree.py` (read by `19_pgls_viridiplantae.py`) |
| `pgls_species_list.txt` | Manually uploaded to timetree.org's "Load a List of Species" tool — not written or read by any script | manual, one-time |
| `flps/` | Raw fLPS 2.0 output, one file per species | `02_run_flps.sh` |

## Known gotchas (do not regress)

- **`PHYLUM_ORDER`/`PHYLUM_COLOURS`/`SUPERGROUP_OF` (`scripts/config.py`) must list every
  phylum/group in the dataset.** Most scripts use `PHYLUM_ORDER` as a filter, not just a sort
  order, so a missing entry is silently dropped from that script's figures/tables — see
  `CLAUDE.md` for the full list of which outputs are and aren't affected by this. A phylum can
  be completely correct in Table 1/2 while silently absent from three supplementary tables and
  three main-text figures, because `14_multiple_testing.py`'s use of `PHYLUM_ORDER` is
  sort-only (with a numeric fallback) rather than a filter.
- Every renamed/regenerated output in this manifest has been cross-checked against the
  manuscript's own cited numbers, not just renamed — treat a filename match as necessary but
  not sufficient; always diff the actual values too.
