# CLAUDE.md — LCR Terminal Enrichment Study

Project-specific context for Claude Code sessions in this repo. Read this before making
changes to scripts, manuscript, or analysis logic.

## Research question

Are low-complexity regions (LCRs) — protein segments dominated by one or a few amino
acid types — enriched at protein **termini** (N- and C-terminal ends), and is that
enrichment a *general architectural property* of proteins rather than a lineage-specific
feature? Prior work (Coletta et al. 2010) showed it in one yeast proteome; Teekas et al.
(2024) showed it across 308 Tetrapoda species. This project asks whether the pattern:

1. holds across all Metazoa, including the most basally branching phyla,
2. extends to all major eukaryotic supergroups, and
3. — the primary novel claim — extends to **prokaryotes** (Bacteria, Archaea), making
   terminal LCR enrichment a domain-independent property of protein architecture.

**Framing (current, as of manuscript_v10):** the paper is now **mechanism-first**, not
evolutionary. The three-domain breadth is presented as evidence of *generality*, and the
amino-acid identity + N/C polarity of terminal LCRs as a *mechanistic fingerprint* (shared
N-terminal-processing/co-translational layer vs. lineage-specific degron biology). The
earlier deep-time framing — LUCA/LECA dating and the Scenario A (ancestral) vs. B
(convergent) origin debate — was **dropped** (it was unresolved and exposed the design's
phylogenetic-pseudoreplication weakness). Do **not** re-introduce LUCA/origin claims. See
[[v7-mechanism-pivot]].

**v8 narrative fix (do not regress):** v7's Introduction framed the prokaryote test as a
clean either/or — "general architecture → appears in prokaryotes; eukaryotic
ubiquitin-proteasome degron biology → does not." That is a **false dichotomy**, and v7's own
Discussion contradicted it by describing bacterial C-terminal ssrA (tmRNA) tags and the
ClpS/ClpAP N-end-rule pathway. Prokaryotes have their own terminally acting targeting and
proteostasis machinery, so terminal enrichment in bacteria excludes only a *strictly
eukaryote-specific ubiquitin-proteasome* explanation, **not** a convergent degradation-driven
one. v8 replaces the dichotomy with two separate evidence lines: taxonomic **breadth**
(generality) and compositional/directional **fingerprint** (mechanism). It also reframes
Ntountoumi et al. (2019)'s ribosomal C-terminal observation as the *motivating precedent* —
is it ribosome-specific, or the visible edge of a proteome-wide pattern? — rather than
pretending the prokaryote positional question was a blank slate. Do not restore the
either/or, and do not state prokaryote novelty without naming the Ntountoumi precedent.

**Ntountoumi 2019 — verified against primary source (2026-07-25):** full text checked at
[PMC6821194](https://pmc.ncbi.nlm.nih.gov/articles/PMC6821194/) / [NAR gkz730](https://doi.org/10.1093/nar/gkz730).
Confirmed: the paper's only original positional claim is in the Results subsection *"LCRs are
frequent in ribosomal proteins"* — "Intriguingly, we observed that the vast majority of the
bacterial and archaeal ribosomal LCRs were located at the C-terminus of their proteins." It
uses **no positional binning scheme**, reports **no terminal-vs-internal fraction for the
proteome**, and applies **no positional null model**. Its N-/C-terminal vs central *functional*
distinction is cited from prior *S. cerevisiae* work (Coletta 2010), not original. Its own
Results subheadings are: prevalence; amino acid content; LCR-containing homologs; ribosomal
proteins; amino acid enrichment by functional category; neural-network web server. Our
"first systematic, proteome-wide positional analysis in prokaryotes" claim is therefore
**sound as worded** — do not re-verify.

**Persi 2023 & Saravanan 2025 — also verified (2026-07-25):** neither contains any positional
analysis. Persi et al. (PNAS 120:e2300154120, [PMC10120016](https://pmc.ncbi.nlm.nih.gov/articles/PMC10120016/))
has no occurrence of N-/C-terminus or "terminal" at all; its Results are *LCRs and gene paralogy
in groups of prokaryotic species* / *in individual genomes* / *Excluding length of genes as a
confounding factor* (that last is a useful precedent for our own length control, currently
uncited). Saravanan et al. ([Front Bioinform 5:1673480](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2025.1673480/full))
analyses three enterobacteria at species/population level with no positional component —
notably it **cites Teekas et al. 2024 in its bibliography** yet still does not test position in
prokaryotes, which strengthens rather than threatens our novelty claim. All three prokaryote
LCR references are therefore confirmed non-positional, so the Discussion sentence asserting
they "did not apply positional binning across the proteome" is verified for all three.

**v9 phantom-citation fix (integrity — do not regress):** v3–v8 cited a non-existent
"companion study" as `(Chan, preceding study)` / `(Chan, in preparation)`, crediting it with
"61 metazoan invertebrate species across 16 phyla". **There is no such study.** Those numbers
are `manuscript_v2.md` — this project's own earlier draft (61 metazoan species, 16 phyla,
4 protist lineages, median 21 LCRs/bacterial proteome; all match verbatim). Neither citation
had a reference-list entry. v9 removes both and reframes the 10 dependent passages so the
paper stands as ONE study of 772 proteomes: gaps are now stated against the *published*
literature (Coletta = 1 yeast proteome; Teekas = Tetrapoda only), and the Methods
"as described previously" pointer — which had no valid referent, a reproducibility hole —
now cites the Teekas et al. (2024) framework explicitly. **Never reintroduce a companion /
preceding / in-preparation self-citation.** The invertebrate result belongs to this paper.

**v9 reference list (verified 2026-07-25):** 28 entries, every one carrying a DOI verified
against the Crossref API. Author lists expanded to full where ≤10 authors; "et al." retained
only for Boija (20), van der Lee (18), Virtanen (35+). One factual error corrected:
Holdsworth et al. 2020 listed **Estavillo GM** as 5th author — the real 5th author is
**Zubrycka A** (doi:10.1111/jipb.12882). Note this list has a history of author fabrication
(SUBMISSION_GUIDE records an earlier fix to Piatkov's "fabricated middle authors"), so treat
author lists as suspect until checked against DOI. Beware: naive Crossref bibliographic
search returns shadow records — "Faculty Opinions recommendation of…", peer-review reports,
decision letters — which yielded 11 wrong DOIs on the first pass; always title-match before
accepting a DOI.

**v10 fern/gymnosperm expansion (2026-07-25):** the 16 species downloaded by `01c` had been
fLPS'd but never analysed — the analysis ran 1 Jul 06:40, the expansion landed 09:08 the same
morning, and `main.py --from-phase 2` was never re-run. Re-running it added them: **756 → 772
proteomes, 962,447 → 1,004,572 LCRs, Viridiplantae 118 → 132, Chlorophyta 2 → 6.** Prokaryote
results are byte-identical (Bacteria 27.56% OR 3.431; Archaea 25.73% OR 3.132; bootstrap and
GEE unchanged) — the primary claim was untouched.

**Two pipeline traps discovered — do not repeat:**
1. `main.py` skips any step whose output exists and still reports "All 19 steps completed
   successfully". Without `--force` a re-run is a silent no-op.
2. `main.py` runs `build_viridiplantae_backbone_tree.py` **before** script 19, but the tree
   builder reads script 19's output. On a first re-run the tree is therefore built from the
   *stale* table, and `load_tree_vcv` silently drops any genus missing from the tree — so the
   PGLS excluded ferns/gymnosperms while `pgls_tier_summary.tsv` displayed them. **Always run
   tree→PGLS twice**: `main.py --from-phase 2 --force`, then
   `build_viridiplantae_backbone_tree.py`, then `19_pgls_viridiplantae.py`.
   The ladder in the tree builder now includes tiers 4/5 (seed-plant crown 350 Ma = midpoint of
   Morris et al. 2018's 330-370; euphyllophyte crown 400 Ma) and uses `join_opt` so an empty
   tier collapses instead of emitting malformed Newick.

**v10 PGLS (complete 10-tier ladder, 100/100 genera matched):** slope 0.253, p=0.185,
R²=0.018. The negative result **holds and is stronger** — it no longer rests on a predictor
with two unsampled tiers. Fern/gymnosperm sampling is narrow (both ferns Pteridaceae, all
three gymnosperms conifers) because Ginkgo/cycads/Gnetales lack protein-level annotation.

**Pipeline IS runnable on the Mac** for phases 2+ (contrary to the old note below): use
`env PATH=/Users/kchan/lcr_study/.venv/bin:$PATH python3 main.py ...` — `main.py` invokes
children as bare `python3`, and system Python 3.9 has no scipy. Only downloads/fLPS need
the network. Backup before any re-run: `results_backup_pre_expansion_20260725.tgz`.

**Length quartiles are GLOBAL, not per-domain** (`07_confound_test.py:54-56`): breaks come
from the pooled host-protein lengths of the whole dataset. v9 and earlier mis-described this
as "within each domain or phylum"; v10 corrects it. This is why bacterial quartile percentages
shifted when plant species were added even though bacterial LCRs did not change.

**Mechanism citations added in v8:** the previously uncited bacterial ssrA and ClpS/ClpAP
claims now cite Keiler et al. 1996 (*Science* 271:990–993) and Erbse et al. 2006 (*Nature*
439:753–756) respectively, in both the Introduction and Discussion.

Current dataset: **772 proteomes, 43 phyla/groups, all 3 domains of life, 1,004,572 LCRs**.
Target venue: molecular / computational-biology journal (e.g. *NAR*, *PLoS Comp Biol*,
*Protein Science*) — retargeted away from *GBE*/*MBE* with the mechanism-first pivot.

## Method (mirrors Teekas et al. 2024)

- **LCR detection:** fLPS 2.0 (`-m 3`, minimum LCR length 3 aa), precompiled binary in
  `flps2/`. Only **SINGLE-type** rows retained (single dominant residue), purity ≥70%.
- **Terminal definition:** protein divided into 20 equal-length positional bins; bins 1
  and 20 (outer 5% at each end) = "terminal." Null expectation = 10% (2/20 bins).
- **Primary statistic:** one-sided Fisher's exact test vs. the 10% null, at species and
  pooled-phylum/domain level.
- **Length-confound control:** proteins stratified into global length quartiles within
  each domain/phylum; enrichment tested per quartile. Enrichment *increasing* Q1→Q3 (as
  observed in Bacteria/Archaea) is evidence against a length artefact.
- **Asymmetry:** `pct_nterm / pct_cterm` per species/phylum — captures N- vs C-terminal
  directional bias (unstable near-zero; treat ratios as indicative when C-term LCR
  count < 50 for a group).
- **Purity gradient:** Mann-Whitney U, terminal vs. internal LCR purity — tests whether
  enrichment is positional-only or also qualitative.
- **Multiple testing:** Holm-Bonferroni across all 43 phylum-level tests.
- **Driver analysis:** singleton-LCR proteins vs. multi-LCR proteins, tested separately.

## Pipeline (orchestrated by `main.py`)

```
Phase 1 — data acquisition
  01a_download_metazoa.py                 Ensembl Metazoa (~metazoan invertebrates)
  01b_download_outgroups.py               UniProt reference proteomes (non-metazoan)
  01c_download_viridiplantae_expansion.py Ensembl Plants (Viridiplantae expansion)
  02_run_flps.sh                          Run fLPS2 on all *.longest.fa (skips existing)

Phase 2 — core analysis
  03_analyse.py           Parse fLPS output → per-LCR table + enrichment stats
  04_visualise.py          Core figures: heatmap, bar chart, U-profile

Phase 3 — extended analyses
  05_asymmetry.py           N vs C terminal split
  06_aa_composition.py      Amino acid identity, terminal vs internal
  07_confound_test.py       Protein-length quartile stratification
  08_purity_gradient.py     Terminal vs internal purity (Mann-Whitney)
  09_phylum_stats.py        Between-phylum KW test + within-phylum CV
  10_driver_analysis.py     Singleton- vs multi-LCR protein enrichment
  11_protein_level_test.py  Protein-level binomial test (sensitivity)

Phase 4 — GO enrichment
  12_go_enrichment.py       GO term enrichment of terminal-LCR proteins (BioMart)

Phase 5 — robustness & mechanism analyses (post v4, addressing peer-review critique)
  14_multiple_testing.py             Holm-Bonferroni across 43 phylum tests
  15_sensitivity_analysis.py         fLPS parameter sensitivity (length/purity combos)
  16_domain_aa_composition.py        Per-domain terminal AA composition (Bact/Arch/Plant/Metazoa)
  fetch_uniprot_annotations.py       PREREQ for 17: fetch UniProtKB SIGNAL annotations
  17_signal_peptide_stratification.py  Bacterial N-term LCRs vs signal-peptide presence
  18_mixed_effects_model.py          Bootstrap CI + cluster-robust GEE for prokaryote pooling
  build_viridiplantae_backbone_tree.py  PREREQ for 19: conservative dated backbone tree
  19_pgls_viridiplantae.py           PGLS test of the N/C asymmetry gradient in plants
  20_llps_analysis.py                PLAAC/LLPS propensity, terminal vs internal LCRs

Archived (mechanism-first pivot, 2026-07-11 — see archive/README.md):
  13_timetree_phylogeny.py           Deep-time dated phylogeny (LUCA/LECA framing, dropped)
Optional/standalone (not in main.py):
  00_setup.sh                        fLPS2 binary + Python deps; run once
  build_timetree_from_api.py         timetree.org species-level tree (future within-tier PGLS)
```

Run via `python3 main.py` (see `README.md` for `--phase`, `--from-phase`, `--dry-run`,
`--skip-download`, `--force` flags). All extended-analysis scripts read from
`results/` tables produced by `03_analyse.py` — nothing upstream needs rerunning
unless the species set or fLPS parameters change.

## Key results (manuscript v9 — numbers from current `results/*.tsv`)

- **Prokaryotes are terminal-enriched** (primary novel finding): Bacteria **27.6%** (92
  spp, pooled OR 3.43), Archaea **25.7%** (27 spp, OR 3.13) — the highest pooled terminal
  fractions of any group. **Novelty wording — do NOT overstate:** Ntountoumi et al. (2019)
  already noted C-terminal LCRs in the *ribosomal-protein* subset of bacteria/archaea, so the
  defensible claim is the first *systematic, proteome-wide* positional analysis in prokaryotes,
  not the "first positional analysis" outright. Always acknowledge the Ntountoumi ribosomal
  precedent when stating novelty. Length-quartile control rules out a length artefact (enrichment
  rises Q1→Q3: Bacteria 21.2%→36.7%, Archaea 20.6%→37.3%). Robust to a cluster-robust GEE
  (Q4 vs Q1 p<0.001) and species bootstrap. High within-domain CV (~50% / ~45%) means the
  pooled value is a domain aggregate driven by LCR-rich taxa (mainly Actinobacteria), not a
  species-typical value — state this honestly.
- **22/24 metazoan phyla** survive Holm-Bonferroni (Acanthocephala and Nematomorpha are
  underpowered single-species groups). Porifera is now **4 spp / 17.7%** — mid-range, *not*
  the exceptional "pre-neural anchor" the old draft claimed.
- **All major eukaryotic supergroups** are pooled-significant; **41/43** phyla/groups
  overall survive Holm-Bonferroni. Only **two** fail: Acanthocephala and Nematomorpha, both
  single-species. Chlorophyta (now 6 spp, p_holm<0.001) and Rhodophyta (p_holm=0.044) cleared
  the threshold once algal sampling was expanded — do not describe them as provisional.
  Heterolobosea (n=1) is the highest eukaryotic value (31.0%).
- **N/C asymmetry is lineage-specific, not a clean N-vs-C dichotomy:** ONLY **Viridiplantae**
  is robustly N-terminal dominant (median ratio 2.76). **Bacteria are BALANCED at the
  proteome level (median 1.00)** — the earlier "bacteria are N-dominant" claim was
  cherry-picked Actinobacteria and is **wrong**; do not repeat it. Archaea are C-dominant
  (0.71); Metazoa/Fungi balanced-to-C; Ciliophora/Metamonada/Platyhelminthes strongly C.
  Signal-peptide-bearing bacterial proteins carry a 3× N-terminal LCR subpopulation (34.1%
  vs 14.6%, OR 3.03, p=7.1e-5) — a secreted minority, not a whole-proteome bias.
- **Composition fingerprint (mechanism keystone):** leucine is the most terminally-biased
  residue in *every* domain (2.8–5.3×); methionine is second in all three eukaryotic groups
  — pointing to N-terminal signal/anchor sequences and initiator-Met retention (a shared
  N-terminal-processing layer), distinct from acidic/basic degron residues.
- **Singleton-LCR proteins** drive the signal in 41/42 groups (multi-LCR in 29/42). Terminal
  purity is elevated in only 9/42 groups (positional, not compositional; prokaryotes not
  purer). Plant mosses→grasses N-polarity gradient is **not robust** to phylogenetic
  correction (PGLS slope 0.253, p=0.185). LLPS propensity is organism-specific, not a
  general terminal signature.
- **Mechanistic model (layered, mechanism-first):** a shared N-terminal-processing /
  co-translational baseline across all cellular life, elaborated by lineage-specific degron
  biology (C-degron in animals, N-degron/PRT6-PCO in land plants). Framed as *mechanism +
  generality*, NOT deep-time origin.

## Known gaps / in-progress work (do not overstate these as resolved)

- Choanoflagellates (*Monosiga*, *Salpingoeca*) and Asgard archaea are still absent. Under
  the mechanism-first framing these are **optional generality extensions**, no longer
  critical — they only mattered for the dropped LUCA-vs-convergence origin question.
- Scripts 15–20 have been run (outputs present in `results/`). Script **13 (deep-time
  phylogeny) is archived** in `archive/` — the mechanism-first pivot dropped it; do not
  re-wire it into `main.py`.
- Single-species groups (report as data points, not phylum confirmations): Perkinsozoa,
  Haptophyta, Brachiopoda, Acanthocephala, Nematomorpha, Priapulida, Chordata, Heterolobosea.
- `manuscript_v6.md` and earlier, `README.md`'s narrative sections, and `findings_summary.md`
  reflect the **older evolutionary framing** — historical, not current. `manuscript_v10.md`
  is the source of truth; v7–v9 are historical. **v9 and earlier describe the pre-expansion
  756-proteome dataset and must not be submitted** — only v10 matches current `results/`.

## Working conventions for this repo

- Manuscript versions live in `manuscript/` (`manuscript_v10.md` is current, mechanism-first);
  don't overwrite prior versions — create a new `_v{n+1}.md` for substantive revisions.
- `findings_summary.md` is a **historical** working synthesis + self-review doc that drove
  the v4→v5 revision under the old evolutionary framing; its "Consolidated Priority List" is
  largely addressed and it is not the current task list.
- Numeric claims in prose (percentages, p-values, species counts) must match the
  corresponding `results/*.tsv` — if editing analysis scripts, regenerate affected
  tables/figures before updating manuscript text.
- This environment cannot run the pipeline (no Python env with fLPS binary guaranteed
  active, no live data fetch) — reason from existing `results/`, `figures/`, and script
  logic; hand off exact commands for the user to run rather than attempting execution.
