# Taxon-sampling expansion & de-biasing — implementation plan

Status: **Workstreams A + B implemented in `scripts/01b_download_outgroups.py`**
(no new scripts). Workstream D (analysis guardrails) still pending.
Scope: `scripts/01b`, `results/species_manifest.tsv`, downstream sensitivity analysis
Goal: reduce taxon-sampling bias in the LCR study without over-claiming where data
genuinely does not exist.

Every accession below was checked live against the NCBI Datasets API before being wired
in: it exists, is the current version, and carries `protein_coding` annotation (a protein
FASTA is downloadable). The download + longest-isoform-per-gene collapse was validated to
reproduce NCBI's own gene counts exactly (e.g. Bolinopsis 21763 proteins → 14118 genes;
Salpingoeca 11725 → 11618). Targets with only a genome and no NCBI protein annotation
(Hormiphora, Pleurobrachia, Xenoturbella bocki, Novocrania, …) were dropped and replaced
with annotated congeners.

---

## 1. Diagnosis (from `results/species_manifest.tsv`, 723 species)

The sample reflects **what has been sequenced**, not an even span of the tree.

| Domain                    | species |
|---------------------------|--------:|
| Metazoa                   | 353 |
| Non-metazoan Eukaryota    | 260 |
| Bacteria                  | 87  |
| Archaea                   | 21  |

**Metazoa is dominated by one class:** Insecta = 217 (**61 % of all animals**).
The rest is a thin tail; 12 phyla rest on a single species, 3 on two.

**The outgroup is dominated by availability bias:** Viridiplantae (118) and
Apicomplexa (43, mostly medically-sequenced parasites).

### 1a. Three distinct problems (each needs a different fix)

- **Recoverable download failures** — species we *intended* to have but lost to dead
  UniProt IDs. 21 outgroup rows failed `UNIPROT_404` (see
  `results/download_status_outgroups.tsv`), **including the animals' closest non-animal
  relatives** `salpingoeca_rosetta`, `capsaspora_owczarzaki`, `sphaeroforma_arctica`.
  Only `monosiga_brevicollis` survives, so the Holozoa outgroup is effectively n=1.
  → **Fix by recovery, not new curation.**
- **Fixable thin clades** — phyla with 1–2 species where more public genomes exist and
  the node matters phylogenetically. → **Fix by adding species from NCBI.**
- **Hard-limited clades** — phyla where 1–2 genomes is all that exists on Earth
  (Priapulida, Nematomorpha, Acanthocephala). → **Cannot fix by sampling; handle in
  analysis + reporting.**

### 1b. Singleton / doubleton triage (Metazoa)

| Phylum (current n)        | Verdict            | Action |
|---------------------------|--------------------|--------|
| Ctenophora (1)            | Fix — top priority | base-of-tree node; add 2–3 |
| Porifera (2)              | Fix — top priority | rival base-of-tree node; add 2–3 |
| Xenacoelomorpha (1)       | Fix — top priority | bilaterian root; add 1–2 |
| Rotifera (1)              | Fix — easy         | add 1–2 |
| Myriapoda (1)             | Fix — easy         | non-insect arthropod; add 1–2 |
| Nemertea (1)              | Fix                | add 1 |
| Hemichordata (1)          | Fix                | add 1 (~2 is ceiling) |
| Brachiopoda (1)           | Fix                | add 1 |
| Placozoa (1)              | Fix (modest)       | add 1 |
| Chordata (1)              | **Choice, not error** | 1 amphioxus is a deliberate scoping call for an invertebrate study |
| Priapulida (1)            | Hard-limited       | analysis-side only |
| Nematomorpha (1)          | Hard-limited       | analysis-side only |
| Acanthocephala (1)        | Hard-limited       | partly covered by Rotifera (nested inside Syndermata) |

Outgroup thin spots to widen: Chlorophyta (2), plus true eukaryote-wide breadth is
narrow behind the plant/parasite bulk (best repaired via EukProt-curated references).

---

## 2. Design principle: one NCBI code path, added to the existing loader

The root cause of both the 404s and the "skipped, download manually" rows in `01c` is
that the loaders **only** understand UniProt reference proteomes, whose IDs rot. The fix
fetches a proteome by **NCBI assembly accession** (stable) and reduces it to the
pipeline's longest-isoform-per-gene contract.

**As built (no new files):** added directly to `scripts/01b_download_outgroups.py` —
- `download_ncbi(accession, out_fa) -> int | None` — pulls PROT_FASTA + GENOME_GFF from
  the NCBI Datasets REST API as a zip (in memory), maps protein→gene from the GFF3 CDS
  attributes, keeps the longest protein per gene, and writes the same 60-col FASTA via
  the existing `_write_fasta` helper. Returns gene count or `None` (callers log and
  continue via the existing `run_pool`).
- `NCBI_PROTEOMES` — a 33-row `(species_key, display, accession, taxon_id)` table:
  21 recoveries + 12 metazoan additions. Domain/phylum resolve from `taxon_id` through
  the existing `resolve_taxonomy`, so metazoan rows automatically get `domain=Metazoa`
  and the manifest's phylum labels.
- `main()` runs a third parallel pass (`process_ncbi`) after UniProt/Ensembl. NCBI-owned
  keys are skipped in the UniProt pass (no more dead 404 requests) and their binomials
  suppress Ensembl duplicates.

`01c`'s manual-download weakness is left as-is for now; the same `download_ncbi` can be
reused there later if the plant gaps need closing automatically.

**Contract preserved everywhere** (unchanged): every species yields
`data/proteomes/{species_key}.longest.fa` (one seq/gene) and one manifest row
`species_key  display_name  phylum  fa_path  domain  taxon_id`. Manifests are upserted by
`species_key`, so re-runs are idempotent and never clobber other scripts' rows.

---

## 3. Workstreams

### Workstream A — recover the 21 `UNIPROT_404` outgroups  *(done)*

Added to `NCBI_PROTEOMES` in `01b`; the UniProt pass now skips these keys. Verification
**also caught three wrong taxon_ids that were latent bugs in `01b`'s `UNIPROT_OUTGROUPS`**
(the mismatched proteome-id/taxon pairs are likely why some 404'd in the first place):

| species_key            | old taxid (wrong)              | corrected taxid | assembly |
|------------------------|--------------------------------|-----------------|----------|
| naegleria_gruberi      | 214684 (= *Cryptococcus* JEC21) | 5762           | GCF_000004985.1 |
| nitrosopumilus_maritimus | 335283 (= *Nitrosomonas eutropha*) | 436308     | GCF_000018465.1 |
| giardia_intestinalis   | 184922 (strain, no annot.)     | 5741           | GCF_000002435.2 |

Verified recoveries (accession → protein count): capsaspora GCF_000151315.2 (8621),
salpingoeca GCF_000188695.1 (11618), sphaeroforma GCF_001186125.1 (18213), entamoeba
GCF_000208925.1, leishmania GCF_000002725.2, naegleria GCF_000004985.1, oryza
GCF_034140825.1, populus GCF_000002775.5, rhizopus GCF_000149305.1, plus 12
bacteria/archaea references. Full list is the recovery block of `NCBI_PROTEOMES`.

**Why first:** restores the single most important reference set (Holozoa = the yardstick
for "is this gene family an animal invention?") plus deep prokaryote references the
study's ancient-vs-eukaryotic boundary test depends on. Pure recovery — no new sampling
decisions.

### Workstream B — metazoan expansion, added to `NCBI_PROTEOMES` in `01b`  *(done)*

No new script. 12 verified, protein-annotated species across 7 fixable phyla (domain
resolves to `Metazoa`; phylum to the existing manifest labels). Original species guesses
that turned out to have **no NCBI protein annotation** were replaced with annotated
congeners (noted below):

| Phylum | added (accession, genes) | note |
|--------|--------------------------|------|
| Porifera (2→**4**)        | Ephydatia muelleri GCA_049114765.1 (30360); Sycon ciliatum GCF_964019385.1 (20530) | base-of-tree node |
| Ctenophora (1→**2**)      | Bolinopsis microptera GCF_026151205.1 (14118) | Hormiphora/Pleurobrachia/Beroe have no NCBI annotation → used Bolinopsis |
| Xenacoelomorpha (1→**3**) | Convolutriloba macropyga GCF_964194025.1 (18007); Symsagittifera roscoffensis GCF_963678635.1 (13514) | *Xenoturbella bocki* unannotated on NCBI |
| Rotifera (1→**3**)        | Brachionus plicatilis GCA_003710015.1 (52286); Rotaria sordida GCA_905250125.1 (59060) | |
| Myriapoda (1→**3**)       | Scutigera coleoptrata GCA_982266805.1 (18157, centipede); Chamberlinius hualienensis GCA_054772095.1 (13817, millipede) | spans two myriapod classes |
| Nemertea (1→**2**)        | Tubulanus polymorphus GCF_964204645.1 (15035) | |
| Hemichordata (1→**2**)    | Ptychodera flava GCF_041260155.1 (31795) | |
| Placozoa (1→**2**)        | Trichoplax sp. H2 GCA_003344405.1 (12174) | |

**Cannot fix via NCBI:** **Brachiopoda** — only *Lingula anatina* is protein-annotated on
NCBI (already the singleton), so it stays n=1. **Chordata** stays n=1 by deliberate
scoping (invertebrate study). These move to Workstream D.

### Workstream C — outgroup breadth via EukProt  *(optional, larger)*

For eukaryote-wide rooting rather than parasite/plant-dominated sampling, add a curated
slice of EukProt proteomes (already one-per-gene) through the same fetcher/manifest path.
Scoped as a follow-up; not required for A/B.

### Workstream D — analysis guardrails for hard-limited clades

These do **not** add species; they stop the sparse tail from driving conclusions.

1. **Rarefaction / robustness.** *(done)* Added `taxon_sampling_robustness()` to
   `scripts/15_sensitivity_analysis.py`; writes `results/sampling_robustness.tsv`
   (Supp Table S7). Re-weights `lcr_positions.tsv` four ways — pooled vs equal-per-species,
   Insecta excluded, Insecta capped to the median-phylum LCR count, singletons excluded —
   plus a standalone test of each singleton phylum. scipy-free (normal-approx binomial).
   **Result: the terminal-enrichment signal is robust** (see §7).
2. **Completeness column.** *(pending — needs the BUSCO tool + the proteomes)* Record
   BUSCO %complete per proteome as a manifest side-car, so absence reads as "not detected
   in an incomplete proteome" vs real loss. Not implemented here because BUSCO can't be
   run in this environment; it is the remaining manual step.
3. **Reporting rule (manuscript).** For any n=1 phylum, report presence only; never infer
   gene-family loss/absence. `sampling_robustness.tsv` backs this: singleton phyla each
   carry <0.5% of all LCRs, so they neither bias the pooled result nor support a
   standalone phylum-level loss claim.

## 8. Viridiplantae PGLS de-biasing (implemented in `01c` + `19`)

The PGLS in `19_pgls_viridiplantae.py` was tilted: after one-per-genus dedup, **90% of
species (76/84) were eudicots/monocots/grasses (tiers 7–9)**, the early-diverging tiers
that define the gradient were n=1–2, and **tiers 4 (ferns) and 5 (gymnosperms) were empty**
because `01c` had never successfully run. Two latent defects compounded it:
- `01c` relied on UniProt proteome ids (mostly blank/unverified for these taxa) and
  *skipped* anything without one → ferns/gymnosperms were undownloadable.
- With no `viridiplantae_timetree.nwk`, `19` silently ran **OLS, not PGLS** — treating 55
  correlated eudicots as independent, so the reported `p≈0` is overconfident. (The
  tier-*mean* trend is clean, r=0.95, so the direction is real; only the resolution and
  confidence at the base were the problem.)

**Fix, implemented:**
- `01c` rewritten to fetch **accession-verified NCBI proteomes** via the same
  `download_ncbi` path as `01b` (16 species; every accession checked live for current +
  protein annotation). Tier lift (genera, after dedup): Chlorophyta 2→6, Charophytes 1→3,
  Bryophytes 2→4, Lycophytes 1→2, Ferns 0→2, Gymnosperms 0→3, Basal angiosperms 2→4 —
  **all 7 early tiers now populated (≥2, and ≥3 for five of them).**
- `19`'s `GENUS_TIER` extended so the new genera (and several original `01c` genera that
  were never listed — klebsormidium, mesotaenium, anthoceros, …) actually receive a tier
  instead of being silently dropped.

**Hard limits (NCBI-unreachable — portal-only gene models):** hornworts (0 NCBI-annotated
genomes), and *Mesotaenium / Anthoceros / Azolla / Salvinia / Isoetes / Ginkgo / Picea*
(genome-only). Annotated congeners were substituted (Closterium, Ceratopteris/Adiantum,
Taxus/Cryptomeria/Pinus longaeva). Tiers 3–4 remain at n=2 for this reason.

**Still required for a *real* PGLS:** supply `results/viridiplantae_timetree.nwk` (timetree.org);
without it `19` stays in OLS fallback. Even then, keep one-per-genus dedup — eudicots still
dominate the raw count (55), so phylogenetic weighting remains essential.

## 7. Fairness re-evaluation (result)

Measured on the current 723-species manifest (`results/lcr_positions.tsv`, 929,802 LCRs):

- **Composition is skewed but not degenerate.** Pielou evenness J′=0.68, Gini=0.77.
  Insecta is 30.2% of all LCRs (64.9% of Metazoan LCRs); plants (Archaeplastida) and
  Metazoa together supply 82% of all LCRs; prokaryotes only 0.3% (few LCRs pass the
  purity filter — biological, not a defect).
- **The headline survives re-weighting** (the key point). For every supergroup, the
  pooled (LCR-weighted) %terminal and the equal-per-species mean agree to within a few
  points (e.g. Metazoa 17.3 vs 17.8; Archaeplastida 24.3 vs 23.5) — the signal is a
  per-species property, not an artefact of LCR-rich species.
- **Not Insecta-driven:** Metazoa %terminal is 17.3 full → **18.3 with Insecta removed**
  (and 18.3 with Insecta capped). Dropping all 14 singleton phyla leaves it unchanged
  (17.3).
- **Pervasive across the thin tail:** all 14 singleton phyla are *individually*
  terminal-enriched (z 2.4–13.2), including the base-of-tree nodes (Ctenophora 20.5%,
  Xenacoelomorpha 21.8%, Placozoa 18.8%).
- **Genuinely fragile bits (unchanged by sampling):** Prokaryote %terminal rests on only
  3,076 LCRs (species sd=13); and any *phylum-level* claim for an n=1 phylum. Neither is
  used to carry the headline.

Verdict: the sample is compositionally biased, but the central terminal-enrichment
conclusion is **robust to that bias** — a quantified answer, not a hope. The expansion
(Workstreams A+B) further hardens the deep-branch coverage the conclusion leans on.

---

## 4. Order of execution

1. ~~Shared fetcher + smoke test~~ — done: `download_ncbi` in `01b`, validated against
   NCBI gene counts.
2. ~~Workstream A + B~~ — done: 33 accessions in `NCBI_PROTEOMES`.
3. **Run it:** `python scripts/01b_download_outgroups.py` (optionally `--limit 1` first
   as a live smoke test). Downloads the 33 proteomes and upserts the manifest.
4. Re-run `bash scripts/02_run_flps.sh` (picks up new `.longest.fa`), then
   `python main.py --from-phase 2` (or `--force`).
5. **Workstream D guardrails** (pending): BUSCO column + Insecta down-sampling in
   `15_sensitivity_analysis.py`; refresh `09_phylum_stats.py` / `phylum_summary*.tsv`.
6. Workstream C (EukProt) if breadth still insufficient.

## 5. Validation / done criteria

- `grep -c UNIPROT_404 results/download_status_outgroups.tsv` drops toward 0
  (or those rows show `OK_NCBI_FALLBACK`).
- Manifest contains ≥3 Holozoa (Choanoflagellata + Filasterea + Ichthyosporea).
- Every previously-singleton **fixable** phylum reaches n≥2 (target n≥3 for the
  base-of-tree nodes Ctenophora + Porifera).
- No file rewritten with CRLF / BOM (synced-drive hygiene); manifest diff is
  add-only rows.
- Sensitivity analysis shows the headline result is stable to Insecta down-sampling.

## 6. Risks & decisions

- **Isoform collapsing without a UniProt `gene:` tag.** *(resolved)* NCBI protein FASTA
  lacks gene IDs in-header; `download_ncbi` maps protein→gene from the GFF3 CDS
  attributes (`gene=` / `locus_tag=` / `Parent=`). Validated: collapsed gene counts match
  NCBI's `protein_coding` totals exactly on GCF and GCA assemblies alike.
- **REST vs CLI.** Uses the Datasets **REST** download endpoint (no `datasets` CLI needed)
  — pulls PROT_FASTA + GENOME_GFF as an in-memory zip; nothing written to `Z:` except the
  final `.longest.fa`.
- **Accession drift.** Every accession was verified live against the API before wiring.
  A future re-verify is a single API pass over `NCBI_PROTEOMES`.
- **Hard limits are real.** Brachiopoda (NCBI-annotated = *Lingula* only), Priapulida,
  Nematomorpha, Acanthocephala stay at n=1 by necessity — Workstream D is the honest
  answer there, not more downloading.
- **`01b` now also fetches metazoan ingroup species.** Slightly outside the script's
  original "outgroups" name, but it already had the generic taxonomy + manifest-upsert
  machinery, so this was the lowest-risk home and avoids a new script.
