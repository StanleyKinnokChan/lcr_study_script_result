# Terminal Low-Complexity Regions as a General Architectural Property of Proteins: A Positional and Compositional Fingerprint Across 756 Proteomes and 43 Lineages

**Stanley Kinnok Chan**¹

¹ *Independent Researcher, London, United Kingdom*

*Correspondence: stanleykinnok.chan@gmail.com · ORCID: 0000-0003-4242-7507*

---

## Abstract

Low-complexity regions (LCRs) — protein segments dominated by one or a few amino acid types — are positionally non-random within eukaryotic proteins, with enrichment at N- and C-termini first described in a single yeast proteome (Coletta et al. 2010, *BMC Systems Biology*) and since confirmed across Tetrapoda (Teekas et al. 2024, *Open Biology*) and all major invertebrate phyla (Chan, preceding study). Whether terminal enrichment is a *general architectural property* of proteins that extends beyond eukaryotes to prokaryotes has not been tested systematically: large-scale prokaryote LCR studies (Ntountoumi et al. 2019; Persi et al. 2023) characterised LCR prevalence and evolution, and terminal (C-terminal) localisation has been noted only for the ribosomal-protein subset (Ntountoumi et al. 2019), but no proteome-wide positional analysis has been reported. Here, applying a single uniform framework (fLPS 2.0, 20-bin positional analysis, SINGLE-type purity ≥70% filter) to 756 proteomes spanning 43 phyla/groups and all three domains of life (962,447 LCR records), we show that terminal LCR enrichment is a domain-independent feature of protein architecture. Pooled bacteria (92 proteomes, 2,906 LCRs, **27.6% terminal**, odds ratio 3.43 relative to the 10% positional null) and archaea (27 proteomes, 513 LCRs, **25.7%**, OR 3.13) are both terminally enriched — to our knowledge the first systematic, proteome-wide analysis of LCR position in prokaryotes. A protein-length-quartile control is decisive: enrichment rises from the shortest to the third length quartile in both domains (Bacteria Q1→Q3: 20.9%→37.0%; Archaea Q1→Q3: 20.7%→36.0%), the opposite of what a short-protein length artefact would produce, and is confirmed by a cluster-robust GEE (Q4 vs Q1 log-odds +1.13, 95% CI [0.90, 1.37]) and a species-level bootstrap. Across eukaryotes, every major supergroup tested is enriched; 39 of 43 phyla/groups survive Holm-Bonferroni correction (Chlorophyta, Rhodophyta, Acanthocephala, and Nematomorpha remain underpowered/provisional). The amino-acid identity and N/C polarity of terminal LCRs form a mechanistic fingerprint. Leucine is the single most terminally biased dominant residue in every domain (2.8–5.3×) and methionine is second in all three eukaryotic groups, implicating N-terminal signal/anchor sequences and initiator-methionine retention — a co-translational / N-terminal-processing layer shared by all cellular life — rather than the acidic/basic degron residues. N/C polarity is lineage-specific, not universal: only Viridiplantae is robustly N-terminal dominant (median asymmetry ratio 2.81), whereas bacteria are balanced at the proteome level (median 1.00), archaea and several protist and lophotrochozoan lineages are C-dominant, and animals and fungi are balanced-to-C-dominant. In bacteria, signal-peptide-bearing proteins carry a threefold-higher N-terminal LCR fraction (34.1% vs 14.6%; OR 3.03, Fisher p = 7.1 × 10⁻⁵), defining a secreted-protein subpopulation rather than a whole-proteome N-terminal bias. Singleton-LCR proteins drive the signal in 41 of 42 groups. We interpret terminal LCR enrichment as a layered architectural property: a shared N-terminal-processing / co-translational baseline present across cellular life, elaborated by lineage-specific degron biology (C-degron in animals, N-degron/PRT6–PCO in plants).

**Keywords:** low-complexity regions, protein termini, protein architecture, compositional bias, co-translational folding, N-terminal processing, signal peptides, intrinsically disordered regions, degron biology, fLPS2

---

## Introduction

Proteins are not compositionally uniform along their length. Low-complexity regions (LCRs) — segments dominated by one or a small number of amino acid types — are distributed non-randomly within protein sequences, clustering in functionally important contexts such as disordered linkers, prion-like domains, and polyamino acid tracts (Wootton and Federhen 1996; Marcotte et al. 1999). The biological significance of LCRs has grown with the recognition that many intrinsically disordered regions (IDRs) are LCR-containing (Romero et al. 2001; van der Lee et al. 2014), and that phase separation — a fundamental mechanism of condensate formation and gene regulation — is often driven by low-complexity IDR sequences (Boija et al. 2018; Shin and Brangwynne 2017; Alberti et al. 2019).

The positional distribution of LCRs *within* proteins — where along the sequence they sit — has received far less attention than their prevalence or functional roles, yet it is a basic feature of protein architecture. Coletta and colleagues (2010) were the first to ask this question: using a single yeast (*Saccharomyces cerevisiae*) proteome, they showed that LCRs are significantly enriched near sequence extremities (p = 7.6 × 10⁻⁶) and that terminal versus central LCR position correlates with distinct functional roles (terminal LCRs: protein connectivity; central LCRs: transcription). This observation was not extended beyond yeast. Teekas and colleagues (2024) revisited the question at scale, reporting that across 308 Tetrapoda species (12 clades) LCRs are significantly enriched in the terminal 5% of protein sequences (the first and last bins of a 20-bin positional map), with 15–25% of all LCRs in terminal positions despite these bins representing only 10% of positional space. This enrichment co-occurred with positively selected sites and was consistent across birds, mammals, reptiles, and amphibians. A companion study (Chan, in preparation) extended the finding to 61 metazoan invertebrate species across 16 phyla and to a small set of non-metazoan eukaryote and prokaryote species, establishing terminal LCR enrichment as pan-metazoan and likely pan-eukaryotic.

In prokaryotes, LCRs are functionally significant and evolutionarily conserved (Ntountoumi et al. 2019), and their formation is compensatorily related to gene duplication — LCR-forming short repeats serve as a transient short-term adaptive mechanism that fades as gene paralogy accumulates (Persi et al. 2023; Saravanan et al. 2025). These studies characterised LCR prevalence, amino acid composition, and gene-level evolutionary context across thousands of bacterial proteomes. The one positional observation to date is that LCRs in bacterial and archaeal *ribosomal* proteins tend to lie at the C-terminus (Ntountoumi et al. 2019); beyond that single protein family, whether LCRs accumulate at protein termini across the proteome in bacteria and archaea, as they do in eukaryotes, has not been tested systematically. If terminal enrichment were a general architectural property of proteins rather than a eukaryote-specific one, it should appear in prokaryotes as well; if it were driven exclusively by eukaryotic ubiquitin-proteasome degron biology, it should not.

Several additional questions of generality motivate the present study. First, the metazoan survey covered only a fraction of animal diversity; Porifera, Xenacoelomorpha, Tardigrada, Rotifera, Nemertea, Nematomorpha, and Collembola were absent. Second, the eukaryotic outgroup analysis covered only four protist lineages, leaving most eukaryotic supergroups — including the entire SAR supergroup — untested. Third, the earlier prokaryote analysis was severely underpowered: with a median of roughly 20 LCRs per bacterial species, individual-species Fisher's tests were uninterpretable.

Here we address these questions by scaling the analysis to 756 proteomes spanning 43 phyla/groups and all three domains of life, encompassing 962,447 LCRs. We ask: (1) For the first time, does terminal LCR enrichment extend to prokaryotes, and does it survive protein-length confound controls? (2) Does the signal span all major eukaryotic supergroups, establishing terminal enrichment as domain-independent? (3) Does terminal enrichment hold across all newly tested metazoan phyla? (4) Is the N/C asymmetry of terminal LCRs a single universal pattern, or is polarity lineage-specific — and, if lineage-specific, what does its variation reveal about the underlying mechanism? Rather than dating the origin of the property, we treat its breadth as evidence of *generality* and use the amino-acid identity and directional polarity of terminal LCRs as a mechanistic fingerprint.

---

## Methods

### Species selection and proteome acquisition

A total of 756 proteomes were assembled from three sources: (1) Ensembl Metazoa release 63 (metazoan invertebrate species); (2) Ensembl Plants and Ensembl release 110 (Viridiplantae, algae, and additional fungi); (3) UniProt reference proteomes (2024 release; UniProt Consortium 2023) for additional eukaryotes and all prokaryotes. Species were selected to maximise taxonomic breadth within each major lineage. The full list is provided in Supplementary Table S1.

Metazoan sampling represents 24 phyla, including newly added Porifera, Xenacoelomorpha, Tardigrada, Rotifera, Collembola, Nematomorpha, and Nemertea, which extend metazoan coverage across the diversity of the animal tree. Phyla represented by a single species (Brachiopoda, Priapulida, Chordata/Cephalochordata, Acanthocephala, Nematomorpha) are flagged throughout as providing single data points rather than phylum-level confirmations; their enrichment values should be interpreted with appropriate caution until additional species are included.

Non-metazoan eukaryote sampling spans all recognised eukaryotic supergroups: Opisthokonta (Fungi: 11 sp.); SAR–Stramenopiles (Oomycota: 26 sp., Bacillariophyta: 5 sp.); SAR–Alveolata (Apicomplexa: 43 sp., Ciliophora: 7 sp., Perkinsozoa: 1 sp.); SAR–Rhizaria (3 sp.); Excavata–Euglenozoa (16 sp.); Excavata–Metamonada (6 sp.); Excavata–Heterolobosea (1 sp.); Amoebozoa (10 sp.); Archaeplastida–Viridiplantae (118 sp.); Archaeplastida–Rhodophyta (3 sp.); Archaeplastida–Chlorophyta (2 sp.); Haptophyta (1 sp.); Cryptophyta (3 sp.); and a residual Protist category (diverse) representing lineages not assigned to the above named groups. The Apicomplexa species include multiple *Plasmodium*, *Eimeria*, *Cryptosporidium*, *Babesia*, *Theileria*, *Hammondia*, and *Toxoplasma* representatives. The Viridiplantae species cover all major angiosperm orders as well as *Selaginella* (lycophyte), *Physcomitrium* and *Marchantia* (bryophytes), and multiple algae (*Chara braunii*). Choanoflagellates, the unicellular sister group to Metazoa, are absent from the current dataset; their inclusion would extend the generality test and is a priority for a subsequent revision (see Discussion).

Prokaryote sampling was expanded substantially: 92 bacterial species spanning Proteobacteria (α, β, γ, δ, ε), Firmicutes, Actinobacteria, Bacteroidetes, Spirochaetes, Cyanobacteria, and additional phyla; 27 archaeal species spanning Euryarchaeota, Crenarchaeota, Thaumarchaeota, and Nanoarchaeota. Species were selected to represent phylum-level diversity. Broader archaeal sampling (including Asgard lineages) would further extend the generality test and is a priority expansion.

Where a species had multiple protein isoforms per gene, only the longest isoform was retained by parsing gene identifiers from FASTA headers. Protein FASTA files were downloaded from Ensembl FTP or UniProt FTP servers.

### LCR detection and filtering

LCRs were detected using fLPS 2.0 (Harrison 2017) with parameters identical to Teekas et al. (2024) and the prior analysis: minimum LCR length = 3 amino acids (-m 3); only SINGLE-type records (single-residue compositional bias) retained in post-processing; purity filter ≥70% (dominant amino acid count / LCR length) applied in post-processing. The pre-compiled macOS/Linux binary was used without recompilation.

The SINGLE-type filter restricts analysis to the most compositionally extreme LCRs (single amino acid dominates ≥70% of the segment). This excludes mixed-composition biases — for example, segments enriched in both G and P typical of disordered linkers — which may have different positional distributions. Sensitivity of the terminal enrichment signal to parameter choice (minimum length ≥6 aa and purity ≥80%; minimum length ≥3 aa and purity ≥60%) and to inclusion of MULTI-type records is reported in Supplementary Table S6.

### Positional binning, enrichment statistics, and multiple testing correction

Positional binning, terminal LCR definition, species-level Fisher's exact test, and protein-level sensitivity analysis were as described previously (one-sided Fisher's exact test; null 10% terminal; terminal = bins 1 or 20 of 20). Phylum-level summaries pool all LCRs from member species before applying Fisher's test. Pooled enrichment is reported as an odds ratio relative to the 10% positional null; p-values are reserved for genuinely bounded per-test comparisons. All statistical tests were computed in Python using SciPy (Virtanen et al. 2020).

Domain-level pooled analysis was added for prokaryotes: all LCRs from Bacteria (n = 92 species) and from Archaea (n = 27 species) were pooled and tested against the 10% null. This approach loses between-species heterogeneity information but provides a powerful aggregate test. The within-domain coefficient of variation of species-level pct_terminal is reported alongside pooled results to characterise heterogeneity (Bacteria CV ≈ 50%; Archaea CV ≈ 45%). To confirm the pooled result is robust to within-species clustering, a cluster-robust GEE logistic regression — terminal (0/1) ~ protein-length-quartile + domain, with species as the cluster (binomial family, exchangeable working correlation) — was fitted, complemented by a species-level bootstrap of the domain means.

To account for multiple comparisons across 43 phylum-level tests, Holm-Bonferroni correction was applied to all phylum-level Fisher's exact p-values. Results are reported as both uncorrected (for comparison with prior work) and corrected. Phyla discussed as significant in the main text have survived Holm-Bonferroni correction unless explicitly noted otherwise.

### Protein-length confound test

The length-stratified analysis was extended to Bacteria, Archaea, and all eukaryotic lineages. Protein-length quartile boundaries were defined globally across all LCRs within each domain or phylum pooled together (not per-species). LCRs were then assigned to quartiles by their host protein's length, and Fisher's exact test for terminal enrichment applied within each quartile. This tests whether terminal enrichment is driven by short proteins mechanically populating terminal bins.

### N-terminal versus C-terminal asymmetry

The asymmetry ratio (pct_nterm / pct_cterm) was calculated per species and summarised per phylum/domain by its median. Because the raw ratio is numerically unstable when pct_cterm approaches zero, all significance testing used a log-odds asymmetry score; the ratio is reported only for interpretability, and is flagged as indicative for phyla with fewer than 50 C-terminal LCRs. A PGLS (phylogenetic generalised least-squares) regression of the Viridiplantae N/C asymmetry ratio on evolutionary tier (an ordinal scale from green algae through grasses) was fitted under a Brownian-motion model to test whether an apparent N-terminal amplification gradient survives correction for phylogenetic non-independence. The variance-covariance matrix was derived from a time-calibrated backbone phylogeny of the sampled Viridiplantae (one representative per genus; deep divergences fixed at established dates, within-tier genera placed as dated polytomies). This backbone is deliberately conservative — it resolves the between-tier structure the gradient concerns but not within-tier topology.

### Signal peptide stratification (Bacteria)

Bacterial LCRs were stratified by protein class (signal peptide present vs. absent) and tested for N-terminal enrichment (bin 1) in each class. Signal-peptide status was taken from UniProtKB SIGNAL sequence-feature annotations, retrieved per organism by taxon identifier from the UniProt REST API. Proteins were matched to annotations by UniProt accession; proteomes not sourced from UniProt (Ensembl Genomes and NCBI Datasets records, which lack UniProt accessions) could not be matched and are reported as an unannotated class rather than misassigned. N-terminal LCR fractions were compared between the with- and without-signal-peptide classes by Fisher's exact test on the pooled bin-1 counts.

### LLPS propensity of terminal versus internal LCRs

Phase-separation propensity of terminal versus internal LCRs was compared in a panel of model organisms using a composition-based proxy applied to the dominant residue of each SINGLE-type LCR. Per-residue weights span the two established phase-separation regimes: aromatic π–π / cation-π stacking (F/Y/W = 1.0, R = 0.5; Vernon et al. 2018) and prion-like poly-Q/poly-N composition (Q/N = 0.8; the signal scored by PLAAC, Lancaster et al. 2014), with flexible spacers (G/S = 0.5) and all other residues scoring 0. Terminal and internal proxy-score distributions were compared per organism by one-sided Mann-Whitney U test (H₁: terminal > internal). The proxy is the offline substitute for PLAAC; PLAAC scores can be supplied in place of the proxy where available.

### Per-domain amino acid composition of terminal LCRs

For each domain group (Bacteria, Archaea, Viridiplantae, Metazoa, and other eukaryotes), the dominant residue of every terminal and internal SINGLE-type LCR was tallied, and a terminal-to-internal enrichment ratio computed per amino acid (frequency among terminal LCRs divided by frequency among internal LCRs). This complements the pooled compositional analysis (which reports the most abundant residues at termini) by identifying which residues are most positionally biased toward termini within each domain, and tests whether the same residues are enriched universally or whether prokaryote terminal LCRs form a distinct compositional class.

### Amino acid identity, driver analysis, purity gradient, GO enrichment

These analyses were completed for all 43 phyla/groups in the dataset. GO term enrichment analysis was performed for species with Ensembl BioMart annotation and sufficient protein identifier overlap; the same five species as the prior analysis reached significance (*Pediculus humanus*, *Tetranychus urticae*, *Caenorhabditis brenneri*, *Lottia gigantea*, *Strongylocentrotus purpuratus*). Systematic GO enrichment across the full 756-species dataset was not completed due to BioMart identifier incompatibilities for non-model organisms. For UniProt-sourced prokaryote proteomes, GO terms are embedded in flat-file annotations and do not require BioMart; GO enrichment for these species is planned as a supplementary extension.

### Data and code availability

All proteomes are publicly available from Ensembl Metazoa, Ensembl Plants, and UniProt
reference proteomes (accessions in Supplementary Table S1). All analysis scripts and the
derived result tables are openly available at [GitHub repository URL — add before
submission] (archived at Zenodo, DOI [add on release]). The pipeline is fully reproducible
from the raw proteome downloads via `main.py`.

---

## Results

### Primary finding: terminal LCR enrichment extends to prokaryotes

We report, for the first time, that LCRs are significantly enriched at the protein termini of both bacteria and archaea — establishing terminal enrichment as a property of protein architecture that is not confined to eukaryotes.

**Bacteria** (92 spp., 2,906 LCRs, 801 terminal, **27.6% terminal**, pooled odds ratio 3.43 relative to the 10% positional null) and **Archaea** (27 spp., 513 LCRs, **25.7% terminal**, pooled OR 3.13) both show significant terminal enrichment when pooled. These terminal fractions exceed those of most eukaryotic lineages (typically 14–24%) and exceed the Tetrapoda range of Teekas et al. (2024).

**The length-confound control is decisive.** Stratifying prokaryote LCRs by host protein length quartile (boundaries defined globally within each domain):

| Domain | Q1 (shortest) | Q2 | Q3 | Q4 (longest) |
|---|---|---|---|---|
| Bacteria | 20.85%*** | 30.50%*** | 37.01%*** | 28.90%*** |
| Archaea | 20.66%** | 26.80%** | 36.00%** | 32.56%** |

All eight quartile tests are significant (Fisher's exact; *** p < 0.001, ** p < 0.01). Enrichment increases monotonically from Q1 to Q3 in both domains — the opposite of what a length-confound artefact would predict (which would inflate terminal fractions disproportionately in short proteins). At Q4 the two domains diverge slightly (Bacteria declines to 28.9%, Archaea remains high at 32.6%), suggesting different length-dependent dynamics worth further investigation, but in neither domain is terminal enrichment an artefact of short median protein lengths.

**Within-domain heterogeneity.** The high CV in bacteria (~50%) and archaea (~45%) — versus ~14% in Viridiplantae across 118 species — reflects that individual prokaryote species are too LCR-poor for reliable individual-species estimates (median ~20 LCRs per bacterial species under our purity filter). The pooled result is therefore a domain-level aggregate driven by LCR-rich taxa (primarily Actinobacteria), not a species-typical value, and we present it as such. Species with sufficient LCRs for individual testing include the Actinobacteria (*Streptomyces coelicolor*, *Streptomyces griseus*, *Mycobacterium tuberculosis*), where individual enrichment is significant (25–36% terminal, p < 0.005). A cluster-robust GEE logistic regression (species as cluster) confirms that terminal enrichment, and its increase across protein-length quartiles, remains highly significant after accounting for species-level clustering (intercept −1.813 in log-odds, above the −2.20 corresponding to the 10% null; Q4 vs Q1 log-odds +1.133, 95% CI [0.899, 1.367], p < 0.001; Bacteria vs Archaea not significantly different, p = 0.49). A species-level bootstrap (5,000 resamples) gives concordant domain means (Bacteria 27.01% [24.33, 29.70]; Archaea 26.63% [22.39, 31.21]).

**Between-domain comparison.** The pooled bacterial (27.6%) and archaeal (25.7%) terminal fractions are the highest of any group in the dataset, exceeding the typical eukaryotic range (14–24%) and the large arthropod samples (Insecta 16.8%, Crustacea 18.1%, Chelicerata 18.1%); they are comparable to the most terminally enriched eukaryotic lineage, Viridiplantae (24.6%).

### Terminal LCR enrichment is general across 43 lineages and all three domains

Across all 756 analysed species, 962,447 LCRs were analysed. The U-shaped LCR positional profile (bins 1 and 20 elevated relative to internal bins) is visually apparent across species in the bin heatmap (Figure 1) and the per-phylum profile overlay (Figure 3).

After Holm-Bonferroni correction across 43 phylum-level tests, pooled terminal enrichment is significant in **39 of 43** phyla/groups (Table 1; Table 2; Figure 2). Four groups do not survive correction and are treated as provisional or underpowered: **Chlorophyta** (p_raw = 0.025 → 0.076), **Rhodophyta** (0.015 → 0.059), **Acanthocephala** (0.071 → 0.091), and **Nematomorpha** (0.046 → 0.091). Within Metazoa, 22 of 24 phyla survive correction (Acanthocephala and Nematomorpha fail, both single-species and small). Several phyla remain represented by a single species; their enrichment values are indicative data points rather than phylum confirmations.

**Table 1. Metazoan phylum-level terminal LCR enrichment.**

| Phylum | N spp. | % Terminal | Holm-corrected sig | Prior study? |
|---|---|---|---|---|
| Porifera | 4 | 17.7 | Yes | NEW |
| Xenacoelomorpha | 3 | 20.2 | Yes | NEW |
| Ctenophora | 2 | 21.1 | Yes | ✓ |
| Placozoa | 2 | 20.2 | Yes | ✓ |
| Cnidaria | 13 | 19.5 | Yes | expanded |
| Tardigrada | 2 | 19.0 | Yes | NEW |
| Rotifera | 3 | 16.4 | Yes | NEW |
| Nematoda | 15 | 19.1 | Yes | expanded |
| Platyhelminthes | 7 | 14.6 | Yes | expanded |
| Annelida | 5 | 16.6 | Yes | expanded |
| Nematomorpha | 1* | 16.6 | No (fails Holm; p_raw 0.046) | NEW |
| Nemertea | 2 | 20.4 | Yes | NEW |
| Brachiopoda | 1* | 21.2 | Yes | ✓ |
| Mollusca | 23 | 18.7 | Yes | expanded |
| Crustacea | 21 | 18.1 | Yes | expanded |
| Chelicerata | 28 | 18.1 | Yes | expanded |
| Myriapoda | 3 | 17.2 | Yes | ✓ |
| Collembola | 2 | 20.3 | Yes | NEW |
| Insecta | 218 | 16.8 | Yes | major expansion |
| Echinodermata | 6† | 19.6 | Yes | expanded |
| Hemichordata | 2 | 17.9 | Yes | ✓ |
| Chordata | 1* | 15.8 | Yes | ✓ |
| Priapulida | 1* | 14.0 | Yes | ✓ |
| Acanthocephala | 1* | 16.4 | No (fails Holm; p_raw 0.071) | NEW |

*Single-species representation: enrichment value is an indicative data point rather than a phylum-level confirmation.
†Echinodermata: *Acanthaster planci* excluded post-hoc (empty fLPS output); 6 analysed.

The within-phylum coefficient of variation (CV) is available for well-sampled phyla (Supplementary Table S3) and is low where sampling is adequate (e.g. Echinodermata, Porifera, Tardigrada, Collembola, Mollusca, Crustacea, Viridiplantae ~14%, Insecta, Oomycota), confirming that the signal is a consistent phylum-level property in well-sampled lineages.

**Terminal enrichment is length-independent.** Length-stratified analysis confirms enrichment in at least two quartiles for all phyla/groups with sufficient LCR counts. Pooled across the entire dataset, the shortest quartile is the *least* enriched (Q1 = 14.0%, Q2 = 22.3%, Q3 = 22.6%, Q4 = 20.0%; all significant), the opposite of what a length-confound artefact would produce.

**Singleton-LCR proteins drive the signal.** Proteins carrying a single LCR show significant terminal enrichment in **41 of 42** phyla/groups; multi-LCR proteins show enrichment in **29 of 42**. Singleton-LCR proteins are the primary driver across the full breadth of the dataset, and the multi-LCR class contributes in a majority of groups.

### Newly sampled basal and understudied phyla: Porifera, Xenacoelomorpha, and Tardigrada

**Porifera (sponges; 4 spp., 17.7%).** Sponges now sampled at four species show terminal enrichment squarely within the mid-range of the metazoan distribution, confirming the property in the phylum without singling it out as exceptional.

**Xenacoelomorpha (3 spp., 20.2%).** With three species now analysed, Xenacoelomorpha shows significant, phylum-consistent terminal enrichment. The phylogenetic position of the clade remains debated; the enrichment value is robust across the sampled species.

**Tardigrada (2 spp., 19.0%).** *Hypsibius exemplaris* and *Paramacrobiotus metropolitanus* are the first tardigrades tested. Tardigrades are known for unusually high proteome disorder; their terminal LCR enrichment within the typical bilaterian range indicates that amplified proteome-wide disorder does not alter the terminal positioning pattern.

### Pan-eukaryotic coverage: all supergroups represented

All major eukaryotic lineages tested show significant pooled terminal enrichment after Holm-Bonferroni correction, with two exceptions among the algae: Chlorophyta and Rhodophyta do not survive correction and are treated as provisional pending additional species. Table 2 lists corrected significance and includes the newly added Heterolobosea (Excavata).

**Table 2. Non-metazoan eukaryote terminal LCR enrichment by lineage.**

| Lineage / Supergroup | N spp. | % Terminal | Holm-corrected sig |
|---|---|---|---|
| Fungi (Opisthokonta) | 11 | 17.4 | Yes |
| Oomycota (SAR–Stramenopiles) | 26 | 21.2 | Yes |
| Bacillariophyta/diatoms (SAR–Stram.) | 5 | 17.5 | Yes |
| Apicomplexa (SAR–Alveolata) | 43 | 12.4 | Yes |
| Ciliophora (SAR–Alveolata) | 7 | 24.2 | Yes |
| Perkinsozoa (SAR–Alveolata) | 1* | 22.0 | Yes |
| Rhizaria (SAR) | 3 | 18.7 | Yes |
| Euglenozoa (Excavata) | 16 | 18.3 | Yes |
| Metamonada (Excavata) | 6 | 18.6 | Yes |
| Heterolobosea (Excavata) | 1* | 31.0 | Yes |
| Amoebozoa | 10 | 16.8 | Yes |
| Viridiplantae (Archaeplastida) | 120 | 24.5 | Yes |
| Rhodophyta (Archaeplastida) | 3 | 13.2 | No (fails Holm; p_raw 0.015) |
| Chlorophyta (Archaeplastida) | 2 | 11.2 | No (fails Holm; p_raw 0.025) |
| Haptophyta | 1* | 13.1 | Yes |
| Cryptophyta | 3 | 19.7 | Yes |
| Protist (diverse) | 12 | 17.1 | Yes |

*Single-species representation.

Euglenozoa (16 spp.) and Metamonada (6 spp.) are among the earliest-diverging eukaryotic clades; their significant enrichment (18.3% and 18.6%) shows the property is not restricted to crown eukaryotes. Heterolobosea, newly added as a single species, shows the highest eukaryotic terminal fraction in the dataset (31.0%) and is reported here as a single data point.

**Apicomplexa** (43 spp., 12.4%): The atypically low terminal fraction reflects the high abundance of asparagine-rich LCRs that are predominantly internal in *Plasmodium* and related species (Muralidharan and Goldberg 2013), functioning in immune evasion and not subject to the positional constraints acting on terminal LCRs in other lineages. Purity gradient analysis confirms that Apicomplexa terminal LCRs are significantly purer than internal LCRs (Δ = +0.007), indicating that LCRs occupying terminal positions in Apicomplexa are subject to qualitative constraints absent from the asparagine-dominated internal pool.

**Chlorophyta and Rhodophyta** (11.2% and 13.2%): the lowest observed enrichment levels, neither significant after multiple-testing correction. In Chlorophyta, *Chlamydomonas reinhardtii* (~10.7% terminal) drives the low phylum-level value while *Ostreococcus lucimarinus* shows individually significant enrichment. Both algal results should be treated as provisional until additional species are added.

### Amino acid composition of terminal LCRs is a mechanistic fingerprint

Terminal amino acid composition by raw abundance (C, E, K most abundant at termini; Q, N not enriched) replicates prior findings across metazoans. The per-domain analysis, which instead ranks residues by their terminal-to-internal enrichment ratio, reveals a strikingly consistent pattern that points to a shared mechanism (Supplementary Table S7). **Leucine is the single most terminally biased dominant residue in every domain group** (2.81× in Bacteria, 5.28× in Archaea, 4.33× in Viridiplantae, 4.81× in Metazoa, 3.45× in other eukaryotes), and **methionine is second across all three eukaryotic groups** (2.35–2.85×). The prokaryotes then diverge: the next most biased residues are arginine (2.05×) and aspartate (1.72×) in Bacteria, and asparagine (2.88×) and arginine/glutamine (2.16×) in Archaea.

The universal leucine bias and the eukaryotic methionine bias implicate hydrophobic signal/anchor sequences and initiator-methionine retention — hallmarks of N-terminal signal-peptide biology and N-terminal processing shared by all cellular life — rather than the acidic/basic degron residues, as the most consistent cross-domain compositional signature of terminal LCRs. This is the compositional half of the mechanistic fingerprint developed below.

### Signal-peptide stratification defines a secreted-protein subpopulation in bacteria

Signal-peptide stratification links the compositional fingerprint to a concrete N-terminal-processing mechanism. Among UniProt-annotated bacterial proteins, N-terminal LCRs (bin 1) are substantially more frequent in signal-peptide-bearing proteins than in those without: **34.1% vs 14.6%** (odds ratio 3.03, Fisher's exact p = 7.1 × 10⁻⁵; 85 vs 425 LCRs across 14 species with matchable annotations). N-terminal signal sequences are therefore a genuine contributor to bacterial N-terminal LCRs — but they define a *subpopulation* of secreted proteins rather than a whole-proteome property: secreted proteins are a minority of the proteome, and non-secreted proteins still show N-terminal LCRs (14.6%) well above the 5% single-bin null, so a signal-peptide-independent N-terminal mechanism also operates. This subpopulation structure is exactly what is expected given that the bacterial proteome is *balanced* in overall N/C polarity (see next section): a strongly N-biased secretory minority is superimposed on an otherwise balanced background. The analysis is limited to the UniProt-sourced subset (Ensembl/NCBI-sourced proteomes lack matchable accessions) and warrants extension as annotation coverage broadens.

### N/C asymmetry is lineage-specific polarity, not a universal pattern

A consistent U-shaped LCR density profile (elevated in bins 1 and 20, depressed in internal bins) is conserved across all domains (Figure 3). The *direction* of the asymmetry — captured by the per-species median asymmetry ratio (pct_nterm / pct_cterm), with significance assessed on a log-odds asymmetry score because the raw ratio is unstable near zero C-terminal counts — is, however, strongly lineage-specific (Figure 4). There is no clean universal dichotomy; the diversity of polarity is itself the central observation.

**Only Viridiplantae is robustly N-terminal dominant** (median ratio 2.81, n = 120). The N-dominance is most extreme in grasses (Poaceae): per-species *Avena*, *Oryza*, *Secale*, and *Triticum* fall in the 3.5–4.5 range. Eudicots are lower but consistently N-leaning; bryophytes are near-balanced. Land plants are the one supergroup with a strong, reproducible N-terminal polarity.

**Bacteria are balanced at the proteome level** (median ratio 1.00, n = 81). We explicitly retract any suggestion of whole-proteome bacterial N-terminal dominance: the strong N-terminal signal seen in bacteria is confined to the secreted, signal-peptide-bearing subpopulation (previous section) and to a few LCR-rich Actinobacteria, and does not characterise the bacterial proteome as a whole.

**Archaea are C-terminal dominant** (median 0.71). Other C-dominant lineages (median ratio < 0.8) include Metamonada (0.46), Platyhelminthes (0.54), Ciliophora (0.66), Rhodophyta (0.71), Mollusca (0.74), Annelida (0.79), Hemichordata (0.79), and Rotifera (0.80). The extreme C-terminal bias in Platyhelminthes, particularly *Schistosoma*, is unexplained and warrants further investigation.

**N-leaning lineages (median ratio > 1.2)** other than Viridiplantae (2.81) are Heterolobosea (1.97, single species), Bacillariophyta (1.66), Apicomplexa (1.30), Rhizaria (1.24), Crustacea (1.22), and Oomycota (1.22).

**Metazoa and Fungi are balanced-to-C-dominant** (Fungi median 0.85; most metazoans 0.8–1.1).

This distribution of polarities — one strongly N-dominant supergroup, a balanced bacterial proteome, and a scatter of C-dominant and mildly N-leaning lineages — is incompatible with a single universal directional mechanism. It is the directional half of the mechanistic fingerprint: a shared, largely symmetric baseline overlaid by lineage-specific processes that tilt polarity in different directions in different clades.

### Protein-level confirmation

Binomial testing confirms terminal LCR enrichment at the protein level across species with ≥10 LCR-containing proteins. The fraction of proteins carrying ≥1 terminal LCR (observed median ~19%, versus 10% null expectation) is significantly elevated in the overwhelming majority of species (binomial test). Terminal LCR enrichment is a per-protein property, not an artefact of a few highly LCR-rich proteins.

### Terminal LCR purity: a minority signature in selected lineages

Purity gradient analysis (one-sided Mann-Whitney U test) was completed for all 43 phyla/groups. The majority show no significant difference between terminal and internal LCR purity, consistent with terminal enrichment being primarily positional rather than qualitative. **Nine phyla** show significantly higher purity in terminal LCRs: **Apicomplexa** (Δ = +0.007), **Rhizaria** (Δ = +0.026, largest effect), **Euglenozoa** (Δ = +0.009), **Viridiplantae** (Δ = +0.002), **Rotifera** (Δ = +0.005), **Annelida** (Δ = +0.013), **Chelicerata** (Δ = +0.006), **Echinodermata** (Δ = +0.014), and **Hemichordata** (Δ = +0.018). Rhizaria shows the largest effect, warranting further investigation of amino acid composition at Rhizaria termini. Terminal LCR purity is *not* significantly elevated in Bacteria (Δ = −0.0003) or Archaea (Δ = −0.010), confirming that prokaryote terminal enrichment is positional rather than qualitative.

### Driver analysis and GO enrichment

Driver analysis is unchanged in direction by the expanded dataset: singleton-LCR proteins drive enrichment in 41 of 42 phyla/groups, and multi-LCR proteins contribute significantly in 29 of 42. GO enrichment was significant in the same five species as previously; the absence of functional concentration in terminal LCRs across well-annotated metazoan species is consistent with a mechanism acting on protein architecture regardless of gene function.

---

## Discussion

### The prokaryote positional finding establishes terminal enrichment as a general architectural property

No prior study has systematically tested whether LCRs are positionally enriched at protein termini across bacterial or archaeal proteomes. The one prior positional observation is Ntountoumi et al.'s (2019) note that LCRs in bacterial and archaeal ribosomal proteins tend to lie at the C-terminus — an incidental, single-protein-family finding rather than a proteome-wide test. The broader prokaryote LCR literature (Ntountoumi et al. 2019; Persi et al. 2023; Saravanan et al. 2025) characterised LCR prevalence, amino acid composition, conservation, and gene-level evolutionary dynamics, but did not apply positional binning across the proteome. Coletta et al. (2010) first showed proteome-wide terminal enrichment, but in a single yeast proteome. A systematic positional test across prokaryotes has therefore remained outstanding.

The pooled bacterial terminal fraction (27.6%, OR 3.43) and archaeal fraction (25.7%, OR 3.13) both exceed the typical eukaryotic range (14–24%). The length-stratified control is decisive: enrichment increases from Q1 to Q3 in both domains (the opposite of a length artefact), and all eight quartile tests are significant. This cannot be dismissed as a sampling or methodological artefact, and it is robust to species-level clustering (GEE Q4 vs Q1 log-odds +1.13, 95% CI [0.90, 1.37]) and to a species bootstrap. High within-domain CV (bacteria ~50%; archaea ~45%) indicates real heterogeneity across species; the pooled result reflects a domain-level aggregate driven by LCR-rich taxa (primarily Actinobacteria) rather than a uniform property of every bacterium, and we present it honestly as such. The demonstration that the property holds in all three domains is best read as evidence of *generality* — terminal LCR enrichment is a feature of protein architecture that does not depend on eukaryote-specific machinery — rather than as a claim about the deep-time origin of the property, which these data cannot address.

### A mechanistic fingerprint: shared N-terminal processing plus lineage-specific degron biology

The amino-acid identity and directional polarity of terminal LCRs together form a fingerprint that separates a shared baseline layer from lineage-specific elaborations.

**Compositional signature — a shared N-terminal-processing layer.** When residues are ranked by terminal-to-internal enrichment rather than raw abundance, leucine is the single most terminally biased residue in every domain, and methionine is second across all three eukaryotic groups. Both point to N-terminal biology shared by all cellular life rather than to eukaryote-specific degrons: leucine-rich hydrophobic runs are the defining feature of N-terminal signal peptides, signal anchors, and transmembrane segments, all of which concentrate near protein termini; and a terminal methionine bias is the expected signature of initiator-methionine retention and N-terminal processing. This compositional evidence dovetails with the signal-peptide stratification result, in which signal-peptide-bearing bacterial proteins carry a threefold-higher N-terminal LCR fraction.

**Directional signature — lineage-specific polarity.** The lineage-specific distribution of N/C polarity shows that this shared baseline is overlaid by different processes in different clades. Only land plants are strongly N-dominant; bacteria are balanced at the proteome level with N-dominance confined to the secretory subpopulation; archaea and several protist and lophotrochozoan lineages are C-dominant; animals and fungi are balanced-to-C. If a single universal mechanism set the polarity, one direction would predominate; instead, polarity varies by lineage, which is direct evidence for lineage-specific amplification on top of a shared, largely symmetric baseline.

Candidate lineage-specific mechanisms are well characterised. In animals, C-terminal degrons recognised by E3 ubiquitin ligases (Koren et al. 2018; Varshavsky 2019), tubulin-code C-terminal tails (Janke and Magiera 2020), and RNA-binding C-terminal regions all concentrate low-complexity sequence at the C-terminus, consistent with the balanced-to-C polarity of Metazoa and Fungi. In land plants, the N-degron pathway involving PRT1/PRT6 E3 ligases and Plant Cysteine Oxidase (PCO) enzymes links N-terminal cysteine oxidation to oxygen-sensing and proteasomal degradation of ERF-VII transcription factors (Gibbs et al. 2011; Holdsworth et al. 2020); this pathway is elaborated in land plants and broadly coincides with their strong N-terminal polarity. In bacteria, N-terminal signal peptides and export sequences drive the secretory N-terminal subpopulation, while C-terminal ssrA degradation tags (tmRNA-mediated) act at the C-terminus. Ribosomal proteins contribute a further C-terminal subpopulation: their LCRs were previously observed to localise to the C-terminus in bacteria and archaea (Ntountoumi et al. 2019). Our proteome-wide analysis places that single-family observation in context — these opposing terminal subpopulations (N-terminal secretory, C-terminal ssrA and ribosomal) coexist, and their net effect is the balanced bacterial proteome we observe rather than a single dominant polarity. N-formylmethionine as an N-degron (Piatkov et al. 2015) is sometimes cited, but fMet is co-translationally cleaved from most bacterial proteins by methionine aminopeptidase, and the bacterial ClpS/ClpAP N-end-rule pathway recognises hydrophobic residues (L, F, W, Y) rather than the E, K, C residues typical of eukaryotic terminal LCRs; the fMet-degron hypothesis therefore requires direct testing.

We synthesise these into a **layered architectural model**: a shared layer — co-translational ribosome kinetics near start and stop codons (Irastortza-Olaziregi and Amster-Choder 2021; Pechmann and Frydman 2013) together with N-terminal processing (signal-peptide cleavage, initiator-methionine excision) — generates a baseline terminal LCR enrichment in all cellular life, with a compositional signature dominated by leucine and (in eukaryotes) methionine; upon this baseline, lineage-specific degron and processing biology superimpose directional polarity (C-degron in animals, N-degron/PRT6–PCO in plants, secretory signal peptides in bacteria). We present this as a mechanistic model of protein architecture, not as a reconstruction of evolutionary history.

### The Viridiplantae N/C polarity gradient does not survive phylogenetic correction

The 118-species plant dataset confirms strong N-terminal polarity in Viridiplantae and shows an apparent gradient in tier means — near-balanced in bryophytes, intermediate in eudicots, extreme in grasses (3.5–4.5). This apparent gradient must, however, be interpreted cautiously: a Brownian-motion PGLS on evolutionary tier does not recover a significant slope (0.296, p = 0.157, R² = 0.024, n = 84 genera) despite a strong non-phylogenetic correlation (Pearson r = 0.59, p < 0.001). The raw correlation is substantially inflated by the phylogenetic clustering of the large, closely related grass and eudicot samples; a per-lineage evolutionary gradient is not robustly supported once shared ancestry is modelled. We report this negative result explicitly as a rigor check: the supergroup-level N-dominance of land plants is robust, but the fine-grained mosses-to-grasses "gradient" is not. If the PRT6/PCO N-degron pathway contributed the plant N-terminal polarity, its substrates would be expected to carry N-terminal LCRs; however, the per-domain composition analysis shows Viridiplantae termini dominated by leucine and methionine rather than the cysteine/arginine residues these pathways recognise, so any PRT6/PCO contribution is not the primary compositional signature of plant terminal LCRs.

### Terminal LCRs and phase separation: an open connection

LCRs are strongly associated with liquid-liquid phase separation (LLPS) and biomolecular condensate formation (Boija et al. 2018; Alberti et al. 2019). We tested whether terminal LCRs have systematically elevated LLPS propensity relative to internal LCRs, using a composition-based proxy across a panel of model organisms. The result is organism-specific rather than universal: of seven organisms, only *Arabidopsis thaliana* shows significantly higher terminal LLPS propensity (Mann-Whitney U, p = 0.033), while in *Drosophila melanogaster* and *Saccharomyces cerevisiae* internal LCRs are the more phase-prone class, and the remainder show no significant difference. Terminal LCRs therefore do not carry a domain-general elevation in phase-separation propensity; where an association exists it is lineage-specific. Confirmatory PLAAC scoring and broader taxonomic sampling would sharpen this picture, but the positional enrichment reported here is not, in itself, driven by a uniform phase-separation signature.

### The Platyhelminthes, algal, and single-species outliers

**Platyhelminthes** (14.6% terminal; C-dominant, median 0.54) shows the lowest metazoan terminal enrichment; the expansion to 7 species confirms this is a genuine phylum-level property. The systematic C-terminal bias, particularly in *Schistosoma*, contrasts with other lophotrochozoans and is mechanistically unexplained.

**Chlorophyta and Rhodophyta** do not survive Holm-Bonferroni correction and are provisional pending additional species. **Acanthocephala** and **Nematomorpha** are single-species phyla that fail correction on power grounds; their effect sizes are within the bilaterian range.

### Limitations

1. **Prokaryote statistical model.** The pooled Fisher's test does not itself correct for between-species clustering, but a cluster-robust GEE logistic regression confirms the result is not an artefact of treating clustered LCRs as independent: the baseline terminal probability exceeds the 10% null (intercept −1.813 in log-odds, above −2.20), the increase across protein-length quartiles is strongly significant (Q4 vs Q1 log-odds +1.133, 95% CI [0.899, 1.367], p < 0.001), and Bacteria and Archaea do not differ (p = 0.49). A species-level bootstrap gives concordant domain means (Bacteria 27.01% [24.33, 29.70]; Archaea 26.63% [22.39, 31.21]). The high bacterial CV (~50%) still indicates the pooled point estimate is influenced by a subset of LCR-rich taxa, and we present it as a domain aggregate rather than a species-typical value.

2. **fLPS parameter sensitivity.** A three-parameter sensitivity analysis (Supplementary Table S6) confirms that SINGLE-type terminal enrichment is robust to the LCR length threshold and purity cutoff. Prokaryote terminal enrichment persists and in fact strengthens under stringent parameters (length ≥6, purity ≥80%: Bacteria 34.14%, Archaea 29.01%) and remains significant under relaxed parameters (length ≥3, purity ≥60%: Bacteria 23.98%, Archaea 22.36%). Terminal enrichment does not extend to MULTI-type LCRs (Bacteria 9.06%, Archaea 11.34%; both non-significant), indicating the signal is specific to the most compositionally extreme, single-residue-dominated regions. Under the stringent setting several already-weak eukaryotic groups (Metamonada, Chlorophyta, Rhodophyta, Haptophyta, Priapulida, Chordata, Acanthocephala, Nematomorpha) lose significance; this caveat qualifies those specific groups but not the prokaryote or well-sampled eukaryote results.

3. **Multiple testing.** Holm-Bonferroni correction was applied to all 43 phylum-level tests; 39 survive. Chlorophyta, Rhodophyta, Acanthocephala, and Nematomorpha do not and are treated as provisional/underpowered.

4. **Single-species groups.** Several phyla and protist lineages (including the newly added Heterolobosea) are represented by one species. These provide initial data points but not phylum-level confirmations, and the domain aggregates should not be over-read as species-typical.

5. **Sampling breadth.** Choanoflagellates and broader archaeal diversity are absent; adding them would extend the generality test to further lineages and is a priority for the next revision (framed as a test of generality, not of origin).

6. **GO enrichment gap.** GO analysis covers only 5 of 756 species (BioMart compatibility); prokaryote GO enrichment via UniProt flat-file annotations is planned.

7. **SINGLE-type filter.** Mixed-composition LCRs are excluded; their positional distribution is unknown.

8. **Asymmetry ratio instability.** The ratio is numerically unstable for phyla with few C-terminal LCRs; all significance testing therefore used a log-odds asymmetry score, with ratios reported only for interpretability.

---

## Conclusions

Applying a single positional framework to 756 proteomes spanning 43 phyla/groups and all three domains of life (962,447 LCRs), we show that terminal LCR enrichment is a general architectural property of proteins rather than a eukaryote-specific one. Pooled bacteria (27.6%, OR 3.43, 92 species) and archaea (25.7%, OR 3.13, 27 species) are terminally enriched — to our knowledge the first systematic, proteome-wide positional analysis of LCRs in prokaryotes — and the enrichment is not explained by protein length: it increases across length quartiles in the direction opposite to a length artefact and is robust to species-level clustering (GEE, bootstrap) and to parameter choice (sensitivity analysis). Across eukaryotes, every major supergroup is represented and 39 of 43 phyla/groups survive Holm-Bonferroni correction (Chlorophyta, Rhodophyta, Acanthocephala, and Nematomorpha remain provisional/underpowered). Singleton-LCR proteins drive the signal in 41 of 42 groups; multi-LCR proteins contribute in 29 of 42.

The amino-acid identity and N/C polarity of terminal LCRs form a mechanistic fingerprint. Leucine is the most terminally biased residue in every domain and methionine second in all eukaryotic groups, implicating an N-terminal-processing / co-translational layer shared by all cellular life; in bacteria, signal-peptide-bearing proteins carry a threefold-higher N-terminal LCR fraction (34.1% vs 14.6%; OR 3.03, p = 7.1 × 10⁻⁵), defining a secreted-protein subpopulation. N/C polarity is lineage-specific rather than universal — only Viridiplantae is robustly N-dominant, bacteria are balanced, archaea and several other lineages are C-dominant — which is incompatible with a single universal directional mechanism. Together the findings support a layered architectural model: a shared N-terminal-processing / co-translational baseline, elaborated by lineage-specific degron biology (C-degron in animals, N-degron/PRT6–PCO in plants) and secretory signal-peptide biology in bacteria.

---

## Acknowledgements

The author thanks the Ensembl Metazoa, Ensembl Plants, and UniProt teams for providing freely downloadable proteomes, and Paul Harrison for making fLPS 2.0 freely available.

---

## Declarations

**Funding.** This research received no external funding.

**Competing interests.** The author declares no competing interests.

**Author contributions.** S.K.C. conceived and designed the study, performed all analyses, and wrote the manuscript.

**Use of AI assistance.** A large-language-model assistant (Anthropic Claude) was used to support manuscript editing, to reconcile numeric claims in the text against the analysis output tables, and to provide simulated editorial review. The author designed the study, implemented and executed all analyses, verified all results, and takes full responsibility for the content of this manuscript.

---

## References

Alberti S, Gladfelter A, Mittag T (2019). Considerations and challenges in studying liquid-liquid phase separation and biomolecular condensates. *Cell* 176:419–434.

Boija A, et al. (2018). Transcription factors activate genes through the phase-separation capacity of their activation domains. *Cell* 175:1842–1855.

Coletta A, Pinney JW, Solís DY, Marsh J, Pettifer SR, Attwood TK (2010). Low-complexity regions within protein sequences have position-dependent roles. *BMC Systems Biology* 4:43.

Gibbs DJ, Lee SC, Isa NM, Gramuglia S, Fukao T, Bassel GW, Correia CS, Corbineau F, Theodoulou FL, Bailey-Serres J, Holdsworth MJ (2011). Homeostatic response to hypoxia is regulated by the N-end rule pathway in plants. *Nature* 479:415–418.

Harrison PM (2017). fLPS: Fast discovery of compositional biases for the protein universe. *BMC Bioinformatics* 18:476.

Holdsworth MJ, Vicente J, Sharma G, Abbas M, Estavillo GM (2020). The plant N-degron pathways of ubiquitin-mediated proteolysis. *Journal of Integrative Plant Biology* 62:70–89.

Irastortza-Olaziregi M, Amster-Choder O (2021). Coupled transcription-translation in prokaryotes: an old couple with new surprises. *Frontiers in Microbiology* 11:624830.

Janke C, Magiera MM (2020). The tubulin code and its role in controlling microtubule properties and functions. *Nature Reviews Molecular Cell Biology* 21:307–326.

Koren I, et al. (2018). The eukaryotic proteome is shaped by E3 ubiquitin ligases targeting C-terminal degrons. *Cell* 173:1622–1635.

Lancaster AK, Nutter-Upham A, Lindquist S, King OD (2014). PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition. *Bioinformatics* 30(17):2501–2502.

Marcotte EM, et al. (1999). A census of protein repeats. *Journal of Molecular Biology* 293:151–160.

Muralidharan V, Goldberg DE (2013). Asparagine repeats in *Plasmodium falciparum* proteins: Good for nothing? *PLoS Pathogens* 9:e1003488.

Ntountoumi C, et al. (2019). Low complexity regions in the proteins of prokaryotes perform important functional roles and are highly conserved. *Nucleic Acids Research* 47:9998–10009.

Pechmann S, Frydman J (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology* 20:237–243.

Persi E, Wolf YI, Karamycheva S, Makarova KS, Koonin EV (2023). Compensatory relationship between low-complexity regions and gene paralogy in the evolution of prokaryotes. *Proceedings of the National Academy of Sciences* 120:e2300154120.

Piatkov KI, Vu TTM, Hwang C-S, Varshavsky A (2015). Formyl-methionine as a degradation signal at the N-termini of bacterial proteins. *Microbial Cell* 2(10):376–393.

Romero P, et al. (2001). Sequence complexity of disordered protein. *Proteins* 42:38–48.

Saravanan V, Kravetz A, Battistuzzi FU (2025). Higher frequency of prokaryotic low complexity regions in core and orthologous genes. *Frontiers in Bioinformatics* 5:1673480.

Shin Y, Brangwynne CP (2017). Liquid phase condensation in cell physiology and disease. *Science* 357:eaaf4382.

Teekas L, Sharma S, Vijay N (2024). Terminal regions of a protein are a hotspot for low complexity regions and selection. *Open Biology* 14:230439.

UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research* 51:D523–D531.

van der Lee R, et al. (2014). Classification of intrinsically disordered regions and proteins. *Chemical Reviews* 114:6589–6631.

Varshavsky A (2019). N-degron and C-degron pathways of protein degradation. *Proceedings of the National Academy of Sciences* 116:358–366.

Vernon RM, et al. (2018). Pi-pi contacts are an overlooked protein feature relevant to phase separation. *eLife* 7:e31486.

Virtanen P, et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17:261–272.

Wootton JC, Federhen S (1996). Analysis of compositionally biased regions in sequence databases. *Methods in Enzymology* 266:554–571.

---

## Figure Legends

**Figure 1. LCR positional distribution across 756 proteomes.** Heatmap showing the fraction of LCRs in each of 20 equally spaced positional bins (bin 1 = N-terminal 5%; bin 20 = C-terminal 5%) for each analysed species. Rows ordered by domain and phylum (Bacteria → Archaea → non-metazoan eukaryotes → Metazoa). Colour scale: fraction of a species' LCRs in that bin (0–0.15). Dashed blue lines mark the terminal bins.

**Figure 2. Terminal LCR enrichment by phylum.** Bar chart showing pooled % terminal LCRs per phylum/group, ordered by lineage. Black dots show individual-species values. Red dotted line: 10% null expectation. Grey shaded band: Tetrapoda range from Teekas et al. (2024; 15–25%). Thirty-nine of 43 groups significantly exceed the null after Holm-Bonferroni correction; the four exceptions (Chlorophyta, Rhodophyta, Acanthocephala, Nematomorpha) are underpowered/provisional. Bacteria (27.6%) and Archaea (25.7%) show among the highest pooled terminal fractions.

**Figure 3. Conserved U-shaped LCR positional profile across all phyla.** Mean fraction of LCRs per positional bin (1–20) per phylum/domain, averaged across member species. Dashed grey vertical lines: terminal bins 1 and 20. Dotted horizontal line: 5% uniform null. All phyla show elevated LCR density at both termini and a depressed internal plateau.

**Figure 4. Lineage-specific N/C terminal polarity.** Per-species asymmetry ratio (pct_nterm / pct_cterm; log-scaled for display; significance assessed on a log-odds asymmetry score; ratio indicative only for phyla with <50 C-terminal LCRs), grouped by phylum/domain. Only Viridiplantae (median 2.81; grasses 3.5–4.5) is robustly N-terminal dominant. Bacteria are balanced at the proteome level (median 1.00). Archaea (0.71) and several protist and lophotrochozoan lineages are C-dominant; Metazoa and Fungi are balanced-to-C. The diversity of polarity across lineages is incompatible with a single universal directional mechanism.

**Figure 5. Protein-length stratified terminal LCR enrichment in prokaryotes.** Grouped bar chart showing % terminal LCRs in four protein-length quartiles (Q1 shortest; Q4 longest; boundaries defined globally within each domain) for Bacteria (blue) and Archaea (red). Error bars: 95% CI from Fisher's exact test. The Q1→Q3 increase in both domains (Bacteria 20.9% → 37.0%; Archaea 20.7% → 36.0%) rules out a length-confound artefact.

---

## Supplementary Tables

**Supplementary Table S1.** Full species list (756 entries): species name, phylum, data source (Ensembl / UniProt), proteome ID, protein count, LCR count, pct_terminal.

**Supplementary Table S2.** Driver analysis results for all 43 phyla/groups: pct_terminal for singleton-LCR proteins, multi-LCR proteins; odds ratios and p-values for each class. Singleton-LCR proteins significant in 41/42; multi-LCR proteins significant in 29/42.

**Supplementary Table S3.** Within-phylum coefficient of variation of pct_terminal: phyla/groups with n ≥ 2 species. Includes phylum, n_species, mean_pct_terminal, std_pct_terminal, CV.

**Supplementary Table S4.** Length-stratified analysis for all 43 phyla/groups: pct_terminal by protein-length quartile; Fisher's exact p-value and significance per quartile.

**Supplementary Table S5.** Purity gradient analysis for all 43 phyla/groups: mean purity of terminal vs. internal LCRs, Δ purity, Mann-Whitney U p-value, significance flag. Nine phyla significant (Apicomplexa, Rhizaria, Euglenozoa, Viridiplantae, Rotifera, Annelida, Chelicerata, Echinodermata, Hemichordata).

**Supplementary Table S6.** fLPS parameter sensitivity analysis: pct_terminal and pooled significance for all 43 phyla/groups under three parameter combinations — (current: length ≥3, purity ≥70%), (stringent: length ≥6, purity ≥80%), (relaxed: length ≥3, purity ≥60%) — plus a MULTI-type LCR positional analysis. SINGLE-type prokaryote terminal enrichment is significant across all settings and strengthens under the stringent setting (Bacteria 34.14%, Archaea 29.01%); MULTI-type LCRs show no terminal enrichment (Bacteria 9.06%, Archaea 11.34%).

**Supplementary Table S7.** Per-domain amino acid composition of terminal LCRs: dominant residue enrichment ratio (terminal / internal) per amino acid, for Bacteria, Archaea, Viridiplantae, Metazoa, and other eukaryotes separately. Leucine is the most terminally enriched residue in every domain (2.8–5.3×); methionine is second across the eukaryotic groups (2.35–2.85×).

**Supplementary Table S8.** Signal peptide stratification of bacterial N-terminal LCRs: per-species N-terminal (bin 1) LCR count and fraction for proteins with a signal peptide, without a signal peptide, and unannotated, based on UniProtKB SIGNAL feature annotations. Pooled across species with matchable UniProt annotations, N-terminal LCRs are 3-fold enriched in signal-peptide-bearing proteins (34.1% of 85 LCRs) relative to non-secreted proteins (14.6% of 425 LCRs; Fisher's exact odds ratio 3.03, p = 7.1 × 10⁻⁵). This defines a secreted-protein subpopulation, not a whole-proteome N-terminal bias.

**Supplementary Table S9.** LLPS propensity (composition-based proxy) of terminal versus internal LCRs per model organism: n terminal/internal LCRs, median proxy score for each, percent aromatic (F/Y/W), and one-sided Mann-Whitney U p-value (terminal > internal). Of seven organisms, only *Arabidopsis thaliana* reaches significance (p = 0.033).

**Supplementary Figures**

**Supplementary Figure 1.** Per-species terminal LCR % distributions for each phylum, violin/strip chart format.

**Supplementary Figure 2.** Full length-stratified enrichment results for all 43 phyla/groups (4 quartiles × 43 groups heatmap).

**Supplementary Figure 3.** Asymmetry ratio distributions for Viridiplantae only, coloured by plant order. Shows the strong supergroup-level N-terminal polarity; note the mosses-to-grasses tier gradient is not robust to phylogenetic correction (Supplementary Figure 7).

**Supplementary Figure 4.** Within-phylum CV versus number of species, all phyla with n ≥ 2.

**Supplementary Figure 5.** Amino acid composition of terminal versus internal LCRs: enrichment ratio (terminal/internal) per amino acid, pooled across all metazoans. C, E, K enriched; Q, N depleted at termini.

**Supplementary Figure 6.** Purity gradient: distributions of purity for terminal (bins 1+20) and internal (bins 2–19) LCRs for the nine significant phyla (Apicomplexa, Rhizaria, Euglenozoa, Viridiplantae, Rotifera, Annelida, Chelicerata, Echinodermata, Hemichordata). Violin plots with Mann-Whitney U p-values.

**Supplementary Figure 7.** PGLS regression of N/C asymmetry ratio on evolutionary tier across 84 Viridiplantae genus representatives, under a Brownian-motion model on a time-calibrated backbone phylogeny. Contrasts the strong non-phylogenetic correlation (OLS/Pearson r = 0.59, p < 0.001) with the non-significant phylogenetically-corrected slope (PGLS slope = 0.296, p = 0.157, R² = 0.024).
