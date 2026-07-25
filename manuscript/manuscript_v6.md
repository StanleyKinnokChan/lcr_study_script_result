# Terminal Low-Complexity Regions Are Enriched Across All Domains of Life: Evidence from 724 Proteomes Spanning 42 Phyla

**Stanley Kinnok Chan**¹

¹ *[Institution, City, Country]*

*Correspondence: stanleykinnok.chan@gmail.com*

---

## Abstract

Low-complexity regions (LCRs) — protein segments dominated by one or a few amino acid types — are known to be positionally non-random within eukaryotic proteins, with enrichment at N- and C-termini first described in a single yeast proteome (Coletta et al. 2010, *BMC Systems Biology*), and since confirmed across Tetrapoda (Teekas et al. 2024, *Open Biology*) and all major invertebrate phyla (Chan, preceding study). Whether this terminal enrichment extends to prokaryotes has never been tested: prior large-scale prokaryote LCR studies (Ntountoumi et al. 2019; Persi et al. 2023) analysed LCR prevalence and evolution but did not apply positional binning. Here we provide the first evidence that terminal LCR enrichment extends to all three domains of life, applying the same computational framework (fLPS 2.0, 20-bin positional analysis, purity ≥70% SINGLE-type filter) to 724 proteomes spanning 42 phyla and 929,801 LCR records. Pooled analysis of 87 bacterial proteomes (27.3% terminal, p~0) and 21 archaeal proteomes (26.2% terminal, p~0) reveals significant terminal enrichment in both prokaryotic domains — enrichment that increases from the shortest to the intermediate protein-length quartile (Bacteria Q1→Q3: 20.6%→37.0%; Archaea Q1→Q3: 19.1%→35.7%), ruling out a length-confound artefact. This prokaryote finding is the primary novel result of this study; all prior positional work was confined to eukaryotes. Among eukaryotes, all 16 major lineages tested show significant pooled enrichment, including all SAR supergroup members, both deep-branching Excavata clades, and all Archaeplastida; 40 of 42 phyla/groups survive Holm-Bonferroni correction across the 42 phylum-level tests (Chlorophyta p=0.025 does not survive correction and is treated as provisional). Among metazoans, 23 of 24 phyla are confirmed including newly tested Porifera (23.1%), Tardigrada (19.0%), and Collembola (20.3%); Xenacoelomorpha is represented by a single species (*Hofstenia miamia*, 21.8%) and should be treated as a data point rather than a phylum confirmation pending additional sampling. The N/C asymmetry reveals a supergroup-specific mechanistic signature not previously described: Viridiplantae and LCR-rich bacteria are N-terminal dominant (grasses: ratios 3.5–4.5; bacteria: ratios 2.5–4.7), whereas most Metazoa and Fungi are C-terminal dominant or balanced, and Ciliophora, Metamonada, and Platyhelminthes are C-terminal dominant. This asymmetry pattern is incompatible with a single universal mechanism and discriminates between candidate explanations. Whether the prokaryote enrichment is ancestral (consistent with a LUCA-level origin, ~3.5–4 Bya) or arose convergently cannot be resolved from the current data; both scenarios are equally supported. Singleton-LCR proteins drive the signal in 41 of 42 phyla/groups. These findings support a layered mechanistic model — universal translational kinetics and N-terminal processing generating baseline terminal enrichment, amplified by lineage-specific degron biology in eukaryotes.

**Keywords:** low-complexity regions, protein termini, compositional bias, prokaryote protein architecture, eukaryotic evolution, intrinsically disordered regions, fLPS2, pan-genomics, LUCA

---

## Introduction

Proteins are not compositionally uniform along their length. Low-complexity regions (LCRs) — segments dominated by one or a small number of amino acid types — are distributed non-randomly within protein sequences, clustering in functionally important contexts such as disordered linkers, prion-like domains, and polyamino acid tracts (Wootton and Federhen 1996; Marcotte et al. 1999). The biological significance of LCRs has grown considerably with the recognition that many intrinsically disordered regions (IDRs) are LCR-containing (Romero et al. 2001; van der Lee et al. 2014), and that phase separation — now understood as a fundamental mechanism of condensate formation and gene regulation — is often driven by low-complexity IDR sequences (Boija et al. 2018; Shin and Brangwynne 2017; Alberti et al. 2019).

The positional distribution of LCRs *within* proteins has received far less attention than their prevalence or functional roles. Coletta and colleagues (2010) were the first to ask where within a protein LCRs sit: using a single yeast (*Saccharomyces cerevisiae*) proteome, they showed that LCRs are significantly enriched near sequence extremities (p=7.6×10⁻⁶) and that terminal versus central LCR position correlates with distinct functional roles (terminal LCRs: protein connectivity; central LCRs: transcription). This observation was not extended beyond yeast. Teekas and colleagues (2024) revisited the question at scale, reporting that across 308 Tetrapoda species (12 clades) LCRs are significantly enriched in the terminal 5% of protein sequences (the first and last bins of a 20-bin positional map), with 15–25% of all LCRs in terminal positions despite these bins representing only 10% of positional space. This enrichment co-occurred with positively selected sites and was consistent across birds, mammals, reptiles, and amphibians. A companion study (Chan, in preparation) extended this finding to 61 metazoan invertebrate species across 16 phyla and to 33 non-metazoan eukaryote and prokaryote species (94 total), establishing terminal LCR enrichment as pan-metazoan and likely pan-eukaryotic.

In prokaryotes, LCRs are known to be functionally significant and evolutionarily conserved (Ntountoumi et al. 2019), and their formation is compensatorily related to gene duplication — LCR-forming short repeats serve as a transient short-term adaptive mechanism that fades as gene paralogy accumulates (Persi et al. 2023; Saravanan et al. 2025). These studies characterised LCR prevalence, amino acid composition, and gene-level evolutionary context across thousands of bacterial proteomes. Critically, none of them asked where within individual proteins the LCRs reside. The prokaryote positional question — whether LCRs accumulate at protein termini in bacteria and archaea as they do in eukaryotes — has remained entirely untested.

Several additional open questions motivate the present study. First, the metazoan survey covered only a fraction of animal phylogenetic diversity; Porifera (sponges), Xenacoelomorpha, Tardigrada, Rotifera, Nemertea, Nematomorpha, and Collembola were absent. Second, the eukaryotic outgroup analysis covered only four protist lineages, leaving most eukaryotic supergroups — including the entire SAR supergroup — untested. Third, the earlier prokaryote analysis was severely underpowered: with a median of 21 LCRs per bacterial species, individual-species Fisher's tests were uninterpretable.

Here we address all four questions by scaling the analysis to 724 proteomes spanning 42 phyla and all three domains of life, encompassing 929,801 LCRs. We ask: (1) For the first time, does terminal LCR enrichment extend to prokaryotes, and does it survive protein-length confound controls? (2) Does the signal span all major eukaryotic supergroups? (3) Does terminal enrichment hold across all newly tested metazoan phyla, including the most basally branching animals? (4) Is the N/C asymmetry pattern consistent across lineages, or do supergroup-specific reversals carry mechanistic information?

---

## Methods

### Species selection and proteome acquisition

A total of 724 proteomes were assembled from three sources: (1) Ensembl Metazoa release 63 (metazoan invertebrate species); (2) Ensembl Plants and Ensembl release 110 (Viridiplantae, algae, and additional fungi); (3) UniProt reference proteomes 2024 (additional eukaryotes and all prokaryotes). Species were selected to maximise taxonomic breadth within each major lineage. The full list is provided in Supplementary Table S1.

Metazoan sampling represents 24 phyla: Acanthocephala (1 sp.), Annelida (5), Brachiopoda (1), Chelicerata (28), Chordata/Cephalochordata (1), Cnidaria (13), Collembola (2), Crustacea (21), Ctenophora (1), Echinodermata (7, of which *Acanthaster planci* was excluded post-hoc due to empty fLPS output, leaving 6 analysed), Hemichordata (1), Insecta (218), Mollusca (23), Myriapoda (1), Nematoda (15), Nematomorpha (1), Nemertea (1), Placozoa (1), Platyhelminthes (7), Porifera (2), Priapulida (1), Rotifera (1), Tardigrada (2), and Xenacoelomorpha (1). The newly added phyla — Porifera, Xenacoelomorpha, Tardigrada, Rotifera, Collembola, Nematomorpha, and Nemertea — extend metazoan coverage to include the most basally branching animals. Phyla represented by a single species (Acanthocephala, Xenacoelomorpha, Ctenophora, Placozoa, Rotifera, Nematomorpha, Nemertea, Brachiopoda, Myriapoda, Hemichordata) are flagged throughout as providing single data points rather than phylum-level confirmations; their enrichment values should be interpreted with appropriate caution until additional species are included.

Non-metazoan eukaryote sampling spans all recognised eukaryotic supergroups: Opisthokonta (Fungi: 10 sp.); SAR–Stramenopiles (Oomycota: 26 sp., Bacillariophyta: 5 sp.); SAR–Alveolata (Apicomplexa: 43 sp., Ciliophora: 7 sp., Perkinsozoa: 1 sp.); SAR–Rhizaria (3 sp.); Excavata–Euglenozoa (15 sp.); Excavata–Metamonada (5 sp.); Amoebozoa (9 sp.); Archaeplastida–Viridiplantae (118 sp.); Archaeplastida–Rhodophyta (3 sp.); Archaeplastida–Chlorophyta (2 sp.); Haptophyta (1 sp.); Cryptophyta (3 sp.); and a residual Protist category (9 sp. including Aureococcus, Blastocystis, Ectocarpus, Fonticula, Hondaea, Nannochloropsis, Symbiodinium, and Thecamonas) representing lineages not assigned to the above named groups. The 43 Apicomplexa species include 6 Plasmodium, 7 Eimeria, 6 Cryptosporidium, and multiple Babesia, Theileria, Hammondia, Toxoplasma, and related species. The 118 Viridiplantae species cover all major angiosperm orders as well as Selaginella (lycophyte), Physcomitrium and Marchantia (bryophytes), and multiple algae (Chara braunii). Note that choanoflagellates — the unicellular sister group to Metazoa and the critical outgroup for dating the origin of animal-type terminal LCR patterns — are absent from the current dataset; their inclusion is a priority for a subsequent revision (see Discussion).

Prokaryote sampling was expanded substantially: 87 bacterial species spanning Proteobacteria (α, β, γ, δ, ε), Firmicutes, Actinobacteria, Bacteroidetes, Spirochaetes, Cyanobacteria, and additional phyla; 21 archaeal species spanning Euryarchaeota, Crenarchaeota, Thaumarchaeota, and Nanoarchaeota. Species were selected to represent phylum-level diversity. Asgard archaea (Heimdallarchaeota, Lokiarchaeota), which represent the closest known prokaryotic relatives of eukaryotes, are not included in the current dataset; their positional LCR properties would be particularly informative for the Scenario A/B question (see Discussion) and constitute a priority expansion.

Where a species had multiple protein isoforms per gene, only the longest isoform was retained by parsing gene identifiers from FASTA headers. Protein FASTA files were downloaded from Ensembl FTP or UniProt FTP servers.

### LCR detection and filtering

LCRs were detected using fLPS 2.0 (Harrison 2017) with parameters identical to Teekas et al. (2024) and the prior 94-species analysis: minimum LCR length = 3 amino acids (-m 3); only SINGLE-type records (single-residue compositional bias) retained in post-processing; purity filter ≥70% (dominant amino acid count / LCR length) applied in post-processing. The pre-compiled macOS/Linux binary was used without recompilation.

The SINGLE-type filter restricts analysis to the most compositionally extreme LCRs (single amino acid dominates ≥70% of the segment). This excludes mixed-composition biases — for example, segments enriched in both G and P typical of disordered linkers — which may have different positional distributions. Sensitivity of the terminal enrichment signal to parameter choice (minimum length ≥6 aa and purity ≥80%; minimum length ≥3 aa and purity ≥60%) and to inclusion of MULTI-type records is reported in Supplementary Table S6.

### Positional binning, enrichment statistics, and multiple testing correction

Positional binning, terminal LCR definition, species-level Fisher's exact test, and protein-level sensitivity analysis were as described previously (one-sided Fisher's exact test; null 10% terminal; terminal = bins 1 or 20 of 20). Phylum-level summaries pool all LCRs from member species before applying Fisher's test.

Domain-level pooled analysis was added for prokaryotes: all LCRs from Bacteria (n=87 species) and from Archaea (n=21 species) were pooled and tested against the 10% null. This approach loses between-species heterogeneity information but provides a powerful aggregate test. The within-domain coefficient of variation of species-level pct_terminal is reported alongside pooled results to characterise heterogeneity (Bacteria CV=50.1%; Archaea CV=44.7%). To confirm the pooled result is robust to within-species clustering, a cluster-robust GEE logistic regression — terminal (0/1) ~ protein-length-quartile + domain, with species as the cluster (binomial family, exchangeable working correlation) — was fitted, complemented by a species-level bootstrap of the domain means; results are reported in the Discussion. A fuller nested random-effects model with an explicit (1|phylum) term remains a follow-up analysis.

To account for multiple comparisons across 42 phylum-level tests, Holm-Bonferroni correction was applied to all phylum-level Fisher's exact p-values. Results are reported as both uncorrected (for comparison with prior work) and corrected. Phyla discussed as significant in the main text have survived Holm-Bonferroni correction unless explicitly noted otherwise.

### Protein-length confound test

The length-stratified analysis was extended to Bacteria, Archaea, and all eukaryotic lineages. Protein-length quartile boundaries were defined globally across all LCRs within each domain or phylum pooled together (not per-species). LCRs were then assigned to quartiles by their host protein's length, and Fisher's exact test for terminal enrichment applied within each quartile. This tests whether terminal enrichment is driven by short proteins mechanically populating terminal bins.

### N-terminal versus C-terminal asymmetry

The asymmetry ratio (pct_nterm / pct_cterm) was calculated per species and summarised per phylum/domain. This ratio is numerically unstable when pct_cterm approaches zero; for phyla with fewer than 50 C-terminal LCRs, the ratio is flagged as indicative only. A PGLS (phylogenetic generalised least-squares) regression of the N/C asymmetry ratio on evolutionary tier (an ordinal scale from green algae through grasses) was fitted under a Brownian-motion model to test whether the progressive N-terminal amplification gradient survives correction for phylogenetic non-independence. The variance-covariance matrix was derived from a time-calibrated backbone phylogeny of the sampled Viridiplantae (one genus representative per genus; deep divergences fixed at established dates, within-tier genera placed as dated polytomies). This backbone is deliberately conservative — it resolves the between-tier structure the gradient concerns but not within-tier topology; a species-level dated tree (timetree.org) would provide additional within-tier resolution.

### Signal peptide stratification (Bacteria)

Bacterial LCRs were stratified by protein class (signal peptide present vs. absent) and tested for N-terminal enrichment (bin 1) in each class. Signal-peptide status was taken from UniProtKB SIGNAL sequence-feature annotations, retrieved per organism by taxon identifier from the UniProt REST API. Proteins were matched to annotations by UniProt accession; proteomes not sourced from UniProt (Ensembl Genomes and NCBI Datasets records, which lack UniProt accessions) could not be matched and are reported as an unannotated class rather than misassigned. N-terminal LCR fractions were compared between the with- and without-signal-peptide classes by Fisher's exact test on the pooled bin-1 counts.

### LLPS propensity of terminal versus internal LCRs

Phase-separation propensity of terminal versus internal LCRs was compared in a panel of model organisms using a composition-based proxy applied to the dominant residue of each SINGLE-type LCR. Per-residue weights span the two established phase-separation regimes: aromatic π–π / cation-π stacking (F/Y/W = 1.0, R = 0.5; Vernon et al. 2018) and prion-like poly-Q/poly-N composition (Q/N = 0.8; the signal scored by PLAAC, Lancaster et al. 2014), with flexible spacers (G/S = 0.5) and all other residues scoring 0. Terminal and internal proxy-score distributions were compared per organism by one-sided Mann-Whitney U test (H₁: terminal > internal). The proxy is the offline substitute for PLAAC; PLAAC scores can be supplied in place of the proxy where available.

### Per-domain amino acid composition of terminal LCRs

For each domain group (Bacteria, Archaea, Viridiplantae, Metazoa, and other eukaryotes), the dominant residue of every terminal and internal SINGLE-type LCR was tallied, and a terminal-to-internal enrichment ratio computed per amino acid (frequency among terminal LCRs divided by frequency among internal LCRs). This complements the pooled compositional analysis (which reports the most abundant residues at termini) by identifying which residues are most positionally biased toward termini within each domain, and tests whether the same residues are enriched universally or whether prokaryote terminal LCRs form a distinct compositional class.

### Amino acid identity, driver analysis, purity gradient, GO enrichment

These analyses were completed for all 42 phyla/groups in the dataset. GO term enrichment analysis was performed for species with Ensembl BioMart annotation and sufficient protein identifier overlap; the same five species as the prior analysis reached significance (Pediculus humanus, Tetranychus urticae, Caenorhabditis brenneri, Lottia gigantea, Strongylocentrotus purpuratus). Systematic GO enrichment across the full 724-species dataset was not completed due to BioMart identifier incompatibilities for non-model organisms. For UniProt-sourced prokaryote proteomes, GO terms are embedded in flat-file annotations and do not require BioMart; GO enrichment for these species is planned as a supplementary extension.

### Code availability

All analysis scripts are available at [repository URL]. The pipeline is fully reproducible from raw proteome downloads.

---

## Results

### Primary finding: terminal LCR enrichment extends to prokaryotes

We report, for the first time, that LCRs are significantly enriched at the protein termini of both bacteria and archaea.

**Bacteria** (87 spp., 2,706 LCRs, 738 terminal, **27.3%**, pooled p~0) and **Archaea** (21 spp., 370 LCRs, 97 terminal, **26.2%**, pooled p~0) both show significant terminal enrichment when pooled. These terminal fractions exceed those of most eukaryotic lineages (typically 14–24%) and exceed the Tetrapoda range of Teekas et al. (2024).

**The length-confound control is decisive.** Stratifying prokaryote LCRs by host protein length quartile (boundaries defined globally within each domain):

| Domain | Q1 (shortest) | Q2 | Q3 | Q4 (longest) |
|---|---|---|---|---|
| Bacteria | 20.6%*** | 29.9%*** | 37.0%*** | 29.3%*** |
| Archaea | 19.1%** | 29.9%** | 35.7%** | 37.9%** |

All eight quartile tests are significant (Fisher's exact; *** p<0.001, ** p<0.01). Enrichment increases monotonically from Q1 to Q3 in both domains — the opposite of what a length-confound artefact would predict (which would inflate terminal fractions disproportionately in short proteins). The Q4 values (bacteria 29.3%, archaea 37.9%) diverge: archaea Q4 continues to increase while bacteria Q4 declines, suggesting potentially different length-dependent dynamics in the two domains that warrant further investigation. Prokaryote terminal enrichment is not an artefact of short median protein lengths.

**Within-domain heterogeneity.** The high CV in bacteria (50.1%) and archaea (44.7%) — versus 13.8% in Viridiplantae across 118 species — reflects that individual prokaryote species are too LCR-poor for reliable individual-species estimates (median ~21 LCRs per bacterial species under our purity filter). The pooled result is robust but should be interpreted as a domain-level aggregate rather than a species-typical value. Species with sufficient LCRs for individual-species testing include the Actinobacteria (Streptomyces coelicolor 125 LCRs, Streptomyces griseus 89 LCRs, Mycobacterium tuberculosis 42 LCRs), where individual enrichment is significant (25–36% terminal, p<0.005). A cluster-robust GEE logistic regression (species as cluster) confirms that terminal enrichment, and its increase across protein-length quartiles, remains highly significant after accounting for species-level clustering (Q4 vs Q1 p<0.001; Bacteria vs Archaea not significantly different, p=0.49).

**Between-domain comparisons.** Pairwise phylum-level tests (after Holm correction) show that bacteria are significantly more terminal-enriched than Insecta (p<0.001), Crustacea (p=0.0015), and Chelicerata (p=0.0006). Bacteria and Viridiplantae are not significantly different from each other in pooled terminal fraction.

### Terminal LCR enrichment across the animal kingdom: 24 phyla

Across all 723 analysed species (724 downloaded; *Acanthaster planci* excluded), 280,687 insect LCRs, 328,474 plant LCRs, and 929,801 LCRs total were analysed. The U-shaped LCR positional profile (bins 1 and 20 elevated relative to internal bins) is visually apparent across species in the bin heatmap (Figure 1) and the per-phylum profile overlay (Figure 3).

After Holm-Bonferroni correction across 42 phylum-level tests, pooled terminal enrichment is significant in 40 of 42 phyla/groups (Table 1; Figure 2). Within Metazoa, 23 of 24 phyla are significant. The single phylum without pooled significance is Acanthocephala (1 species, n=128 LCRs, 16.4% terminal, p=0.07; uncorrected); this is consistent with insufficient statistical power from a single small-proteome species. Twelve metazoan phyla are currently represented by a single species; their enrichment values are indicative data points rather than phylum confirmations.

**Table 1. Metazoan phylum-level terminal LCR enrichment.**

| Phylum | N spp. | Total LCRs | % Terminal | Pooled sig (corrected) | Prior study? |
|---|---|---|---|---|---|
| Porifera | 2 | 1,549 | 23.1 | Yes | NEW |
| Xenacoelomorpha | 1* | 412 | 21.8 | Yes | NEW |
| Ctenophora | 1* | 542 | 20.5 | Yes | ✓ |
| Placozoa | 1* | 234 | 18.8 | Yes | ✓ |
| Cnidaria | 13 | 11,181 | 19.5 | Yes | expanded |
| Tardigrada | 2 | 1,843 | 19.0 | Yes | NEW |
| Rotifera | 1* | 3,887 | 15.9 | Yes | NEW |
| Nematoda | 15 | 12,578 | 19.1 | Yes | expanded |
| Platyhelminthes | 7 | 5,425 | 14.6 | Yes | expanded |
| Annelida | 5 | 6,872 | 16.6 | Yes | expanded |
| Nematomorpha | 1* | 181 | 16.6 | Yes | NEW |
| Nemertea | 1* | 741 | 20.2 | Yes | NEW |
| Brachiopoda | 1* | 1,250 | 21.2 | Yes | ✓ |
| Mollusca | 23 | 24,710 | 18.7 | Yes | expanded |
| Crustacea | 21 | 39,253 | 18.1 | Yes | expanded |
| Chelicerata | 28 | 28,033 | 18.1 | Yes | expanded |
| Myriapoda | 1* | 699 | 14.0 | Yes | ✓ |
| Collembola | 2 | 2,918 | 20.3 | Yes | NEW |
| Insecta | 218 | 280,687 | 16.8 | Yes | major expansion |
| Echinodermata | 6† | 6,356 | 19.6 | Yes | expanded |
| Hemichordata | 1* | 797 | 18.8 | Yes | ✓ |
| Chordata | 1* | 1,121 | 15.8 | Yes | ✓ |
| Priapulida | 1* | 1,273 | 14.0 | Yes | ✓ |
| Acanthocephala | 1* | 128 | 16.4 | No (p=0.07, uncorrected) | NEW |

*Single-species representation: enrichment value is an indicative data point rather than a phylum-level confirmation.
†Echinodermata: 7 species downloaded; *Acanthaster planci* excluded post-hoc (empty fLPS output); 6 analysed.

The within-phylum coefficient of variation (CV) is now available for 29 phyla/groups (Supplementary Table S3). Low CV is observed for well-sampled phyla: Echinodermata CV=4.6%, Porifera CV=6.0%, Tardigrada CV=7.2%, Collembola CV=8.5%, Mollusca CV=11.6%, Crustacea CV=11.6%, Viridiplantae CV=13.8%, Insecta CV=15.6%, Oomycota CV=16.2%. This confirms the signal is a consistent phylum-level property in well-sampled lineages.

**Terminal enrichment is length-independent.** Length-stratified analysis confirms enrichment in at least two quartiles for all 42 phyla/groups with sufficient LCR counts. Pooled across all metazoans: Q1 (shortest proteins) = 13.9%, Q2 = 22.2%, Q3 = 22.7%, Q4 = 20.0%, all significant (p~0). The modest reduction in Q1 relative to Q2–Q4 is inconsistent with a length-confound artefact.

**Singleton-LCR proteins drive the signal.** Proteins carrying a single LCR show significant terminal enrichment in 41 of 42 phyla/groups; multi-LCR proteins show enrichment in 23/42. Singleton-LCR proteins are confirmed as the primary driver across the full phylogenetic breadth of the dataset.

### New phylogenetic anchors: Porifera, Xenacoelomorpha, and Tardigrada

**Porifera (sponges; 2 spp., 23.1%).** Amphimedon queenslandica (22.6%) and Halichondria panicea (23.1%) both show significant terminal enrichment (p<0.001). Sponges diverged from the animal stem ≥600–650 million years ago, have no neurons, and no muscles. Their terminal LCR enrichment establishes this property as pre-neural and ancestrally metazoan.

**Xenacoelomorpha (1 sp., 21.8%).** *Hofstenia miamia* shows 21.8% terminal LCRs (p=2×10⁻⁶). The phylogenetic position of Xenacoelomorpha remains debated; current consensus increasingly places them as sister to Nephrozoa. This is a single species: the result provides an initial data point for the clade but does not confirm phylum-level enrichment. The addition of *Symsagittifera roscoffensis* (transcriptome available) is a priority for the next revision.

**Tardigrada (2 spp., 19.0%).** *Hypsibius exemplaris* (18.5%) and *Paramacrobiotus metropolitanus* (19.5%) are the first tardigrades tested. Tardigrades are known for unusually high proteome disorder; their terminal LCR enrichment within the typical bilaterian range indicates that amplified proteome-wide disorder does not alter the terminal positioning pattern.

### Pan-eukaryotic coverage: all supergroups confirmed

All 16 major eukaryotic lineages tested show significant pooled terminal enrichment after Holm-Bonferroni correction, with one exception: Chlorophyta (p=0.025 uncorrected, does not survive correction; treated as provisional pending additional species). Table 2 lists corrected significance.

**Table 2. Non-metazoan eukaryote terminal LCR enrichment by lineage.**

| Lineage / Supergroup | N spp. | % Terminal | Pooled sig (corrected) |
|---|---|---|---|
| Fungi (Opisthokonta) | 10 | 17.8 | Yes |
| Oomycota (SAR–Stramenopiles) | 26 | 21.2 | Yes |
| Bacillariophyta/diatoms (SAR–Stram.) | 5 | 17.5 | Yes |
| Apicomplexa (SAR–Alveolata) | 43 | 12.4 | Yes |
| Ciliophora (SAR–Alveolata) | 7 | 24.2 | Yes |
| Perkinsozoa (SAR–Alveolata) | 1* | 22.0 | Yes |
| Rhizaria (SAR) | 3 | 18.7 | Yes |
| Euglenozoa (Excavata) | 15 | 18.3 | Yes |
| Metamonada (Excavata) | 5 | 18.6 | Yes |
| Amoebozoa | 9 | 16.8 | Yes |
| Viridiplantae (Archaeplastida) | 118 | 24.5 | Yes |
| Rhodophyta (Archaeplastida) | 3 | 13.2 | Yes |
| Chlorophyta (Archaeplastida) | 2 | 11.2 | Provisional (p=0.025 uncorrected; does not survive Holm-Bonferroni) |
| Haptophyta | 1* | 13.1 | Yes |
| Cryptophyta | 3 | 19.7 | Yes |
| Protist (diverse) | 9 | 17.4 | Yes |

*Single-species representation.

Euglenozoa (15 spp.) and Metamonada (5 spp.) are among the earliest-diverging eukaryotic clades. Their significant enrichment at 18.3% and 18.6% places the origin of terminal LCR bias at or before the last eukaryotic common ancestor (LECA, ~1.5–2 Bya).

**Apicomplexa** (43 spp., 12.4%): The atypically low terminal fraction reflects the high abundance of asparagine-rich LCRs that are predominantly internal in *Plasmodium* and related species (Muralidharan and Goldberg 2013), functioning in immune evasion and not subject to the positional selection acting on terminal LCRs in other lineages. Purity gradient analysis confirms that Apicomplexa terminal LCRs are significantly purer than internal LCRs (Δ=+0.007, p~0), indicating that LCRs occupying terminal positions in Apicomplexa are subject to qualitative constraints absent from the asparagine-dominated internal pool.

**Chlorophyta** (2 spp., 11.2%, p=0.025 uncorrected): The lowest observed enrichment level, and not significant after multiple-testing correction. *Chlamydomonas reinhardtii* (~10.7% terminal) drives the low phylum-level value; *Ostreococcus lucimarinus* shows individually significant N- and C-terminal enrichment. The Chlorophyta result should be treated as provisional until additional species are added.

### The U-shaped positional profile and N/C asymmetry: supergroup-specific signatures

A consistent U-shaped LCR density profile (elevated in bins 1 and 20, depressed in internal bins) is conserved across all domains (Figure 3). However, the relative magnitude of N-terminal versus C-terminal enrichment — the asymmetry ratio (pct_nterm / pct_cterm); note this ratio is indicative only for phyla with fewer than 50 C-terminal LCRs — reveals a striking supergroup-specific pattern not previously reported (Figure 4):

**Viridiplantae — strongly N-terminal dominant.** Across 118 plant species, mean asymmetry ratio ≈ 3.0. Most extreme in grasses (Poaceae): Avena species (3.9–4.5), Oryza species (~3.7–4.1), Secale cereale (4.3), Triticum aestivum (3.8). Eudicots show lower but consistent N-terminal dominance (Arabidopsis thaliana: 1.93; Glycine max: 3.3). Bryophytes are near-balanced: Physcomitrium patens (moss, ~1.0), Marchantia polymorpha (~2.3). This gradient — lowest in mosses, intermediate in eudicots, highest in grasses — is suggestive of progressive amplification of N-terminal LCR enrichment during land plant evolution. However, the pattern does not survive correction for phylogenetic non-independence: although the non-phylogenetic correlation of asymmetry ratio with evolutionary tier is strong (Pearson r=0.59, p<0.001, n=84 genera), a Brownian-motion PGLS on the same data finds no significant gradient (slope=0.30, p=0.16, R²=0.02). The strong raw correlation is therefore substantially inflated by the phylogenetic clustering of the large, closely related samples of grasses and eudicots; the trend in tier means is real, but a per-lineage evolutionary gradient is not robustly supported once shared ancestry is modelled (see Discussion).

**Metazoa and Fungi — C-terminal dominant or balanced.** Mean asymmetry ratio across metazoans ≈ 0.85–1.1. Mollusca and Platyhelminthes are systematically C-dominant (ratios ~0.68 and ~0.55 respectively); the extreme C-terminal bias in Platyhelminthes, particularly *Schistosoma* (ratio 0.47–0.61), is unexplained and warrants further investigation. Fungi are approximately balanced to slightly C-dominant (Saccharomyces cerevisiae: 1.2; Schizosaccharomyces pombe: 0.58; Aspergillus fumigatus: 1.0).

**SAR–Ciliophora — C-terminal dominant.** Stentor coeruleus (ratio 0.64), Paramecium tetraurelia (0.58), Ichthyophthirius multifiliis (0.43), Stylonychia lemnae (0.76).

**Excavata–Metamonada — strongly C-terminal dominant.** Giardia intestinalis (0.54), Spironucleus salmonicida (0.35), Tritrichomonas foetus (0.38).

**Bacteria with sufficient LCRs — N-terminal dominant.** Streptomyces coelicolor (2.5), Streptomyces griseus (3.7), Deinococcus radiodurans (3.3), Geobacter sulfurreducens (3.5), Chloroflexus aurantiacus (2.8), Bacteroides thetaiotaomicron (4.7). The high N-terminal fraction in bacteria parallels the Viridiplantae pattern and contrasts with Metazoa and Fungi.

### Protein-level confirmation

Binomial testing confirms terminal LCR enrichment at the protein level across species with ≥10 LCR-containing proteins. The fraction of proteins carrying ≥1 terminal LCR (observed median ~19%, versus 10% null expectation) is significantly elevated in the overwhelming majority of species (p<0.05, binomial test). Terminal LCR enrichment is a per-protein property, not an artefact of a few highly LCR-rich proteins.

### Terminal LCR purity: a minority signature in selected lineages

Purity gradient analysis (one-sided Mann-Whitney U test) was completed for all 42 phyla/groups. The majority show no significant difference between terminal and internal LCR purity, consistent with terminal enrichment being primarily positional rather than qualitative. Seven phyla show significantly higher purity in terminal LCRs: **Apicomplexa** (Δ=+0.007, p~0), **Rhizaria** (Δ=+0.026, p~0), **Euglenozoa** (Δ=+0.009, p=7.6×10⁻⁵), **Viridiplantae** (Δ=+0.001, p=0.003), **Annelida** (Δ=+0.013, p=0.001), **Chelicerata** (Δ=+0.006, p=0.001), and **Echinodermata** (Δ=+0.014, p=0.002). Rhizaria shows the largest effect (Δ=+0.026), warranting further investigation of amino acid composition at Rhizaria termini. Terminal LCR purity is not significantly elevated in Bacteria (Δ=+0.0002, p=0.81) or Archaea (Δ=−0.004, p=0.33), confirming that prokaryote terminal enrichment is positional rather than qualitative.

### Amino acid composition, driver analysis, and GO enrichment

Terminal amino acid composition (C, E, K most abundant at termini; Q, N not enriched) replicates prior findings across metazoans. The per-domain analysis, which instead ranks residues by their terminal-to-internal enrichment ratio, reveals a different and strikingly consistent pattern: leucine is the most terminally-biased dominant residue in all five domain groups (ratio 2.8× in Bacteria, 5.3× in Archaea, 4.3× in Viridiplantae, 4.8× in Metazoa, 3.5× in other eukaryotes), with methionine the second most terminally-biased residue in Metazoa, Viridiplantae, and other eukaryotes (2.3–2.8×). The prokaryotes then diverge: the next most biased residues are arginine and aspartate in Bacteria (2.1×, 1.7×) and asparagine and glutamine in Archaea (2.9×, 2.2×). The universal leucine bias and the eukaryotic methionine bias implicate hydrophobic signal/anchor sequences and initiator-methionine retention — hallmarks of N-terminal signal-peptide biology and N-terminal processing — rather than the acidic/basic degron residues, as the most consistent cross-domain compositional signature of terminal LCRs (Supplementary Table S7). Driver analysis (singleton-LCR proteins drive enrichment across all 42 phyla/groups) is not qualitatively changed by the expanded dataset. GO enrichment was significant in the same five species as previously; the absence of functional concentration in terminal LCRs across well-annotated metazoan species remains consistent with a mechanism acting on protein architecture regardless of gene function.

---

## Discussion

### The prokaryote positional finding is the primary novel result

No prior study has tested whether LCRs are positionally enriched at protein termini in bacteria or archaea. Prior large-scale prokaryote LCR work — Ntountoumi et al. (2019) in >1,500 proteomes, Persi et al. (2023) in genome-level evolutionary analysis, and Saravanan et al. (2025) in pan-genome core vs. accessory comparisons — characterised LCR prevalence, amino acid composition, conservation, and gene-level evolutionary dynamics, but none applied positional binning within the protein sequence. Coletta et al. (2010) first showed positional enrichment at protein termini but in a single yeast proteome. The prokaryote positional question has therefore been open since 2010.

The pooled bacterial terminal fraction (27.3%) and archaeal fraction (26.2%) both exceed the typical eukaryotic range (14–24%) and are significant at p~0. The length-stratified control is decisive: enrichment increases from Q1 to Q3 in both domains (the opposite of a length artefact), and all eight quartile tests are significant. This cannot be dismissed as a sampling or methodological artefact. High within-domain CV (bacteria 50.1%; archaea 44.7%) indicates heterogeneity across species; the pooled result reflects a domain-level aggregate driven by LCR-rich taxa (primarily Actinobacteria) rather than a uniform property of all bacteria. Nonetheless, the enrichment in individually testable species (Streptomyces, Deinococcus, Geobacter) is consistent with the pooled result.

### Two equally viable mechanistic scenarios

Two interpretations of the prokaryote finding are consistent with the data. We present them as equally supported; the current dataset does not resolve between them.

**Scenario A (ancestral origin):** Terminal LCR enrichment reflects a mechanism present in the Last Universal Common Ancestor (LUCA, ~3.5–4 Bya). Under this scenario, the mechanism must operate in all cellular organisms independently of eukaryote-specific features. Translation kinetics — ribosome pausing near start and stop codons, generating compositionally simple sequences at protein termini (Irastortza-Olaziregi and Amster-Choder 2021; Pechmann and Frydman 2013) — is a candidate universal mechanism consistent with ribosome biology across all domains.

**Scenario B (convergent):** Prokaryote terminal LCR enrichment arose independently through lineage-specific mechanisms distinct from eukaryotic degron biology. Candidate prokaryote-specific mechanisms include N-terminal signal peptides and export sequences (abundant at N-termini of secreted proteins), C-terminal ssrA degradation tags (tmRNA-mediated; absent in eukaryotes), and general N-terminal processing. Note that N-formylmethionine (fMet) as an N-degron (Piatkov et al. 2015) is a sometimes-cited candidate, but fMet is co-translationally cleaved from the majority of bacterial proteins by methionine aminopeptidase; the Piatkov et al. pathway specifically targets Met-retained proteins, and the ClpS/ClpAP N-end rule pathway in bacteria recognises hydrophobic residues (L, F, W, Y), not the E, K, C residues typically enriched in eukaryotic terminal LCRs. The fMet-degron hypothesis therefore requires direct testing before it can be considered a plausible mechanism.

The N/C asymmetry data provides partial discrimination: bacteria (when LCR-rich) are predominantly N-terminal enriched, as are Viridiplantae, whereas most Metazoa and Fungi are C-terminal enriched or balanced. If prokaryote terminal enrichment arose through the same mechanism as metazoan C-degron-driven enrichment, a similar C-terminal bias would be expected. The N-terminal dominance in bacteria is more consistent with N-terminal processing biology (signal peptide cleavage) or with Scenario B.

Signal peptide stratification provides direct evidence on this point. Among UniProt-annotated bacterial proteins, N-terminal LCRs (bin 1) are substantially more frequent in signal-peptide-bearing proteins than in those without (34.1% vs. 14.6%; odds ratio 3.03, Fisher's exact p = 7.1 × 10⁻⁵; 85 vs. 425 LCRs across 14 species with matchable annotations). This confirms that N-terminal signal sequences are a genuine contributor to bacterial N-terminal LCR enrichment — consistent with Scenario B. However, signal peptides do not account for the pattern in full: non-secreted proteins still show N-terminal enrichment (14.6%) well above the 5% single-bin null, so a signal-peptide-independent mechanism also operates. The result therefore favours a mixed interpretation — signal-peptide biology amplifies an underlying, more general N-terminal signal — rather than a purely secretory explanation. This analysis is limited to the UniProt-sourced subset (Ensembl/NCBI-sourced proteomes lack matchable accessions), and warrants extension as annotation coverage broadens.

### Critical missing taxa: choanoflagellates and Asgard archaea

Two groups are conspicuously absent from the current dataset and are the highest priority for the next sampling revision.

**Choanoflagellates** (*Monosiga brevicollis*, *Salpingoeca rosetta*) are the unicellular sister group to all Metazoa and have sequenced UniProt reference proteomes. Their N/C asymmetry direction would directly answer whether the C-terminal dominance seen in animals predates multicellularity or arose at the animal stem. If choanoflagellates are C-terminal dominant like most animals, the switch from prokaryote-like N-terminal dominance to C-terminal dominance occurred before animal multicellularity. If N-terminal dominant, the switch occurred in parallel with or after the origin of animals.

**Asgard archaea** (Heimdallarchaeota, Lokiarchaeota, Thorarchaeota) are the closest known prokaryotic relatives of eukaryotes under current phylogenies and are available in NCBI. Their positional LCR properties would directly test whether the eukaryotic terminal enrichment pattern is inherited from the archaeal ancestor (Scenario A expectation: Asgard archaea should match eukaryote direction) or arose de novo in eukaryotes (Scenario B expectation: Asgard archaea should match other bacteria/archaea).

### The Viridiplantae N/C asymmetry gradient and a mechanistic proposal

The 118-species plant dataset confirms and substantially extends the N-terminal dominance first noted in Arabidopsis. The trend in tier means — near-balanced in bryophytes (Physcomitrium: ~1.0), intermediate in eudicots (Arabidopsis: 1.93), and extreme in grasses (3.5–4.5) — is visually consistent across angiosperm orders. This apparent gradient must, however, be interpreted cautiously: a Brownian-motion PGLS on evolutionary tier does not recover a significant slope (0.30, p=0.16) despite the strong non-phylogenetic correlation (r=0.59), indicating that the pattern is driven largely by the many closely related grass and eudicot genera rather than by an independent per-lineage trend. We therefore present progressive N-terminal amplification as a suggestive, unconfirmed hypothesis. Because our backbone tree models within-tier genera as polytomies — a conservative structure that maximally down-weights within-tier variation — a species-level dated phylogeny may recover additional signal; resolving this is a priority for the next revision.

A mechanistic basis is plausible. Land plants have a well-characterised N-degron pathway involving PRT1 and PRT6 E3 ligases and Plant Cysteine Oxidase (PCO) enzymes that link N-terminal cysteine oxidation to oxygen-sensing and proteasomal degradation of ERF-VII transcription factors (Gibbs et al. 2014; Holdsworth et al. 2020). This pathway is absent or simplified in animals and has become progressively more elaborate during land plant evolution, broadly in step with the higher N-terminal asymmetry of angiosperms relative to bryophytes — though, as noted above, the fine-grained mosses-to-grasses gradient is not statistically robust to phylogenetic correction. If PRT6/PCO-pathway substrates were enriched for N-terminal LCRs, this would drive N-terminal LCR accumulation in plants; however, the per-domain composition analysis (Supplementary Table S7) shows that Viridiplantae termini are dominated by leucine and methionine rather than the cysteine/arginine residues these N-degron pathways recognise, so any PRT6/PCO contribution is not the primary compositional signature of plant terminal LCRs.

### Mechanistic framework: layered model

The original mechanistic proposal for terminal LCR enrichment invoked C-degron and N-degron pathways of the eukaryotic ubiquitin-proteasome system (Koren et al. 2018; Varshavsky 2019). This explanation remains valid for the eukaryotic data, supported by the abundance of C, E, and K — known degron residues — among metazoan terminal LCRs and by the absence of GO functional enrichment (consistent with architecture-wide selection irrespective of gene function).

The per-domain composition analysis, however, adds a second, more universal signal that the degron account alone does not capture. When residues are ranked by terminal-to-internal enrichment rather than raw abundance, leucine is the single most terminally-biased residue in every domain, and methionine is second across all three eukaryotic groups. Both point to N-terminal biology that is shared by all cellular life rather than to eukaryote-specific degrons: leucine-rich hydrophobic runs are the defining feature of N-terminal signal peptides, signal anchors, and transmembrane segments, all of which concentrate near protein termini; and a terminal methionine bias is the expected signature of initiator-methionine retention and N-terminal processing. This compositional evidence dovetails directly with the signal-peptide stratification result above — in which signal-peptide-bearing bacterial proteins carry a threefold-higher N-terminal LCR fraction — and locates the shared, ancestral layer of the model in N-terminal signal-sequence and processing biology, upon which the lineage-specific degron elaborations (C-degron in animals, N-degron/PRT6–PCO in plants) are superimposed.

The prokaryote finding requires broadening the framework. Three mechanisms may operate universally or semi-universally:

1. **Ribosome pausing at start and stop codons** — characterised in all cellular life (Irastortza-Olaziregi and Amster-Choder 2021) — would predict symmetric N- and C-terminal enrichment. The observed C-dominance in most eukaryotes argues this is a contributing but not sole mechanism.

2. **N-terminal processing** — signal peptide cleavage, N-terminal methionine excision — is universal and preferentially generates N-terminal sequence constraints, consistent with N-terminal dominance in bacteria and plants.

3. **C-terminal disordered tails** — functional in eukaryotic degron recognition (Koren et al. 2018), tubulin code biology (Janke and Magiera 2020), and RNA-binding — are the most eukaryote-specific mechanism and explain C-terminal dominance in animals and fungi.

We propose a layered model: a universal layer (ribosome kinetics, N-terminal processing) generates baseline terminal LCR enrichment in all cellular life, superimposed by lineage-specific amplification — C-degron biology in animals, N-degron pathway elaboration in plants, and N-terminal processing biology in bacteria. The mechanistic separation of the N-terminal and C-terminal enrichment signals across domains is the key predictive framework, and the per-domain amino acid composition analysis provides direct support: the universal leucine bias and eukaryotic methionine bias identify N-terminal signal-sequence and processing biology as the shared baseline layer, distinct from the acidic/basic degron residues that mark the lineage-specific eukaryotic amplifications.

### Terminal LCRs and phase separation: an open connection

LCRs are strongly associated with liquid-liquid phase separation (LLPS) and biomolecular condensate formation (Boija et al. 2018; Alberti et al. 2019). Terminal disordered tails are known nucleators of condensate assembly in specific protein families (e.g., RNA-binding proteins, transcription factors). We tested whether terminal LCRs have systematically elevated LLPS propensity relative to internal LCRs, using a composition-based proxy (aromatic and prion-like poly-Q/poly-N residues; a substitute for PLAAC, Lancaster et al. 2014) across a panel of model organisms. The result is organism-specific rather than universal: of seven organisms, only *Arabidopsis thaliana* shows significantly higher terminal LLPS propensity (Mann-Whitney U, p = 0.033), while in *Drosophila melanogaster* and *Saccharomyces cerevisiae* internal LCRs are the more phase-prone class, and the remainder show no significant difference. Terminal LCRs therefore do not carry a domain-general elevation in phase-separation propensity; where an association exists it is lineage-specific. Confirmatory PLAAC scoring and broader taxonomic sampling would sharpen this picture, but the positional enrichment reported here is not, in itself, driven by a uniform phase-separation signature.

### The Platyhelminthes, Chlorophyta, and Acanthocephala outliers

**Platyhelminthes** (14.6% terminal; C-dominant) shows the lowest metazoan terminal enrichment. The expansion to 7 species confirms this is a genuine phylum-level property. The systematic C-terminal bias, particularly in *Schistosoma* (ratio 0.47–0.61), contrasts with all other lophotrochozoans and is mechanistically unexplained.

**Chlorophyta** (11.2%; p=0.025 uncorrected; does not survive Holm-Bonferroni): *Chlamydomonas reinhardtii* drives the low value; *Ostreococcus lucimarinus* shows individually significant enrichment. The Chlorophyta conclusion is provisional pending additional species.

**Acanthocephala** (p=0.07) reflects power limitation from a single species with 128 LCRs. The observed effect size (16.4%) is within the bilaterian range.

### Limitations

1. **Prokaryote statistical model:** The pooled Fisher's test does not correct for between-species clustering, but a cluster-robust GEE logistic regression (binomial family, exchangeable working correlation, species as the cluster) confirms that the result is not an artefact of treating clustered LCRs as independent: the baseline terminal probability exceeds the 10% null (intercept 95% CI [−2.07, −1.56] in log-odds, excluding the null value of −2.20), the increase across protein-length quartiles is strongly significant (Q4 vs Q1 log-odds +1.13, 95% CI [0.90, 1.37], p<0.001), and Bacteria and Archaea do not differ (p=0.49). A species-level bootstrap (5,000 resamples) gives concordant domain means (Bacteria 27.0% [24.3–29.7], Archaea 26.6% [22.4–31.2]). The 50.1% bacterial CV still indicates the pooled point estimate is influenced by a subset of LCR-rich taxa, but the enrichment itself is robust to species-level clustering.

2. **fLPS parameter sensitivity:** A three-parameter-combination sensitivity analysis (Supplementary Table S6) confirms that the SINGLE-type terminal enrichment is robust to the LCR length threshold and purity cutoff — prokaryote terminal enrichment persists under stringent (length ≥6, purity ≥80%: Bacteria 34.1%, Archaea 29.0%) and relaxed (length ≥3, purity ≥60%: Bacteria 24.0%, Archaea 22.4%) parameters alike. Terminal enrichment does not extend to MULTI-type LCRs (Bacteria 9.1%, Archaea 11.3%; both non-significant), indicating the signal is specific to the most compositionally extreme, single-residue-dominated regions.

3. **Multiple testing:** Holm-Bonferroni correction has been applied to all 42 phylum-level tests. Chlorophyta (p=0.025) does not survive correction and is treated as provisional.

4. **Single-species phyla:** Twelve metazoan phyla and several protist lineages are represented by one species. These provide initial data points but not phylum-level confirmations.

5. **Missing key outgroups:** Choanoflagellates and Asgard archaea are absent. These are the two highest-priority additions for the next revision.

6. **GO enrichment gap:** GO analysis covers only 5 of 723 species (BioMart compatibility); prokaryote GO enrichment via UniProt flat-file annotations is planned.

7. **SINGLE-type filter:** Mixed-composition LCRs are excluded; their positional distribution is unknown.

8. **Asymmetry ratio instability:** The ratio is numerically unstable for phyla with few C-terminal LCRs; log-odds asymmetry scoring is planned.

9. **LUCA/LECA temporal inference:** The data are consistent with a LUCA-level origin but cannot distinguish this from convergent evolution. Scenarios A and B remain equally viable.

---

## Conclusions

We report the first evidence that terminal LCR enrichment extends to prokaryotes: pooled bacteria (27.3%, 87 species) and archaea (26.2%, 21 species) show significant terminal enrichment that is not explained by protein length and increases across protein-length quartiles in the direction opposite to a length-confound artefact. This is the primary novel finding; all prior positional work was confined to eukaryotes. Signal peptide stratification shows that N-terminal signal sequences contribute to — but do not fully explain — the bacterial N-terminal signal (signal-peptide proteins 34.1% vs. 14.6% N-terminal LCRs, p = 7.1 × 10⁻⁵, yet non-secreted proteins remain enriched above null). Whether the residual, signal-peptide-independent enrichment reflects a mechanism present since LUCA (~3.5–4 Bya) or convergent evolution in prokaryotes cannot yet be determined; both scenarios remain viable pending per-domain amino acid composition analysis and Asgard archaea sampling.

Among eukaryotes, 40 of 42 phyla/groups survive Holm-Bonferroni correction, establishing terminal LCR bias as a LECA-level property (≥1.5 Bya); Chlorophyta (p=0.025 uncorrected) is provisional. Among metazoans, 23 of 24 phyla are confirmed including the basally-branching Porifera and Tardigrada; Xenacoelomorpha is a single-species data point pending expansion. Singleton-LCR proteins drive the signal in 41 of 42 phyla/groups.

The N/C asymmetry is a supergroup-specific mechanistic signature: Viridiplantae and LCR-rich bacteria are N-terminal dominant; most Metazoa and Fungi are C-terminal dominant or balanced; Ciliophora, Metamonada, and Platyhelminthes are C-terminal dominant. This pattern is incompatible with a single universal mechanism and discriminates between candidate explanations. Together, the findings support a layered mechanistic model in which universal translational kinetics and N-terminal processing generate baseline terminal LCR enrichment in all cellular life, amplified by lineage-specific degron and processing biology.

---

## Acknowledgements

The author thanks the Ensembl Metazoa, Ensembl Plants, and UniProt teams for providing freely downloadable proteomes, and Paul Harrison for making fLPS 2.0 freely available.

---

## References

Alberti S, Gladfelter A, Mittag T (2019). Considerations and challenges in studying liquid-liquid phase separation and biomolecular condensates. *Cell* 176:419–434.

Berriman M, et al. (2009). The genome of the blood fluke *Schistosoma mansoni*. *Nature* 460:352–358.

Boija A, et al. (2018). Transcription factors activate genes through the phase-separation capacity of their activation domains. *Cell* 175:1842–1855.

Coletta A, Pinney JW, Solís DY, Marsh J, Pettifer SR, Attwood TK (2010). Low-complexity regions within protein sequences have position-dependent roles. *BMC Systems Biology* 4:43.

Eme L, Sharpe SC, Brown MW, Roger AJ (2014). On the age of eukaryotes: evaluating evidence from fossils and molecular clocks. *Cold Spring Harbor Perspectives in Biology* 6:a016139.

Gibbs DJ, et al. (2014). Homeostatic response to hypoxia is regulated by the N-end rule pathway in plants. *Nature* 479:415–418.

Harrison PM (2017). fLPS: Fast discovery of compositional biases for the protein universe. *BMC Bioinformatics* 18:476.

Holdsworth MJ, Vicente J, Sharma G, Abbas M, Estavillo GM (2020). The plant N-degron pathways of ubiquitin-mediated proteolysis. *Journal of Integrative Plant Biology* 62:70–89.

Irastortza-Olaziregi M, Amster-Choder O (2021). Coupled transcription-translation in prokaryotes: an old couple with new surprises. *Frontiers in Microbiology* 11:619430.

Janke C, Magiera MM (2020). The tubulin code and its role in controlling microtubule properties and functions. *Nature Reviews Molecular Cell Biology* 21:307–326.

Koren I, et al. (2018). The eukaryotic proteome is shaped by E3 ubiquitin ligases targeting C-terminal degrons. *Cell* 173:1622–1635.

Lancaster AK, et al. (2014). PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition. *Bioinformatics* 30:2–3.

Marcotte EM, et al. (1999). A census of protein repeats. *Journal of Molecular Biology* 293:151–160.

Muralidharan V, Goldberg DE (2013). Asparagine repeats in *Plasmodium falciparum* proteins: Good for nothing? *PLoS Pathogens* 9:e1003488.

Ntountoumi C, et al. (2019). Low complexity regions in the proteins of prokaryotes perform important functional roles and are highly conserved. *Nucleic Acids Research* 47:9998–10009.

Pechmann S, Frydman J (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology* 20:237–243.

Persi E, Wolf YI, Karamycheva S, Makarova KS, Koonin EV (2023). Compensatory relationship between low-complexity regions and gene paralogy in the evolution of prokaryotes. *Proceedings of the National Academy of Sciences* 120:e2300154120.

Piatkov KI, Oh J-H, Liu Y, Bhatt DL, Varshavsky A (2015). Formyl-methionine as a degradation signal at the N-termini of bacterial proteins. *Microbial Cell* 2:376–393. [PMC4745127 — verify full author list before submission]

Romero P, et al. (2001). Sequence complexity of disordered protein. *Proteins* 42:38–48.

Saravanan S, Kravetz AN, Battistuzzi FU (2025). Higher frequency of prokaryotic low complexity regions in core and orthologous genes. *Frontiers in Bioinformatics* 5:1673480.

Shin Y, Brangwynne CP (2017). Liquid phase condensation in cell physiology and disease. *Science* 357:eaaf4382.

Teekas L, Sharma S, Vijay N (2024). Terminal regions of a protein are a hotspot for low complexity regions and selection. *Open Biology* 14:230439.

UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research* 51:D523–D531.

van der Lee R, et al. (2014). Classification of intrinsically disordered regions and proteins. *Chemical Reviews* 114:6589–6631.

Varshavsky A (2019). N-degron and C-degron pathways of protein degradation. *Proceedings of the National Academy of Sciences* 116:358–366.

Virtanen P, et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17:261–272.

Wootton JC, Federhen S (1996). Analysis of compositionally biased regions in sequence databases. *Methods in Enzymology* 266:554–571.

---

## Figure Legends

**Figure 1. LCR positional distribution across 723 proteomes.** Heatmap showing the fraction of LCRs in each of 20 equally spaced positional bins (bin 1 = N-terminal 5%; bin 20 = C-terminal 5%) for each analysed species. Rows ordered by domain and phylum (Bacteria → Archaea → non-metazoan eukaryotes → Metazoa). Colour scale: fraction of a species' LCRs in that bin (0–0.15). Dashed blue lines mark the terminal bins.

**Figure 2. Terminal LCR enrichment by phylum.** Bar chart showing pooled % terminal LCRs per phylum/group, ordered phylogenetically. Black dots show individual-species values. Red dotted line: 10% null expectation. Grey shaded band: Tetrapoda range from Teekas et al. (2024; 15–25%). All groups except Acanthocephala (1 sp., insufficient power) and Chlorophyta (provisional) significantly exceed the null after Holm-Bonferroni correction. Bacteria (27.3%) and Archaea (26.2%) show the highest pooled terminal fractions.

**Figure 3. Conserved U-shaped LCR positional profile across all phyla.** Mean fraction of LCRs per positional bin (1–20) per phylum/domain, averaged across member species. Dashed grey vertical lines: terminal bins 1 and 20. Dotted horizontal line: 5% uniform null. All phyla show elevated LCR density at both termini and a depressed internal plateau.

**Figure 4. Supergroup-specific N/C terminal asymmetry.** Per-species asymmetry ratio (pct_nterm / pct_cterm; log-scaled for display; indicative only for phyla with <50 C-terminal LCRs), grouped by phylum/domain. Viridiplantae (118 spp.) are strongly N-terminal dominant (median ratio ~3.0; grasses 3.5–4.5); bryophytes near-balanced (~1.0). Metazoa and Fungi broadly balanced to C-dominant. Ciliophora, Metamonada, and Platyhelminthes are systematically C-dominant. Bacteria with sufficient LCRs are N-terminal dominant.

**Figure 5. Protein-length stratified terminal LCR enrichment in prokaryotes.** Grouped bar chart showing % terminal LCRs in four protein-length quartiles (Q1 shortest; Q4 longest; boundaries defined globally within each domain) for Bacteria (blue) and Archaea (red). Error bars: 95% CI from Fisher's exact test. The Q1→Q3 increase in bacteria (20.6% → 37.0%) rules out a length-confound artefact.

---

## Supplementary Tables

**Supplementary Table S1.** Full species list (724 entries): species name, phylum, data source (Ensembl / UniProt), proteome ID, protein count, LCR count, pct_terminal.

**Supplementary Table S2.** Driver analysis results for all 42 phyla/groups: pct_terminal for singleton-LCR proteins, multi-LCR proteins; odds ratios and p-values for each class.

**Supplementary Table S3.** Within-phylum coefficient of variation of pct_terminal: 29 phyla/groups with n ≥ 2 species. Includes phylum, n_species, mean_pct_terminal, std_pct_terminal, CV.

**Supplementary Table S4.** Length-stratified analysis for all 42 phyla/groups: pct_terminal by protein-length quartile; Fisher's exact p-value and significance per quartile.

**Supplementary Table S5.** Purity gradient analysis for all 42 phyla/groups: mean purity of terminal vs. internal LCRs, Δ purity, Mann-Whitney U p-value, significance flag.

**Supplementary Table S6.** fLPS parameter sensitivity analysis: pct_terminal and pooled significance for all 42 phyla/groups under three parameter combinations — (current: length ≥3, purity ≥70%), (stringent: length ≥6, purity ≥80%), (relaxed: length ≥3, purity ≥60%) — plus a MULTI-type LCR positional analysis. SINGLE-type terminal enrichment is significant across all parameter settings; MULTI-type LCRs show no terminal enrichment.

**Supplementary Table S7.** Per-domain amino acid composition of terminal LCRs: dominant residue enrichment ratio (terminal / internal) per amino acid, for Bacteria, Archaea, Viridiplantae, Metazoa, and other eukaryotes separately. Leucine is the most terminally-enriched residue in every domain (2.8–5.3×); methionine is second across the eukaryotic groups.

**Supplementary Table S8.** Signal peptide stratification of bacterial N-terminal LCRs: per-species N-terminal (bin 1) LCR count and fraction for proteins with a signal peptide, without a signal peptide, and unannotated, based on UniProtKB SIGNAL feature annotations. Pooled across species with matchable UniProt annotations, N-terminal LCRs are 3-fold enriched in signal-peptide-bearing proteins (34.1% of 85 LCRs) relative to non-secreted proteins (14.6% of 425 LCRs; Fisher's exact odds ratio 3.03, p = 7.1 × 10⁻⁵).

**Supplementary Table S9.** LLPS propensity (composition-based proxy) of terminal versus internal LCRs per model organism: n terminal/internal LCRs, median proxy score for each, percent aromatic (F/Y/W), and one-sided Mann-Whitney U p-value (terminal > internal). Of seven organisms, only *Arabidopsis thaliana* reaches significance (p = 0.033).

**Supplementary Figures**

**Supplementary Figure 1.** Per-species terminal LCR % distributions for each phylum, violin/strip chart format.

**Supplementary Figure 2.** Full length-stratified enrichment results for all 42 phyla/groups (4 quartiles × 42 groups heatmap).

**Supplementary Figure 3.** Asymmetry ratio distributions for Viridiplantae only, coloured by plant order. Shows progressive grass amplification relative to bryophytes and basal eudicots.

**Supplementary Figure 4.** Within-phylum CV versus number of species, all 29 phyla with n ≥ 2.

**Supplementary Figure 5.** Amino acid composition of terminal versus internal LCRs: enrichment ratio (terminal/internal) per amino acid, pooled across all metazoans. C, E, K enriched; Q, N depleted at termini.

**Supplementary Figure 6.** Purity gradient: distributions of purity for terminal (bins 1+20) and internal (bins 2–19) LCRs for the seven significant phyla (Apicomplexa, Rhizaria, Euglenozoa, Viridiplantae, Annelida, Chelicerata, Echinodermata). Violin plots with Mann-Whitney U p-values.

**Supplementary Figure 7.** Phylogenetic tree (TimeTree-derived backbone) annotated with pooled pct_terminal per phylum/group. Branch colours indicate supergroup. Node annotations mark LUCA (~3.8 Bya), LECA (~1.8 Bya), and stem Metazoa (~700 Mya).

**Supplementary Figure 8.** PGLS regression of N/C asymmetry ratio on evolutionary tier across 84 Viridiplantae genus representatives, under a Brownian-motion model on a time-calibrated backbone phylogeny. Contrasts the strong non-phylogenetic correlation (OLS/Pearson r=0.59, p<0.001) with the non-significant phylogenetically-corrected slope (PGLS slope=0.30, p=0.16, R²=0.02).
