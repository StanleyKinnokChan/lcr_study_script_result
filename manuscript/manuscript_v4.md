# Terminal Low-Complexity Regions Are Enriched Across All Domains of Life: Evidence from 724 Proteomes Spanning 42 Phyla

**Stanley Kinnok Chan**¹

¹ *[Institution, City, Country]*

*Correspondence: stanleykinnok.chan@gmail.com*

---

## Abstract

Low-complexity regions (LCRs) — protein segments dominated by one or a few amino acid types — are concentrated at the N- and C-termini of proteins across Tetrapoda (Teekas et al. 2024, *Open Biology*) and all major invertebrate phyla (Chan, preceding study). Whether this terminal enrichment extends beyond Metazoa to all eukaryotes, and whether it is present in prokaryotes, remains untested at scale. Here we apply the same computational framework (fLPS 2.0, 20-bin positional analysis, purity ≥70% SINGLE-type filter) to 724 proteomes spanning 42 phyla and all three domains of life, encompassing 929,801 LCR records. Among metazoans, terminal LCR enrichment is significant in 24 phyla including newly tested Porifera (23.1%), Xenacoelomorpha (21.8%), Tardigrada (19.0%), Rotifera (15.9%), Nemertea (20.2%), and Collembola (20.3%); only Acanthocephala (1 species; 16.4%, p=0.07) falls short of significance, likely due to limited power. Among non-metazoan eukaryotes, all 16 major eukaryotic lineages tested show significant pooled enrichment, including all SAR supergroup members (Stramenopiles: Oomycota 21.2%, Bacillariophyta 17.5%; Alveolata: Apicomplexa 12.4%, Ciliophora 24.2%; Rhizaria 18.7%), Excavata (Euglenozoa 18.3%, Metamonada 18.6%), Amoebozoa (16.8%), and all tested Archaeplastida (Viridiplantae 24.5%, Rhodophyta 13.2%, Chlorophyta 11.2%). Critically, pooled analysis of 87 bacterial proteomes (27.3%, p~0) and 21 archaeal proteomes (26.2%, p~0) reveals significant terminal LCR enrichment across both prokaryotic domains — enrichment persisting across all four protein-length quartiles, ruling out a length-confound artefact. The N/C terminal asymmetry shows a supergroup-specific pattern: Metazoa and Fungi are predominantly C-terminal enriched; Viridiplantae are strongly N-terminal enriched (mean pct_nterm ~17% vs. pct_cterm ~6%; grasses extreme at ratios 3.5–4.5); bacteria with sufficient LCRs are N-terminal enriched; while Ciliophora, Metamonada, and Platyhelminthes are C-terminal enriched. These findings establish terminal LCR enrichment as a near-universal feature of protein architecture across all domains of life, constraining mechanistic explanations to processes operating in all cellular organisms — including ribosomal translation kinetics at start/stop codons and universal N-terminal processing pathways — rather than eukaryote-specific degron systems alone. Within-domain variability for bacteria (CV=50.1%) is high, and the pooled prokaryote figure reflects a small number of LCR-rich species; the result is interpreted as a domain-level aggregate rather than a species-typical value.

**Keywords:** low-complexity regions, protein termini, compositional bias, eukaryotic evolution, prokaryote protein architecture, intrinsically disordered regions, fLPS2, pan-genomics

---

## Introduction

Proteins are not compositionally uniform along their length. Low-complexity regions (LCRs) — segments dominated by one or a small number of amino acid types — are distributed non-randomly within protein sequences, clustering in functionally important contexts such as disordered linkers, prion-like domains, and polyamino acid tracts (Wootton and Federhen 1996; Marcotte et al. 1999). The biological significance of LCRs has grown considerably with the recognition that many intrinsically disordered regions (IDRs) are LCR-containing (Romero et al. 2001; van der Lee et al. 2014), and that phase separation is often driven by low-complexity IDR sequences (Boija et al. 2018; Shin and Brangwynne 2017).

The positional distribution of LCRs within proteins has received less attention than their prevalence or functional roles. Teekas and colleagues (2024) reported that across all major Tetrapoda clades, LCRs are significantly enriched in the terminal 5% of protein sequences (the first and last bins of a 20-bin positional map), with 15–25% of all LCRs in terminal positions despite these bins representing only 10% of positional space. This enrichment was observed consistently across birds, mammals, reptiles, and amphibians. A companion study (Chan, in preparation) extended this finding to 61 metazoan invertebrate species across 16 phyla and to 33 non-metazoan eukaryote and prokaryote species (94 total), establishing terminal LCR enrichment as a pan-metazoan and likely pan-eukaryotic property approximately 1 billion years old.

LCRs in prokaryotes are functional and evolutionarily conserved (Ntountoumi et al. 2019), and recent work suggests that LCR content and gene paralogy are compensatorily related in prokaryotic evolution (Luo et al. 2023) — providing context for interpreting variation in bacterial LCR abundance across lineages. Several open questions remain. First, the metazoan survey covered only a fraction of animal phylogenetic diversity; many phyla now have sequenced genomes that were absent. Porifera (sponges), Xenacoelomorpha, Tardigrada, Rotifera, Nemertea, Nematomorpha, and Collembola were not represented. Second, the eukaryotic outgroup analysis covered only four protist lineages, leaving most eukaryotic supergroups untested; in particular, the SAR supergroup (Stramenopiles, Alveolata, Rhizaria) and Excavata had minimal representation. Third, the prokaryote analysis was severely underpowered: with median 21 LCRs per bacterial species, individual-species Fisher's tests were not interpretable. Whether bacteria and archaea show terminal LCR enrichment at all remained unresolved.

Here we address all three gaps by scaling the analysis to 724 proteomes spanning 42 phyla and all three domains of life, encompassing 929,801 LCRs. We ask: (1) Does terminal LCR enrichment hold across all newly tested metazoan phyla, including the most basally branching animals? (2) Does the signal span all eukaryotic supergroups? (3) With 87 bacterial and 21 archaeal proteomes, is terminal LCR enrichment detectable in prokaryotes? (4) Is the N/C asymmetry pattern consistent across lineages, or do lineage-specific reversals exist?

---

## Methods

### Species selection and proteome acquisition

A total of 724 proteomes were assembled from three sources: (1) Ensembl Metazoa release 63 (metazoan invertebrate species); (2) Ensembl Plants and Ensembl release 110 (Viridiplantae, algae, and additional fungi); (3) UniProt reference proteomes 2024 (additional eukaryotes and all prokaryotes). Species were selected to maximise taxonomic breadth within each major lineage. The full list is provided in Supplementary Table S1.

Metazoan sampling represents 24 phyla: Acanthocephala (1 sp.), Annelida (5), Brachiopoda (1), Chelicerata (28), Chordata/Cephalochordata (1), Cnidaria (13), Collembola (2), Crustacea (21), Ctenophora (1), Echinodermata (7, of which Acanthaster planci was excluded post-hoc due to empty fLPS output), Hemichordata (1), Insecta (218), Mollusca (23), Myriapoda (1), Nematoda (15), Nematomorpha (1), Nemertea (1), Placozoa (1), Platyhelminthes (7), Porifera (2), Priapulida (1), Rotifera (1), Tardigrada (2), and Xenacoelomorpha (1). The newly added phyla — Porifera, Xenacoelomorpha, Tardigrada, Rotifera, Collembola, Nematomorpha, and Nemertea — extend metazoan coverage to include the most basally branching animals and several ecdysozoan and lophotrochozoan phyla absent from prior analyses.

Non-metazoan eukaryote sampling spans all recognised eukaryotic supergroups: Opisthokonta (Fungi: 10 sp.); SAR–Stramenopiles (Oomycota: 26 sp., Bacillariophyta: 5 sp.); SAR–Alveolata (Apicomplexa: 43 sp., Ciliophora: 7 sp., Perkinsozoa: 1 sp.); SAR–Rhizaria (3 sp.); Excavata–Euglenozoa (15 sp.); Excavata–Metamonada (5 sp.); Amoebozoa (9 sp.); Archaeplastida–Viridiplantae (118 sp.); Archaeplastida–Rhodophyta (3 sp.); Archaeplastida–Chlorophyta (2 sp.); Haptophyta (1 sp.); Cryptophyta (3 sp.); and a residual Protist category (9 sp. including Aureococcus, Blastocystis, Ectocarpus, Fonticula, Hondaea, Nannochloropsis, Symbiodinium, and Thecamonas) representing lineages not assigned to the above named groups. The 43 Apicomplexa species include 6 Plasmodium, 7 Eimeria, 6 Cryptosporidium, and multiple Babesia, Theileria, Hammondia, Toxoplasma, and related species. The 118 Viridiplantae species cover all major angiosperm orders as well as Selaginella (lycophyte), Physcomitrium and Marchantia (bryophytes), and multiple algae (Chara braunii).

Prokaryote sampling was expanded substantially: 87 bacterial species spanning Proteobacteria (α, β, γ, δ, ε), Firmicutes, Actinobacteria, Bacteroidetes, Spirochaetes, Cyanobacteria, and additional phyla; 21 archaeal species spanning Euryarchaeota, Crenarchaeota, Thaumarchaeota, and Nanoarchaeota. Species were selected to represent phylum-level diversity rather than to enrich for high-LCR taxa.

Where a species had multiple protein isoforms per gene, only the longest isoform was retained by parsing gene identifiers from FASTA headers. Protein FASTA files were downloaded from Ensembl FTP or UniProt FTP servers.

### LCR detection and filtering

LCRs were detected using fLPS 2.0 (Harrison 2017) with parameters identical to Teekas et al. (2024) and the prior 94-species analysis: minimum LCR length = 3 amino acids (-m 3); only SINGLE-type records (single-residue compositional bias) retained in post-processing; purity filter ≥70% (dominant amino acid count / LCR length) applied in post-processing. The pre-compiled macOS/Linux binary was used without recompilation.

### Positional binning, enrichment statistics, and domain-level pooling

Positional binning, terminal LCR definition, species-level Fisher's exact test, and protein-level sensitivity analysis were as described previously (one-sided Fisher's exact test; null 10% terminal; terminal = bins 1 or 20 of 20). Phylum-level summaries pool all LCRs from member species before applying Fisher's test.

Domain-level pooled analysis was added for prokaryotes: all LCRs from Bacteria (n=87 species) and from Archaea (n=21 species) were pooled and tested against the 10% null. This approach loses between-species heterogeneity information but provides a powerful aggregate test. The within-domain coefficient of variation of species-level pct_terminal is reported alongside pooled results to characterise heterogeneity (Bacteria CV=50.1%; Archaea CV=44.7%).

### Protein-length confound test (updated)

The length-stratified analysis was extended to Bacteria, Archaea, and all eukaryotic lineages. LCRs from each phylum/domain were stratified into four quartiles by host protein length (global distribution within group) and Fisher's exact test for terminal enrichment applied within each quartile. This directly tests whether the terminal enrichment signal is driven by short proteins mechanically populating terminal bins.

### N-terminal versus C-terminal asymmetry, amino acid identity, driver analysis, purity gradient, GO enrichment

These analyses were completed for all 42 phyla/groups in the dataset. GO term enrichment analysis was performed for species with Ensembl BioMart annotation and sufficient protein identifier overlap; the same five species as the prior analysis reached significance (Pediculus humanus, Tetranychus urticae, Caenorhabditis brenneri, Lottia gigantea, Strongylocentrotus purpuratus). Systematic GO enrichment across the full 724-species dataset was not completed due to BioMart identifier incompatibilities for non-model organisms; this remains a limitation.

### Code availability

All analysis scripts are available at [repository URL]. The pipeline is fully reproducible from raw proteome downloads.

---

## Results

### Terminal LCR enrichment across the animal kingdom: 24 phyla

Across all 723 analysed species (724 downloaded; Acanthaster planci excluded), 280,687 insect LCRs, 328,474 plant LCRs, and 929,801 LCRs total were analysed. The metazoan dataset covers 24 phyla. The U-shaped LCR positional profile (bins 1 and 20 elevated relative to internal bins) is visually apparent across species in the bin heatmap (Figure 1) and the per-phylum profile overlay (Figure 3).

Pooled terminal enrichment is significant in 23 of 24 tested metazoan phyla (Table 1; Figure 2). The single phylum without pooled significance is Acanthocephala (1 species, n=128 LCRs, 16.4% terminal, p=0.07); this is consistent with insufficient statistical power from a single small-proteome species rather than absence of enrichment, as the observed effect size (16.4%) is within the typical bilaterian range. All remaining 23 phyla show pooled terminal enrichment at 14.0–23.1%, overlapping the Tetrapoda baseline of 15–25% reported by Teekas et al. (2024).

**Table 1. Metazoan phylum-level terminal LCR enrichment.**

| Phylum | N spp. | Total LCRs | % Terminal | Pooled sig | Prior study? |
|---|---|---|---|---|---|
| Porifera | 2 | 1,549 | 23.1 | Yes | NEW |
| Xenacoelomorpha | 1 | 412 | 21.8 | Yes | NEW |
| Ctenophora | 1 | 542 | 20.5 | Yes | ✓ |
| Placozoa | 1 | 234 | 18.8 | Yes | ✓ |
| Cnidaria | 13 | 11,181 | 19.5 | Yes | expanded |
| Tardigrada | 2 | 1,843 | 19.0 | Yes | NEW |
| Rotifera | 1 | 3,887 | 15.9 | Yes | NEW |
| Nematoda | 15 | 12,578 | 19.1 | Yes | expanded |
| Platyhelminthes | 7 | 5,425 | 14.6 | Yes | expanded |
| Annelida | 5 | 6,872 | 16.6 | Yes | expanded |
| Nematomorpha | 1 | 181 | 16.6 | Yes | NEW |
| Nemertea | 1 | 741 | 20.2 | Yes | NEW |
| Brachiopoda | 1 | 1,250 | 21.2 | Yes | ✓ |
| Mollusca | 23 | 24,710 | 18.7 | Yes | expanded |
| Crustacea | 21 | 39,253 | 18.1 | Yes | expanded |
| Chelicerata | 28 | 28,033 | 18.1 | Yes | expanded |
| Myriapoda | 1 | 699 | 14.0 | Yes | ✓ |
| Collembola | 2 | 2,918 | 20.3 | Yes | NEW |
| Insecta | 218 | 280,687 | 16.8 | Yes | major expansion |
| Echinodermata | 6† | 6,356 | 19.6 | Yes | expanded |
| Hemichordata | 1 | 797 | 18.8 | Yes | ✓ |
| Chordata | 1 | 1,121 | 15.8 | Yes | ✓ |
| Priapulida | 1 | 1,273 | 14.0 | Yes | ✓ |
| Acanthocephala | 1 | 128 | 16.4 | No (p=0.07) | NEW |

†Echinodermata: 7 species downloaded; Acanthaster planci excluded post-hoc (empty fLPS output); 6 analysed.

The within-phylum coefficient of variation (CV) is now available for 29 phyla/groups (Supplementary Table S3). Low CV is observed for well-sampled phyla: Echinodermata CV=4.6%, Porifera CV=6.0%, Tardigrada CV=7.2%, Collembola CV=8.5%, Mollusca CV=11.6%, Crustacea CV=11.6%, Viridiplantae CV=13.8%, Insecta CV=15.6%, Oomycota CV=16.2%. Single-phylum CVs across the broader protist dataset are higher: Rhizaria CV=44.2%, Cryptophyta CV=72.7%, reflecting high between-species variability in these small-n groups. This confirms the signal is a consistent phylum-level property in well-sampled lineages rather than driven by outlier species.

**Terminal enrichment is length-independent.** Length-stratified analysis (four quartiles by protein length) confirms enrichment in at least two quartiles for all 42 phyla/groups with sufficient LCR counts. Pooled across all metazoans: Q1 (shortest proteins) = 13.9%, Q2 = 22.2%, Q3 = 22.7%, Q4 = 20.0%, all significant (p~0). The modest reduction in Q1 relative to Q2–Q4 is inconsistent with a length-confound artefact, which would predict the highest terminal fraction in the shortest proteins.

**Singleton-LCR proteins drive the signal.** Driver analysis was completed for all 42 phyla/groups. Proteins carrying a single LCR show significant terminal enrichment in 41 of 42 groups; the single exception is Acanthocephala (p=0.07), consistent with insufficient power from one small-proteome species. Multi-LCR proteins show significant enrichment in fewer groups (23/42) and at lower effect sizes. The 19 groups lacking significant multi-LCR enrichment are predominantly small-n phyla (Acanthocephala, Xenacoelomorpha, Nematomorpha, Placozoa, Myriapoda, Hemichordata) and unicellular lineages (Archaea, Chlorophyta, Rhodophyta, Haptophyta, Ctenophora), consistent with either insufficient power or a genuine mechanistic difference whereby multi-LCR accumulation at termini is less prevalent in prokaryotes and basal eukaryotes. Singleton-LCR proteins are confirmed as the primary driver across the full phylogenetic breadth of the dataset.

### New phylogenetic anchors: Porifera, Xenacoelomorpha, and Tardigrada

Three newly included metazoan phyla provide critical phylogenetic anchors.

**Porifera (sponges; 2 spp., 23.1%).** Amphimedon queenslandica (22.6%) and Halichondria panicea (23.1%) both show significant terminal enrichment (p<0.001). Sponges have no neurons, no muscles, and diverged from the animal stem ≥600–650 million years ago. Their terminal LCR enrichment, indistinguishable in magnitude from bilaterian averages, establishes this property as pre-neural and ancestrally metazoan.

**Xenacoelomorpha (1 sp., 21.8%).** Hofstenia miamia shows 21.8% terminal LCRs (p=2×10⁻⁶). The phylogenetic position of Xenacoelomorpha remains debated, though current consensus increasingly places them as sister to Nephrozoa (the bilaterian crown), making them the deepest-branching bilaterian clade. In either placement, their enrichment extends the pattern to the base of Bilateria.

**Tardigrada (2 spp., 19.0%).** Hypsibius exemplaris (18.5%) and Paramacrobiotus metropolitanus (19.5%) are the first tardigrades tested. Tardigrades are known for extreme stress tolerance and unusually high proteome disorder; their terminal LCR enrichment is within the typical bilaterian range, indicating that amplified proteome-wide disorder does not alter the terminal positioning pattern.

### Pan-eukaryotic coverage: all supergroups confirmed

All 16 major eukaryotic lineages tested show significant pooled terminal enrichment (Table 2). This includes the complete SAR supergroup and both deep-branching Excavata clades.

**Table 2. Non-metazoan eukaryote terminal LCR enrichment by lineage.**

| Lineage / Supergroup | N spp. | % Terminal | Pooled sig |
|---|---|---|---|
| Fungi (Opisthokonta) | 10 | 17.8 | Yes |
| Oomycota (SAR–Stramenopiles) | 26 | 21.2 | Yes |
| Bacillariophyta/diatoms (SAR–Stram.) | 5 | 17.5 | Yes |
| Apicomplexa (SAR–Alveolata) | 43 | 12.4 | Yes |
| Ciliophora (SAR–Alveolata) | 7 | 24.2 | Yes |
| Perkinsozoa (SAR–Alveolata) | 1 | 22.0 | Yes |
| Rhizaria (SAR) | 3 | 18.7 | Yes |
| Euglenozoa (Excavata) | 15 | 18.3 | Yes |
| Metamonada (Excavata) | 5 | 18.6 | Yes |
| Amoebozoa | 9 | 16.8 | Yes |
| Viridiplantae (Archaeplastida) | 118 | 24.5 | Yes |
| Rhodophyta (Archaeplastida) | 3 | 13.2 | Yes |
| Chlorophyta (Archaeplastida) | 2 | 11.2 | Yes (p=0.025) |
| Haptophyta | 1 | 13.1 | Yes |
| Cryptophyta | 3 | 19.7 | Yes |
| Protist (diverse) | 9 | 17.4 | Yes |

Euglenozoa (15 spp.) and Metamonada (5 spp.) are among the earliest-diverging eukaryotic clades under most current phylogenies. Their significant enrichment at 18.3% and 18.6% places the origin of terminal LCR bias at or before the last eukaryotic common ancestor (LECA, ~1.5–2 Bya), corroborating the prior study's conclusion.

**Apicomplexa** (43 spp., 12.4%): The phylum is now pooled-significant despite most individual species not reaching significance. The low per-species LCR counts under our purity filter (5–5,600 per species; most <300) mean that pooled significance reflects a real but small-effect enrichment. Importantly, Apicomplexa enrichment is predominantly N-terminal: Eimeria species and Cryptosporidium show significant N-terminal but not C-terminal enrichment. Plasmodium falciparum specifically remains non-significant (10.3% terminal, p=0.46), consistent with prior results. The atypically low Apicomplexa terminal fraction (12.4%) is mechanistically informative: *Plasmodium* and related apicomplexan proteomes contain an unusually high abundance of asparagine-rich LCRs that are predominantly internal rather than terminal (Muralidharan and Goldberg 2013). These asparagine repeats, found in ~30% of the *P. falciparum* proteome, appear to function in immune evasion and are not subject to the same positional selection acting on terminal LCRs in other lineages. Their internal dominance dilutes the phylum-level terminal fraction without implying an absence of terminal enrichment in non-asparagine LCR classes. Consistently, purity gradient analysis shows that Apicomplexa terminal LCRs are significantly purer than internal LCRs (Δ=+0.007, p~0; see *Terminal LCR purity* section), suggesting that even within the low-terminal-fraction Apicomplexa pool, the LCRs that do occupy terminal positions are subject to a qualitative constraint absent from the asparagine-dominated internal pool.

**Chlorophyta** (2 spp., 11.2%, p=0.025): The lowest eukaryotic enrichment level. Chlamydomonas reinhardtii (5,150 LCRs, pct_terminal ~10.7%) drives this result. Notably, Ostreococcus lucimarinus (a second chlorophyte) shows 15.75% N-terminal and 10.96% C-terminal enrichment — both individually significant — suggesting the low phylum-level value is substantially driven by Chlamydomonas rather than being a Chlorophyta-wide phenomenon.

### Terminal LCR enrichment extends to prokaryotes

**Bacteria** (87 spp., 2,706 LCRs, 738 terminal, **27.3%**, pooled p~0) and **Archaea** (21 spp., 370 LCRs, 97 terminal, **26.2%**, pooled p~0) both show significant terminal enrichment when pooled. These are among the highest terminal fractions of any domain analysed, exceeding most eukaryotic lineages.

**Length-confound control is decisive.** Stratifying prokaryote LCRs by protein length quartile (Figure 5):

| Domain | Q1 (shortest) | Q2 | Q3 | Q4 (longest) |
|---|---|---|---|---|
| Bacteria | 20.6%*** | 29.9%*** | 37.0%*** | 29.3%*** |
| Archaea | 19.1%** | 29.9%** | 35.7%** | 37.9%** |

All eight quartile tests are significant (Fisher's exact; *** p<0.001, ** p<0.01). Crucially, the enrichment increases from Q1 to Q3 in both domains — the opposite of what a length-confound artefact would predict (which would inflate terminal fractions disproportionately in short proteins). The Q4 reversal (bacteria 29.3%, archaea 37.9%) may reflect structural differences in very large bacterial proteins, but the overall pattern rules out a simple length artefact as the explanation. Prokaryote terminal enrichment is not an artefact of short median protein lengths.

**Within-domain heterogeneity.** The high CV in bacteria (50.1%) and archaea (44.7%) — versus 13.8% in Viridiplantae across 118 species — reflects that individual prokaryote species are too LCR-poor for reliable individual-species estimates (median ~21 LCRs per bacterial species under our purity filter). The pooled result is robust but should not be interpreted as uniform across all bacterial lineages. Species with sufficient LCRs for individual-species testing include the Actinobacteria (Streptomyces coelicolor 125 LCRs, Streptomyces griseus 89 LCRs, Mycobacterium tuberculosis 42 LCRs), where individual enrichment is significant (25–36% terminal, p<0.005).

**Between-domain comparisons.** Pairwise phylum-level tests (after Holm correction) show that bacteria are significantly more terminal-enriched than Insecta (p<0.001), Crustacea (p=0.0015), and Chelicerata (p=0.0006), and lower than Apicomplexa (p=0.000001). Bacteria and Viridiplantae are not significantly different from each other in pooled terminal fraction.

### The U-shaped positional profile and C/N asymmetry: supergroup-specific signatures

A consistent U-shaped LCR density profile (elevated in bins 1 and 20, depressed in internal bins) is conserved across all domains, replicating the metazoan finding at full phylogenetic breadth (Figure 3).

However, the relative magnitude of N-terminal versus C-terminal enrichment (the asymmetry ratio = pct_nterm / pct_cterm) reveals a striking supergroup-specific pattern not previously reported at scale (Figure 4):

**Viridiplantae — strongly N-terminal dominant.** Across 118 plant species, N-terminal enrichment substantially exceeds C-terminal enrichment in virtually all species (mean asymmetry ratio ≈ 3.0). This pattern is most extreme in grasses (Poaceae): Avena species (ratio 3.9–4.5), Oryza species (~3.7–4.1), Secale cereale (4.3), Triticum aestivum (3.8), Panicum hallii (4.0). Eudicots show lower but consistent N-terminal dominance (Arabidopsis thaliana: ratio 1.93; Glycine max: 3.3; most Brassicaceae 1.8–2.1). Importantly, the pattern is *reduced* in the bryophytes: Physcomitrium patens (moss) shows a balanced ratio (~1.0), and Marchantia polymorpha (liverwort) a moderate ratio (~2.3). This phylogenetic gradient — lowest in mosses, intermediate in lycophytes, highest in grasses — suggests progressive amplification of N-terminal LCR enrichment during land plant evolution.

**Metazoa and Fungi — C-terminal dominant or balanced.** Most metazoan phyla show C-terminal enrichment equal to or greater than N-terminal enrichment (mean asymmetry ratio across all metazoan species ≈ 0.85–1.1). Mollusca and Platyhelminthes are systematically C-dominant (mean ratios ~0.68 and ~0.55 respectively); Platyhelminthes show the most extreme C-terminal bias, particularly in Schistosoma (ratio 0.47–0.61) and Schmidtea species. Fungi are approximately balanced to slightly C-dominant (Saccharomyces cerevisiae ratio 1.2; Schizosaccharomyces pombe 0.58; Aspergillus fumigatus 1.0).

**SAR–Ciliophora — C-terminal dominant.** Stentor coeruleus (ratio 0.64), Paramecium tetraurelia (0.58), Ichthyophthirius multifiliis (0.43), Stylonychia lemnae (0.76) — ciliates consistently show more C-terminal than N-terminal enrichment, despite being Alveolata like Apicomplexa.

**Excavata–Metamonada — strongly C-terminal dominant.** Giardia intestinalis (ratio 0.54), Spironucleus salmonicida (0.35), Tritrichomonas foetus (0.38). Among the deepest-branching eukaryotes, terminal enrichment is almost entirely C-terminal.

**Bacteria with sufficient LCRs — N-terminal dominant.** Streptomyces coelicolor (ratio 2.5), Streptomyces griseus (3.7), Deinococcus radiodurans (3.3), Geobacter sulfurreducens (3.5), Chloroflexus aurantiacus (2.8), Bacteroides thetaiotaomicron (4.7). The high N-terminal fraction in bacteria parallels the Viridiplantae pattern.

### Protein-level confirmation

At the protein level, binomial testing confirms the LCR-positional Fisher's test across all species. Across species with ≥10 LCR-containing proteins, the fraction of proteins carrying ≥1 terminal LCR (observed median ~19%, versus 10% null expectation) is significantly elevated in the overwhelming majority of species (p<0.05, binomial test). This demonstrates that terminal LCR enrichment is a protein-level property — proteins are disproportionately more likely to carry a terminal LCR than expected — rather than arising solely from a small number of proteins with many terminal LCRs.

### Terminal LCR purity: a minority signature in selected lineages

Purity gradient analysis (one-sided Mann-Whitney U test; alternative: terminal LCRs have greater purity than internal LCRs) was completed for all 42 phyla/groups. The majority of groups show no significant difference in purity between terminal and internal LCRs, consistent with terminal enrichment being a positional rather than qualitative phenomenon. Seven phyla show significantly higher purity in terminal LCRs: **Apicomplexa** (Δ=+0.007, p~0), **Rhizaria** (Δ=+0.026, p~0), **Euglenozoa** (Δ=+0.009, p=7.6×10⁻⁵), **Viridiplantae** (Δ=+0.001, p=0.003), **Annelida** (Δ=+0.013, p=0.001), **Chelicerata** (Δ=+0.006, p=0.001), and **Echinodermata** (Δ=+0.014, p=0.002). Rhizaria shows the largest effect (Δ=+0.026), indicating that terminal LCRs in this SAR clade (foraminifera, radiolaria, cercozoans) are substantially purer than internal ones — a signal deserving further investigation of amino acid composition at Rhizaria termini. For the remaining six phyla, effect sizes are small (Δ ≤ 0.014), indicating that while terminal LCRs are marginally purer on average, the qualitative difference is small relative to the strong positional enrichment signal. Terminal LCR purity is not significantly elevated relative to internal LCRs in Bacteria (Δ=+0.0002, p=0.81), Archaea (Δ=−0.004, p=0.33), or the vast majority of metazoan phyla including Insecta (p=1.0), confirming that terminal enrichment in prokaryotes and most metazoans is positional rather than qualitative.

### Amino acid composition, driver analysis, and GO enrichment

Terminal amino acid composition (C, E, K enriched at termini; Q, N not enriched) and driver analysis (singleton-LCR proteins drive enrichment across all 42 phyla/groups) replicate the prior findings and are not qualitatively changed by the expanded dataset. GO enrichment was significant in the same five species as previously; no new GO signals were identified in the additional species. The absence of functional concentration in terminal LCRs across well-annotated metazoan species remains consistent with a mechanism acting on protein architecture regardless of function.

---

## Discussion

### Terminal LCR enrichment is a near-universal property of protein architecture

The prior study established terminal LCR enrichment as pan-metazoan and likely pan-eukaryotic. The present analysis, with 8× more species spanning all domains of life, reaches a stronger conclusion: terminal LCR enrichment is statistically detectable in pooled bacteria, archaea, and all eukaryotic supergroups. The phenomenon is near-universal.

For eukaryotes, the new data provide strong confirmatory evidence. The complete coverage of all SAR supergroup members — including Oomycota (water moulds), diatoms, ciliates, and Rhizaria — fills major gaps. Metamonada (Giardia, Spironucleus, Tritrichomonas) are among the most deeply branching eukaryotes; their enrichment (18.6%) corroborates the LECA origin independently of the fungi/plant data.

For prokaryotes, the results are genuinely unexpected. The pooled bacterial terminal fraction (27.3%) and archaeal fraction (26.2%) exceed those of most eukaryotic lineages (typically 14–24%). The length-stratified control is decisive: enrichment increases from Q1 to Q3 in bacteria (20.6% → 37.0%) with a partial reversal at Q4 (29.3%), ruling out a length-confound artefact and implying that intermediate-to-long bacterial proteins accumulate terminal LCRs at disproportionately high rates.

Two interpretations are possible. **Scenario A (ancestral/universal):** Terminal LCR enrichment reflects a mechanism present in the Last Universal Common Ancestor (LUCA, ~3.5–4 Bya). Under this scenario, the mechanism must operate in all cellular organisms, independently of eukaryote-specific features such as the ubiquitin-proteasome system. Translation kinetics — ribosome pausing near start and stop codons, promoting compositionally simple sequences at protein termini (Irastortza-Olaziregi and Amster-Choder 2021; Pechmann and Frydman 2013) — is a candidate universal mechanism consistent with ribosome biology in all domains. **Scenario B (convergent):** Bacteria and archaea independently evolved terminal LCR enrichment through lineage-specific mechanisms. One concrete candidate is formyl-methionine (fMet) as an N-degron: Piatkov et al. (2015) demonstrated that N-terminal formyl-methionine marks bacterial proteins for proteasome-independent degradation via the N-end rule pathway in bacteria, and related fMet-degron logic extends to eukaryotic cytoplasmic N-terminal processing. If fMet retention or downstream N-terminal sequences are LCR-prone in specific bacterial protein classes, this would specifically enrich N-terminal LCRs, consistent with the N-terminal dominance observed in high-LCR bacteria. Under Scenario B, the bacterial terminal LCR signal has a different mechanistic basis than the eukaryotic C-degron-driven enrichment. Additional Scenario B candidates include N-terminal signal peptides and export sequences (abundant at N-termini of bacterial secreted proteins) and C-terminal ssrA degradation tags (tmRNA-mediated, eukaryote-absent).

The N/C asymmetry data partially discriminates between these scenarios. Bacteria (when LCR-rich) are predominantly N-terminal enriched — as are Viridiplantae — whereas most metazoans and Fungi are C-terminal enriched or balanced. If prokaryote terminal enrichment shared the same mechanism as metazoan C-degron-driven enrichment, one would predict a similar C-terminal bias. The N-terminal dominance in bacteria is more consistent with N-terminal processing biology (signal peptide cleavage, N-formylmethionine removal; Piatkov et al. 2015) or with Scenario B. However, the length-confound control argues against a purely artefactual explanation, so Scenario B would still require a genuine biological mechanism in prokaryotes.

We conclude that the prokaryote result is real but mechanistically uncertain. The most informative next step is a prokaryote-specific analysis of which protein classes carry terminal LCRs (e.g., secreted proteins, membrane proteins, ribosomal proteins) to test whether signal peptides or known prokaryotic degradation tags explain the signal.

### The Viridiplantae N/C asymmetry reversal at scale

The 118-species plant dataset confirms and extends the Arabidopsis N-terminal dominance seen in the prior study. The pattern is consistent across all angiosperm orders tested and is most extreme in Poaceae (grasses), where N-terminal LCR fractions of 20–25% with C-terminal fractions of 5–6% give asymmetry ratios of 3.5–4.5. The phylogenetic gradient — lowest in bryophytes (Physcomitrium: ratio ~1.0), intermediate in eudicots (Arabidopsis: 1.93), highest in grasses — suggests that N-terminal LCR amplification is a derived feature of vascular plant and particularly grass proteomes.

A mechanistic basis is plausible. Land plants have a well-characterised N-degron pathway involving PRT1 and PRT6 E3 ligases that recognise N-terminal basic and oxidised-Cys residues for proteasomal targeting (Gibbs et al. 2014; Holdsworth et al. 2020). This pathway is more elaborate in land plants than in animals, encompassing oxygen-sensing through Plant Cysteine Oxidase (PCO) enzymes that oxidise N-terminal cysteine on Group VII ETHYLENE RESPONSE FACTOR (ERF-VII) transcription factors, triggering PRT6-mediated arginylation and subsequent degradation (Holdsworth et al. 2020). If this pathway preferentially targets LCR-containing proteins with N-terminal degron residues, it would drive accumulation of N-terminal LCRs in proteins subject to N-terminal processing. The progression from mosses (low N-terminal bias) to angiosperms (high) would reflect the progressive elaboration of N-degron pathway complexity during land plant evolution — a testable prediction.

The contrasting C-terminal dominance in Mollusca and especially Platyhelminthes is unexplained. Platyhelminthes (flatworms) have highly reorganised genomes with unusual codon usage, and parasitic species (Schistosoma) undergo extensive host adaptation. Whether the systematic C-terminal bias reflects a genuine proteome-wide shift in protein terminal composition or a lineage-specific loss of N-terminal LCR biology remains to be determined.

### Mechanistic framework: from eukaryotic degrons to universal mechanisms

The original mechanistic proposal for terminal LCR enrichment invoked C-degron and N-degron pathways of the eukaryotic ubiquitin-proteasome system (Koren et al. 2018; Varshavsky 2019). This explanation remains valid for the eukaryotic data and is specifically supported by the amino acid composition of terminal LCRs (cysteine, lysine, glutamate enriched — all known degron residues) and the absence of functional enrichment in GO analysis (consistent with a mechanism targeting protein architecture regardless of function).

However, the prokaryote finding requires broadening the framework. Three mechanisms could operate universally or semi-universally:

1. **Ribosome pausing at start and stop codons** — well-characterised in all cellular life (Irastortza-Olaziregi and Amster-Choder 2021) — would predict N-terminal and C-terminal LCR enrichment symmetrically. The observed C-dominance in most eukaryotes argues this is a contributing but not sole mechanism; pausing-associated slow codons at 5′ mRNA ends are enriched in computationally simple (low-complexity) sequences that may directly predispose translated termini to LCR formation.

2. **N-terminal processing** — signal peptide cleavage, N-terminal methionine excision, N-formylmethionine removal (Piatkov et al. 2015) — is universal and preferentially generates N-terminal LCR-like sequences. This is consistent with N-terminal dominance in bacteria and plants, and less prominent in Metazoa where C-degron biology is amplified.

3. **C-terminal disordered tails** — functional in eukaryotic degron recognition (Koren et al. 2018), tubulin code biology (Janke and Magiera 2020), and RNA-binding — are the most eukaryote-specific mechanism and explain C-terminal dominance in animals and fungi.

We propose a layered model: a universal layer (ribosome kinetics, N-terminal processing) generates baseline terminal LCR enrichment in all cellular life, superimposed by lineage-specific amplification through degron biology (eukaryotes), N-degron pathway elaboration (plants; Holdsworth et al. 2020), or export/degradation sequence biology (bacteria; Piatkov et al. 2015).

### The Platyhelminthes, Chlorophyta, and Acanthocephala outliers

**Platyhelminthes** (14.6% terminal; C-dominant) shows the lowest metazoan terminal enrichment despite pooled significance. The expansion to 7 species confirms this is a genuine phylum-level property, not driven by Schistosoma alone. The pattern (low terminal fraction, C-dominated) contrasts with all other lophotrochozoans.

**Chlorophyta** (11.2%; p=0.025): Chlamydomonas reinhardtii remains the only non-animal, non-apicomplexan eukaryote with weak overall enrichment. Its N-terminal LCR fraction is 3.4% (not significant) — the single clear instance of a eukaryote without N-terminal enrichment. The recently added Ostreococcus lucimarinus, by contrast, shows both N- and C-terminal significance. Chlamydomonas is a flagellate with extensive cilia-based biology; its unusual LCR composition (polyglutamate-enriched tubulin tails) may mask the general terminal enrichment signal.

**Acanthocephala** (p=0.07) is an anomaly of power rather than biology. With one species (Pomphorhynchus laevis) and 128 LCRs, any enrichment below ~19% cannot reach p<0.05; the observed 16.4% at borderline significance is the expected result from a phylum-representative species.

### Limitations

1. Prokaryote mechanistic interpretation remains uncertain. The C/N asymmetry of bacterial terminal enrichment (N-dominant) differs from most eukaryotes (C-dominant), suggesting distinct mechanisms.
2. GO enrichment was completed for the same 5 species as the prior study; systematic coverage of the 724-species dataset is pending.
3. The residual "Protist" category (9 species) aggregates diverse, phylogenetically unrelated lineages; phylum-level enrichment statistics for this group are not biologically meaningful.
4. Single-species phyla (Acanthocephala, Xenacoelomorpha, Ctenophora, Placozoa, Rotifera, Nematomorpha, Nemertea, Brachiopoda, Myriapoda, Hemichordata, Perkinsozoa, Haptophyta) carry high uncertainty at the phylum level.

---

## Conclusions

Terminal LCR enrichment is confirmed across 23 of 24 metazoan phyla (24th non-significant due to power limitations), including the basally-branching Porifera, Xenacoelomorpha, and Tardigrada. All 16 major eukaryotic lineages tested show significant pooled enrichment, establishing terminal LCR bias as a LECA-level property (≥1.5 Bya). Pooled analysis of 87 bacterial and 21 archaeal proteomes reveals significant terminal enrichment (27.3% and 26.2% respectively) that is not explained by protein length, though the high within-domain variability (bacteria CV=50.1%) indicates this represents a domain-level aggregate rather than a species-typical value. The pattern suggests the phenomenon may extend to all cellular life, consistent with a mechanism present since LUCA (~3.5–4 Bya). Singleton-LCR proteins drive the signal in 41 of 42 phyla/groups (Acanthocephala not significant due to power limitations from a single small-proteome species). The N/C asymmetry is a supergroup-specific signature: Viridiplantae are strongly N-terminal enriched (most extreme in grasses), most Metazoa and Fungi are C-terminal enriched or balanced, and Ciliophora, Metamonada, Platyhelminthes, and LCR-rich bacteria are C-terminal enriched. Terminal LCR purity is significantly elevated relative to internal LCRs in seven phyla (Apicomplexa, Rhizaria, Euglenozoa, Viridiplantae, Annelida, Chelicerata, Echinodermata) but not in the majority of groups, confirming that enrichment is primarily positional. Together, these findings support a layered mechanistic model in which universal translational kinetics and N-terminal processing generate a baseline terminal enrichment in all life, amplified by lineage-specific degron and processing pathways — N-degron elaboration in plants, C-degron biology in animals, and N-formylmethionine-dependent processing in bacteria.

---

## Acknowledgements

The author thanks the Ensembl Metazoa, Ensembl Plants, and UniProt teams for providing freely downloadable proteomes, and Paul Harrison for making fLPS 2.0 freely available.

---

## References

Alberti S, Gladfelter A, Mittag T (2019). Considerations and challenges in studying liquid-liquid phase separation and biomolecular condensates. *Cell* 176:419–434.

Berriman M, et al. (2009). The genome of the blood fluke *Schistosoma mansoni*. *Nature* 460:352–358.

Boija A, et al. (2018). Transcription factors activate genes through the phase-separation capacity of their activation domains. *Cell* 175:1842–1855.

Eme L, Sharpe SC, Brown MW, Roger AJ (2014). On the age of eukaryotes: evaluating evidence from fossils and molecular clocks. *Cold Spring Harbor Perspectives in Biology* 6:a016139.

Gibbs DJ, et al. (2014). Homeostatic response to hypoxia is regulated by the N-end rule pathway in plants. *Nature* 479:415–418.

Harrison PM (2017). fLPS: Fast discovery of compositional biases for the protein universe. *BMC Bioinformatics* 18:476.

Holdsworth MJ, Vicente J, Sharma G, Abbas M, Estavillo GM (2020). The plant N-degron pathways of ubiquitin-mediated proteolysis. *Journal of Integrative Plant Biology* 62:70–89.

Irastortza-Olaziregi M, Amster-Choder O (2021). Coupled transcription-translation in prokaryotes: an old couple with new surprises. *Frontiers in Microbiology* 11:619430.

Janke C, Magiera MM (2020). The tubulin code and its role in controlling microtubule properties and functions. *Nature Reviews Molecular Cell Biology* 21:307–326.

Koren I, et al. (2018). The eukaryotic proteome is shaped by E3 ubiquitin ligases targeting C-terminal degrons. *Cell* 173:1622–1635.

Lancaster AK, et al. (2014). PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition. *Bioinformatics* 30:2–3.

Luo H, Gao F, Lin Y (2023). Compensatory relationship between low-complexity regions and gene paralogy in the evolution of prokaryotes. *Proceedings of the National Academy of Sciences* 120:e2215514120. [PMC10120016 — verify author list and DOI before submission]

Marcotte EM, et al. (1999). A census of protein repeats. *Journal of Molecular Biology* 293:151–160.

Muralidharan V, Goldberg DE (2013). Asparagine repeats in *Plasmodium falciparum* proteins: Good for nothing? *PLoS Pathogens* 9:e1003488.

Ntountoumi C, et al. (2019). Low complexity regions in the proteins of prokaryotes perform important functional roles and are highly conserved. *Nucleic Acids Research* 47:9998–10009.

Pechmann S, Frydman J (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology* 20:237–243.

Piatkov KI, Oh J-H, Liu Y, Bhatt DL, Varshavsky A (2015). Formyl-methionine as a degradation signal at the N-termini of bacterial proteins. *Microbial Cell* 2:376–393. [PMC4745127 — verify full author list before submission]

Romero P, et al. (2001). Sequence complexity of disordered protein. *Proteins* 42:38–48.

Shin Y, Brangwynne CP (2017). Liquid phase condensation in cell physiology and disease. *Science* 357:eaaf4382.

Teekas L, Sharma S, Vijay N (2024). Terminal regions of a protein are a hotspot for low complexity regions and selection. *Open Biology* 14:230439.

UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research* 51:D523–D531.

van der Lee R, et al. (2014). Classification of intrinsically disordered regions and proteins. *Chemical Reviews* 114:6589–6631.

Varshavsky A (2019). N-degron and C-degron pathways of protein degradation. *Proceedings of the National Academy of Sciences* 116:358–366.

Virtanen P, et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17:261–272.

Wootton JC, Federhen S (1996). Analysis of compositionally biased regions in sequence databases. *Methods in Enzymology* 266:554–571.

---

## Figure Legends

**Figure 1. LCR positional distribution across 723 proteomes.** Heatmap showing the fraction of LCRs in each of 20 equally spaced positional bins (bin 1 = N-terminal 5%; bin 20 = C-terminal 5%) for each analysed species. Rows are species ordered by domain and phylum (Bacteria → Archaea → non-metazoan eukaryotes → Metazoa, phylogenetically arranged within each group). Colour scale: fraction of a species' LCRs in that bin (0–0.15). Dashed blue lines mark the terminal bins. The elevated signal at bins 1 and 20 relative to internal bins, visible across virtually all species regardless of domain, is the core finding of this study.

**Figure 2. Terminal LCR enrichment by phylum.** Bar chart showing pooled % terminal LCRs (bins 1 + 20) per phylum/group, ordered phylogenetically. Black dots show individual-species values. Red dotted line: 10% null expectation (2/20 bins). Grey shaded band: Tetrapoda range from Teekas et al. (2024; 15–25%). All groups except Acanthocephala (1 sp., insufficient power) significantly exceed the null. Bacteria and Archaea show the highest pooled terminal fractions (27.3% and 26.2% respectively). Viridiplantae (24.5%) exceeds most metazoan phyla.

**Figure 3. Conserved U-shaped LCR positional profile across all phyla.** Each line shows the mean fraction of LCRs per positional bin (1–20) for one phylum/domain, averaged across member species. Dashed grey vertical lines mark terminal bins 1 and 20. Dotted horizontal line: 5% uniform null expectation. All phyla show elevated LCR density at both termini and a depressed internal plateau, demonstrating that the enrichment is a shape-conserved feature of protein architecture, not merely an edge artefact.

**Figure 4. Supergroup-specific N/C terminal asymmetry.** Per-species asymmetry ratio (pct_nterm / pct_cterm), grouped by phylum/domain. A ratio >1 indicates N-terminal dominance; <1 indicates C-terminal dominance. Viridiplantae (118 spp.) are consistently N-terminal dominant (median ratio ~3.0; grasses most extreme, 3.5–4.5); bryophytes are near-balanced (~1.0). Metazoa and Fungi are broadly balanced to C-dominant. Ciliophora, Metamonada, and Platyhelminthes are systematically C-dominant. Bacteria with sufficient LCRs are N-terminal dominant. The pattern reveals that while terminal enrichment is near-universal, the directional bias is supergroup-specific.

**Figure 5. Protein-length stratified terminal LCR enrichment in prokaryotes.** Grouped bar chart showing % terminal LCRs in four protein-length quartiles (Q1 shortest; Q4 longest) for Bacteria (blue) and Archaea (red). Error bars: 95% CI from Fisher's exact test. The monotonic increase with protein length in bacteria (20.6% → 37.0% for Q1–Q3) rules out a length-confound artefact; a confound would produce highest terminal % in Q1.

---

## Supplementary Tables

**Supplementary Table S1.** Full species list (724 entries): species name, phylum, data source (Ensembl / UniProt), proteome ID, protein count, LCR count, pct_terminal.

**Supplementary Table S2.** Driver analysis results for all 42 phyla/groups: pct_terminal for singleton-LCR proteins, multi-LCR proteins; odds ratios and p-values for each class.

**Supplementary Table S3.** Within-phylum coefficient of variation of pct_terminal: 29 phyla/groups with n ≥ 2 species. Includes phylum, n_species, mean_pct_terminal, std_pct_terminal, CV.

**Supplementary Table S4.** Length-stratified analysis for all 42 phyla/groups: pct_terminal by protein-length quartile; Fisher's exact p-value and significance per quartile.

**Supplementary Table S5.** Purity gradient analysis for all 42 phyla/groups: mean purity of terminal vs. internal LCRs, Δ purity, Mann-Whitney U p-value, significance flag.

**Supplementary Figures**

**Supplementary Figure 1.** Per-species terminal LCR % distributions for each phylum, violin/strip chart format. Consistent with Table 1 and Table 2 pooled values but shows individual-species spread.

**Supplementary Figure 2.** Full length-stratified enrichment results for all 42 phyla/groups (4 quartiles × 42 groups heatmap).

**Supplementary Figure 3.** Asymmetry ratio (pct_nterm / pct_cterm) distributions for Viridiplantae only, coloured by plant order. Shows progressive grass amplification relative to bryophytes and basal eudicots.

**Supplementary Figure 4.** Within-phylum CV versus number of species, all 29 phyla with n ≥ 2. Shows convergence of CV to low values as species count increases.

**Supplementary Figure 5.** Amino acid composition of terminal versus internal LCRs: enrichment ratio (terminal/internal) per amino acid, pooled across all metazoans. C, E, K enriched; Q, N depleted at termini.

**Supplementary Figure 6.** Purity gradient: distributions of purity for terminal (bins 1+20) and internal (bins 2–19) LCRs for the seven significant phyla (Apicomplexa, Rhizaria, Euglenozoa, Viridiplantae, Annelida, Chelicerata, Echinodermata). Violin plots with Mann-Whitney U p-values.

**Supplementary Figure 7.** Phylogenetic tree (TimeTree-derived backbone) annotated with pooled pct_terminal per phylum/group. Branch colours indicate supergroup. Node annotations mark LUCA (~3.8 Bya), LECA (~1.8 Bya), and stem Metazoa (~700 Mya).
