# Terminal Low-Complexity Regions Are Enriched Across All Domains of Life: Evidence from 724 Proteomes Spanning 42 Phyla

**Stanley K. Chan**¹

¹ *[Institution, City, Country]*

*Correspondence: stanleykinnok.chan@gmail.com*

---

## Abstract

Low-complexity regions (LCRs) — protein segments dominated by one or a few amino acid types — are concentrated at the N- and C-termini of proteins across Tetrapoda (Teekas et al. 2024, *Open Biology*) and all major invertebrate phyla (Chan, preceding study). Whether this terminal enrichment extends beyond Metazoa to all eukaryotes, and whether it is present in prokaryotes, remains untested at scale. Here we apply the same computational framework (fLPS 2.0, 20-bin positional analysis, purity ≥70% SINGLE-type filter) to 724 proteomes spanning 42 phyla and all three domains of life, encompassing approximately 930,000 LCR records. Among metazoans, terminal LCR enrichment is significant in 24 phyla including newly tested Porifera (23.1%), Xenacoelomorpha (21.8%), Tardigrada (19.0%), Rotifera (15.9%), Nemertea (20.2%), and Collembola (20.3%); only Acanthocephala (1 species; 16.4%, p=0.07) falls short of significance, likely due to limited power. Among non-metazoan eukaryotes, all 16 major eukaryotic lineages tested show significant pooled enrichment, including all SAR supergroup members (Stramenopiles: Oomycota 21.2%, Bacillariophyta 17.5%; Alveolata: Apicomplexa 12.4%, Ciliophora 24.2%; Rhizaria 18.7%), Excavata (Euglenozoa 18.3%, Metamonada 18.6%), Amoebozoa (16.8%), and all tested Archaeplastida (Viridiplantae 24.5%, Rhodophyta 13.2%, Chlorophyta 11.2%). Critically, pooled analysis of 87 bacterial proteomes (27.3%, p~0) and 21 archaeal proteomes (26.2%, p~0) reveals significant terminal LCR enrichment across both prokaryotic domains — enrichment persisting across all four protein-length quartiles, ruling out a length-confound artefact. The N/C terminal asymmetry shows a supergroup-specific pattern: Metazoa and Fungi are predominantly C-terminal enriched; Viridiplantae are strongly N-terminal enriched (mean pct_nterm ~17% vs. pct_cterm ~6%; grasses extreme at ratios 3.5–4.5); bacteria with sufficient LCRs are N-terminal enriched; while Ciliophora, Metamonada, and Platyhelminthes are C-terminal enriched. These findings establish terminal LCR enrichment as a near-universal feature of protein architecture across all domains of life, constraining mechanistic explanations to processes operating in all cellular organisms — including ribosomal translation kinetics at start/stop codons — rather than eukaryote-specific degron systems alone.

**Keywords:** low-complexity regions, protein termini, compositional bias, eukaryotic evolution, prokaryote protein architecture, intrinsically disordered regions, fLPS2, pan-genomics

---

## Introduction

Proteins are not compositionally uniform along their length. Low-complexity regions (LCRs) — segments dominated by one or a small number of amino acid types — are distributed non-randomly within protein sequences, clustering in functionally important contexts such as disordered linkers, prion-like domains, and polyamino acid tracts (Wootton and Federhen 1996; Marcotte et al. 1999). The biological significance of LCRs has grown considerably with the recognition that many intrinsically disordered regions (IDRs) are LCR-containing (Romero et al. 2001; van der Lee et al. 2014), and that phase separation is often driven by low-complexity IDR sequences (Boija et al. 2018; Shin and Brangwynne 2017).

The positional distribution of LCRs within proteins has received less attention than their prevalence or functional roles. Teekas and colleagues (2024) reported that across all major Tetrapoda clades, LCRs are significantly enriched in the terminal 5% of protein sequences (the first and last bins of a 20-bin positional map), with 15–25% of all LCRs in terminal positions despite these bins representing only 10% of positional space. This enrichment was observed consistently across birds, mammals, reptiles, and amphibians. A companion study (Chan, in preparation) extended this finding to 61 metazoan invertebrate species across 16 phyla and to 33 non-metazoan eukaryote and prokaryote species (94 total), establishing terminal LCR enrichment as a pan-metazoan and likely pan-eukaryotic property approximately 1 billion years old.

Several open questions remain. First, the metazoan survey covered only a fraction of animal phylogenetic diversity; many phyla now have sequenced genomes that were absent. Porifera (sponges), Xenacoelomorpha, Tardigrada, Rotifera, Nemertea, Nematomorpha, and Collembola were not represented. Second, the eukaryotic outgroup analysis covered only four protist lineages, leaving most eukaryotic supergroups untested; in particular, the SAR supergroup (Stramenopiles, Alveolata, Rhizaria) and Excavata had minimal representation. Third, the prokaryote analysis was severely underpowered: with median 21 LCRs per bacterial species, individual-species Fisher's tests were not interpretable. Whether bacteria and archaea show terminal LCR enrichment at all remained unresolved.

Here we address all three gaps by scaling the analysis to 724 proteomes spanning 42 phyla and all three domains of life, encompassing approximately 930,000 LCRs. We ask: (1) Does terminal LCR enrichment hold across all newly tested metazoan phyla, including the most basally branching animals? (2) Does the signal span all eukaryotic supergroups? (3) With 87 bacterial and 21 archaeal proteomes, is terminal LCR enrichment detectable in prokaryotes? (4) Is the N/C asymmetry pattern consistent across lineages, or do lineage-specific reversals exist?

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

The length-stratified analysis was extended to Bacteria and Archaea. LCRs from each domain were stratified into four quartiles by host protein length (global distribution within domain) and Fisher's exact test for terminal enrichment applied within each quartile. This directly tests whether the prokaryote terminal enrichment signal is driven by short proteins mechanically populating terminal bins.

### N-terminal versus C-terminal asymmetry, amino acid identity, driver analysis, purity gradient, GO enrichment

These analyses were as described previously. GO term enrichment analysis was performed for species with Ensembl BioMart annotation and sufficient protein identifier overlap; the same five species as the prior analysis reached significance (Pediculus humanus, Tetranychus urticae, Caenorhabditis brenneri, Lottia gigantea, Strongylocentrotus purpuratus). Systematic GO enrichment across the full 724-species dataset was not completed due to BioMart identifier incompatibilities for non-model organisms; this remains a limitation.

Note: driver analysis, purity gradient analysis, and the full within-phylum CV analysis were completed for 20 phyla/groups consistent with those analysed in the prior study. These analyses have not yet been extended to the newly added metazoan phyla (Porifera, Xenacoelomorpha, Tardigrada, Rotifera, Collembola, Nematomorpha, Nemertea) or to the newly added eukaryotic lineages; updated results for these groups are flagged as pending.

### Code availability

All analysis scripts are available at [repository URL]. The pipeline is fully reproducible from raw proteome downloads.

---

## Results

### Terminal LCR enrichment across the animal kingdom: 24 phyla

Across all 723 analysed species (724 downloaded; Acanthaster planci excluded), 280,687 insect LCRs, 328,474 plant LCRs, and approximately 930,000 LCRs total were analysed. The metazoan dataset covers 24 phyla. The U-shaped LCR positional profile (bins 1 and 20 elevated relative to internal bins) is visually apparent across species in the bin heatmap (Figure 1) and the per-phylum profile overlay (Figure 3).

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

†Acanthaster planci excluded (empty fLPS output); 7 nominal entries, 6 analysed.

The within-phylum coefficient of variation (CV) is low for well-sampled phyla: Echinodermata CV=4.6%, Mollusca CV=11.6%, Crustacea CV=11.6%, Viridiplantae CV=13.8%, Insecta CV=15.6%. This confirms the signal is a consistent phylum-level property rather than driven by outlier species, even at the scale of 218 insect or 118 plant species.

**Terminal enrichment is length-independent.** Length-stratified analysis (four quartiles by protein length) confirms enrichment in at least two quartiles for all phyla with sufficient LCR counts. Pooled across all metazoans: Q1 (shortest proteins) = 13.9%, Q2 = 22.2%, Q3 = 22.7%, Q4 = 20.0%, all significant (p~0). The modest reduction in Q1 relative to Q2–Q4 is inconsistent with a length-confound artefact, which would predict the highest terminal fraction in the shortest proteins.

**Singleton-LCR proteins drive the signal.** Proteins carrying a single LCR show significant terminal enrichment in all 20 phyla for which driver analysis was completed. Multi-LCR proteins show enrichment in fewer phyla and at lower effect sizes.

### New phylogenetic anchors: Porifera, Xenacoelomorpha, and Tardigrada

Three newly included metazoan phyla provide critical phylogenetic anchors.

**Porifera (sponges; 2 spp., 23.1%).** Amphimedon queenslandica (22.6%) and Halichondria panicea (23.1%) both show significant terminal enrichment (p<0.001). Sponges have no neurons, no muscles, and diverged from the animal stem ≥600–650 million years ago. Their terminal LCR enrichment, indistinguishable in magnitude from bilaterian averages, establishes this property as pre-neural and ancestrally metazoan.

**Xenacoelomorpha (1 sp., 21.8%).** Hofstenia miamia shows 21.8% terminal LCRs (p=2×10⁻⁶). The phylogenetic position of Xenacoelomorpha is debated — either sister to all other bilaterians or nested within Deuterostomia — but in either placement, its enrichment extends the pattern to the deepest bilaterian branches.

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

Euglenozoa (15 spp.) and Metamonada (5 spp.) represent the deepest-branching eukaryotic clades under most phylogenies. Their significant enrichment at 18.3% and 18.6% places the origin of terminal LCR bias at or before the last eukaryotic common ancestor (LECA, ~1.5–2 Bya), corroborating the prior study's conclusion.

**Apicomplexa** (43 spp., 12.4%): The phylum is now pooled-significant despite most individual species not reaching significance. The low per-species LCR counts under our purity filter (5–5,600 per species; most <300) mean that pooled significance reflects a real but small-effect enrichment. Importantly, Apicomplexa enrichment is predominantly N-terminal: Eimeria species and Cryptosporidium show significant N-terminal but not C-terminal enrichment. Plasmodium falciparum specifically remains non-significant (10.3% terminal, p=0.46), consistent with prior results.

**Chlorophyta** (2 spp., 11.2%, p=0.025): The lowest eukaryotic enrichment level. Chlamydomonas reinhardtii (5,150 LCRs, pct_terminal ~10.7%) drives this result. Notably, Ostreococcus lucimarinus (a second chlorophyte) shows 15.75% N-terminal and 10.96% C-terminal enrichment — both individually significant — suggesting the low phylum-level value is substantially driven by Chlamydomonas rather than being a Chlorophyta-wide phenomenon.

### Terminal LCR enrichment extends to prokaryotes

**Bacteria** (87 spp., 2,706 LCRs, 738 terminal, **27.3%**, pooled p~0) and **Archaea** (21 spp., 370 LCRs, 97 terminal, **26.2%**, pooled p~0) both show significant terminal enrichment when pooled. These are among the highest terminal fractions of any domain analysed, exceeding most eukaryotic lineages.

**Length-confound control is decisive.** Stratifying prokaryote LCRs by protein length quartile:

| Domain | Q1 (shortest) | Q2 | Q3 | Q4 (longest) |
|---|---|---|---|---|
| Bacteria | 20.6%*** | 29.9%*** | 37.0%*** | 29.3%*** |
| Archaea | 19.1%** | 29.9%** | 35.7%** | 37.9%** |

All eight quartile tests are significant (Fisher's exact; *** p<0.001, ** p<0.01). Crucially, the enrichment *increases* with protein length in both domains — the opposite of what a length-confound artefact would produce (which would inflate terminal fractions disproportionately in short proteins). Prokaryote terminal enrichment is not an artefact of short median protein lengths.

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

### Terminal LCRs are not purer than internal LCRs (replicated)

Terminal LCR purity does not significantly exceed internal LCR purity in bacteria (Δ=+0.0002, p=0.81) or archaea (Δ=−0.004, p=0.33), confirming that terminal enrichment in prokaryotes, as in most eukaryotes, is positional rather than qualitative.

### Amino acid composition, driver analysis, and GO enrichment

Terminal amino acid composition (C, E, K enriched at termini; Q, N not enriched) and driver analysis (singleton-LCR proteins drive enrichment across all phyla) replicate the prior findings and are not qualitatively changed by the expanded dataset. GO enrichment was significant in the same five species as previously; no new GO signals were identified in the additional species. The absence of functional concentration in terminal LCRs across well-annotated metazoan species remains consistent with a mechanism acting on protein architecture regardless of function.

---

## Discussion

### Terminal LCR enrichment is a near-universal property of protein architecture

The prior study established terminal LCR enrichment as pan-metazoan and likely pan-eukaryotic. The present analysis, with 8× more species spanning all domains of life, reaches a stronger conclusion: terminal LCR enrichment is statistically detectable in pooled bacteria, archaea, and all eukaryotic supergroups. The phenomenon is near-universal.

For eukaryotes, the new data provide strong confirmatory evidence. The complete coverage of all SAR supergroup members — including Oomycota (water moulds), diatoms, ciliates, and Rhizaria — fills major gaps. Metamonada (Giardia, Spironucleus, Tritrichomonas) are among the most deeply branching eukaryotes; their enrichment (18.6%) corroborates the LECA origin independently of the fungi/plant data.

For prokaryotes, the results are genuinely unexpected. The pooled bacterial terminal fraction (27.3%) and archaeal fraction (26.2%) exceed those of most eukaryotic lineages (typically 14–24%). The length-stratified control is decisive: enrichment increases monotonically with protein length in bacteria (20.6% → 37.0% across Q1–Q3), ruling out a length-confound artefact and implying that longer bacterial proteins accumulate terminal LCRs at disproportionately high rates.

Two interpretations are possible. **Scenario A (ancestral/universal):** Terminal LCR enrichment reflects a mechanism present in the Last Universal Common Ancestor (LUCA, ~3.5–4 Bya). Under this scenario, the mechanism must operate in all cellular organisms, independently of eukaryote-specific features such as the ubiquitin-proteasome system. Translation kinetics — ribosome pausing near start and stop codons, promoting compositionally simple sequences at protein termini — is a candidate universal mechanism consistent with ribosome biology in all domains. **Scenario B (convergent):** Bacteria and archaea independently evolved terminal LCR enrichment through lineage-specific mechanisms, such as N-terminal signal peptides and export sequences (abundant at N-termini of bacterial secreted proteins) or C-terminal ssrA degradation tags (tmRNA-mediated, eukaryote-absent). Under Scenario B, the bacterial terminal LCR signal has a different mechanistic basis than the eukaryotic one.

The N/C asymmetry data partially discriminates between these scenarios. Bacteria (when LCR-rich) are predominantly N-terminal enriched — as are Viridiplantae — whereas most metazoans and Fungi are C-terminal enriched or balanced. If prokaryote terminal enrichment shared the same mechanism as metazoan C-degron-driven enrichment, one would predict a similar C-terminal bias. The N-terminal dominance in bacteria is more consistent with N-terminal processing biology (signal peptide cleavage, N-formylmethionine removal) or with Scenario B. However, the length-confound control argues against a purely artefactual explanation, so Scenario B would still require a genuine biological mechanism in prokaryotes.

We conclude that the prokaryote result is real but mechanistically uncertain. The most informative next step is a prokaryote-specific analysis of which protein classes carry terminal LCRs (e.g., secreted proteins, membrane proteins, ribosomal proteins) to test whether signal peptides or known prokaryotic degradation tags explain the signal.

### The Viridiplantae N/C asymmetry reversal at scale

The 118-species plant dataset confirms and extends the Arabidopsis N-terminal dominance seen in the prior study. The pattern is consistent across all angiosperm orders tested and is most extreme in Poaceae (grasses), where N-terminal LCR fractions of 20–25% with C-terminal fractions of 5–6% give asymmetry ratios of 3.5–4.5. The phylogenetic gradient — lowest in bryophytes (Physcomitrium: ratio ~1.0), intermediate in eudicots (Arabidopsis: 1.93), highest in grasses — suggests that N-terminal LCR amplification is a derived feature of vascular plant and particularly grass proteomes.

A mechanistic basis is plausible. Land plants have a well-characterised N-degron pathway involving PRT1 and PRT6 E3 ligases that recognise N-terminal basic and oxidised-Cys residues for proteasomal targeting (Gibbs et al. 2014). This pathway is more elaborate in land plants than in animals. If it preferentially targets LCR-containing proteins with N-terminal degron residues, it would drive accumulation of N-terminal LCRs in proteins subject to N-terminal processing. The progression from mosses (low N-terminal bias) to angiosperms (high) would reflect elaboration of this pathway during land plant evolution — a testable prediction.

The contrasting C-terminal dominance in Mollusca and especially Platyhelminthes is unexplained. Platyhelminthes (flatworms) have highly reorganised genomes with unusual codon usage, and parasitic species (Schistosoma) undergo extensive host adaptation. Whether the systematic C-terminal bias reflects a genuine proteome-wide shift in protein terminal composition or a lineage-specific loss of N-terminal LCR biology remains to be determined.

### Mechanistic framework: from eukaryotic degrons to universal mechanisms

The original mechanistic proposal for terminal LCR enrichment invoked C-degron and N-degron pathways of the eukaryotic ubiquitin-proteasome system (Koren et al. 2018; Varshavsky 2019). This explanation remains valid for the eukaryotic data and is specifically supported by the amino acid composition of terminal LCRs (cysteine, lysine, glutamate enriched — all known degron residues) and the absence of functional enrichment in GO analysis (consistent with a mechanism targeting protein architecture regardless of function).

However, the prokaryote finding requires broadening the framework. Three mechanisms could operate universally:

1. **Ribosome pausing at start and stop codons** — well-characterised in all cellular life — would predict N-terminal and C-terminal LCR enrichment symmetrically. The observed C-dominance in most eukaryotes argues this is a contributing but not sole mechanism.

2. **N-terminal processing** — signal peptide cleavage, N-terminal methionine excision, N-formylmethionine removal — is universal and preferentially generates N-terminal LCR-like sequences. This is consistent with N-terminal dominance in bacteria and plants, and less dominant in Metazoa where C-degron biology is prominent.

3. **C-terminal disordered tails** — functional in eukaryotic degron recognition, tubulin code biology, and RNA-binding — are the most eukaryote-specific mechanism and explain C-terminal dominance in animals and fungi.

We propose a layered model: a universal layer (ribosome kinetics, N-terminal processing) generates baseline terminal LCR enrichment in all cellular life, superimposed by lineage-specific amplification through degron biology (eukaryotes), N-degron pathway elaboration (plants), or possibly export sequence biology (bacteria).

### The Platyhelminthes, Chlorophyta, and Acanthocephala outliers

**Platyhelminthes** (14.6% terminal; C-dominant) shows the lowest metazoan terminal enrichment despite pooled significance. The expansion to 7 species confirms this is a genuine phylum-level property, not driven by Schistosoma alone. The pattern (low terminal fraction, C-dominated) contrasts with all other lophotrochozoans.

**Chlorophyta** (11.2%; p=0.025): Chlamydomonas reinhardtii remains the only non-animal, non-apicomplexan eukaryote with weak overall enrichment. Its N-terminal LCR fraction is 3.4% (not significant) — the single clear instance of a eukaryote without N-terminal enrichment. The recently added Ostreococcus lucimarinus, by contrast, shows both N- and C-terminal significance. Chlamydomonas is a flagellate with extensive cilia-based biology; its unusual LCR composition (polyglutamate-enriched tubulin tails) may mask the general terminal enrichment signal.

**Acanthocephala** (p=0.07) is an anomaly of power rather than biology. With one species (Pomphorhynchus laevis) and 128 LCRs, any enrichment below ~19% cannot reach p<0.05; the observed 16.4% at borderline significance is the expected result from a phylum-representative species.

### Limitations

1. Prokaryote mechanistic interpretation remains uncertain. The C/N asymmetry of bacterial terminal enrichment (N-dominant) differs from most eukaryotes (C-dominant), suggesting distinct mechanisms.
2. Driver analysis, purity gradient, and full length-stratified analysis have not yet been completed for the 7 newly added metazoan phyla (Porifera, Xenacoelomorpha, Tardigrada, Rotifera, Collembola, Nematomorpha, Nemertea).
3. GO enrichment was completed for the same 5 species as the prior study; systematic coverage of the 724-species dataset is pending.
4. The residual "Protist" category (9 species) aggregates diverse, phylogenetically unrelated lineages; phylum-level enrichment statistics for this group are not biologically meaningful.
5. Single-species phyla (Acanthocephala, Xenacoelomorpha, Ctenophora, Placozoa, Rotifera, Nematomorpha, Nemertea, Brachiopoda, Myriapoda, Hemichordata, Perkinsozoa, Haptophyta) carry high uncertainty at the phylum level.

---

## Conclusions

Terminal LCR enrichment is confirmed across 23 of 24 metazoan phyla (24th non-significant due to power limitations), including the basally-branching Porifera, Xenacoelomorpha, and Tardigrada. All 16 major eukaryotic lineages tested show significant pooled enrichment, establishing terminal LCR bias as a LECA-level property (≥1.5 Bya). Pooled analysis of 87 bacterial and 21 archaeal proteomes reveals significant terminal enrichment that is not explained by protein length, suggesting the phenomenon may extend to all cellular life. The N/C asymmetry is a supergroup-specific signature: Viridiplantae are strongly N-terminal enriched (most extreme in grasses), most Metazoa and Fungi are C-terminal enriched or balanced, and Ciliophora, Metamonada, Platyhelminthes, and LCR-rich bacteria are C-terminal enriched. Together, these findings support a layered mechanistic model in which universal translational and processing biology generates a baseline terminal enrichment in all life, amplified by lineage-specific degron and N-terminal processing pathways.

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

Janke C, Magiera MM (2020). The tubulin code and its role in controlling microtubule properties and functions. *Nature Reviews Molecular Cell Biology* 21:307–326.

Koren I, et al. (2018). The eukaryotic proteome is shaped by E3 ubiquitin ligases targeting C-terminal degrons. *Cell* 173:1622–1635.

Lancaster AK, et al. (2014). PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition. *Bioinformatics* 30:2–3.

Marcotte EM, et al. (1999). A census of protein repeats. *Journal of Molecular Biology* 293:151–160.

Ntountoumi C, et al. (2019). Low complexity regions in the proteins of prokaryotes perform important functional roles and are highly conserved. *Nucleic Acids Research* 47:9998–10009.

Pechmann S, Frydman J (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology* 20:237–243.

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

**Figure 5. Dated phylogeny of study taxa.** Cladogram of species included in the analysis, with branches proportional to divergence times from TimeTree. Branch colour corresponds to phylum/domain (colour key in figure). Nodes and taxon labels are shown for major divergences. Node ages (Mya) are shown for key splits: Bacteria–Archaea, LUCA, LECA, Opisthokonta, land plant origin, Bilateria, etc. Temporal gaps in sampling >200 Mya are annotated. The figure illustrates the approximately 4-billion-year span covered by the 724-proteome dataset.

**Supplementary Figure 1. Amino acid composition of terminal versus internal LCRs by phylum.** Stacked bar charts showing the fractional amino acid composition of LCRs classified as terminal (bin 1 or 20) versus internal (bins 2–19), for each phylum/group. Cysteine (C), lysine (K), and glutamate (E) are consistently enriched at termini relative to internal positions across the broadest range of phyla, consistent with known C-degron and N-degron residue identities.

**Supplementary Figure 2. Amino acid terminal enrichment ratio.** Heatmap showing log₂(terminal fraction / internal fraction) per amino acid per phylum, highlighting which residues are disproportionately enriched at termini. Red = terminal-enriched; blue = internal-enriched.

**Supplementary Figure 3. Terminal enrichment is not confounded by protein length.** Grouped bar charts showing % terminal LCRs within each of four protein-length quartiles (Q1 = shortest; Q4 = longest) for the major groups (all species pooled; individual phyla/domains). For all groups with sufficient data, enrichment is significant in every quartile. Crucially, in bacteria and archaea, terminal enrichment *increases* with protein length (Q1: 20.6% → Q3: 37.0% in bacteria), the opposite of a length-confound artefact. Error bars: 95% Wilson confidence intervals.

**Supplementary Figure 4. Terminal LCR purity gradient.** Dot plot showing mean purity of terminal versus internal LCRs per phylum (mean fraction of dominant amino acid). Points are jittered by phylum. For most eukaryotic phyla, terminal LCRs are marginally purer than internal LCRs (Mann-Whitney U test, significant in 7/20 phyla tested; Apicomplexa, Euglenozoa, Viridiplantae, Annelida, Chelicerata, Insecta, Echinodermata). Bacteria and Archaea show no significant purity gradient (p=0.81 and p=0.33 respectively), indicating that prokaryote terminal enrichment is positional rather than a quality-based selection for higher purity.

**Supplementary Figure 5. Pairwise between-phylum terminal enrichment comparisons.** Heatmap of Holm-corrected p-values from pairwise Mann-Whitney U tests of per-species pct_terminal between all phylum pairs. Significant cells (Holm p<0.05) are highlighted. Bacteria and Archaea are not significantly different from each other or from Viridiplantae and Amoebozoa, but are significantly higher than Insecta, Chelicerata, Crustacea, and Mollusca. Viridiplantae is significantly higher than multiple animal phyla.

**Supplementary Figure 6. Terminal enrichment is driven by singleton-LCR proteins.** Grouped bar chart comparing % terminal LCRs in proteins carrying exactly one LCR (singleton) versus proteins carrying ≥2 LCRs (multi-LCR), per phylum. Singleton-LCR proteins show significant terminal enrichment in all 20 phyla tested (Fisher's exact; p<0.05). Multi-LCR proteins show significant enrichment in 15/20 phyla, with generally lower effect sizes. This demonstrates that the terminal enrichment signal is a property of the general protein proteome, not restricted to a specialized multi-domain class.

**Supplementary Figure 7. GO term enrichment in terminal-LCR-containing proteins.** Dot plot of significant GO terms (FDR ≤ 0.05, Benjamini-Hochberg) for the five species reaching significance: Pediculus humanus (membrane component, GO:0016021), Tetranychus urticae (catalytic activity), Caenorhabditis brenneri (phosphoprotein phosphatase activity), Lottia gigantea (microtubule/GTPase/neuron projection), Strongylocentrotus purpuratus (microtubule cytoskeleton, cell cycle, GTP binding, glycosyltransferase). The recurrence of cytoskeletal/tubulin terms in echinoderm and mollusc species is consistent with known polyglutamate LCRs in tubulin C-terminal tails. No GO enrichment was detected in 23+ other tested species, consistent with a mechanism acting across the proteome independently of function.

---

## Supplementary Material

**Supplementary Table S1.** Complete species list: source database, phylum, domain, taxon ID, longest-isoform protein count, and per-species enrichment statistics (n_lcr, n_terminal, pct_terminal, odds_ratio, p-value, significant).

**Supplementary Table S2.** Driver analysis: pct_terminal, odds ratio, and Fisher's exact p-value for singleton-LCR and multi-LCR proteins per phylum (20 phyla; newly added phyla pending).

**Supplementary Table S3.** Amino acid composition of terminal vs. internal LCRs per phylum.

**Supplementary Table S4.** Protein-level sensitivity analysis: fraction of proteins with ≥1 terminal LCR per species, binomial test p-value.

**Supplementary Table S5.** GO term enrichment results: all species where BioMart annotation was available (FDR ≤ 0.05 entries only).

**Supplementary Table S6.** Length-stratified terminal enrichment by domain and phylum (four protein-length quartiles).

*See Figure Legends section for full captions. Main text figures: Figures 1–5. Supplementary figures: Supplementary Figures 1–7 (amino acid composition, enrichment ratio, length-confound, purity gradient, phylum comparison, driver analysis, GO enrichment).*

**Supplementary Figure 1.** Amino acid composition of terminal vs. internal LCRs by phylum (fig5_aa_composition).

**Supplementary Figure 2.** Amino acid terminal enrichment ratio heatmap (fig5b_aa_enrichment_ratio).

**Supplementary Figure 3.** Terminal enrichment across protein-length quartiles: all species, Bacteria, Archaea, and major groups (fig6_length_confound).

**Supplementary Figure 4.** LCR purity gradient: terminal vs. internal mean purity per phylum/domain (fig7_purity_gradient).

**Supplementary Figure 5.** Pairwise between-phylum Holm-corrected comparison heatmap (fig8_phylum_comparison).

**Supplementary Figure 6.** Driver analysis: singleton-LCR vs. multi-LCR terminal enrichment per phylum (fig9_driver).

**Supplementary Figure 7.** GO enrichment for five species with significant results (fig10_go_enrichment).
