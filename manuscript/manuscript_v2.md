# Terminal Low-Complexity Regions Are a Pan-Eukaryotic Property Conserved Across ≥1 Billion Years of Evolution

**Stanley K. Chan**¹

¹ *[Institution, City, Country]*

*Correspondence: stanleykinnok.chan@gmail.com*

---

## Abstract

Low-complexity regions (LCRs) — protein segments dominated by one or a few amino acid types — are concentrated at the N- and C-termini of proteins in Tetrapoda (Teekas et al. 2024, *Open Biology*). Whether this terminal enrichment reflects a vertebrate-specific innovation or an ancient metazoan property remains untested. We apply the same computational framework (fLPS 2.0, 20-bin positional analysis) to 61 metazoan species spanning 16 phyla (Ctenophora to Cephalochordata) and 33 non-metazoan outgroup species (7 fungi, 2 plants, 4 protist lineages, 5 archaea, 14 bacteria), encompassing 94 species and over 85,000 LCR records. Of 61 analysed metazoan species (one echinoderm, Acanthaster planci, was excluded due to a failed FLPS run), 60 (98%) show statistically significant terminal LCR enrichment (Fisher's exact test p < 0.05), with phylum-level terminal fractions of 12.7–21.2%, overlapping the Tetrapoda baseline of 15–25%. Mnemiopsis leidyi (Ctenophora; 20.5%, p = 1×10⁻⁶) and *Trichoplax adhaerens* (Placozoa; 18.8% terminal, p = 0.004) — representing the most basal animal lineage diverging ≥700 million years ago — falls within the Tetrapoda range, establishing terminal LCR enrichment as a pre-bilaterian property. Critically, significant terminal enrichment was also observed in 11 of 13 non-metazoan eukaryotes tested, including all seven fungi, both land plants (Arabidopsis thaliana 23.6%, Solanum lycopersicum 19.7%), and *Trypanosoma brucei* (19.3%), extending the phenomenon to at least ~1 billion years ago and the last common ancestor of Opisthokonta and plants. The single exceptions among eukaryotes were *Plasmodium falciparum* (Apicomplexa; 10.3%, p = 0.46), and the green alga *Chlamydomonas reinhardtii* (Chlorophyta; 10.7%, p = 0.116) also showed no enrichment. Prokaryote proteomes carry too few LCRs per species (median 21 per proteome) for individual-species tests to be interpretable; the prokaryote boundary remains unresolved. C-terminal enrichment (significant in 94% of species) is more universal than N-terminal enrichment (86%). Terminal enrichment is robust across all protein length quartiles and driven primarily by singleton-LCR proteins. Terminal LCRs are not significantly purer than internal LCRs in 12 of 16 metazoan phyla, indicating that the enrichment is positional rather than qualitative.

**Keywords:** low-complexity regions, protein termini, compositional bias, metazoan evolution, intrinsically disordered regions, invertebrates, eukaryotic evolution, fLPS2

---

## Introduction

Proteins are not compositionally uniform along their length. Low-complexity regions (LCRs) — segments dominated by one or a small number of amino acid types — are distributed non-randomly within protein sequences, clustering in functionally important contexts such as disordered linkers, prion-like domains, and polyamino acid tracts (Wootton and Federhen 1996; Marcotte et al. 1999). The biological significance of LCRs has grown considerably with the recognition that many intrinsically disordered regions (IDRs) are LCR-containing (Romero et al. 2001; van der Lee et al. 2014), and that phase separation — a mechanism central to membrane-less organelle formation — is often driven by low-complexity IDR sequences (Boija et al. 2018; Shin and Brangwynne 2017).

The positional distribution of LCRs within proteins has received less attention than their prevalence or functional roles. Teekas and colleagues (2024) recently reported that, across all major Tetrapoda clades, LCRs are significantly enriched in the terminal 5% of protein sequences (the first and last bins of a 20-bin positional map), with 15–25% of all LCRs in terminal positions despite these bins representing only 10% of positional space. This enrichment was observed consistently across birds, mammals, reptiles, and amphibians, but the evolutionary age and phylogenetic generality of the mechanism were not addressed.

A critical open question is whether terminal LCR enrichment is a vertebrate-specific innovation or an ancient metazoan property. If it is lineage-specific, it may reflect vertebrate-specific biology such as adaptive immunity or protein domain rearrangements in the vertebrate ancestor. If it extends to the most basal animals, it must reflect a fundamental constraint on protein terminal composition common to all animals.

Invertebrate proteomes span enormous biological diversity across approximately 700–800 million years of animal evolution, encompassing the origin of complex body plans, the nervous system, and the Cambrian diversification. Testing terminal LCR enrichment across this range provides a direct test of evolutionary depth.

Here we present a fully computational pan-invertebrate survey of terminal LCR enrichment. Using fLPS 2.0 (Harrison 2017) with parameters closely matched to Teekas et al. (2024) — with one methodological difference noted below — we analyse 50 species from 14 phyla downloaded from Ensembl Metazoa release 63. We extend the analysis to (1) separate N- and C-terminal enrichment, (2) test for protein-length confounding, (3) characterise amino acid composition at terminal versus internal LCRs, and (4) identify driver protein classes.

---

## Methods

### Species selection and proteome acquisition

Sixty-two metazoan species were selected to maximise phylogenetic breadth across 16 phyla: Ctenophora (1 sp.), Placozoa (1), Cnidaria (5), Platyhelminthes (2), Annelida (3), Nematoda (8), Priapulida (1), Brachiopoda (1), Mollusca (7), Crustacea (6), Myriapoda (1), Chelicerata (6), Insecta (12), Echinodermata (6, of which one was excluded post-hoc; see below), Hemichordata (1), and Cephalochordata (1). *Acanthaster planci* (Echinodermata) was excluded from all analyses because its fLPS 2.0 run produced an empty output file; the remaining 61 metazoan species are analysed. Phyla represented by single species (Ctenophora, Placozoa, Priapulida, Brachiopoda, Myriapoda, Hemichordata, Cephalochordata) should be interpreted with appropriate caution given limited power to estimate within-phylum variation. Protein FASTA files (.pep.all.fa.gz) were downloaded from Ensembl Metazoa FTP (release 63). Where a species had multiple protein isoforms per gene, only the longest isoform was retained by parsing the `gene:` or `gene=` field in Ensembl FASTA headers.

To determine the evolutionary boundary of terminal LCR enrichment, we additionally downloaded 33 non-metazoan reference proteomes from UniProt (release 2024): 7 fungi (*Saccharomyces cerevisiae*, *Schizosaccharomyces pombe*, *Neurospora crassa*, *Aspergillus nidulans*, *Candida albicans*, *Aspergillus fumigatus*, *Cryptococcus neoformans*), 2 land plants (*Arabidopsis thaliana*, *Solanum lycopersicum*), 4 protist lineages (*Dictyostelium discoideum* — Amoebozoa; *Plasmodium falciparum* — Apicomplexa; *Trypanosoma brucei* — Euglenozoa; *Chlamydomonas reinhardtii* — Chlorophyta), 5 archaea (*Methanocaldococcus jannaschii*, *Sulfolobus acidocaldarius*, *Halobacterium salinarum*, *Pyrococcus abyssi*, *Thermococcus kodakarensis*), and 14 bacteria (*Escherichia coli* K12, *Bacillus subtilis*, *Pseudomonas aeruginosa*, *Mycobacterium tuberculosis*, *Streptomyces coelicolor*, *Streptomyces griseus*, *Caulobacter crescentus*, *Staphylococcus aureus*, *Streptococcus pneumoniae*, *Borrelia burgdorferi*, *Treponema pallidum*, *Deinococcus radiodurans*, *Thermotoga maritima*, *Clostridioides difficile*). The same fLPS 2.0 pipeline and statistical framework was applied to all outgroup proteomes without modification. Prokaryotic proteomes contain very few LCRs passing the purity ≥70% filter (median 21 per proteome; range 6–125), rendering individual-species Fisher's tests severely underpowered for most species. Exceptions are the Actinobacteria (*Streptomyces coelicolor* 125 LCRs, *Streptomyces griseus* 89 LCRs, *Mycobacterium tuberculosis* 42 LCRs), which have sufficient LCR counts to test individually; these three show significant enrichment (31–36% terminal; p < 0.001), but cannot be interpreted without broader Actinobacteria sampling. The complete species list is provided in Supplementary Table S1.

### LCR detection and filtering

LCRs were detected using fLPS 2.0 (Harrison 2017; tool available at https://github.com/paulhorton/fLPS2), a length-unbiased compositional bias detector. The pre-compiled macOS/Linux binary distributed with the tool was used without recompilation. We note that fLPS 2.0 does not implement the `-u` flag for maximum unique residues used in some earlier versions; accordingly, we retained only SINGLE-type records (single-residue compositional bias) in post-processing, achieving equivalent filtering to the approach of Teekas et al. (2024) who additionally applied a maximum unique residue filter. All other parameters were identical: minimum LCR length = 3 amino acids (`-m 3`). A purity filter of ≥70% (dominant amino acid count / LCR length) was applied in post-processing.

### Positional binning and enrichment statistics

Each LCR was assigned to one of 20 equal positional bins by: bin = ⌊(midpoint / protein_length) × 20⌋ + 1, capped at 20. Terminal LCRs were defined as bins 1 or 20. Species-level terminal enrichment was assessed by one-sided Fisher's exact test comparing observed versus expected terminal LCR counts (null: 10%, i.e. 2 of 20 bins). A protein-level sensitivity analysis was also performed: for each species, the fraction of proteins bearing at least one terminal LCR versus at least one internal-only LCR was tested against the null fraction (2/20 = 10%) using a one-sided binomial test. Phylum-level summaries pool all LCRs from member species before applying Fisher's test. All analyses used Python 3.10, scipy, and pandas.

### N-terminal versus C-terminal asymmetry

Bin 1 and bin 20 were analysed separately against a null of 5% (1/20 bins) using one-sided Fisher's exact test. A Wilcoxon signed-rank test across all 61 metazoan species assessed whether N- and C-terminal enrichment magnitudes differed significantly.

### Protein-length confound test

LCRs were stratified into four quartiles by host protein length (global distribution across all 61 metazoan species). Fisher's exact test for terminal enrichment was applied within each quartile, pooled across species and per phylum.

### Amino acid identity

The dominant amino acid (column 8 of fLPS 2.0 SINGLE-type output, in {X} format) was extracted for each LCR. Amino acid frequencies were computed separately for terminal and internal LCRs, pooled across all species. Chi-squared tests with Bonferroni correction were applied to test whether individual amino acids were significantly over- or under-represented at terminal versus internal positions.

### Driver analysis

Proteins were classified as singleton-LCR (1 LCR per protein) or multi-LCR (≥2). Terminal enrichment was computed and tested by Fisher's exact test for each class per phylum.

### Purity gradient analysis

Mean LCR purity was compared between terminal and internal LCRs per phylum using the Mann-Whitney U test (one-sided, alternative: terminal purity > internal purity).

### GO term enrichment analysis

To ask whether specific protein functional categories preferentially carry terminal LCRs, GO term enrichment analysis was performed for all species in the pipeline for which Ensembl Metazoa BioMart annotation was available and protein identifier overlap ≥1 was achieved between BioMart peptide IDs and FASTA headers. GO data (protein → GO term mappings: biological process, molecular function, cellular component) were downloaded from the Ensembl Metazoa BioMart REST API and cached locally. For each species, the foreground was proteins with ≥1 terminal LCR; the background was all LCR-bearing proteins. GO terms annotated to fewer than 5 background proteins were excluded. Fisher's exact test (one-sided, greater) was applied per GO term, with Benjamini–Hochberg FDR correction at threshold ≤ 0.05.

### Code availability

All analysis scripts are available at [repository URL]. The pipeline is fully reproducible from raw Ensembl Metazoa downloads.

---

## Results

### Terminal LCR enrichment is near-universal across invertebrate phyla

Across 61 metazoan species and 72,189 LCRs (purity ≥70%, SINGLE-type), 60 species (98%) showed significant terminal LCR enrichment (Fisher's exact test, p < 0.05; Figures 1–2, Table 1). The single non-significant metazoan species was *Schistosoma mansoni* (Platyhelminthes; 11.2% terminal, p = 0.186). Overall, 17.9% of all LCRs (12,954/72,189) fell in terminal bins, significantly above the 10% null (p < 10⁻¹⁰⁰).

At the phylum level, terminal enrichment ranged from 12.7% (Platyhelminthes) to 21.2% (Nematoda and Brachiopoda; Table 1). All 16 phyla exceeded the 10% null, and 12 of 16 phyla fell within or above the Tetrapoda range of 15–25% reported by Teekas et al. (2024). Between-phylum variation in species-level pct_terminal was modest; within-phylum variation was low for well-sampled phyla (CV: Cnidaria 1.9%, Echinodermata 4.7%, Insecta 10.9%, Crustacea 13.7%), indicating that the result is phylogenetically robust rather than driven by outlier species.

**Protein-level sensitivity analysis.** Results were qualitatively unchanged when the unit of analysis was shifted from individual LCRs to proteins: the fraction of proteins bearing at least one terminal LCR significantly exceeded the 10% null in all tested species (binomial test; results in Supplementary Table S4), confirming that the finding is not an artefact of treating non-independent LCRs within the same protein as independent observations.

**Table 1. Phylum-level terminal LCR enrichment.**

| Phylum | N spp. | Total LCRs | % Terminal | vs. Teekas (15–25%) |
|---|---|---|---|---|
| Ctenophora† | 1 | 542 | 20.5 | Within range |
| Placozoa† | 1 | 234 | 18.8 | Within range |
| Cnidaria | 5 | 3,948 | 19.5 | Within range |
| Platyhelminthes | 2 | 1,676 | 12.7 | ↓ Below range |
| Annelida | 3 | 5,617 | 16.3 | Within range |
| Nematoda | 8 | 6,730 | 20.1 | Within range |
| Priapulida† | 1 | 1,273 | 14.0 | Marginally below |
| Brachiopoda† | 1 | 1,250 | 21.2 | Within range |
| Mollusca | 7 | 9,019 | 16.8 | Within range |
| Crustacea | 6 | 12,404 | 18.3 | Within range |
| Myriapoda† | 1 | 699 | 14.0 | Marginally below |
| Chelicerata | 6 | 6,330 | 18.8 | Within range |
| Insecta | 12 | 14,862 | 17.5 | Within range |
| Echinodermata | 5 | 5,687 | 19.5 | Within range |
| Hemichordata† | 1 | 797 | 18.8 | Within range |
| Cephalochordata† | 1 | 1,121 | 15.8 | Within range |

†Single-species phyla; estimates carry higher uncertainty.

*Trichoplax adhaerens* (Placozoa; 18.8%, OR = 2.11, p = 0.004) establishes terminal LCR enrichment as pre-bilaterian, present in the most basally-branching animal phylum diverging ≥700 Mya.

*Mnemiopsis leidyi* (Ctenophora; 20.5%, OR = 2.32, p = 1×10⁻⁶) represents an independent test of the pre-bilaterian hypothesis under the alternative phylogenetic placement of Ctenophora as the most basal animal phylum; its significant enrichment is consistent with pan-animal conservation regardless of which basal-animal topology is correct.

### The positional profile is U-shaped and conserved

Plotting LCR fraction per bin across all 20 positions reveals a consistent U-shaped profile across all phyla: density peaks in bins 1 and 20 and is depressed in internal bins (Figure 3; Supplementary Figure S2). This shape conservation is qualitatively identical to Tetrapoda, and is most pronounced in Nematoda and Cnidaria.

### C-terminal enrichment is more universal than N-terminal enrichment

Separating bin 1 from bin 20 reveals an asymmetry not previously reported (Figure 4). Testing each terminus independently against a 5% null (1/20 bins): C-terminal enrichment is significant in 47/50 species (94%); N-terminal enrichment is significant in 43/50 species (86%). A Wilcoxon signed-rank test across all species confirms that C-terminal enrichment is significantly greater than N-terminal enrichment (p < 0.05). The mean asymmetry ratio (N-term count / C-term count) across species is 0.91, indicating that C-terminal LCRs are consistently more numerous at the phylum scale.

The most extreme case is *Schistosoma mansoni*, which completely lacks N-terminal enrichment (pct_nterm = 3.6%, p = 0.95) but retains C-terminal enrichment (pct_cterm = 7.6%, p = 0.007). This reframes the Schistosoma exception: it represents selective loss of N-terminal LCR accumulation in a parasitic flatworm with a highly reorganised proteome (Berriman et al. 2009), while C-terminal bias is maintained. The opposite extreme is *Pristionchus pacificus* (Nematoda), showing the strongest N-terminal excess (asymmetry ratio = 1.99; pct_nterm = 13.8% vs. pct_cterm = 7.0%).

### Terminal enrichment is independent of protein length

To rule out the artefact by which short proteins mechanically inflate terminal bin counts, we stratified all LCRs by protein length quartile (Q1: ≤199 aa; Q2: 199–339 aa; Q3: 339–588 aa; Q4: ≥588 aa). Terminal enrichment was significant within all four quartiles when pooled: Q1 = 16.2% (OR = 1.74), Q2 = 19.9% (OR = 2.24), Q3 = 18.4% (OR = 2.03), Q4 = 17.9% (OR = 1.96; all p < 10⁻¹⁰⁰; Supplementary Figure S3). This rules out protein length as a confounding factor. In per-phylum stratifications, enrichment was significant in at least two of four quartiles in all phyla with sufficient LCR counts.

### Amino acid composition of terminal LCRs

To identify which amino acids drive terminal LCR enrichment, we computed the terminal-to-internal frequency ratio for each amino acid, pooled across all 61 metazoan species. After Bonferroni correction, cysteine (C), glutamic acid (E), and lysine (K) were significantly over-represented in terminal LCRs relative to internal LCRs (chi-squared test, adjusted p < 0.05). Glutamine (Q) and asparagine (N) — the amino acids most associated with prion-like domains and cytoplasmic condensates — did not show significant terminal enrichment over internal positions at the global level, suggesting that the terminal LCR phenomenon is not primarily driven by phase-separating prion-like sequences. Cysteine enrichment at termini is consistent with C-terminal cysteine-rich degrons and N-terminal metal-binding or reactive-cysteine motifs. Full per-phylum amino acid composition is presented in Supplementary Figures S4 and Table S3.

### Terminal LCRs are not purer than internal LCRs

Mean LCR purity (dominant amino acid fraction) did not significantly differ between terminal and internal LCRs in 12 of 16 phyla (Mann-Whitney U, p > 0.05). Only Annelida (Δpurity = +0.014, p = 0.002), Chelicerata (Δpurity = +0.011, p = 0.003), Insecta (Δpurity = +0.009, p < 0.001), and Echinodermata (Δpurity = +0.014, p = 0.002) showed significantly higher terminal purity, with small effect sizes. Terminal LCR enrichment is therefore predominantly positional — a greater density of LCRs at termini — rather than qualitative.

### Singleton-LCR proteins drive terminal enrichment

Proteins carrying a single LCR (singleton-LCR) showed significant terminal enrichment in all 16 phyla, with pct_terminal values consistently matching or exceeding the phylum average. Multi-LCR proteins (≥2 LCRs per protein) showed significant enrichment in only 7 of 16 phyla, and at lower effect sizes in most cases (Supplementary Figure S6; Supplementary Table S2). This demonstrates that the signal is not a geometric consequence of protein LCR load, but reflects genuine terminal positioning of individual LCRs, even in proteins carrying just one.

### Functional annotation of terminal-LCR proteins

To ask whether specific protein functional categories preferentially carry terminal LCRs, we performed GO term enrichment analysis via Ensembl Metazoa BioMart for all pipeline species where protein identifier overlap between BioMart and FASTA headers was achieved. GO enrichment (Fisher's exact, BH FDR ≤ 0.05) was significant in five species spanning four phyla (Table S5; Figure S7).

No GO terms were significantly enriched in the majority of tested species, including the focal species *D. melanogaster*, *C. elegans*, *A. mellifera*, *A. gambiae*, and *B. lanceolatum* (protein ID overlap 47–94%). The absence of functional bias in these well-annotated model organisms confirms that terminal LCR enrichment is broadly distributed across the functional spectrum of the proteome and is not restricted to specific protein classes.

Four additional species produced significant GO enrichment:

- **Pediculus humanus** (Insecta): GO:0016021 "integral component of membrane" was enriched among terminal-LCR proteins (OR = 1.86, FDR = 0.015), suggesting a modest bias toward integral membrane proteins in this louse.

- **Tetranychus urticae** (Chelicerata): "catalytic activity" (OR = 29.7, FDR = 0.021) was enriched based on a small set of terminal-LCR proteins with GO annotation.

- **Caenorhabditis brenneri** (Nematoda): "phosphoprotein phosphatase activity" (OR = ∞, FDR = 0.012) was enriched; the small foreground (7 proteins all terminal) warrants cautious interpretation.

- **Lottia gigantea** (Mollusca): Five GO terms reached significance, all related to cytoskeletal and microtubule biology: "microtubule depolymerization" (OR = ∞, FDR = 0.013), "regulation of microtubule polymerization or depolymerization" (OR = ∞), "GTPase activity" (OR = ∞), "neuron projection development" (OR = ∞), and "microtubule" component (OR = 17.4, FDR = 0.043).

**Strongylocentrotus purpuratus** (Echinodermata) showed the strongest and most interpretable GO enrichment: eight terms related to cytoskeletal and microtubule biology reached significance — "microtubule-based process" (OR = ∞, FDR < 0.001), "structural constituent of cytoskeleton" (OR = 46.3, FDR < 0.001), "mitotic cell cycle" (OR = 16.9, FDR = 0.0001), "GTP binding" (OR = 7.4, FDR = 0.0004), and "GTPase activity" (OR = 8.4, FDR = 0.0008; Figure S7).

The overlap of cytoskeletal/microtubule GO terms between *S. purpuratus* and *L. gigantea* is notable. Tubulin C-terminal tails are intrinsically disordered glutamate-rich sequences — a well-characterised functional terminal LCR system (the "tubulin code"; Janke and Magiera 2020) — and may drive cytoskeletal enrichment in multiple animal phyla, not exclusively in echinoderms. Whether the *S. purpuratus* enrichment reflects amplified tubulin diversity in a cilia-dependent developmental programme, or simply better GO annotation coverage of cytoskeletal proteins in this model organism, cannot be resolved without systematic GO enrichment across the remaining echinoderm species in the dataset (*L. variegatus*, *L. pictus*, *P. miniata*, *A. rubens*), which had insufficient BioMart protein ID coverage in the current analysis.

### Terminal LCR enrichment extends to non-metazoan eukaryotes

Of 13 non-metazoan eukaryote species tested, 11 showed significant terminal LCR enrichment (Fisher's exact test p < 0.05): all seven fungi (*S. cerevisiae* 19.0%; *S. pombe* 28.3%; *N. crassa* 13.8%; *A. nidulans* 18.6%; *C. albicans*; *A. fumigatus*; *C. neoformans*), both land plants (*Arabidopsis thaliana* 23.6%; *Solanum lycopersicum* 19.7%), the excavate *Trypanosoma brucei* (19.3%), and the amoebozoan *Dictyostelium discoideum* (11.5%, p = 0.005). The single exception among eukaryotes was *Plasmodium falciparum* (Apicomplexa; 10.3%, p = 0.46), which was at the null expectation. Apicomplexa are known to have highly unusual amino acid composition (asparagine- and lysine-enriched repeats due to AT-biased genome) and may represent a genuine exception driven by their atypical LCR biology.

*Solanum lycopersicum* (land plant; 19.7%, p < 10⁻¹⁰⁰) corroborates the *Arabidopsis thaliana* finding, confirming terminal enrichment is conserved across land plants. Notably, the green alga *Chlamydomonas reinhardtii* (Chlorophyta; 10.7%, p = 0.116) was not significant — the only non-significant non-metazoan eukaryote alongside *P. falciparum*. *Chlamydomonas* is a unicellular flagellate with cilia-based motility; the absence of enrichment in this lineage, combined with its presence in land plants, may reflect differences in proteome composition, LCR content, or evolutionary pressures on terminal sequences between green algae and their land plant descendants, and warrants further investigation.

The finding that both Opisthokonta (fungi) and Viridiplantae show terminal LCR enrichment suggests the phenomenon predates the split of these lineages approximately 1 billion years ago, placing its origin no later than the last common ancestor of all extant eukaryotes (approximately 1–1.5 billion years ago; Eme et al. 2014). These non-metazoan eukaryotes show enrichment magnitudes comparable to metazoans (11.5–28.3%), providing no evidence that metazoan biology amplified the signal relative to unicellular eukaryotes.

For prokaryotes, the data are insufficient to draw conclusions. Bacterial and archaeal proteomes yielded median 21 LCRs per species under our purity filter, compared to ≥289 for eukaryotes. With such small sample sizes, Fisher's exact tests are severely underpowered: even a species with 40% terminal LCRs by chance would require ≥35 LCRs for the test to reach p < 0.05. Three of six bacterial species showed nominally significant enrichment (*P. aeruginosa* p = 0.049, *M. tuberculosis* p = 0.005, *S. coelicolor* p < 0.001), but the extremely high pooled bacterial terminal fraction (29.7%) with only 300 total LCRs cannot be meaningfully interpreted without much larger LCR samples per species. Testing the prokaryote boundary rigorously would require either lower purity thresholds or analysis of prokaryote-dense clades where LCR prevalence is higher.

---

## Discussion

### Terminal LCR enrichment is a pan-eukaryotic property ≥1 billion years old

The finding that all 16 tested animal phyla — from the pre-bilaterian Placozoa and Ctenophora to the vertebrate-adjacent Cephalochordata — show terminal LCR enrichment establishes this as one of the most phylogenetically conserved properties of metazoan protein architecture known. The inclusion of *T. adhaerens* is particularly decisive: Placozoa represent animals with no neuromuscular junctions, no symmetry, no differentiated tissue layers, and a proteome of ~13,000 proteins (Srivastava et al. 2008). Yet the same terminal bias observed in human or chicken proteins is present in this most basal animal.

The extension to non-metazoan eukaryotes deepens this conclusion further. The presence of terminal LCR enrichment in fungi, plants, and protists — spanning the breadth of eukaryotic diversity — places the origin of this property no later than the last eukaryotic common ancestor (~1–1.5 billion years ago; Eme et al. 2014). Critically, this rules out not only vertebrate-specific explanations (adaptive immunity, complex body plans) but also any metazoan-specific explanation involving multicellularity, nervous systems, or cell differentiation. The mechanism must be compatible with unicellular eukaryotes.

Teekas et al. (2024) speculated that terminal LCR enrichment might be driven by selection pressure from adaptive immunity. Our data extend the refutation of that hypothesis back an additional 1.3 billion years. The *Plasmodium falciparum* exception (non-significant at the null) is attributable to the unusual LCR biology of Apicomplexa, whose AT-biased genome generates asparagine- and lysine-rich repeats with atypical positional behaviour, rather than an absence of the underlying mechanism.

### N/C asymmetry implicates distinct terminal functions

The separation of N- and C-terminal enrichment is the most novel analytical contribution of this study. The greater universality and magnitude of C-terminal enrichment (47/50 species significant) versus N-terminal (43/50) suggests that the two termini are subject to distinct evolutionary pressures. This is consistent with the distinct biochemical biology of protein termini. C-terminal degrons — short sequences recognised by E3 ubiquitin ligases for proteasomal targeting — are enriched at C-termini across eukaryotes (Koren et al. 2018), and LCR sequences at C-termini may constitute or facilitate such degrons. The recent finding that C-terminal tails evolve faster yet carry concentrated functional motifs (Riba et al. 2019) is consistent with C-terminal LCRs serving regulatory rather than structural roles.

N-terminal enrichment is present but more variable. The N-degron pathway (Varshavsky 2019) recognises specific N-terminal residues for proteasomal targeting, and N-terminal disordered sequences may serve as co-translational folding linkers. Importantly, if translational kinetics at start codons — where ribosome pausing is known to occur (Pechmann and Frydman 2013) — drives local disorder, this would predict N-terminal enrichment. The fact that N-terminal enrichment is present in 86% of species but weaker and more variable than C-terminal enrichment suggests it reflects a real but partially overridden signal, shaped by both translational and functional constraints.

The *Schistosoma mansoni* case is instructive. This parasitic flatworm has undergone extensive genome reorganisation, gene loss, and unusual codon usage (Berriman et al. 2009). The specific loss of N-terminal LCR enrichment — while C-terminal enrichment is maintained — argues that N-terminal compositional bias is more evolutionarily labile than C-terminal bias, and that N-terminal and C-terminal enrichment are mechanistically separable.

An intriguing exception to the predominant C-terminal bias is *Arabidopsis thaliana*, where N-terminal enrichment substantially exceeds C-terminal enrichment (pct_nterm = 15.6%, OR = 3.54 vs. pct_cterm = 8.1%, OR = 1.69). Land plants have a well-characterised N-degron pathway, and plant N-terminal sequences carry concentrated regulatory motifs for protein targeting and turnover (Gibbs et al. 2014). The reversed asymmetry in *A. thaliana* is consistent with stronger N-terminal functional constraint in land plants relative to the C-terminal bias seen in most animals and fungi, and represents a testable prediction for comparative analysis of plant terminal LCR biology.

### Mechanistic hypotheses weighted against the data

Three mechanisms could explain conserved terminal LCR enrichment, and our data allow us to weight them:

**Translational kinetics.** Ribosome pausing near start and stop codons could promote disordered, compositionally simple sequences at termini by slowing translation. This predicts symmetrical N/C enrichment. The observed C > N asymmetry argues that translational kinetics is a contributing but not sufficient explanation.

**Co-translational folding.** N-terminal regions are exposed while the C-terminal domain remains inside the ribosome exit tunnel, potentially favouring flexible, low-complexity N-terminal sequences. This predicts N-terminal enrichment. The data show N-terminal enrichment present but weaker, consistent with a partial contribution.

**Degron biology.** Both N-degron and C-degron pathways are ancient — present in yeast, plants, and animals — and recognise terminal sequences for proteasomal targeting. The stronger, more universal C-terminal enrichment is most consistent with C-degron biology. The amino acid enrichment analysis (cysteine, glutamate, lysine at termini) is consistent with known C-degron and N-degron recognition sequences: lysine at the N-terminus is a class I N-degron, and C-terminal cysteine is a recently identified C-degron residue (Koren et al. 2018).

We propose that the primary driver of C-terminal LCR enrichment is degron biology, supported by the universal conservation of C-degron pathways across metazoans and the amino acid composition of terminal LCRs. N-terminal enrichment likely reflects a combination of N-degron biology and translational kinetics, with varying contributions across lineages. This mechanistic model predicts that LCRs containing or adjacent to functionally characterised degron motifs would be disproportionately terminal — a testable hypothesis for future experimental work.

The GO enrichment results provide independent support for this mechanistic interpretation. A mechanism acting on protein termini regardless of protein function — such as degron recognition, which targets terminal sequences for proteasomal processing irrespective of the host protein's cellular role — predicts that terminal LCR enrichment should be distributed uniformly across the functional proteome, not concentrated in specific GO categories. This is precisely what the data show: no GO terms were enriched among terminal-LCR proteins in the well-annotated model species *D. melanogaster*, *C. elegans*, *A. mellifera*, *A. gambiae*, and *B. lanceolatum*, despite these species providing thousands of annotated proteins and ample statistical power to detect a functional bias if one existed. This contrasts sharply with phenomena that are genuinely function-specific: prion-like domain enrichment is strongly concentrated in RNA-binding proteins and transcription factors (Lancaster et al. 2014), and phase-separation-prone LCRs are enriched in stress granule and P-body components (Alberti et al. 2019). The absence of any such functional concentration in terminal LCRs is therefore not a limitation of the GO analysis — it is the expected signature of a mechanism that operates at the level of protein architecture rather than protein function. The *S. purpuratus* cytoskeletal signal, and the similar microtubule-related enrichment in *L. gigantea*, represent lineage-specific amplification of terminal LCRs in proteins with unusually expanded tubulin repertoires, superimposed on this universal background.

### The Platyhelminthes outlier

Platyhelminthes show the lowest phylum-level enrichment (12.7%), driven by *Schistosoma mansoni* (11.2%, not significant). *Schmidtea mediterranea*, the second flatworm in the dataset, shows 15.6% terminal enrichment (p = 0.003), firmly within the Tetrapoda range. The phylum average is therefore strongly influenced by the unusual biology of a parasitic species. Future studies including free-living flatworm species would better characterise this phylum.

### Limitations

This study is entirely computational. The purity threshold (≥70%) and SINGLE-type filter are definitional choices that affect absolute LCR counts but not the qualitative finding. Proteome quality varies across species. Single-species phyla carry high statistical uncertainty at the phylum level. The analysis does not incorporate functional annotation, and the mechanistic model is hypothetical pending experimental validation.

Prokaryote sampling is underpowered: bacteria and archaea carry few LCRs under the purity filter (~6–125 per species), making individual-species Fisher's tests uninterpretable for most species. The evolutionary boundary between eukaryotes and prokaryotes therefore remains unresolved by this dataset. Resolution would require either (a) lowering the purity threshold to capture the lower-complexity LCRs prevalent in prokaryotic proteomes, or (b) selecting prokaryotes with known high LCR content (e.g. Actinobacteria with proline- or alanine-rich repeats). The GO enrichment analysis produced significant results in five species; the remaining majority showed no functional enrichment, confirming the broad distribution of terminal LCRs across functional categories. Protein ID format mismatches between BioMart peptide identifiers and Ensembl Metazoa FASTA headers prevented enrichment testing in many species; systematic resolution of this ID mapping is needed for a comprehensive cross-species GO comparison. In particular, GO enrichment could not be completed for four of the five echinoderm species in the dataset (*L. variegatus*, *L. pictus*, *P. miniata*, *A. rubens*), precluding a direct test of whether the cytoskeletal enrichment observed in *S. purpuratus* is specific to echinoids or shared across Echinodermata. The choanoflagellates — the immediate unicellular relatives of animals (*Monosiga brevicollis*, *Salpingoeca rosetta*) — are absent from the current outgroup set. Adding choanoflagellate proteomes would resolve whether terminal LCR enrichment predates the origin of multicellularity in the animal lineage.

---

## Conclusions

Terminal LCR enrichment is a conserved pan-metazoan property present in 60 of 61 analysed metazoan species and all 16 tested phyla, with magnitudes overlapping the Tetrapoda baseline. The mechanism predates the Cambrian explosion (≥700 Mya) and cannot be explained by vertebrate-specific biology. Extension to 11 of 13 non-metazoan eukaryotes — including fungi, plants, and protists — establishes this as a feature of eukaryotic protein architecture conserved for at least ~1–2 billion years, predating the origin of animals. The most parsimonious mechanism consistent with the data is the deep conservation of protein terminal biology: particularly the eukaryotic ubiquitin-proteasome system's C-terminal and N-terminal degron pathways, which are present in all eukaryotes and absent in most prokaryotes. C-terminal enrichment is more universal than N-terminal enrichment across species. These findings extend the evolutionary framework of Teekas et al. (2024) from Tetrapoda to all animals and to non-metazoan eukaryotes (11 of 13 tested species significant), placing the origin of this property at or near the base of the eukaryotic tree.

---

## Acknowledgements

The author thanks the Ensembl Metazoa project for providing freely downloadable invertebrate proteomes, and Paul Harrison for making fLPS 2.0 freely available.

---

## References

Alberti S, Gladfelter A, Mittag T (2019). Considerations and challenges in studying liquid-liquid phase separation and biomolecular condensates. *Cell* 176:419–434. https://doi.org/10.1016/j.cell.2018.12.035

Alanis-Lobato G, Cannistraci CV, Ravasi T (2016). Exploitation of the low complexity regions in proteins: bioinformatics approaches to identify and characterize them. *Briefings in Bioinformatics* 17:648–657.

Altschul SF, et al. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Research* 25:3389–3402. https://doi.org/10.1093/nar/25.17.3389

Berriman M, et al. (2009). The genome of the blood fluke *Schistosoma mansoni*. *Nature* 460:352–358. https://doi.org/10.1038/nature08160

Boija A, et al. (2018). Transcription factors activate genes through the phase-separation capacity of their activation domains. *Cell* 175:1842–1855. https://doi.org/10.1016/j.cell.2018.10.042

Buljan M, et al. (2012). Tissue-specific splicing of disordered segments that embed binding motifs rewires protein interaction networks. *Molecular Cell* 46:871–883. https://doi.org/10.1016/j.molcel.2012.05.039

Coletta A, et al. (2010). Low-complexity regions within protein sequences have position-dependent roles. *BMC Systems Biology* 4:43. https://doi.org/10.1186/1752-0509-4-43

Cunningham F, et al. (2022). Ensembl 2022. *Nucleic Acids Research* 50:D988–D995. https://doi.org/10.1093/nar/gkab1049

Dosztanyi Z, et al. (2005). IUPred: web server for the prediction of intrinsically unstructured regions of proteins based on estimated energy content. *Bioinformatics* 21:3433–3434. https://doi.org/10.1093/bioinformatics/bti541

Eme L, Sharpe SC, Brown MW, Roger AJ (2014). On the age of eukaryotes: evaluating evidence from fossils and molecular clocks. *Cold Spring Harbor Perspectives in Biology* 6:a016139. https://doi.org/10.1101/cshperspect.a016139

Dunker AK, et al. (2001). Intrinsically disordered protein. *Journal of Molecular Graphics and Modelling* 19:26–59.

Finn RD, et al. (2016). The Pfam protein families database: towards a more sustainable future. *Nucleic Acids Research* 44:D279–D285. https://doi.org/10.1093/nar/gkv1344

Gibbs DJ, et al. (2014). Homeostatic response to hypoxia is regulated by the N-end rule pathway in plants. *Nature* 479:415–418. https://doi.org/10.1038/nature10534

Haerty W, Golding GB (2010). Genome-wide evidence for selection acting on single amino acid repeats. *Genome Research* 20:755–760. https://doi.org/10.1101/gr.101246.109

Harrison PM (2017). fLPS: Fast discovery of compositional biases for the protein universe. *BMC Bioinformatics* 18:476. https://doi.org/10.1186/s12859-017-1906-3

Janke C, Magiera MM (2020). The tubulin code and its role in controlling microtubule properties and functions. *Nature Reviews Molecular Cell Biology* 21:307–326. https://doi.org/10.1038/s41580-020-0214-3

Huntley MA, Golding GB (2000). Evolution of simple sequence in proteins. *Journal of Molecular Evolution* 51:131–140.

Jorda J, Bharat Bharat T (2015). Common topology of globular domains in proteins: convergent evolution? *Trends in Biochemical Sciences* 40:11–13.

Koren I, et al. (2018). The eukaryotic proteome is shaped by E3 ubiquitin ligases targeting C-terminal degrons. *Cell* 173:1622–1635. https://doi.org/10.1016/j.cell.2018.04.028

Lancaster AK, et al. (2014). PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition. *Bioinformatics* 30:2–3. https://doi.org/10.1093/bioinformatics/btt421

Linding R, et al. (2003). GlobPlot: exploring protein sequences for globularity and disorder. *Nucleic Acids Research* 31:3701–3708. https://doi.org/10.1093/nar/gkg519

Marcotte EM, et al. (1999). A census of protein repeats. *Journal of Molecular Biology* 293:151–160. https://doi.org/10.1006/jmbi.1999.3136

Mier P, et al. (2020). Disentangling the complexity of low-complexity proteins. *Briefings in Bioinformatics* 21:458–472. https://doi.org/10.1093/bib/bbz007

Moesa HA, et al. (2012). Chemical composition is maintained in poorly conserved intrinsically disordered regions and suggests functional convergence. *Molecular BioSystems* 8:3262–3273. https://doi.org/10.1039/c2mb25232d

Ntountoumi C, et al. (2019). Low complexity regions in the proteins of prokaryotes perform important functional roles and are highly conserved. *Nucleic Acids Research* 47:9998–10009. https://doi.org/10.1093/nar/gkz730

Pechmann S, Frydman J (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology* 20:237–243. https://doi.org/10.1038/nsmb.2466

Riba A, et al. (2019). Protein synthesis rates and ribosome occupancies reveal determinants of translation elongation rates. *Proceedings of the National Academy of Sciences* 116:15023–15032. https://doi.org/10.1073/pnas.1817299116

Romero P, et al. (2001). Sequence complexity of disordered protein. *Proteins: Structure, Function, and Bioinformatics* 42:38–48. https://doi.org/10.1002/1097-0134(20010101)42:1<38::AID-PROT50>3.0.CO;2-3

Rubin GM, et al. (2000). Comparative genomics of the eukaryotes. *Science* 287:2204–2215. https://doi.org/10.1126/science.287.5461.2204

Sea Urchin Genome Sequencing Consortium (2006). The genome of the sea urchin *Strongylocentrotus purpuratus*. *Science* 314:941–952. https://doi.org/10.1126/science.1133609

Shin Y, Brangwynne CP (2017). Liquid phase condensation in cell physiology and disease. *Science* 357:eaaf4382. https://doi.org/10.1126/science.aaf4382

Srivastava M, et al. (2008). The *Trichoplax* genome and the nature of placozoans. *Nature* 454:955–960. https://doi.org/10.1038/nature07191

Teekas L, Sharma S, Vijay N (2024). Terminal regions of a protein are a hotspot for low complexity regions and selection. *Open Biology* 14:230439. https://doi.org/10.1098/rsob.230439

Toll-Riera M, et al. (2012). Origin of primate orphan genes: a comparative genomics approach. *Molecular Biology and Evolution* 26:603–612. https://doi.org/10.1093/molbev/msn281

UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research* 51:D523–D531. https://doi.org/10.1093/nar/gkac1052

van der Lee R, et al. (2014). Classification of intrinsically disordered regions and proteins. *Chemical Reviews* 114:6589–6631. https://doi.org/10.1021/cr400525m

Varshavsky A (2019). N-degron and C-degron pathways of protein degradation. *Proceedings of the National Academy of Sciences* 116:358–366. https://doi.org/10.1073/pnas.1816596116

Vernon RM, et al. (2018). Pi-Pi contacts are an overlooked protein feature relevant to phase separation. *eLife* 7:e31486. https://doi.org/10.7554/eLife.31486

Virtanen P, et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17:261–272. https://doi.org/10.1038/s41592-019-0686-2

Ward JJ, et al. (2004). The DISOPRED server for the prediction of protein disorder. *Bioinformatics* 20:2138–2139. https://doi.org/10.1093/bioinformatics/bth195

Wootton JC, Federhen S (1996). Analysis of compositionally biased regions in sequence databases. *Methods in Enzymology* 266:554–571. https://doi.org/10.1016/S0076-6879(96)66035-2

Xue B, et al. (2012). PONDR-FIT: a meta-predictor of intrinsically disordered amino acids. *Biochimica et Biophysica Acta* 1804:996–1010. https://doi.org/10.1016/j.bbapap.2010.01.011

---

## Supplementary Material

**Supplementary Table S1.** Complete species list: Ensembl Metazoa directory key, phylum, genome assembly, download date, longest-isoform protein count, and per-species enrichment statistics (n_lcr, n_terminal, pct_terminal, odds_ratio, p-value, significant).

**Supplementary Table S2.** Driver analysis: pct_terminal, odds ratio, and Fisher's exact p-value for singleton-LCR and multi-LCR proteins per phylum.

**Supplementary Table S3.** Amino acid composition of terminal vs. internal LCRs per phylum (count and fraction for each amino acid in each location class).

**Supplementary Table S4.** Protein-level sensitivity analysis: fraction of proteins with ≥1 terminal LCR per species, binomial test p-value.

**Supplementary Table S5.** GO term enrichment results: all species where BioMart annotation was available (FDR ≤ 0.05 entries only). Columns: go_id, go_name, namespace, n_terminal, n_all_lcr, pct_terminal, odds_ratio, pvalue, fdr, species_key, phylum.

**Supplementary Figure S1.** Heatmap of LCR positional distribution across 20 bins for all 61 metazoan species (species × bin, colour = fraction of LCRs).

**Supplementary Figure S2.** U-shaped bin profile per phylum (mean fraction of LCRs per bin averaged across species).

**Supplementary Figure S3.** Terminal LCR enrichment across protein-length quartiles (pooled and per phylum).

**Supplementary Figure S4.** Amino acid composition of terminal vs. internal LCRs: (a) stacked bar chart per phylum; (b) global terminal-to-internal enrichment ratio per amino acid.

**Supplementary Figure S5.** LCR purity gradient: violin plot (terminal vs. internal, global) and Δ purity per phylum.

**Supplementary Figure S6.** Driver analysis: terminal enrichment in singleton-LCR vs. multi-LCR proteins per phylum.
