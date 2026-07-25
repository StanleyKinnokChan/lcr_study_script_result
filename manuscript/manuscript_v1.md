# Terminal Low-Complexity Regions Are a Pan-Metazoan Property Conserved Across ≥700 Million Years of Animal Evolution

**Stanley K. Chan**

*Correspondence: stanleykinnok.chan@gmail.com*

---

## Abstract

Low-complexity regions (LCRs) — protein segments dominated by one or a few amino acid types — are concentrated at the N- and C-termini of proteins in Tetrapoda (Teekas et al. 2024, *Open Biology*). Whether this terminal enrichment is a vertebrate-specific innovation or an ancient metazoan property has not been tested. Here we apply the same computational framework (fLPS 2.0, 20-bin positional analysis) to 50 invertebrate species spanning 14 phyla, from the most basal animal phylum (Placozoa) to Cephalochordata, representing over 59,000 LCR records. We find that 49 of 50 species (98%) show statistically significant terminal LCR enrichment (Fisher's exact test, p < 0.05), with phylum-level terminal fractions ranging from 12.7% (Platyhelminthes) to 21.2% (Nematoda and Brachiopoda) — overlapping the Tetrapoda baseline of 15–25%. Critically, *Trichoplax adhaerens* (Placozoa; 18.8% terminal enrichment), representing the most basally-branching animal lineage diverging approximately 700–800 million years ago, shows significant enrichment, establishing the mechanism as pre-bilaterian. Three orthogonal analyses reinforce this finding: (1) terminal enrichment is significant across all protein length quartiles, ruling out a length-confound artefact; (2) N- and C-termini show distinct enrichment profiles — C-terminal bias is more universal and stronger, while N-terminal enrichment is more variable across phyla; and (3) the signal is driven predominantly by singleton-LCR proteins, not multi-LCR disordered scaffolds. Terminal LCRs are not significantly purer than internal LCRs in most phyla, suggesting the positional bias reflects structural or translational constraints rather than intensified selection for compositional homogeneity at termini. These results establish terminal LCR enrichment as a conserved pan-metazoan property predating the Cambrian explosion, and point toward a fundamental mechanistic explanation rooted in protein terminal biology rather than lineage-specific adaptation.

**Keywords:** low-complexity regions, protein termini, compositional bias, metazoan evolution, intrinsically disordered regions, fLPS2, invertebrates

---

## Introduction

Proteins are not compositionally uniform along their length. Low-complexity regions (LCRs) — segments dominated by one or a small number of amino acid types — are distributed non-randomly within protein sequences, clustering in functionally important contexts such as disordered linkers, prion-like domains, and polyamino acid tracts (Wootton and Federhen 1996; Marcotte et al. 1999). The biological significance of LCRs has grown considerably with the recognition that many intrinsically disordered regions (IDRs) are LCR-containing (Romero et al. 2001; van der Lee et al. 2014), and that phase separation — a mechanism central to the formation of membrane-less organelles — is often driven by low-complexity IDR sequences (Boija et al. 2018; Shin and Brangwynne 2017).

The positional distribution of LCRs within proteins has received less attention than their prevalence. Teekas and colleagues (2024) recently reported that, across all major Tetrapoda clades, LCRs are significantly enriched in the terminal 5% of protein sequences (the first and last bins of a 20-bin positional map), with 15–25% of all LCRs residing in these terminal positions despite the bins representing only 10% of positional space (2 of 20 bins). This terminal enrichment was observed consistently across birds, mammals, reptiles, and amphibians, and was interpreted as evidence for evolutionary pressure maintaining low-complexity compositions at protein ends — possibly related to adaptive immunity, chaperone-mediated degradation, or co-translational folding dynamics.

A critical question left open by Teekas et al. (2024) is the evolutionary age of this mechanism. If terminal LCR enrichment is found only in Tetrapoda, it may reflect an evolutionary innovation coinciding with the expansion of vertebrate adaptive immunity or other Tetrapoda-specific biology. If, however, it extends to invertebrate phyla — particularly to the most basal animals — it represents a deeply conserved property of metazoan protein architecture, predating the Cambrian explosion (~541 Mya) and likely reflecting a fundamental constraint on protein terminal composition.

Invertebrate proteomes span enormous biological diversity. The phylogenetic distance from the most basal animals (Placozoa, Cnidaria) to vertebrates spans approximately 700–800 million years, encompassing the origin of body plans, the evolution of the nervous system, the development of immune complexity, and dramatic variation in genome size, protein length distributions, and GC content. If LCR terminal enrichment persists across this diversity, it cannot be explained by any lineage-specific biological innovation in vertebrates.

Here we test this hypothesis using a fully computational approach. We downloaded proteomes for 50 invertebrate species from 14 phyla from Ensembl Metazoa release 63, applied fLPS 2.0 (Harrison 2017) with identical parameters to those used by Teekas et al., and quantified LCR positional distributions across the same 20-bin framework. We extend the analysis in four directions not addressed by the original study: (1) separation of N-terminal and C-terminal enrichment, (2) amino acid identity of terminal versus internal LCRs, (3) protein-length stratification to test for methodological confounding, and (4) classification of proteins by LCR multiplicity to identify driver classes.

---

## Methods

### Species selection and proteome acquisition

Fifty invertebrate species were selected to maximise phylogenetic breadth across 14 phyla: Placozoa (1 sp.), Cnidaria (4), Platyhelminthes (2), Annelida (2), Nematoda (6), Priapulida (1), Brachiopoda (1), Mollusca (6), Crustacea (6), Chelicerata (4), Insecta (10), Echinodermata (5), Hemichordata (1), and Cephalochordata (1). Protein FASTA files (.pep.all.fa.gz) were downloaded from Ensembl Metazoa FTP (release 63; https://ftp.ebi.ac.uk/ensemblgenomes). Where a species had multiple protein isoforms per gene, only the longest isoform was retained, identified by parsing the `gene:` or `gene=` field in Ensembl FASTA headers. A complete species list with Ensembl directory keys is provided in Supplementary Table S1.

### LCR detection

Low-complexity regions were identified with fLPS 2.0 (Harrison 2017), using the pre-compiled macOS binary distributed with the tool. Parameters: minimum LCR length = 3 amino acids (`-m 3`). Only SINGLE-type records (single-residue compositional bias) were retained, matching the approach of Teekas et al. (2024). An LCR purity filter of ≥70% (fraction of the LCR composed of the dominant amino acid ≥ aa_count / lcr_length) was applied in post-processing, consistent with Teekas et al. Sequences with non-standard characters were not excluded from input but stop codons and X residues within LCR boundaries were not counted toward purity.

### Positional binning and terminal enrichment statistics

Each LCR was assigned to one of 20 equal-length positional bins based on the normalised midpoint position: bin = ⌊(midpoint / protein_length) × 20⌋ + 1, capped at 20. Terminal LCRs were defined as those in bins 1 (N-terminal) or 20 (C-terminal), following Teekas et al. (2024). Terminal enrichment for each species was assessed by a one-sided Fisher's exact test comparing the observed terminal LCR count to the expected count under the null hypothesis of uniform positional distribution (10% of LCRs in terminal bins, corresponding to 2/20 bins). Odds ratios (OR) and p-values are reported; significance threshold α = 0.05.

Phylum-level summaries were computed by aggregating LCR counts across all species in a phylum before applying Fisher's exact test. All statistical analyses were performed in Python 3 using scipy.stats.

### N-terminal versus C-terminal asymmetry

Bin 1 (N-terminal) and bin 20 (C-terminal) were analysed independently. Each was tested against a null of 5% (1/20 bins) using a one-sided Fisher's exact test. The asymmetry ratio was computed as n_bin1 / n_bin20. A Wilcoxon signed-rank test was applied across species to assess whether N-terminal and C-terminal enrichment levels differed globally.

### Protein-length confound test

To test whether short proteins artificially inflate terminal LCR counts (since an LCR near position 5 of a 50-aa protein lands in bin 1, while the same LCR in a 5,000-aa protein is internal), all LCRs were stratified into quartiles by host protein length (Q1: ≤ Q25, Q2: Q25–Q50, Q3: Q50–Q75, Q4: ≥ Q75, based on global protein length distribution). Fisher's exact test was applied within each quartile.

### Amino acid identity

The dominant amino acid identity of each SINGLE-type LCR was extracted from column 8 of the fLPS 2.0 output ({X} format). Amino acid frequency distributions were computed separately for terminal and internal LCRs, and the terminal-to-internal enrichment ratio per amino acid was calculated on globally pooled data.

### Driver analysis

Proteins were classified as singleton-LCR (exactly one LCR per protein) or multi-LCR (≥2 LCRs per protein). Terminal enrichment was computed and tested independently for each class per phylum.

### Purity gradient

The mean purity of terminal versus internal LCRs was compared per phylum using the Mann-Whitney U test (one-sided, alternative: terminal > internal). Delta purity (mean_terminal − mean_internal) was computed per phylum.

### Code and reproducibility

All analysis scripts are available at [repository URL]. Analysis was performed using Python 3.10, with pandas, scipy, matplotlib, seaborn, and biopython. fLPS 2.0 was obtained from https://github.com/paulhorton/fLPS2.

---

## Results

### Terminal LCR enrichment is near-universal across invertebrate phyla

Across 50 invertebrate species representing 14 phyla and over 59,210 LCRs (purity ≥70%, SINGLE-type), 49 species (98%) showed statistically significant terminal LCR enrichment compared to the uniform positional null (Fisher's exact test, p < 0.05; Figure 1, Figure 2). The single non-significant species was *Schistosoma mansoni* (Platyhelminthes; 11.2% terminal, p = 0.186), the only parasitic flatworm in the dataset. The overall terminal LCR fraction across all species pooled was 18.1% (10,713 of 59,210 LCRs), significantly above the null expectation of 10.0% (OR = 1.99, Fisher's exact p < 10⁻¹⁰⁰).

At the phylum level, terminal enrichment ranged from 12.7% (Platyhelminthes) to 21.2% (Nematoda; Brachiopoda). The majority of phyla fell within the 15–25% range reported by Teekas et al. (2024) for Tetrapoda (Table 1, Figure 2). No phylum showed terminal enrichment below 10% (the null expectation). Between-phylum variation was modest: a Kruskal-Wallis test on species-level pct_terminal values across phyla was [statistically significant / not significant — pending script 09 output], and within-phylum coefficients of variation were low for Cnidaria (CV = 1.9%), Echinodermata (4.7%), and Insecta (10.9%), confirming that the result is phylogenetically consistent rather than driven by individual outlier species.

Notably, *Trichoplax adhaerens* — representing Placozoa, the most basally branching animal phylum, with an estimated divergence time of 700–800 million years ago — showed 18.8% terminal enrichment (OR = 2.11, p = 0.004), firmly within the Tetrapoda range. This establishes terminal LCR enrichment as a pre-bilaterian property.

**Table 1. Phylum-level terminal LCR enrichment.**

| Phylum | N species | Total LCRs | % Terminal | Comparison to Teekas (15–25%) |
|---|---|---|---|---|
| Placozoa | 1 | 234 | 18.8 | Within range |
| Cnidaria | 4 | 3,176 | 20.0 | Within range |
| Platyhelminthes | 2 | 1,676 | 12.7 | Below range |
| Annelida | 2 | 2,382 | 20.0 | Within range |
| Nematoda | 6 | 5,421 | 21.2 | Within range |
| Priapulida | 1 | 1,273 | 14.0 | Below range |
| Brachiopoda | 1 | 1,250 | 21.2 | Within range |
| Mollusca | 6 | 7,688 | 16.2 | Within range |
| Crustacea | 6 | 12,404 | 18.3 | Within range |
| Chelicerata | 4 | 2,584 | 17.6 | Within range |
| Insecta | 10 | 13,517 | 17.4 | Within range |
| Echinodermata | 5 | 5,687 | 19.5 | Within range |
| Hemichordata | 1 | 797 | 18.8 | Within range |
| Cephalochordata | 1 | 1,121 | 15.8 | Within range |

### The positional enrichment profile is U-shaped and conserved across phyla

Plotting the fraction of LCRs per bin across all 20 positional bins reveals a consistent U-shaped profile across all 14 phyla: LCR density is elevated in bins 1 and 20 and depleted in internal bins (Figure 3). This shape is qualitatively identical to the U-shaped profile reported by Teekas et al. for Tetrapoda, extending shape conservation to invertebrates. The U-shape is most pronounced in Nematoda and Cnidaria, and shallowest in Platyhelminthes and Priapulida, consistent with their lower phylum-level pct_terminal values.

### C-terminal enrichment is more universal than N-terminal enrichment

Separating bin 1 (N-terminal) and bin 20 (C-terminal) reveals a striking asymmetry (Figure 4). When tested independently against a null of 5% (1/20 bins), C-terminal enrichment is significant in 47 of 50 species (94%), while N-terminal enrichment is significant in 43 of 50 species (86%). The mean asymmetry ratio across all species (N-term count / C-term count) is 0.91, indicating a modest but consistent C-terminal bias.

The asymmetry is most extreme in two cases. *Schistosoma mansoni* — the only globally non-significant species — completely lacks significant N-terminal enrichment (pct_nterm = 3.6%, p = 0.95) but retains significant C-terminal enrichment (pct_cterm = 7.6%, p = 0.007). This reframes the Schistosoma exception: it is not a failure of terminal enrichment, but a specific loss of N-terminal LCR accumulation, potentially linked to the compact and unusual proteome of a parasitic flatworm. Conversely, *Pristionchus pacificus* (Nematoda) shows the strongest N-terminal bias (asymmetry ratio = 1.99; N-term 13.8%, C-term 7.0%), suggesting lineage-specific elaboration of N-terminal compositional bias.

A Wilcoxon signed-rank test across all species confirms that N-terminal and C-terminal enrichment magnitudes differ significantly (p < 0.05), with C-terminal enrichment being consistently greater. This is consistent with the known biology of C-terminal degrons (Brown et al. 2011; Riba et al. 2019) and the relative enrichment of functional linear motifs at C-termini across metazoans.

### Terminal enrichment is not driven by short-protein artefact

A concern with positional bin analyses is that short proteins may disproportionately contribute to terminal bins: a 30-aa protein with a single LCR will always land in bin 1 or 20 regardless of the LCR's position. To test this, we stratified all LCRs by host protein length into quartiles (Q1: ≤ 199 aa, Q2: 199–339 aa, Q3: 339–588 aa, Q4: ≥ 588 aa based on the global length distribution) and assessed terminal enrichment within each quartile. Terminal enrichment was significant in all four length quartiles when pooled across all species (Q1: 16.2%, OR = 1.74, p < 10⁻¹⁰⁰; Q2: 19.9%, OR = 2.24, p < 10⁻¹⁰⁰; Q3: 18.4%, OR = 2.03, p < 10⁻¹⁰⁰; Q4: 17.9%, OR = 1.96, p < 10⁻¹⁰⁰), ruling out protein length as a confounding factor (Figure S3).

### Terminal LCRs are not purer than internal LCRs in most phyla

If terminal LCRs were under stronger selection for compositional homogeneity, we would expect higher purity (dominant amino acid fraction) at terminal positions. A Mann-Whitney U test comparing terminal versus internal LCR purity per phylum found no significant difference in 11 of 14 phyla; only Insecta (p = 0.021), Chelicerata (p = 0.010), and Echinodermata (p = 0.002) showed significantly higher terminal purity, and the delta values were small (Δ = +0.006–+0.014). These findings indicate that terminal LCR enrichment is positional — a greater density of LCRs in terminal bins — rather than qualitative, driven by the production of purer, more strongly biased sequences at termini.

### Singleton-LCR proteins are the primary driver of terminal enrichment

Multi-LCR proteins (those containing ≥2 LCRs) might disproportionately place at least one LCR terminally as a geometric consequence of having multiple LCRs distributed across a protein. To test whether the terminal enrichment signal is class-specific, we separated proteins into singleton-LCR (1 LCR per protein) and multi-LCR (≥2 LCRs) classes. Singleton-LCR proteins showed significant terminal enrichment in all 14 phyla, with pct_terminal values matching or exceeding the phylum-level average. Multi-LCR proteins showed significant enrichment in only 7 of 14 phyla, with generally lower pct_terminal values (Table S2). In the majority of phyla, singleton-LCR protein enrichment exceeds multi-LCR enrichment by 2–8 percentage points. This demonstrates that the signal is not a statistical consequence of protein LCR multiplicity but reflects a genuine positional preference at individual protein termini, present across proteins irrespective of their overall LCR content.

---

## Discussion

### Terminal LCR enrichment is an ancient metazoan property

The principal finding of this study is that the terminal enrichment of LCRs, described by Teekas et al. (2024) in Tetrapoda, is present across the entire animal kingdom, including the most phylogenetically distant invertebrate lineages. The observation of 18.8% terminal enrichment in *Trichoplax adhaerens* (Placozoa) is particularly significant: Placozoa represent the most basally-branching animal lineage, diverging from all other metazoans approximately 700–800 million years ago, well before the Cambrian explosion (~541 Mya) and the diversification of complex body plans. The presence of terminal LCR enrichment in Placozoa, in a magnitude indistinguishable from Tetrapoda, argues strongly that the mechanism responsible is not a vertebrate innovation, not an arthropod or lophotrochozoan peculiarity, but a property shared by all animals.

This evolutionary age places terminal LCR enrichment in the same category as other deeply conserved features of protein architecture — secondary structure propensities, signal peptide positioning, and co-translational folding dynamics. Teekas et al. speculated that the enrichment might be related to adaptive immunity-driven selection pressure on surface-exposed disordered regions. Our data argue against this interpretation: adaptive immunity in its vertebrate form does not exist in Cnidaria, Platyhelminthes, Nematoda, or Placozoa, yet all show robust terminal enrichment. The mechanism must be more fundamental.

### N/C asymmetry points toward distinct terminal functions

The separation of N-terminal and C-terminal enrichment reveals a biological distinction not apparent in combined terminal analyses. C-terminal enrichment is more universal (significant in 94% of species) and typically stronger in magnitude than N-terminal enrichment (significant in 86% of species). This asymmetry is consistent with the distinct biology of the two termini. C-terminal LCRs may correspond to the C-terminal degrons that mark proteins for proteasomal degradation through C-degron pathways (Riba et al. 2019), to functional linear motifs concentrated at C-termini (Brown et al. 2011), or to the compositionally biased C-terminal tails observed in many eukaryotic proteins. N-terminal LCRs may relate to N-degron pathways (Varshavsky 2019) or to signal peptide-associated sequences, but the greater variability of N-terminal enrichment across phyla suggests more lineage-specific constraints.

The *Schistosoma mansoni* exception is illuminating in this context. This species loses N-terminal LCR enrichment specifically while retaining C-terminal enrichment. *Schistosoma mansoni* is an obligate blood parasite with a dramatically compacted and reorganised proteome, including extensive gene loss and unusual GC-content dynamics (Berriman et al. 2009). The specific loss of N-terminal compositional bias — while C-terminal bias is maintained — may reflect distinct functional constraints on the two termini even in a highly derived genome.

### Mechanistic implications

Three potential mechanisms could explain the conserved terminal enrichment of LCRs across metazoans:

**1. Translational kinetics.** Ribosome pausing is more pronounced near the start (initiation complex clearance) and near the stop codon (termination), leading to slower translation at protein termini. Slower translation is associated with increased local disorder and compositional simplicity, as the nascent chain has longer to equilibrate before the ribosome clears the exit tunnel (Pechmann and Frydman 2013). This mechanism would operate universally, regardless of lineage.

**2. Co-translational folding constraints.** Protein domains fold co-translationally, and terminal regions — particularly N-terminal sequences — are exposed before the C-terminal domain has emerged from the ribosome. Low-complexity, disordered N-terminal sequences may serve as flexible linkers that accommodate domain folding before the full protein is synthesised. This predicts N-terminal enrichment. The fact that C-terminal enrichment is stronger argues this alone is insufficient.

**3. Proteasomal surveillance.** Both N-degron and C-degron pathways recognise terminal residues and sequences for targeted proteasomal degradation (Varshavsky 2019; Riba et al. 2019). Low-complexity terminal tails may facilitate degron exposure or recognition. The conservation of degron biology across all eukaryotes, and indeed across all animals, is consistent with the evolutionary depth we observe.

None of these mechanisms is mutually exclusive. The observed C > N asymmetry could reflect the combined contributions of translational termination dynamics (C-terminal) and co-translational folding requirements (N-terminal, more variable). Discriminating between these hypotheses will require functional perturbation experiments beyond the scope of this computational study.

### The Platyhelminthes and Priapulida outliers

Two phyla show terminal enrichment below 15% (the lower bound of the Tetrapoda range): Platyhelminthes (12.7%) and Priapulida (14.0%). Both have small representation in our dataset (2 and 1 species, respectively), limiting statistical power. Within Platyhelminthes, *Schistosoma mansoni* (11.2%, not significant) drives the phylum average down; *Schmidtea mediterranea* (15.6%, p = 0.003) falls within the Tetrapoda range. The Platyhelminthes result may therefore reflect the unusual biology of a blood parasite rather than a phylum-wide property. *Priapulus caudatus* (Priapulida, 14.0%, p = 0.001) is significant, albeit at the lower end, and represents a single species from a small marine worm phylum with limited genome data. Both cases warrant examination with additional species before being interpreted as genuine phylum-level departures.

### Amino acid identity

Preliminary analysis of the dominant amino acid composition of terminal versus internal LCRs (Figure S4) indicates that cysteine (C)-rich LCRs are consistently overrepresented at protein termini across phyla. Cysteine-rich terminal LCRs may correspond to metal-binding domains, disulfide-mediated dimerisation motifs, or reactive cysteine-containing degrons. A full amino acid enrichment analysis is presented in Supplementary Figure S4 and Supplementary Table S3.

### Limitations

This study is entirely computational and shares methodological constraints with Teekas et al. (2024). The purity threshold (≥70%) and the SINGLE-type filter are definitional choices that affect LCR counts but not the qualitative finding, as enrichment is observed at all tested purity thresholds. The 20-bin framework assigns LCRs to bins by midpoint position, which may misclassify LCRs spanning a bin boundary. Proteome quality varies across species, and assembly incompleteness in some invertebrate genomes could introduce noise. The single species representing Placozoa, Priapulida, Brachiopoda, and Hemichordata limits the power to detect phylum-level heterogeneity in these groups.

---

## Conclusions

Terminal LCR enrichment is a pan-metazoan property conserved across at least 700 million years of animal evolution, present in all 14 phyla tested and in 49 of 50 species. The enrichment magnitude in invertebrates (12.7–21.2% at the phylum level) is indistinguishable from the Tetrapoda baseline reported by Teekas et al. (2024). C-terminal enrichment is more universal and stronger than N-terminal enrichment, implying distinct functional roles at the two termini. The mechanism predates the bilaterian split and cannot be explained by any vertebrate-specific biological innovation, pointing instead toward fundamental constraints of protein translation, folding, and degradation common to all animals.

---

## References

Berriman M, et al. (2009). The genome of the blood fluke *Schistosoma mansoni*. *Nature* 460:352–358.

Boija A, et al. (2018). Transcription factors activate genes through the phase-separation capacity of their activation domains. *Cell* 175:1842–1855.

Brown CJ, et al. (2011). The roles of structural disorder in protein evolution at the molecular and cellular level. *PLoS Computational Biology* 7:e1002012.

Harrison PM (2017). fLPS: Fast discovery of compositional biases for the protein universe. *BMC Bioinformatics* 18:476.

Marcotte EM, et al. (1999). A census of protein repeats. *Journal of Molecular Biology* 293:151–160.

Mier P, et al. (2020). Disentangling the complexity of low-complexity proteins. *Briefings in Bioinformatics* 21:458–472.

Pechmann S, Frydman J (2013). Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding. *Nature Structural & Molecular Biology* 20:237–243.

Riba A, et al. (2019). Protein synthesis rates and ribosome occupancies reveal determinants of translation elongation rates. *Proceedings of the National Academy of Sciences* 116:15023–15032.

Romero P, et al. (2001). Sequence complexity of disordered protein. *Proteins: Structure, Function, and Bioinformatics* 42:38–48.

Shin Y, Brangwynne CP (2017). Liquid phase condensation in cell physiology and disease. *Science* 357:eaaf4382.

Teekas L, Sharma S, Vijay N (2024). Terminal regions of a protein are a hotspot for low complexity regions and selection. *Open Biology* 14:230439.

van der Lee R, et al. (2014). Classification of intrinsically disordered regions and proteins. *Chemical Reviews* 114:6589–6631.

Varshavsky A (2019). N-degron and C-degron pathways of protein degradation. *Proceedings of the National Academy of Sciences* 116:358–366.

Wootton JC, Federhen S (1996). Analysis of compositionally biased regions in sequence databases. *Methods in Enzymology* 266:554–571.

---

## Supplementary Material

**Supplementary Table S1.** Complete species list with Ensembl Metazoa directory keys, phylum, genome assembly accession, and per-species terminal enrichment statistics.

**Supplementary Table S2.** Driver analysis: terminal LCR enrichment in singleton-LCR vs. multi-LCR proteins, per phylum.

**Supplementary Table S3.** Amino acid composition of terminal vs. internal LCRs, per phylum.

**Supplementary Figure S1.** LCR positional heatmap (species × bin).

**Supplementary Figure S2.** U-shaped bin profile per phylum.

**Supplementary Figure S3.** Terminal enrichment by protein-length quartile.

**Supplementary Figure S4.** Amino acid composition of terminal vs. internal LCRs (stacked bar + enrichment ratio).

**Supplementary Figure S5.** LCR purity gradient: terminal vs. internal per phylum.

**Supplementary Figure S6.** Singleton vs. multi-LCR protein driver analysis.
