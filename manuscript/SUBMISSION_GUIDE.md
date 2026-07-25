# Submission guide — $0 route: bioRxiv → Peer Community In Genomics

Goal: publish `manuscript_v10.md` at **zero cost**, open-access, peer-reviewed, indexed —
without any institutional affiliation. Do the steps in order. Never enter payment details;
both routes below are free by design.

---

## 0. Before you upload — fill these in (only you can)

Placeholders left in `manuscript_v10.md`, plus prep tasks:

- [x] **City, Country** — done; affiliation reads "Independent Researcher, London, United Kingdom".
- [x] **ORCID iD** — done; 0000-0003-4242-7507 is in the correspondence line.
- [ ] **GitHub repo** — push this project (scripts + `results/` tables) to a public GitHub
      repo; paste its URL into the "Data and code availability" section. Then create a free
      **Zenodo** release from that repo to get a permanent DOI (GitHub→Zenodo integration is
      one click) and add that DOI too. Reviewers will expect working links.
- [x] **Reference list verified** — 28 references, **every one carrying a DOI** verified against
      the Crossref API (v9). Author lists expanded to full where ≤10 authors; "et al." retained
      only where the list exceeds 10 (Boija 20, van der Lee 18, Virtanen 35+), which every major
      style permits. v9 also corrected a wrong author: Holdsworth et al. 2020 listed
      "Estavillo GM" as 5th author; the real 5th author is **Zubrycka A**. Zero citation orphans,
      zero uncited references, zero self-citations. Earlier round: Six
      errors were corrected: Piatkov (fabricated middle authors), Gibbs (year 2014→2011),
      Lancaster/PLAAC (pages), Irastortza-Olaziregi (article number), Saravanan (first initial);
      two uncited references (Eme, Berriman) removed. The "et al." style task is now closed —
      only the three >10-author entries retain it, which every major style permits.
- [x] **Convert to PDF** — bioRxiv and PCI want a PDF (or Word) with figures. On the Mac,
      from the repo root:
      ```bash
      pandoc manuscript/manuscript_v10.md -o manuscript_v10.pdf \
        --pdf-engine=xelatex -V geometry:margin=1in \
        --include-in-header=manuscript/pdf-header.tex
      ```
      **`--include-in-header` is mandatory.** Without it pandoc still exits 0, but xelatex
      silently drops `≥`, `⁻`, `⁵`, `⁶` as missing glyphs — "purity ≥70%" prints as
      "purity 70%" and "7.1 × 10⁻⁵" loses its exponent. `pdf-header.tex` maps those
      characters to math mode. Always check the build log is empty; any
      `Missing character` line means the PDF is corrupt regardless of exit code.
      Over ssh, PATH is minimal — prefix with
      `export PATH=/Library/TeX/texbin:/opt/homebrew/bin:$PATH`.
      Figures already exist in `figures/` (fig1–fig10, suppfig_*). For bioRxiv you can embed
      them in the PDF or upload the PDFs as separate figure files — either is accepted.

---

## 1. Post the preprint to bioRxiv (free, do this first)

This makes the paper openly readable worldwide and timestamps your priority, regardless of
which journal it later lands in.

1. Go to https://www.biorxiv.org → "Submit a Manuscript." Create a free account.
2. **Corresponding author:** yourself; affiliation field → type "Independent Researcher".
   A personal email is fine.
3. **Subject category:** *Bioinformatics* (primary). Reasonable alternatives: *Genomics* or
   *Molecular Biology*.
4. **License:** choose **CC-BY 4.0** (required later by PCI; maximises reuse). If you want to
   be more conservative you may pick CC-BY-NC-ND, but CC-BY keeps every downstream option open.
5. **Title / Abstract:** paste from the manuscript (see block below).
6. Upload the PDF (+ figure files if separate). Declare no competing interests; funding "None".
7. Submit. bioRxiv screens (a few days), then posts with a citable DOI.

---

## 2. Submit the same preprint to Peer Community In Genomics (free peer review)

This is your guaranteed-free, peer-reviewed, indexed publication.

1. Go to https://genomics.peercommunityin.org → "Submit your preprint."
2. Provide the **bioRxiv DOI** from step 1 (PCI reviews the preprint in place).
3. Write a short cover/submission comment: one paragraph on the finding (first positional
   LCR analysis across all three domains; mechanistic N/C-polarity + composition fingerprint)
   and why it fits a genomics/computational audience.
4. Suggest 3–5 potential recommenders/reviewers if asked (people who work on low-complexity
   regions, protein disorder, or comparative proteomics — e.g. authors you cite).
5. If a recommender takes it on, you get real peer review. On a positive recommendation you
   can either (a) keep it as a **PCI-recommended preprint** (citable, with a DOI), or
   (b) publish the final version free in **Peer Community Journal** (diamond OA, DOAJ/Scopus
   indexed). Both cost nothing.

**If PCI declines or is too slow (fallback, still $0):** submit the same manuscript to a
subscription/hybrid journal and choose the **non-open-access** option (free to publish):
- *Bioinformatics* (Oxford) — "Genome analysis" section.
- *Proteins: Structure, Function, and Bioinformatics* (Wiley).
Your bioRxiv preprint keeps the content free to read even though the journal version is paywalled.
When the form offers open access, decline — do not pay an APC.

---

## Title and abstract to paste

**Title:**
Terminal Low-Complexity Regions as a General Architectural Property of Proteins: A Positional
and Compositional Fingerprint Across 772 Proteomes and 43 Lineages

**Abstract:** use the Abstract paragraph verbatim from `manuscript_v10.md`.

**Keywords:** low-complexity regions, protein termini, protein architecture, compositional
bias, co-translational folding, N-terminal processing, signal peptides, intrinsically
disordered regions, degron biology, fLPS2

---

## The one rule
Both routes are free by design. If any submission form presents a **mandatory** fee, stop and
switch to the next venue. You never need to pay to publish this paper.
