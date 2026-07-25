#!/usr/bin/env python3
"""
Download expansion proteomes for under-represented Viridiplantae tiers.

EnsemblPlants r63 is agronomy-focused: ~90% of its green-plant species are eudicots,
non-grass monocots and grasses (PGLS tiers 7-9), while the early-diverging lineages that
actually define the algae→land-plant gradient in 19_pgls_viridiplantae.py are sparse or
absent:

  Tier 0 Chlorophyta:          n=2  (Chlamydomonas, Ostreococcus) → target n ≥ 4
  Tier 1 Charophytes:          n=1  (Chara)                       → target n ≥ 3
  Tier 2 Bryophytes:           n=2  (Physcomitrium, Marchantia)   → target n ≥ 4
  Tier 3 Lycophytes:           n=1  (Selaginella)                 → target n ≥ 2
  Tier 4 Ferns (absent):       n=0                                → target n ≥ 2
  Tier 5 Gymnosperms (absent): n=0                                → target n ≥ 3
  Tier 6 Basal angiosperms:    n=2  (Amborella, Nymphaea)         → target n ≥ 4

Source: NCBI Datasets, pinned by assembly accession
  Every accession below was verified live against the NCBI Datasets API: it exists, is
  the current version, and carries protein_coding annotation (a protein FASTA is
  downloadable). This replaces the old UniProt-proteome-id approach, whose ids drift/404
  and which left ferns + gymnosperms un-downloadable. Proteomes are reduced to the
  longest isoform per gene via the assembly GFF3 — identical contract to 01a/01b.

Hard limits (NOT reachable via NCBI — gene models are portal-only, e.g. Phytozome /
ConGenIE / FernBase, in incompatible formats):
  - Hornworts (Anthocerotophyta): zero NCBI-annotated genomes.
  - Mesotaenium, Zygnema, Penium (charophytes); Azolla, Salvinia, Isoetes; Ginkgo,
    Picea: genome-only on NCBI. Congeners with NCBI annotation are used instead
    (Closterium for the Zygnematophyceae sister; Taxus/Cryptomeria/Pinus longaeva for
    gymnosperms; Ceratopteris/Adiantum for ferns).

Running this script:
  python3 scripts/01c_download_viridiplantae_expansion.py
  bash scripts/02_run_flps.sh            # fLPS on the new proteomes
  python3 main.py --from-phase 2         # re-run analysis incl. 19_pgls_viridiplantae.py

Outputs (appended to existing files):
  data/proteomes/<species_key>.longest.fa
  results/species_manifest.tsv
"""

import io
import sys
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import NamedTuple

import re
import requests

from config import PROJECT_DIR, RESULTS_DIR

OUT_DIR  = PROJECT_DIR / "data" / "proteomes"
MANIFEST = RESULTS_DIR / "species_manifest.tsv"
MANIFEST_HEADER = ["species_key", "display_name", "phylum", "fa_path", "domain", "taxon_id"]

NCBI_DATASETS_DL = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{}"
                    "/download?include_annotation_type=PROT_FASTA"
                    "&include_annotation_type=GENOME_GFF")

OUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


class Species(NamedTuple):
    key:      str    # species_key in manifest (snake_case)
    display:  str    # display_name
    phylum:   str    # "Viridiplantae" (charophytes + land plants) or "Chlorophyta" (green algae)
    accession: str   # NCBI assembly accession (verified: current + protein-annotated)
    taxon_id: str    # NCBI taxon ID
    tier:     int    # PGLS tier (documentation only; 19_pgls_viridiplantae.py re-derives it by genus)


# ── Expansion species (accession-verified NCBI Datasets) ────────────────────────
# phylum="Chlorophyta" for green algae (tier 0), "Viridiplantae" for charophytes and
# land plants (tiers 1-6), matching the labels 01b writes and 19's VIRIDIPLANTAE_PHYLA.
EXPANSION_SPECIES: list[Species] = [
    # Tier 0 — Chlorophyta (green algae)
    Species("volvox_carteri",           "Volvox carteri f. nagariensis", "Chlorophyta",  "GCF_000143455.1", "3068",   0),
    Species("micromonas_commoda",       "Micromonas commoda",            "Chlorophyta",  "GCF_000090985.2", "296587", 0),
    Species("chlorella_variabilis",     "Chlorella variabilis",          "Chlorophyta",  "GCF_000147415.1", "554065", 0),
    Species("coccomyxa_subellipsoidea", "Coccomyxa subellipsoidea C-169","Chlorophyta",  "GCF_000258705.1", "574566", 0),

    # Tier 1 — Charophytes (streptophyte algae; sister to land plants)
    Species("klebsormidium_nitens",     "Klebsormidium nitens",          "Viridiplantae","GCA_000708835.1", "105231", 1),
    Species("closterium_nies54",        "Closterium sp. NIES-54",        "Viridiplantae","GCA_949281255.1", "3014311",1),  # Zygnematophyceae, closest land-plant sister on NCBI

    # Tier 2 — Bryophytes
    Species("sphagnum_fallax",          "Sphagnum fallax",               "Viridiplantae","GCA_021442195.1", "53036",  2),
    Species("ceratodon_purpureus",      "Ceratodon purpureus",           "Viridiplantae","GCA_014871385.1", "3225",   2),

    # Tier 3 — Lycophytes
    Species("diphasiastrum_complanatum","Diphasiastrum complanatum",     "Viridiplantae","GCA_029204225.1", "34168",  3),

    # Tier 4 — Ferns (Polypodiopsida)
    Species("ceratopteris_richardii",   "Ceratopteris richardii",        "Viridiplantae","GCA_020310875.1", "49495",  4),
    Species("adiantum_capillus_veneris","Adiantum capillus-veneris",     "Viridiplantae","GCA_014529385.2", "13818",  4),

    # Tier 5 — Gymnosperms (three conifer families; the only NCBI-annotated gymnosperms)
    Species("taxus_chinensis",          "Taxus chinensis",               "Viridiplantae","GCA_019776745.2", "29808",  5),
    Species("cryptomeria_japonica",     "Cryptomeria japonica",          "Viridiplantae","GCF_030272615.1", "3369",   5),
    Species("pinus_longaeva",           "Pinus longaeva",                "Viridiplantae","GCA_056508025.1", "3344",   5),

    # Tier 6 — Basal angiosperms / magnoliids (ANA + magnoliid grade)
    Species("cinnamomum_kanehirae",     "Cinnamomum micranthum f. kanehirae","Viridiplantae","GCA_003546025.1","337451",6),
    Species("aristolochia_fimbriata",   "Aristolochia fimbriata",        "Viridiplantae","GCA_019845555.1", "158543", 6),
]


# ── NCBI Datasets download (longest isoform per gene, via GFF3) ──────────────────

def _protein_to_gene(gff_text: str) -> dict[str, str]:
    """Map protein accession -> gene key from GFF3 CDS attributes (isoform collapse)."""
    p2g: dict[str, str] = {}
    for line in gff_text.splitlines():
        if not line or line[0] == "#":
            continue
        c = line.split("\t")
        if len(c) < 9 or c[2] != "CDS":
            continue
        attr = c[8]
        pid = re.search(r"protein_id=([^;]+)", attr)
        if not pid:
            continue
        gene = (re.search(r"\bgene=([^;]+)", attr)
                or re.search(r"locus_tag=([^;]+)", attr)
                or re.search(r"Parent=([^;,]+)", attr))
        p2g[pid.group(1)] = gene.group(1) if gene else pid.group(1)
    return p2g


def download_ncbi(accession: str, out_fa: Path) -> int | None:
    """Fetch a protein set by NCBI assembly accession, keep the longest isoform per gene
    (grouped through the assembly GFF3), and write out_fa (60-col, LF). Returns the gene
    count, or None on failure — same one-seq-per-gene contract as 01a/01b.
    Retries transient network/DNS errors (macOS resolver can drop parallel lookups)."""
    content = None
    for attempt in range(5):
        try:
            r = requests.get(NCBI_DATASETS_DL.format(accession), timeout=600)
            r.raise_for_status()
            content = r.content
            break
        except requests.exceptions.RequestException as e:
            if attempt == 4:
                print(f"    [!] NCBI download error after retries: {e}")
                return None
            time.sleep(2 * (attempt + 1))
    try:
        zf = zipfile.ZipFile(io.BytesIO(content))
    except zipfile.BadZipFile as e:
        print(f"    [!] bad zip for {accession}: {e}")
        return None
    faa = gff = None
    for name in zf.namelist():
        if name.endswith("protein.faa"):
            faa = zf.read(name).decode("utf-8", "replace")
        elif name.endswith("genomic.gff"):
            gff = zf.read(name).decode("utf-8", "replace")
    if not faa:
        print(f"    [!] no protein.faa in NCBI package for {accession}")
        return None

    p2g = _protein_to_gene(gff) if gff else {}
    longest: dict[str, tuple[str, str]] = {}
    order: list[str] = []
    header, parts = None, []

    def flush():
        if header is None:
            return
        seq = "".join(parts)
        gene = p2g.get(header.split()[0], header.split()[0])
        cur = longest.get(gene)
        if cur is None:
            longest[gene] = (header, seq)
            order.append(gene)
        elif len(seq) > len(cur[1]):
            longest[gene] = (header, seq)

    for line in faa.splitlines():
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            parts = []
        else:
            parts.append(line.strip())
    flush()

    with open(out_fa, "w", encoding="utf-8", newline="\n") as fh:
        for gene in order:
            hdr, seq = longest[gene]
            fh.write(f">{hdr}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    return len(order)


def upsert_manifest(rows: list[dict]) -> None:
    existing: dict[str, dict] = {}
    if MANIFEST.exists() and MANIFEST.stat().st_size > 0:
        lines = MANIFEST.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        for line in lines[1:]:
            rec = dict(zip(header, line.split("\t")))
            existing[rec.get("species_key", "")] = rec
    for r in rows:
        existing[r["species_key"]] = r
    with open(MANIFEST, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(MANIFEST_HEADER) + "\n")
        for key in sorted(existing):
            rec = existing[key]
            f.write("\t".join(str(rec.get(col, "")) for col in MANIFEST_HEADER) + "\n")


# ── Per-species worker (thread-pool safe) ───────────────────────────────────────

def process(sp: Species) -> tuple[dict | None, str]:
    out_fa = OUT_DIR / f"{sp.key}.longest.fa"
    row = {"species_key": sp.key, "display_name": sp.display, "phylum": sp.phylum,
           "fa_path": str(out_fa), "domain": "Non-metazoan Eukaryota", "taxon_id": sp.taxon_id}
    if out_fa.exists() and out_fa.stat().st_size > 0:
        print(f"  [EXISTS] tier {sp.tier}  {sp.display}")
        return row, "EXISTS"
    print(f"  [DOWNLOAD] tier {sp.tier}  {sp.display}  ({sp.accession})")
    n = download_ncbi(sp.accession, out_fa)
    if not n:
        print(f"    [!] failed: {sp.display} ({sp.accession})")
        return None, "FAILED"
    print(f"    [+] {n} genes -> {out_fa.name}")
    return row, "OK"


def main(workers: int = 3):
    print(f"Viridiplantae tier expansion: {len(EXPANSION_SPECIES)} target species "
          f"(NCBI Datasets, accession-verified)\n")

    rows: list[dict] = []
    failed: list[str] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(process, sp): sp for sp in EXPANSION_SPECIES}
        for fut in as_completed(futures):
            try:
                row, status = fut.result()
            except Exception as e:
                print(f"    [!] worker error: {e}")
                continue
            if row is not None:
                rows.append(row)
            elif status == "FAILED":
                failed.append(futures[fut].display)

    if rows:
        upsert_manifest(rows)

    print(f"\n{'='*60}")
    print(f"Downloaded / present : {len(rows)}/{len(EXPANSION_SPECIES)}")
    if failed:
        print(f"Failed               : {len(failed)}")
        for name in failed:
            print(f"  {name}")

    # Projected tier composition (one representative per genus, as 19 de-duplicates)
    from collections import defaultdict
    per_tier = defaultdict(set)
    for sp in EXPANSION_SPECIES:
        per_tier[sp.tier].add(sp.display.split()[0])
    print("\nGenera added per PGLS tier (before merging with existing EnsemblPlants):")
    for tier in sorted(per_tier):
        print(f"  Tier {tier}: {', '.join(sorted(per_tier[tier]))}")

    print(f"\nNext: bash scripts/02_run_flps.sh  &&  python3 main.py --from-phase 2")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Viridiplantae tier expansion (NCBI Datasets).")
    ap.add_argument("--workers", type=int, default=3,
                    help="parallel downloads (default 3; use 1 if DNS/network is flaky)")
    main(workers=ap.parse_args().workers)
