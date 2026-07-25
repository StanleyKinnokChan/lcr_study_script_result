#!/usr/bin/env python3
"""
Download non-metazoan outgroup proteomes from two complementary sources and register
them in results/species_manifest.tsv, alongside the metazoan rows written by 01a.

Why three sources
  1. UniProt reference proteomes — a hand-picked, taxonomically balanced backbone of
     prokaryotes (Bacteria + Archaea) and key reference eukaryotes. Canonical FASTA is
     already one sequence per gene.
  2. Ensembl Genomes release 63 — guaranteed to have FASTA on the FTP (UniProt proteome
     ids drift and 404 over time). Used to broaden the eukaryotic outgroups from the
     release-63 lookup tables in the project root.
  3. NCBI Datasets, accession-pinned — used to (a) recover reference proteomes whose
     UniProt proteome id has since 404'd (incl. the animal-sister Holozoa: Salpingoeca,
     Capsaspora, Sphaeroforma), and (b) add specific under-sampled invertebrate phyla
     (base-of-tree ctenophores/sponges, xenacoelomorphs, rotifers, myriapods, …) that
     Ensembl Metazoa r63 covers with only a single species. Pinned by stable GCA/GCF
     accession and reduced to longest-isoform-per-gene via the assembly GFF3, so the
     one-seq-per-gene contract is identical to the other two sources.

Outgroup boundary being tested
  - enrichment in Bacteria/Archaea  -> mechanism is ancient (translational, >3.5 Bya)
  - enrichment in Eukaryota only     -> arose ~2 Bya (eukaryotic)
  - enrichment in Metazoa only       -> animal-specific (~700 Mya)

Selection (lean profile)
  Plants   : all unique species         (~125)
  Protists : all unique species         (~134)
  Fungi    : pan_compara references only (small, balanced — avoids Ascomycota swamp)
  Bacteria : pan_compara references only (Ensembl Bacteria has ~10k species; the
             curated UniProt set plus pan-compara references is plenty)
  Flip any division's "select" below between "all" / "pan_compara" / "peptide_compara".

Output & contract (identical to 01a)
  data/proteomes/{species_key}.longest.fa  (one seq per gene)
  results/species_manifest.tsv  columns: species_key display_name phylum fa_path domain taxon_id
  Phylum/domain resolved from NCBI Taxonomy, cached in results/taxonomy_cache.tsv
  (shared with 01a). UniProt rows take precedence over Ensembl rows for the same species.
"""

import argparse
import gzip
import io
import re
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

# ── Paths & config ──────────────────────────────────────────────────────────────

RELEASE = 63
PROJECT_DIR = Path(__file__).parent.parent
OUT_DIR     = PROJECT_DIR / "data" / "proteomes"
RESULTS_DIR = PROJECT_DIR / "results"
MANIFEST    = RESULTS_DIR / "species_manifest.tsv"
TAX_CACHE   = RESULTS_DIR / "taxonomy_cache.tsv"
STATUS_FILE = RESULTS_DIR / "download_status_outgroups.tsv"

ENSEMBL_BASE  = f"https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-{RELEASE}"
UNIPROT_FTP   = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes"
EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_DATASETS_DL = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{}"
                    "/download?include_annotation_type=PROT_FASTA"
                    "&include_annotation_type=GENOME_GFF")

MANIFEST_HEADER = ["species_key", "display_name", "phylum", "fa_path", "domain", "taxon_id"]

OUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Ensembl release-63 divisions to mine for outgroups.
#   file    : release-63 lookup table in the project root
#   ftp_dir : division directory name on the Ensembl Genomes FTP
#   select  : "all" | "pan_compara" | "peptide_compara"
ENSEMBL_DIVISIONS = {
    "plants":   {"file": "species_EnsemblPlant.txt",    "ftp_dir": "plants",   "select": "all"},
    "protists": {"file": "species_EnsemblProtists.txt", "ftp_dir": "protists", "select": "all"},
    "fungi":    {"file": "species_EnsemblFungi.txt",     "ftp_dir": "fungi",    "select": "pan_compara"},
    "bacteria": {"file": "species_EnsemblBacteria.txt",  "ftp_dir": "bacteria", "select": "pan_compara"},
}

# ── UniProt curated backbone ────────────────────────────────────────────────────
# species_key -> (display_name, kingdom_dir, proteome_id, taxon_id)
UNIPROT_OUTGROUPS = {
    # Fungi — Ascomycota
    "saccharomyces_cerevisiae":      ("Saccharomyces cerevisiae",      "Eukaryota", "UP000002311", "559292"),
    "schizosaccharomyces_pombe":     ("Schizosaccharomyces pombe",     "Eukaryota", "UP000002485", "284812"),
    "neurospora_crassa":             ("Neurospora crassa",             "Eukaryota", "UP000001805", "367110"),
    "aspergillus_nidulans":          ("Aspergillus nidulans",          "Eukaryota", "UP000000560", "227321"),
    "candida_albicans":              ("Candida albicans",              "Eukaryota", "UP000000559", "237561"),
    "aspergillus_fumigatus":         ("Aspergillus fumigatus",         "Eukaryota", "UP000002530", "330879"),
    # Fungi — Basidiomycota / Mucoromycota
    "cryptococcus_neoformans":       ("Cryptococcus neoformans",       "Eukaryota", "UP000010091", "235443"),
    "ustilago_maydis":               ("Ustilago maydis",               "Eukaryota", "UP000000561", "237631"),
    "rhizopus_delemar":              ("Rhizopus delemar",              "Eukaryota", "UP000019020", "246409"),
    # Plants — Embryophyta
    "arabidopsis_thaliana":          ("Arabidopsis thaliana",          "Eukaryota", "UP000006548", "3702"),
    "oryza_sativa":                  ("Oryza sativa",                  "Eukaryota", "UP000007015", "39947"),
    "populus_trichocarpa":           ("Populus trichocarpa",           "Eukaryota", "UP000006728", "3694"),
    "solanum_lycopersicum":          ("Solanum lycopersicum",          "Eukaryota", "UP000004994", "4081"),
    # Protists — diverse non-animal eukaryote lineages
    "dictyostelium_discoideum":      ("Dictyostelium discoideum",      "Eukaryota", "UP000002195", "44689"),
    "entamoeba_histolytica":         ("Entamoeba histolytica",         "Eukaryota", "UP000000673", "294381"),
    "plasmodium_falciparum":         ("Plasmodium falciparum",         "Eukaryota", "UP000001450", "36329"),
    "trypanosoma_brucei":            ("Trypanosoma brucei",            "Eukaryota", "UP000008524", "185431"),
    "leishmania_major":              ("Leishmania major",              "Eukaryota", "UP000000542", "347515"),
    "giardia_intestinalis":          ("Giardia intestinalis",          "Eukaryota", "UP000001057", "184922"),
    "chlamydomonas_reinhardtii":     ("Chlamydomonas reinhardtii",     "Eukaryota", "UP000006906", "3055"),
    "naegleria_gruberi":             ("Naegleria gruberi",             "Eukaryota", "UP000009136", "214684"),
    # Opisthokont protists — immediate relatives of Metazoa
    "monosiga_brevicollis":          ("Monosiga brevicollis",          "Eukaryota", "UP000001357", "81824"),
    "salpingoeca_rosetta":           ("Salpingoeca rosetta",           "Eukaryota", "UP000011781", "946362"),
    "sphaeroforma_arctica":          ("Sphaeroforma arctica",          "Eukaryota", "UP000027617", "72019"),
    "capsaspora_owczarzaki":         ("Capsaspora owczarzaki",         "Eukaryota", "UP000001498", "595528"),
    # Archaea
    "methanocaldococcus_jannaschii": ("Methanocaldococcus jannaschii", "Archaea",   "UP000000805", "243232"),
    "halobacterium_salinarum":       ("Halobacterium salinarum",       "Archaea",   "UP000000554", "64091"),
    "methanobacterium_thermoautotrophicum": ("Methanobacterium thermoautotrophicum", "Archaea", "UP000000853", "187420"),
    "archaeoglobus_fulgidus":        ("Archaeoglobus fulgidus",        "Archaea",   "UP000001013", "224325"),
    "pyrococcus_abyssi":             ("Pyrococcus abyssi",             "Archaea",   "UP000000810", "272844"),
    "methanosarcina_acetivorans":    ("Methanosarcina acetivorans",    "Archaea",   "UP000002571", "188937"),
    "thermococcus_kodakarensis":     ("Thermococcus kodakarensis",     "Archaea",   "UP000000536", "69014"),
    "haloarcula_marismortui":        ("Haloarcula marismortui",        "Archaea",   "UP000001426", "272569"),
    "sulfolobus_acidocaldarius":     ("Sulfolobus acidocaldarius",     "Archaea",   "UP000001018", "330779"),
    "thermoplasma_acidophilum":      ("Thermoplasma acidophilum",      "Archaea",   "UP000000578", "273075"),
    "nitrosopumilus_maritimus":      ("Nitrosopumilus maritimus",      "Archaea",   "UP000000748", "335283"),
    # Bacteria — broad phylum sampling
    "escherichia_coli_k12":          ("Escherichia coli K12",          "Bacteria",  "UP000000625", "83333"),
    "pseudomonas_aeruginosa":        ("Pseudomonas aeruginosa",        "Bacteria",  "UP000002438", "208964"),
    "salmonella_typhimurium":        ("Salmonella typhimurium",        "Bacteria",  "UP000008962", "99287"),
    "caulobacter_crescentus":        ("Caulobacter crescentus",        "Bacteria",  "UP000001816", "190650"),
    "bacillus_subtilis":             ("Bacillus subtilis",             "Bacteria",  "UP000001570", "224308"),
    "staphylococcus_aureus":         ("Staphylococcus aureus",         "Bacteria",  "UP000008816", "93061"),
    "streptococcus_pneumoniae":      ("Streptococcus pneumoniae",      "Bacteria",  "UP000000586", "171101"),
    "mycobacterium_tuberculosis":    ("Mycobacterium tuberculosis",    "Bacteria",  "UP000001584", "83332"),
    "streptomyces_coelicolor":       ("Streptomyces coelicolor",       "Bacteria",  "UP000001973", "100226"),
    "streptomyces_avermitilis":      ("Streptomyces avermitilis",      "Bacteria",  "UP000000708", "227882"),
    "streptomyces_griseus":          ("Streptomyces griseus",          "Bacteria",  "UP000001685", "455632"),
    "synechocystis_pcc6803":         ("Synechocystis sp. PCC 6803",    "Bacteria",  "UP000001425", "1148"),
    "borrelia_burgdorferi":          ("Borrelia burgdorferi",          "Bacteria",  "UP000001807", "224326"),
    "treponema_pallidum":            ("Treponema pallidum",            "Bacteria",  "UP000000811", "243276"),
    "deinococcus_radiodurans":       ("Deinococcus radiodurans",       "Bacteria",  "UP000002524", "243230"),
    "thermus_thermophilus":          ("Thermus thermophilus HB8",      "Bacteria",  "UP000000592", "300852"),
    "aquifex_aeolicus":              ("Aquifex aeolicus",              "Bacteria",  "UP000000949", "224324"),
    "thermotoga_maritima":           ("Thermotoga maritima",           "Bacteria",  "UP000008183", "243274"),
    "clostridioides_difficile":      ("Clostridioides difficile",      "Bacteria",  "UP000001978", "272563"),
}

# ── NCBI Datasets accession-pinned proteomes ────────────────────────────────────
# (species_key, display_name, assembly_accession, taxon_id)
# Every accession was checked live against the NCBI Datasets API: it exists, is the
# current version, and carries protein_coding annotation (a protein FASTA is
# downloadable). Domain + phylum are resolved from taxon_id like every other row.
# Keys shared with UNIPROT_OUTGROUPS take this route instead (their UniProt ids 404).
NCBI_PROTEOMES: list[tuple[str, str, str, str]] = [
    # ── Recover UNIPROT_404 reference proteomes ─────────────────────────────────
    # Animal-sister Holozoa — the yardstick for "is this gene family animal-specific?"
    ("capsaspora_owczarzaki",        "Capsaspora owczarzaki",        "GCF_000151315.2", "595528"),
    ("salpingoeca_rosetta",          "Salpingoeca rosetta",          "GCF_000188695.1", "946362"),
    ("sphaeroforma_arctica",         "Sphaeroforma arctica",         "GCF_001186125.1", "72019"),
    # Other early-diverging eukaryote references
    ("entamoeba_histolytica",        "Entamoeba histolytica",        "GCF_000208925.1", "294381"),
    ("giardia_intestinalis",         "Giardia intestinalis",         "GCF_000002435.2", "5741"),
    ("leishmania_major",             "Leishmania major",             "GCF_000002725.2", "347515"),
    ("naegleria_gruberi",            "Naegleria gruberi",            "GCF_000004985.1", "5762"),
    ("oryza_sativa",                 "Oryza sativa",                 "GCF_034140825.1", "39947"),
    ("populus_trichocarpa",          "Populus trichocarpa",          "GCF_000002775.5", "3694"),
    ("rhizopus_delemar",             "Rhizopus delemar",             "GCF_000149305.1", "246409"),
    # Deep prokaryote references (the ancient-vs-eukaryotic boundary the study tests)
    ("salmonella_typhimurium",       "Salmonella typhimurium",       "GCF_000006945.2", "99287"),
    ("streptomyces_avermitilis",     "Streptomyces avermitilis",     "GCF_000009765.2", "227882"),
    ("synechocystis_pcc6803",        "Synechocystis sp. PCC 6803",   "GCF_000009725.1", "1148"),
    ("thermoplasma_acidophilum",     "Thermoplasma acidophilum",     "GCF_000195915.1", "273075"),
    ("thermus_thermophilus",         "Thermus thermophilus HB8",     "GCF_000091545.1", "300852"),
    ("aquifex_aeolicus",             "Aquifex aeolicus",             "GCF_000008625.1", "224324"),
    ("archaeoglobus_fulgidus",       "Archaeoglobus fulgidus",       "GCF_000008665.1", "224325"),
    ("haloarcula_marismortui",       "Haloarcula marismortui",       "GCF_000011085.1", "272569"),
    ("methanobacterium_thermoautotrophicum", "Methanobacterium thermoautotrophicum", "GCF_000008645.1", "187420"),
    ("methanosarcina_acetivorans",   "Methanosarcina acetivorans",   "GCF_000007345.1", "188937"),
    ("nitrosopumilus_maritimus",     "Nitrosopumilus maritimus SCM1","GCF_000018465.1", "436308"),

    # ── Metazoan expansion: fixable singleton/doubleton phyla ────────────────────
    # Ensembl Metazoa r63 holds only one species for each of these; add a verified
    # second/third so no clade-level signal rests on a single genome. Domain resolves
    # to "Metazoa" and phylum to the labels already used in the manifest.
    ("ephydatia_muelleri",           "Ephydatia muelleri",           "GCA_049114765.1", "6052"),     # Porifera
    ("sycon_ciliatum",               "Sycon ciliatum",               "GCF_964019385.1", "27933"),    # Porifera
    ("bolinopsis_microptera",        "Bolinopsis microptera",        "GCF_026151205.1", "2820187"),  # Ctenophora
    ("convolutriloba_macropyga",     "Convolutriloba macropyga",     "GCF_964194025.1", "536237"),   # Xenacoelomorpha
    ("symsagittifera_roscoffensis",  "Symsagittifera roscoffensis",  "GCF_963678635.1", "84072"),    # Xenacoelomorpha
    ("brachionus_plicatilis",        "Brachionus plicatilis",        "GCA_003710015.1", "10195"),    # Rotifera
    ("rotaria_sordida",              "Rotaria sordida",              "GCA_905250125.1", "392033"),   # Rotifera
    ("scutigera_coleoptrata",        "Scutigera coleoptrata",        "GCA_982266805.1", "29022"),    # Myriapoda (centipede)
    ("chamberlinius_hualienensis",   "Chamberlinius hualienensis",   "GCA_054772095.1", "1551368"),  # Myriapoda (millipede)
    ("tubulanus_polymorphus",        "Tubulanus polymorphus",        "GCF_964204645.1", "672921"),   # Nemertea
    ("ptychodera_flava",             "Ptychodera flava",             "GCF_041260155.1", "63121"),    # Hemichordata
    ("trichoplax_sp_h2",             "Trichoplax sp. H2",            "GCA_003344405.1", "287889"),   # Placozoa
]

# ── Lookup-table parsing (shared layout with 01a) ───────────────────────────────

def read_lookup(path: Path, url: str) -> list[dict]:
    if path.exists() and path.stat().st_size > 0 and not _looks_like_html(path):
        text = path.read_text(encoding="utf-8", errors="replace")
    else:
        print(f"   Local lookup missing/invalid — fetching {url}")
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        text = r.text
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) < 14:
            continue
        rows.append({
            "name":            c[0].strip(),
            "species":         c[1].strip(),
            "taxon_id":        c[3].strip(),
            "pan_compara":     c[9].strip(),
            "peptide_compara": c[10].strip(),
            "core_db":         c[13].strip(),
        })
    return rows


def _looks_like_html(path: Path) -> bool:
    head = path.read_text(encoding="utf-8", errors="replace")[:200].lower()
    return "<!doctype html" in head or "<html" in head


def clean_name(name: str) -> str:
    return re.split(r"\s+-\s+|\s*\(", name.strip())[0].strip()


def binomial_key(name: str) -> str:
    toks = clean_name(name).split()
    if len(toks) >= 2:
        return f"{toks[0]}_{toks[1]}".lower()
    return (toks[0].lower() if toks else clean_name(name).lower())


def assembly_score(row: dict) -> tuple:
    sp = row["species"]
    return (
        1 if row["peptide_compara"] == "Y" else 0,
        1 if row["pan_compara"] == "Y" else 0,
        1 if sp.endswith("rs") else 0,
        sp,
    )


def select_rows(rows: list[dict], mode: str) -> list[dict]:
    if mode == "pan_compara":
        rows = [r for r in rows if r["pan_compara"] == "Y"]
    elif mode == "peptide_compara":
        rows = [r for r in rows if r["peptide_compara"] == "Y"]
    best: dict[str, dict] = {}
    for row in rows:
        key = binomial_key(row["name"])
        if key not in best or assembly_score(row) > assembly_score(best[key]):
            best[key] = row
    return sorted(best.values(), key=lambda r: r["species"])


# ── FTP path construction & download (Ensembl) ──────────────────────────────────

def collection_of(core_db: str) -> str | None:
    m = re.match(r"(.+_collection)_core_\d", core_db)
    return m.group(1) if m else None


def pep_dir_url(ftp_dir: str, species: str, core_db: str) -> str:
    coll = collection_of(core_db)
    if coll:
        return f"{ENSEMBL_BASE}/{ftp_dir}/fasta/{coll}/{species}/pep/"
    return f"{ENSEMBL_BASE}/{ftp_dir}/fasta/{species}/pep/"


def find_pep_filename(pep_url: str) -> str | None:
    r = requests.get(pep_url, timeout=120)
    if r.status_code != 200:
        return None
    m = re.search(r'href="([^"]+\.pep\.all\.fa\.gz)"', r.text)
    return m.group(1) if m else None


def download(url: str, dest: Path) -> bool:
    try:
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            tmp = dest.with_suffix(dest.suffix + ".part")
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
            tmp.replace(dest)
        return True
    except requests.exceptions.RequestException as e:
        print(f"      [!] download error: {e}")
        return False


def extract_longest(gz_path: Path, out_fa: Path) -> int:
    """One sequence per Ensembl 'gene:' (the longest isoform)."""
    longest: dict[str, tuple[str, str]] = {}
    order: list[str] = []
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as fh:
        header, gene_id, parts = None, None, []

        def flush():
            if header is None:
                return
            seq = "".join(parts)
            gid = gene_id or header.split()[0]
            cur = longest.get(gid)
            if cur is None:
                longest[gid] = (header, seq)
                order.append(gid)
            elif len(seq) > len(cur[1]):
                longest[gid] = (header, seq)

        for line in fh:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                m = re.search(r"\bgene:(\S+)", header)
                gene_id = m.group(1) if m else None
                parts = []
            else:
                parts.append(line.strip())
        flush()
    _write_fasta(longest, order, out_fa)
    return len(order)


def _write_fasta(records: dict[str, tuple[str, str]], order: list[str], out_fa: Path) -> None:
    with open(out_fa, "w", encoding="utf-8") as f:
        for key in order:
            hdr, seq = records[key]
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


# ── UniProt download ────────────────────────────────────────────────────────────

def download_uniprot(proteome_id: str, taxon_id: str, kingdom_dir: str, out_fa: Path) -> int | None:
    """Fetch a UniProt reference proteome (canonical, one seq per gene) -> out_fa."""
    urls = [
        f"{UNIPROT_FTP}/{kingdom_dir}/{proteome_id}/{proteome_id}_{taxon_id}.fasta.gz",
        f"{UNIPROT_FTP}/{kingdom_dir}/{proteome_id}/{proteome_id}.fasta.gz",
    ]
    gz = None
    for url in urls:
        try:
            r = requests.get(url, timeout=180)
        except requests.exceptions.RequestException as e:
            print(f"      [!] {e}")
            continue
        if r.status_code == 200:
            gz = r.content
            break
    if gz is None:
        print("      [!] UniProt proteome not found (404)")
        return None
    content = gzip.decompress(gz).decode("utf-8", errors="replace")
    records: dict[str, tuple[str, str]] = {}
    order: list[str] = []
    header, parts = None, []
    for line in content.splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = (header, "".join(parts))
                order.append(header)
            header, parts = line[1:].strip(), []
        else:
            parts.append(line.strip())
    if header is not None:
        records[header] = (header, "".join(parts))
        order.append(header)
    _write_fasta(records, order, out_fa)
    return len(order)


# ── NCBI Datasets accession-pinned download ─────────────────────────────────────

def _ncbi_protein_to_gene(gff_text: str) -> dict[str, str]:
    """Map protein accession -> gene key from GFF3 CDS attributes, so isoforms of the
    same gene collapse to one sequence (matching NCBI's protein_coding gene count)."""
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
    """Fetch a protein set by NCBI assembly accession via the Datasets API, keep the
    longest isoform per gene (grouped through the assembly GFF3), and write out_fa.
    Returns the gene count, or None on failure — same one-seq-per-gene contract as the
    UniProt/Ensembl paths. Retries transient network/DNS errors (macOS resolver can
    drop parallel lookups)."""
    content = None
    for attempt in range(5):
        try:
            r = requests.get(NCBI_DATASETS_DL.format(accession), timeout=600)
            r.raise_for_status()
            content = r.content
            break
        except requests.exceptions.RequestException as e:
            if attempt == 4:
                print(f"      [!] NCBI download error after retries: {e}")
                return None
            time.sleep(2 * (attempt + 1))
    try:
        zf = zipfile.ZipFile(io.BytesIO(content))
    except zipfile.BadZipFile as e:
        print(f"      [!] bad zip for {accession}: {e}")
        return None
    faa = gff = None
    for name in zf.namelist():
        if name.endswith("protein.faa"):
            faa = zf.read(name).decode("utf-8", "replace")
        elif name.endswith("genomic.gff"):
            gff = zf.read(name).decode("utf-8", "replace")
    if not faa:
        print(f"      [!] no protein.faa in NCBI package for {accession}")
        return None

    p2g = _ncbi_protein_to_gene(gff) if gff else {}
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
    _write_fasta(longest, order, out_fa)
    return len(order)


# ── NCBI taxonomy resolution (phylum + domain) ──────────────────────────────────

def load_tax_cache() -> dict[str, tuple[str, str]]:
    cache = {}
    if TAX_CACHE.exists():
        for line in TAX_CACHE.read_text(encoding="utf-8").splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 3:
                cache[parts[0]] = (parts[1], parts[2])
    return cache


def save_tax_cache(cache: dict[str, tuple[str, str]]) -> None:
    with open(TAX_CACHE, "w", encoding="utf-8") as f:
        f.write("taxon_id\tdomain\tphylum\n")
        for tid, (domain, phylum) in sorted(cache.items()):
            f.write(f"{tid}\t{domain}\t{phylum}\n")


def _label_from_lineage(lineage: list[tuple[str, str]]) -> tuple[str, str]:
    names = {n for n, _ in lineage}
    ranks = {r: n for n, r in lineage if r and r != "no rank"}

    def has(x):
        return x in names

    if has("Bacteria"):
        return ("Bacteria", "Bacteria")
    if has("Archaea"):
        return ("Archaea", "Archaea")
    if has("Metazoa"):
        phy = ranks.get("phylum")
        if phy == "Arthropoda":
            cls, sub = ranks.get("class"), ranks.get("subphylum")
            if cls == "Insecta":
                return ("Metazoa", "Insecta")
            if sub == "Chelicerata" or cls in {"Arachnida", "Merostomata", "Pycnogonida"} or has("Xiphosura"):
                return ("Metazoa", "Chelicerata")
            if sub == "Myriapoda" or cls in {"Chilopoda", "Diplopoda", "Symphyla", "Pauropoda"}:
                return ("Metazoa", "Myriapoda")
            if sub == "Crustacea" or has("Crustacea") or cls in {
                "Branchiopoda", "Malacostraca", "Copepoda", "Hexanauplia",
                "Thecostraca", "Ostracoda", "Cephalocarida", "Remipedia", "Ichthyostraca",
            }:
                return ("Metazoa", "Crustacea")
            if cls in {"Collembola", "Entognatha", "Diplura", "Protura"} or has("Collembola"):
                return ("Metazoa", "Collembola")
            return ("Metazoa", sub or cls or "Arthropoda")
        return ("Metazoa", phy or "Metazoa_incertae_sedis")

    domain = "Non-metazoan Eukaryota"
    if has("Fungi"):
        return (domain, "Fungi")
    if has("Viridiplantae"):
        return (domain, "Chlorophyta" if has("Chlorophyta") else "Viridiplantae")
    if has("Amoebozoa"):
        return (domain, "Amoebozoa")
    for clade, label in [
        ("Apicomplexa", "Apicomplexa"), ("Ciliophora", "Ciliophora"),
        ("Euglenozoa", "Euglenozoa"), ("Heterolobosea", "Heterolobosea"),
        ("Fornicata", "Metamonada"), ("Metamonada", "Metamonada"),
        ("Parabasalia", "Metamonada"), ("Bacillariophyta", "Bacillariophyta"),
        ("Oomycota", "Oomycota"), ("Rhodophyta", "Rhodophyta"),
        ("Cryptophyceae", "Cryptophyta"), ("Haptophyta", "Haptophyta"),
        ("Rhizaria", "Rhizaria"), ("Cercozoa", "Rhizaria"),
    ]:
        if has(clade):
            return (domain, label)
    return (domain, ranks.get("phylum") or "Protist")


def _efetch_taxonomy(ids: list[str]) -> dict[str, tuple[str, str]]:
    out: dict[str, tuple[str, str]] = {}
    data = {"db": "taxonomy", "id": ",".join(ids), "retmode": "xml"}
    for attempt in range(4):
        try:
            r = requests.post(EUTILS_EFETCH, data=data, timeout=120)
            r.raise_for_status()
            root = ET.fromstring(r.text)
            break
        except (requests.exceptions.RequestException, ET.ParseError) as e:
            if attempt == 3:
                print(f"   [!] NCBI taxonomy fetch failed: {e}")
                return out
            time.sleep(2 * (attempt + 1))
    for taxon in root.findall("Taxon"):
        tid = taxon.findtext("TaxId")
        lineage = []
        lex = taxon.find("LineageEx")
        if lex is not None:
            for t in lex.findall("Taxon"):
                lineage.append((t.findtext("ScientificName"), t.findtext("Rank")))
        lineage.append((taxon.findtext("ScientificName"), taxon.findtext("Rank")))
        if tid:
            out[tid] = _label_from_lineage(lineage)
    return out


def resolve_taxonomy(taxon_ids: list[str], cache: dict[str, tuple[str, str]]) -> None:
    todo = sorted({t for t in taxon_ids if t and t not in cache})
    if not todo:
        return
    print(f"Resolving {len(todo)} taxon ids from NCBI Taxonomy...")
    for i in range(0, len(todo), 200):
        cache.update(_efetch_taxonomy(todo[i:i + 200]))
        save_tax_cache(cache)
        print(f"   {min(i + 200, len(todo))}/{len(todo)} resolved")
        time.sleep(0.34)


# ── Manifest upsert (preserves rows owned by 01a) ───────────────────────────────

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
    with open(MANIFEST, "w", encoding="utf-8") as f:
        f.write("\t".join(MANIFEST_HEADER) + "\n")
        for key in sorted(existing):
            rec = existing[key]
            f.write("\t".join(str(rec.get(col, "")) for col in MANIFEST_HEADER) + "\n")


def _row(species_key, display, phylum, out_fa, domain, taxon_id) -> dict:
    return {
        "species_key":  species_key,
        "display_name": display,
        "phylum":       phylum,
        "fa_path":      str(out_fa),
        "domain":       domain,
        "taxon_id":     taxon_id,
    }


# ── Per-species workers (thread-pool safe: own files, own HTTP requests) ────────

def process_uniprot(n: int, total: int, task: tuple, cache: dict) -> tuple:
    key, display, kingdom_dir, pid, taxon_id = task
    domain, phylum = cache.get(taxon_id, ("Non-metazoan Eukaryota", "Unknown"))
    log = [f"[U {n}/{total}] {display}  ({phylum})"]
    out_fa = OUT_DIR / f"{key}.longest.fa"
    if out_fa.exists():
        log.append(f"      [-] {out_fa.name} already present")
        return _row(key, display, phylum, out_fa, domain, taxon_id), (key, "EXISTS"), "\n".join(log)
    n_seq = download_uniprot(pid, taxon_id, kingdom_dir, out_fa)
    if n_seq is None:
        return None, (key, "UNIPROT_404"), "\n".join(log)
    log.append(f"      [+] {n_seq} genes -> {out_fa.name}")
    return _row(key, display, phylum, out_fa, domain, taxon_id), (key, "OK"), "\n".join(log)


def process_ensembl(n: int, total: int, task: tuple, cache: dict) -> tuple:
    div, ftp_dir, species, display, core_db, taxon_id = task
    domain, phylum = cache.get(taxon_id, ("Non-metazoan Eukaryota", "Unknown"))
    log = [f"[E {n}/{total}] [{div}] {display}  ({phylum})"]
    out_fa = OUT_DIR / f"{species}.longest.fa"
    if out_fa.exists():
        log.append(f"      [-] {out_fa.name} already present")
        return _row(species, display, phylum, out_fa, domain, taxon_id), (species, "EXISTS"), "\n".join(log)
    pep_url  = pep_dir_url(ftp_dir, species, core_db)
    filename = find_pep_filename(pep_url)
    if not filename:
        log.append(f"      [!] no .pep.all.fa.gz under {pep_url}")
        return None, (species, "NO_PEP_FILE"), "\n".join(log)
    gz_path = OUT_DIR / filename
    if not gz_path.exists() and not download(pep_url + filename, gz_path):
        return None, (species, "DOWNLOAD_FAILED"), "\n".join(log)
    n_seq = extract_longest(gz_path, out_fa)
    log.append(f"      [+] {n_seq} genes -> {out_fa.name}")
    return _row(species, display, phylum, out_fa, domain, taxon_id), (species, "OK"), "\n".join(log)


def process_ncbi(n: int, total: int, task: tuple, cache: dict) -> tuple:
    key, display, accession, taxon_id = task
    domain, phylum = cache.get(taxon_id, ("Non-metazoan Eukaryota", "Unknown"))
    log = [f"[N {n}/{total}] {display}  ({phylum})  {accession}"]
    out_fa = OUT_DIR / f"{key}.longest.fa"
    if out_fa.exists():
        log.append(f"      [-] {out_fa.name} already present")
        return _row(key, display, phylum, out_fa, domain, taxon_id), (key, "EXISTS"), "\n".join(log)
    n_seq = download_ncbi(accession, out_fa)
    if n_seq is None:
        return None, (key, "NCBI_FAILED"), "\n".join(log)
    log.append(f"      [+] {n_seq} genes -> {out_fa.name}")
    return _row(key, display, phylum, out_fa, domain, taxon_id), (key, "OK"), "\n".join(log)


def run_pool(worker, tasks: list, cache: dict, workers: int) -> tuple[list, list]:
    """Run `worker(n, total, task, cache)` across a thread pool; collect rows+status."""
    manifest_rows, status = [], []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(worker, n, len(tasks), task, cache)
                   for n, task in enumerate(tasks, 1)]
        for fut in as_completed(futures):
            try:
                mrow, st, log = fut.result()
            except Exception as e:                       # one species must not sink the run
                print(f"      [!] worker error: {e}")
                continue
            print(log)
            if mrow is not None:
                manifest_rows.append(mrow)
            status.append(st)
    return manifest_rows, status


# ── Main ────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Download non-metazoan outgroup proteomes.")
    ap.add_argument("--limit", type=int, default=0, help="cap species per source (testing)")
    ap.add_argument("--workers", type=int, default=8, help="parallel download workers (default 8)")
    args = ap.parse_args()

    print(f"Outgroup proteomes — UniProt + Ensembl Genomes r{RELEASE} + NCBI Datasets\n")

    # 1. Accession-pinned NCBI task list (recoveries + invertebrate expansion). Its keys
    #    own those species: the UniProt pass skips them (their proteome ids 404), and its
    #    binomials make the Ensembl divisions skip the same species too.
    ncbi_tasks = list(NCBI_PROTEOMES)          # (key, display, accession, taxon_id)
    ncbi_keys = {t[0] for t in ncbi_tasks}

    # 2. Build the UniProt task list, skipping NCBI-owned keys; record every binomial so
    #    the Ensembl divisions skip duplicates (of both UniProt and NCBI species).
    seen_binomials: set[str] = set()
    uniprot_tasks = []
    for key, (display, kingdom_dir, pid, taxon_id) in UNIPROT_OUTGROUPS.items():
        seen_binomials.add(binomial_key(display))
        if key in ncbi_keys:
            continue
        uniprot_tasks.append((key, display, kingdom_dir, pid, taxon_id))
    for _, display, _, _ in ncbi_tasks:
        seen_binomials.add(binomial_key(display))
    if args.limit:
        uniprot_tasks = uniprot_tasks[:args.limit]
        ncbi_tasks = ncbi_tasks[:args.limit]

    # 3. Build the Ensembl task list per division (filtered, deduped, non-overlapping).
    ensembl_tasks = []   # (division, ftp_dir, species, display, core_db, taxon_id)
    for div, cfg in ENSEMBL_DIVISIONS.items():
        path = PROJECT_DIR / cfg["file"]
        url  = f"{ENSEMBL_BASE}/{cfg['ftp_dir']}/species_Ensembl{div.capitalize()}.txt"
        rows = select_rows(read_lookup(path, url), cfg["select"])
        kept = 0
        for row in rows:
            bk = binomial_key(row["name"])
            if bk in seen_binomials:
                continue
            seen_binomials.add(bk)
            ensembl_tasks.append((div, cfg["ftp_dir"], row["species"],
                                  clean_name(row["name"]), row["core_db"], row["taxon_id"]))
            kept += 1
        print(f"   {div:<9} select={cfg['select']:<15} -> {kept} new species")
        if args.limit:
            ensembl_tasks = ensembl_tasks[:args.limit]
    print()

    # 4. Resolve all taxonomy up front (shared cache with 01a).
    cache = load_tax_cache()
    all_ids = ([t[4] for t in uniprot_tasks] + [t[5] for t in ensembl_tasks]
               + [t[3] for t in ncbi_tasks])
    resolve_taxonomy(all_ids, cache)

    # 5. UniProt downloads (parallel).
    print(f"\n== UniProt reference proteomes ({args.workers} workers) ==")
    u_rows, u_status = run_pool(process_uniprot, uniprot_tasks, cache, args.workers)

    # 6. Ensembl downloads (parallel).
    print(f"\n== Ensembl Genomes release-63 ({args.workers} workers) ==")
    e_rows, e_status = run_pool(process_ensembl, ensembl_tasks, cache, args.workers)

    # 7. NCBI Datasets downloads (parallel) — recoveries + invertebrate expansion.
    print(f"\n== NCBI Datasets accession-pinned ({args.workers} workers) ==")
    n_rows, n_status = run_pool(process_ncbi, ncbi_tasks, cache, args.workers)

    manifest_rows = u_rows + e_rows + n_rows
    status = u_status + e_status + n_status

    # 8. Register everything; metazoan rows from 01a are preserved by species_key.
    upsert_manifest(manifest_rows)
    with open(STATUS_FILE, "w", encoding="utf-8") as f:
        f.write("species_key\tstatus\n")
        for sp, st in status:
            f.write(f"{sp}\t{st}\n")

    ok = sum(1 for _, st in status if st in ("OK", "EXISTS"))
    total = len(uniprot_tasks) + len(ensembl_tasks) + len(ncbi_tasks)
    print(f"\nDone: {ok}/{total} outgroup species ready. Manifest: {MANIFEST}")
    print(f"Status log: {STATUS_FILE}")
    print("Next: bash scripts/02_run_flps.sh  (picks up new .longest.fa files)")


if __name__ == "__main__":
    sys.exit(main())
