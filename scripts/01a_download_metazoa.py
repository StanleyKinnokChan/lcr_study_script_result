#!/usr/bin/env python3
"""
Download longest-isoform protein FASTAs for the full Ensembl Metazoa release-63
species set, and register them in results/supp_table_S1_species_list.tsv for the LCR pipeline.

Pipeline contract (must hold for 02_run_flps.sh and 03_analyse.py downstream):
  - Every species ends up as data/proteomes/{species_key}.longest.fa   (one seq per gene)
  - Every species has a row in results/supp_table_S1_species_list.tsv with columns:
      species_key  display_name  phylum  fa_path  domain  taxon_id

Species source
  The release-63 lookup table species_EnsemblMetazoa.txt (project root) is the source
  of truth: it lists every species in the release with its exact assembly-versioned
  `species` id and its `core_db` (which encodes collection membership, needed for the
  FTP path). Release 63 is guaranteed to have FASTA on the FTP, so we pin to it.
  If the local file is absent, the same table is fetched from the FTP as a fallback.

  All unique binomial species are taken (one best assembly per binomial), giving the
  broadest possible metazoan ingroup. Redundant strains / alternate assemblies of the
  same species collapse to the best-annotated assembly.

Phylum / domain
  The lookup table has no phylum column, only taxonomy_id. Phylum and domain are
  resolved from NCBI Taxonomy (batched E-utilities, cached to results/taxonomy_cache.tsv
  so reruns are offline). Arthropoda is split into Insecta / Crustacea / Chelicerata /
  Myriapoda to match the manifest convention used by the existing dataset.

FTP layout (release 63), confirmed live:
  standalone species : {division}/fasta/{species}/pep/{...}.pep.all.fa.gz
  collection species : {division}/fasta/{collection}/{species}/pep/{...}.pep.all.fa.gz
  where {collection} is the part of core_db before "_core_<release>_..." (it always
  ends in "_collection").
"""

import argparse
import gzip
import re
import sys
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

# ── Paths & config ──────────────────────────────────────────────────────────────

RELEASE      = 63
DIVISION_FTP = "metazoa"                       # FTP directory name for this division
DOMAIN_LABEL = "Metazoa"                       # all rows written by this script

PROJECT_DIR = Path(__file__).parent.parent
LOOKUP_FILE = PROJECT_DIR / "species_EnsemblMetazoa.txt"
OUT_DIR     = PROJECT_DIR / "data" / "proteomes"
RESULTS_DIR = PROJECT_DIR / "results"
MANIFEST    = RESULTS_DIR / "supp_table_S1_species_list.tsv"
TAX_CACHE   = RESULTS_DIR / "taxonomy_cache.tsv"
STATUS_FILE = RESULTS_DIR / "download_status_metazoa.tsv"

BASE_URL      = f"https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-{RELEASE}"
LOOKUP_URL    = f"{BASE_URL}/{DIVISION_FTP}/species_EnsemblMetazoa.txt"
EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

MANIFEST_HEADER = ["species_key", "display_name", "phylum", "fa_path", "domain", "taxon_id"]

OUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Lookup-table parsing ────────────────────────────────────────────────────────
# Columns (tab-separated): 0 name, 1 species, 2 division, 3 taxonomy_id, 4 assembly,
# 5 assembly_accession, 6 genebuild, 7 variation, 8 microarray, 9 pan_compara,
# 10 peptide_compara, 11 genome_alignments, 12 other_alignments, 13 core_db, 14 species_id

def read_lookup(path: Path, url: str) -> list[dict]:
    """Read the release lookup table from the local file, or fetch it from the FTP."""
    if path.exists() and path.stat().st_size > 0 and not _looks_like_html(path):
        text = path.read_text(encoding="utf-8", errors="replace")
        print(f"Using local lookup table: {path.name}")
    else:
        print(f"Local lookup missing/invalid — fetching {url}")
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
    """'Aedes albopictus - (Asian tiger mosquito, Foshan)' -> 'Aedes albopictus'."""
    return re.split(r"\s+-\s+|\s*\(", name.strip())[0].strip()


def binomial_key(name: str) -> str:
    """Genus+species, lowercased ('genus_species'); collapses strains/subspecies."""
    toks = clean_name(name).split()
    if len(toks) >= 2:
        return f"{toks[0]}_{toks[1]}".lower()
    return (toks[0].lower() if toks else clean_name(name).lower())


def assembly_score(row: dict) -> tuple:
    """Higher is better: prefer gene-tree annotation, pan-compara, RefSeq builds."""
    sp = row["species"]
    return (
        1 if row["peptide_compara"] == "Y" else 0,
        1 if row["pan_compara"] == "Y" else 0,
        1 if sp.endswith("rs") else 0,           # RefSeq-annotated assemblies
        sp,                                      # stable tie-break
    )


def dedup_best_per_species(rows: list[dict]) -> list[dict]:
    best: dict[str, dict] = {}
    for row in rows:
        key = binomial_key(row["name"])
        if key not in best or assembly_score(row) > assembly_score(best[key]):
            best[key] = row
    return sorted(best.values(), key=lambda r: r["species"])


# ── FTP path construction & download ────────────────────────────────────────────

def collection_of(core_db: str) -> str | None:
    """'fungi_ascomycota5_collection_core_63_116_1' -> 'fungi_ascomycota5_collection'."""
    m = re.match(r"(.+_collection)_core_\d", core_db)
    return m.group(1) if m else None


def pep_dir_url(species: str, core_db: str) -> str:
    coll = collection_of(core_db)
    if coll:
        return f"{BASE_URL}/{DIVISION_FTP}/fasta/{coll}/{species}/pep/"
    return f"{BASE_URL}/{DIVISION_FTP}/fasta/{species}/pep/"


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
    """Write one sequence per gene (the longest isoform) to out_fa.

    Ensembl pep headers carry 'gene:<id>'. Records sharing a gene id collapse to their
    longest sequence; records without a gene id are each kept as their own gene.
    """
    longest: dict[str, tuple[str, str]] = {}   # gene_id -> (header, sequence)
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

    with open(out_fa, "w", encoding="utf-8") as f:
        for gid in order:
            hdr, seq = longest[gid]
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    return len(order)


# ── NCBI taxonomy resolution (phylum + domain) ──────────────────────────────────

METAZOAN_FALLBACK_PHYLA = {
    "Porifera", "Ctenophora", "Placozoa", "Cnidaria", "Platyhelminthes", "Nematoda",
    "Annelida", "Mollusca", "Arthropoda", "Echinodermata", "Hemichordata", "Chordata",
    "Brachiopoda", "Priapulida", "Rotifera", "Tardigrada", "Nemertea", "Nematomorpha",
    "Acanthocephala", "Xenacoelomorpha", "Bryozoa", "Entoprocta", "Gastrotricha",
    "Chaetognatha", "Kinorhyncha", "Onychophora", "Phoronida", "Loricifera",
}


def load_tax_cache() -> dict[str, tuple[str, str]]:
    cache = {}
    if TAX_CACHE.exists():
        for line in TAX_CACHE.read_text(encoding="utf-8").splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 3:
                cache[parts[0]] = (parts[1], parts[2])   # taxon_id -> (domain, phylum)
    return cache


def save_tax_cache(cache: dict[str, tuple[str, str]]) -> None:
    with open(TAX_CACHE, "w", encoding="utf-8") as f:
        f.write("taxon_id\tdomain\tphylum\n")
        for tid, (domain, phylum) in sorted(cache.items()):
            f.write(f"{tid}\t{domain}\t{phylum}\n")


def _label_from_lineage(lineage: list[tuple[str, str]]) -> tuple[str, str]:
    """lineage: ordered [(scientific_name, rank), ...] incl. the organism itself.
    Returns (domain, phylum_label)."""
    names = {n for n, _ in lineage}
    ranks = {r: n for n, r in lineage if r and r != "no rank"}

    def has(x):
        return x in names

    if has("Bacteria"):
        return ("Bacteria", "Bacteria")
    if has("Archaea"):
        return ("Archaea", "Archaea")

    # Eukaryota below this point
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
    """One batched E-utilities call -> {taxon_id: (domain, phylum)}."""
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
        batch = todo[i:i + 200]
        cache.update(_efetch_taxonomy(batch))
        save_tax_cache(cache)
        print(f"   {min(i + 200, len(todo))}/{len(todo)} resolved")
        time.sleep(0.34)   # stay under NCBI's 3 req/s anonymous limit


# ── Manifest upsert ─────────────────────────────────────────────────────────────

def upsert_manifest(rows: list[dict]) -> None:
    """Insert/replace rows by species_key, preserving rows owned by other scripts."""
    existing: dict[str, dict] = {}
    if MANIFEST.exists() and MANIFEST.stat().st_size > 0:
        lines = MANIFEST.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        for line in lines[1:]:
            vals = line.split("\t")
            rec = dict(zip(header, vals))
            existing[rec.get("species_key", "")] = rec
    for r in rows:
        existing[r["species_key"]] = r
    with open(MANIFEST, "w", encoding="utf-8") as f:
        f.write("\t".join(MANIFEST_HEADER) + "\n")
        for key in sorted(existing):
            rec = existing[key]
            f.write("\t".join(str(rec.get(col, "")) for col in MANIFEST_HEADER) + "\n")


# ── Main ────────────────────────────────────────────────────────────────────────

def process_species(n: int, total: int, row: dict, cache: dict) -> tuple:
    """Download + extract one species. Pure per-species work (own files, own HTTP
    requests) so it is safe to run many at once in a thread pool. Returns
    (manifest_row | None, (species_key, status), log_text)."""
    species = row["species"]
    display = clean_name(row["name"])
    domain, phylum = cache.get(row["taxon_id"], (DOMAIN_LABEL, "Unknown"))
    log = [f"[{n}/{total}] {display}  ({phylum})"]

    out_fa = OUT_DIR / f"{species}.longest.fa"
    if out_fa.exists():
        log.append(f"      [-] {out_fa.name} already present")
        return (_row(species, display, phylum, out_fa, domain, row["taxon_id"]),
                (species, "EXISTS"), "\n".join(log))

    pep_url  = pep_dir_url(species, row["core_db"])
    filename = find_pep_filename(pep_url)
    if not filename:
        log.append(f"      [!] no .pep.all.fa.gz under {pep_url}")
        return None, (species, "NO_PEP_FILE"), "\n".join(log)

    gz_path = OUT_DIR / filename
    if not gz_path.exists() and not download(pep_url + filename, gz_path):
        return None, (species, "DOWNLOAD_FAILED"), "\n".join(log)

    n_seq = extract_longest(gz_path, out_fa)
    log.append(f"      [+] {n_seq} genes -> {out_fa.name}")
    return (_row(species, display, phylum, out_fa, domain, row["taxon_id"]),
            (species, "OK"), "\n".join(log))


def main():
    ap = argparse.ArgumentParser(description="Download Ensembl Metazoa release-63 proteomes.")
    ap.add_argument("--limit", type=int, default=0, help="process only the first N species (testing)")
    ap.add_argument("--workers", type=int, default=8, help="parallel download workers (default 8)")
    args = ap.parse_args()

    print(f"Ensembl Metazoa release {RELEASE} — full-species proteome download\n")
    rows = dedup_best_per_species(read_lookup(LOOKUP_FILE, LOOKUP_URL))
    if args.limit:
        rows = rows[:args.limit]
    print(f"{len(rows)} unique species selected.\n")

    cache = load_tax_cache()
    resolve_taxonomy([r["taxon_id"] for r in rows], cache)

    print(f"Downloading with {args.workers} parallel workers...\n")
    manifest_rows, status = [], []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(process_species, n, len(rows), row, cache)
                   for n, row in enumerate(rows, 1)]
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

    upsert_manifest(manifest_rows)

    with open(STATUS_FILE, "w", encoding="utf-8") as f:
        f.write("species_key\tstatus\n")
        for sp, st in status:
            f.write(f"{sp}\t{st}\n")

    ok = sum(1 for _, st in status if st in ("OK", "EXISTS"))
    print(f"\nDone: {ok}/{len(rows)} species ready. Manifest: {MANIFEST}")
    print(f"Status log: {STATUS_FILE}")
    print("Next: python scripts/01b_download_outgroups.py, then bash scripts/02_run_flps.sh")


def _row(species, display, phylum, out_fa, domain, taxon_id) -> dict:
    return {
        "species_key":  species,
        "display_name": display,
        "phylum":       phylum,
        "fa_path":      str(out_fa),
        "domain":       domain,
        "taxon_id":     taxon_id,
    }


if __name__ == "__main__":
    sys.exit(main())
