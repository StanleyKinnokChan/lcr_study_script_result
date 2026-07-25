#!/usr/bin/env python3
"""
Fetch UniProt signal-peptide + subcellular-location annotations for the
prokaryotic species in the study, producing the annotation TSV consumed by
17_signal_peptide_stratification.py (MODE A).

For every Bacteria/Archaea species in results/species_manifest.tsv the script
queries the UniProtKB REST API by organism taxon id and records, per protein:
  - has_signal_peptide : True if a SIGNAL sequence feature is annotated
  - subcellular_location : the reported subcellular location term(s), if any

Output (appended/updated incrementally, one row per protein):
  results/uniprot_signal_peptides.tsv
      protein_id  has_signal_peptide  subcellular_location  species_key

`protein_id` is the UniProt accession (e.g. C0H3V3) so that it matches the
accession extracted from the fLPS `sp|ACC|NAME` identifiers in script 17.

Notes
-----
* Species whose proteomes did NOT come from UniProt (Ensembl Genomes / NCBI
  Datasets sources in 01b) will usually return no hits by strain taxon id; they
  are simply left unannotated and reported as such — the stratification in
  script 17 handles missing annotations gracefully.
* The run is idempotent: species already present in the output are skipped,
  so the fetch can be resumed after an interruption.

Usage:
  python3 scripts/fetch_uniprot_annotations.py            # Bacteria + Archaea
  python3 scripts/fetch_uniprot_annotations.py --domain Bacteria
  python3 scripts/fetch_uniprot_annotations.py --force    # ignore existing rows
"""

import argparse
import io
import sys
import time
import pandas as pd
import requests
from pathlib import Path

from config import RESULTS_DIR

OUT_FILE = RESULTS_DIR / "uniprot_signal_peptides.tsv"
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,ft_signal,cc_subcellular_location"
REQUEST_TIMEOUT = 120
MAX_RETRIES = 4
RETRY_BACKOFF = 5   # seconds, multiplied by attempt number


def fetch_species(taxon_id: int) -> pd.DataFrame:
    """
    Query UniProtKB for all entries of one organism taxon id.
    Returns a DataFrame with columns [accession, ft_signal, subcellular].
    Empty DataFrame if the organism has no UniProt entries.
    """
    params = {
        "query": f"organism_id:{taxon_id}",
        "fields": FIELDS,
        "format": "tsv",
    }
    last_err = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(UNIPROT_STREAM, params=params, timeout=REQUEST_TIMEOUT)
            if r.status_code == 200:
                if not r.text.strip():
                    return pd.DataFrame(columns=["accession", "ft_signal", "subcellular"])
                df = pd.read_csv(io.StringIO(r.text), sep="\t")
                df.columns = ["accession", "ft_signal", "subcellular"][: len(df.columns)]
                return df
            # 429 / 5xx → back off and retry
            if r.status_code in (429, 500, 502, 503, 504):
                last_err = f"HTTP {r.status_code}"
                time.sleep(RETRY_BACKOFF * attempt)
                continue
            last_err = f"HTTP {r.status_code}: {r.text[:200]}"
            break
        except requests.RequestException as e:
            last_err = str(e)
            time.sleep(RETRY_BACKOFF * attempt)
    print(f"    WARNING: fetch failed (taxon {taxon_id}): {last_err}")
    return pd.DataFrame(columns=["accession", "ft_signal", "subcellular"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", nargs="+", default=["Bacteria", "Archaea"],
                    help="Domains to fetch (default: Bacteria Archaea)")
    ap.add_argument("--force", action="store_true",
                    help="Re-fetch species already present in the output file")
    args = ap.parse_args()

    manifest_path = RESULTS_DIR / "species_manifest.tsv"
    if not manifest_path.exists():
        sys.exit(f"ERROR: {manifest_path} missing — run the download phase first.")

    manifest = pd.read_csv(manifest_path, sep="\t")
    domain_col = "domain" if "domain" in manifest.columns else "phylum"
    targets = manifest[manifest[domain_col].isin(args.domain)].copy()
    if "taxon_id" not in targets.columns:
        sys.exit("ERROR: manifest has no taxon_id column — cannot query UniProt.")
    targets = targets.dropna(subset=["taxon_id"])

    if targets.empty:
        sys.exit(f"No species found for domains {args.domain}.")

    # Resume support: skip species already written
    done_species: set[str] = set()
    existing = None
    if OUT_FILE.exists() and not args.force:
        existing = pd.read_csv(OUT_FILE, sep="\t")
        if "species_key" in existing.columns:
            done_species = set(existing["species_key"].unique())
        print(f"Resuming: {len(done_species)} species already in {OUT_FILE.name}")

    print(f"Fetching UniProt annotations for {len(targets)} species "
          f"({', '.join(args.domain)})")

    collected: list[pd.DataFrame] = []
    for i, (_, row) in enumerate(targets.iterrows(), 1):
        sp_key = row["species_key"]
        if sp_key in done_species:
            continue
        taxon_id = int(row["taxon_id"])
        df = fetch_species(taxon_id)

        if df.empty:
            n_prot = n_sp = 0
            out = pd.DataFrame(columns=["protein_id", "has_signal_peptide",
                                        "subcellular_location", "species_key"])
        else:
            has_sp = df["ft_signal"].fillna("").astype(str).str.contains("SIGNAL")
            out = pd.DataFrame({
                "protein_id":          df["accession"],
                "has_signal_peptide":  has_sp,
                "subcellular_location": (df["subcellular"].fillna("")
                                         if "subcellular" in df.columns else ""),
                "species_key":         sp_key,
            })
            n_prot = len(out)
            n_sp = int(has_sp.sum())

        collected.append(out)
        print(f"  [{i}/{len(targets)}] {row['display_name']}: "
              f"{n_prot:,} proteins, {n_sp:,} with signal peptide")

        # Flush incrementally so long runs are crash-safe
        if collected:
            batch = pd.concat(collected, ignore_index=True)
            if existing is not None:
                batch = pd.concat([existing, batch], ignore_index=True)
            batch.drop_duplicates(subset=["protein_id", "species_key"],
                                  keep="last", inplace=True)
            batch.to_csv(OUT_FILE, sep="\t", index=False)

    if not collected:
        print("Nothing new to fetch — all target species already annotated.")
        return

    final = pd.read_csv(OUT_FILE, sep="\t")
    n_total = len(final)
    n_with_sp = int(final["has_signal_peptide"].astype(str).isin(["True", "true"]).sum())
    n_species = final["species_key"].nunique()
    print(f"\nWrote {OUT_FILE}")
    print(f"  {n_total:,} proteins across {n_species} species; "
          f"{n_with_sp:,} annotated with a signal peptide")


if __name__ == "__main__":
    main()
