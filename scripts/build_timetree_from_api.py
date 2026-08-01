#!/usr/bin/env python3
"""
Build a RESOLVED, dated Viridiplantae phylogeny from the PGLS species list using
TimeTree divergence data, writing results/viridiplantae_timetree.nwk.

This is the higher-resolution alternative to build_viridiplantae_backbone_tree.py.
It must be run in an environment with network access (it calls the TimeTree API).

Strategy (hybrid)
-----------------
The deep land-plant splits (green-algae / streptophyte / land-plant / angiosperm /
monocot-eudicot / grass) are already well established, so those are fixed from the
backbone calibration. What the backbone lacks is WITHIN-tier resolution (the 55
eudicot and 17 grass genera sit in flat polytomies). This script therefore queries
TimeTree's pairwise API — using the NCBI taxon IDs already in supp_table_S1_species_list.tsv —
only for WITHIN-tier pairs, then builds a UPGMA tree from the combined distance
matrix (API divergence times within tiers; backbone node ages between tiers).
UPGMA yields an ultrametric tree, which is what the PGLS VCV in
19_pgls_viridiplantae.py expects.

Missing pairs (a taxon TimeTree cannot resolve) fall back to the backbone age for
that tier, so the script always produces a valid tree; it reports API coverage so
you can judge how much real resolution was obtained. If coverage is poor, prefer
the timetree.org "Load a List of Species" web upload of results/pgls_species_list.txt.

Caching: pairwise results are cached in results/timetree_pairwise_cache.json, so the
run is resumable and re-runs are fast.

Output:
  results/viridiplantae_timetree.nwk   (leaf names = study species_key)
"""

import json
import time
import itertools
import pandas as pd
import requests
from pathlib import Path

from config import RESULTS_DIR

PGLS_TABLE  = RESULTS_DIR / "pgls_viridiplantae.tsv"
MANIFEST    = RESULTS_DIR / "supp_table_S1_species_list.tsv"
CACHE_FILE  = RESULTS_DIR / "timetree_pairwise_cache.json"
OUT_TREE    = RESULTS_DIR / "viridiplantae_timetree.nwk"

API_URL     = "http://timetree.org/api/pairwise/{a}/{b}"
REQUEST_TIMEOUT = 30
SLEEP_BETWEEN   = 0.12          # be polite to the API
MAX_RETRIES     = 3

# Tier crown ages (Ma) and the between-tier split ages, mirroring the backbone
# builder. Between-tier divergence = age of the shallowest clade containing both.
TIER_CROWN = {0: 750.0, 1: 400.0, 2: 430.0, 3: 350.0,
              6: 135.0, 7: 120.0, 8: 110.0, 9: 65.0}
# (age, set_of_tiers) sorted shallow → deep
NESTED_CLADES = [
    (130.0, {8, 9}),
    (160.0, {7, 8, 9}),
    (175.0, {6, 7, 8, 9}),
    (430.0, {3, 6, 7, 8, 9}),
    (480.0, {2, 3, 6, 7, 8, 9}),
    (850.0, {1, 2, 3, 6, 7, 8, 9}),
    (1100.0, {0, 1, 2, 3, 6, 7, 8, 9}),
]


def divergence_age(tier_a: int, tier_b: int) -> float:
    """Backbone divergence age (Ma) between two tiers."""
    if tier_a == tier_b:
        return TIER_CROWN[tier_a]
    for age, tiers in NESTED_CLADES:
        if tier_a in tiers and tier_b in tiers:
            return age
    return NESTED_CLADES[-1][0]


def load_cache() -> dict:
    if CACHE_FILE.exists():
        try:
            return json.loads(CACHE_FILE.read_text())
        except Exception:
            return {}
    return {}


def parse_mya(resp: requests.Response) -> float | None:
    """Defensively extract a divergence time (Ma) from a TimeTree API response."""
    try:
        data = resp.json()
    except ValueError:
        import re
        m = re.search(r'(\d+(?:\.\d+)?)', resp.text)
        return float(m.group(1)) if m else None
    if isinstance(data, dict):
        for key in ("median_time", "mya", "time", "adjusted_age",
                    "precomputed_age", "median"):
            if key in data and data[key] not in (None, "", "NA"):
                try:
                    return float(data[key])
                except (TypeError, ValueError):
                    pass
    return None


def query_pairwise(taxid_a: str, taxid_b: str, cache: dict) -> float | None:
    key = f"{taxid_a}:{taxid_b}"
    rkey = f"{taxid_b}:{taxid_a}"
    if key in cache:
        return cache[key]
    if rkey in cache:
        return cache[rkey]
    mya = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(API_URL.format(a=taxid_a, b=taxid_b),
                             timeout=REQUEST_TIMEOUT)
            if r.status_code == 200:
                mya = parse_mya(r)
                break
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(SLEEP_BETWEEN * attempt * 4)
                continue
            break
        except requests.RequestException:
            time.sleep(SLEEP_BETWEEN * attempt * 4)
    cache[key] = mya
    time.sleep(SLEEP_BETWEEN)
    return mya


def main():
    if not PGLS_TABLE.exists():
        raise SystemExit(f"ERROR: {PGLS_TABLE} missing — run 19_pgls_viridiplantae.py first.")

    pgls = pd.read_csv(PGLS_TABLE, sep="\t")
    manifest = pd.read_csv(MANIFEST, sep="\t")[["species_key", "taxon_id"]]
    df = pgls.merge(manifest, on="species_key", how="left")
    df = df.dropna(subset=["taxon_id"]).reset_index(drop=True)
    df["taxon_id"] = df["taxon_id"].astype(int).astype(str)
    df["tier"] = df["tier"].astype(int)

    species = df["species_key"].tolist()
    tier    = dict(zip(df["species_key"], df["tier"]))
    taxid   = dict(zip(df["species_key"], df["taxon_id"]))
    n = len(species)
    print(f"Building resolved tree for {n} species.")

    cache = load_cache()

    # Distances: API within-tier (multi-species tiers), backbone between-tier.
    within_pairs = [
        (a, b) for a, b in itertools.combinations(species, 2)
        if tier[a] == tier[b]
    ]
    # only query tiers with >1 species (others have no within-tier pairs anyway)
    print(f"Within-tier pairs to resolve via TimeTree API: {len(within_pairs)}")

    api_hits = 0
    dist = {}   # (i,j) -> distance (= 2 * divergence age)
    for idx, (a, b) in enumerate(within_pairs, 1):
        mya = query_pairwise(taxid[a], taxid[b], cache)
        if mya is None:
            mya = divergence_age(tier[a], tier[b])   # fall back to tier crown
        else:
            api_hits += 1
        dist[(a, b)] = 2.0 * mya
        if idx % 100 == 0:
            print(f"  {idx}/{len(within_pairs)} pairs  (API hits so far: {api_hits})")
            CACHE_FILE.write_text(json.dumps(cache))
    CACHE_FILE.write_text(json.dumps(cache))

    coverage = (api_hits / len(within_pairs) * 100) if within_pairs else 0.0
    print(f"\nTimeTree API coverage of within-tier pairs: {api_hits}/{len(within_pairs)} "
          f"({coverage:.1f}%)")
    if coverage < 50:
        print("  WARNING: low API coverage — the tree is close to the backbone. "
              "Consider the timetree.org web upload of results/pgls_species_list.txt.")

    # Between-tier distances from the backbone calibration.
    for a, b in itertools.combinations(species, 2):
        if (a, b) not in dist:
            dist[(a, b)] = 2.0 * divergence_age(tier[a], tier[b])

    # ── UPGMA (ultrametric) via Biopython ─────────────────────────────────────
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from Bio import Phylo

    # lower-triangular matrix (row i has i+1 entries; diagonal = 0)
    matrix = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                a, b = species[j], species[i]
                row.append(dist.get((a, b)) or dist.get((b, a)))
        matrix.append(row)

    dm = DistanceMatrix(names=species, matrix=matrix)
    tree = DistanceTreeConstructor().upgma(dm)
    # UPGMA labels internal nodes; clear them so downstream parsers see only leaves
    for clade in tree.get_nonterminals():
        clade.name = None
    Phylo.write(tree, str(OUT_TREE), "newick")

    print(f"\nWrote {OUT_TREE}  ({n} leaves, UPGMA / ultrametric).")
    print("Re-run: uv run python scripts/19_pgls_viridiplantae.py")


if __name__ == "__main__":
    main()
