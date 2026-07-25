#!/usr/bin/env python3
"""
Build a dated phylogeny for all species in the study and identify large temporal gaps
in coverage. Divergence times are sourced from the TimeTree database (timetree.org)
via their public REST API, supplemented by hard-coded consensus dates for key nodes
where the API is unavailable or returns no estimate.

Outputs
-------
results/timetree_divergences.tsv  — pairwise or node-level divergence estimates
results/timetree_gaps.tsv         — temporal gaps >200 Mya with suggested fill species
figures/fig_timetree.pdf/.png     — dated cladogram coloured by phylum

Usage
-----
    python scripts/13_timetree_phylogeny.py

No credentials needed.  The TimeTree API is public but rate-limited; the script
caches responses in results/timetree_cache.json so re-runs are fast.

Dependencies
------------
    pip install requests ete3 matplotlib pandas
    (ete3 requires PyQt5 or PySide2 for display but renders to file without a GUI)
"""

import json
import time
import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from config import PROJECT_DIR, RESULTS_DIR, FIGURES_DIR, PHYLUM_ORDER, PHYLUM_COLOURS

CACHE_FILE   = RESULTS_DIR / "timetree_cache.json"
OUT_DIVERG   = RESULTS_DIR / "timetree_divergences.tsv"
OUT_GAPS     = RESULTS_DIR / "timetree_gaps.tsv"

TIMETREE_API = "https://timetree.org/api/pairwise"
GAP_THRESHOLD_MYA = 200   # flag gaps wider than this
MAX_WORKERS = 10          # Number of parallel threads

# ── Hard-coded consensus divergence times (Mya) for key backbone nodes ────────
BACKBONE_TIMES: dict[tuple[str, str], float] = {
    ("Bacteria",     "Eukaryota"):       2100.0,
    ("Archaea",      "Eukaryota"):       2100.0,
    ("Bacteria",     "Archaea"):         3500.0,
    ("Amoebozoa",    "Opisthokonta"):    1500.0,
    ("Amoebozoa",    "Viridiplantae"):   1500.0,
    ("Fungi",        "Metazoa"):          900.0,
    ("Choanoflagellatea", "Metazoa"):     800.0,
    ("Porifera",     "Eumetazoa"):        700.0,
    ("Ctenophora",   "Bilateria"):        700.0,
    ("Placozoa",     "Bilateria"):        680.0,
    ("Cnidaria",     "Bilateria"):        650.0,
    ("Deuterostomia","Protostomia"):      600.0,
    ("Lophotrochozoa","Ecdysozoa"):       550.0,
    ("Platyhelminthes","Annelida"):       550.0,
    ("Nematoda",     "Arthropoda"):       550.0,
    ("Nematoda",     "Priapulida"):       540.0,
    ("Crustacea",    "Insecta"):          450.0,
    ("Crustacea",    "Chelicerata"):      500.0,
    ("Crustacea",    "Myriapoda"):        450.0,
    ("Chelicerata",  "Insecta"):          490.0,
    ("Myriapoda",    "Insecta"):          440.0,
    ("Echinodermata","Hemichordata"):     500.0,
    ("Echinodermata","Chordata"):          530.0,
    ("Hemichordata", "Chordata"):         520.0,
    ("Echinoidea",   "Asteroidea"):       480.0,
    ("Echinoidea",   "Holothuroidea"):    480.0,
    ("Echinoidea",   "Ophiuroidea"):      470.0,
}

PHYLUM_REPRESENTATIVE: dict[str, str] = {
    # Prokaryotes
    "Bacteria":          "Escherichia coli",
    "Archaea":           "Methanocaldococcus jannaschii",
    # SAR
    "Oomycota":          "Phytophthora infestans",
    "Bacillariophyta":   "Thalassiosira pseudonana",
    "Apicomplexa":       "Plasmodium falciparum",
    "Ciliophora":        "Tetrahymena thermophila",
    "Perkinsozoa":       "Perkinsus marinus",
    "Rhizaria":          "Bigelowiella natans",
    # Excavata
    "Euglenozoa":        "Trypanosoma brucei",
    "Metamonada":        "Giardia intestinalis",
    # Amoebozoa
    "Amoebozoa":         "Dictyostelium discoideum",
    # Archaeplastida
    "Viridiplantae":     "Arabidopsis thaliana",
    "Chlorophyta":       "Chlamydomonas reinhardtii",
    "Rhodophyta":        "Cyanidioschyzon merolae",
    # Other eukaryotes
    "Haptophyta":        "Emiliania huxleyi",
    "Cryptophyta":       "Guillardia theta",
    "Protist":           "Aureococcus anophagefferens",
    # Opisthokonta (non-metazoan)
    "Fungi":             "Saccharomyces cerevisiae",
    # Basal Metazoa
    "Porifera":          "Amphimedon queenslandica",
    "Ctenophora":        "Mnemiopsis leidyi",
    "Placozoa":          "Trichoplax adhaerens",
    "Xenacoelomorpha":   "Hofstenia miamia",
    "Cnidaria":          "Nematostella vectensis",
    # Lophotrochozoa
    "Platyhelminthes":   "Schmidtea mediterranea",
    "Rotifera":          "Adineta vaga",
    "Nemertea":          "Lineus longissimus",
    "Annelida":          "Capitella teleta",
    "Brachiopoda":       "Lingula anatina",
    "Mollusca":          "Lottia gigantea",
    "Acanthocephala":    "Pomphorhynchus laevis",
    # Ecdysozoa
    "Nematomorpha":      "Gordionus sp.",
    "Nematoda":          "Caenorhabditis elegans",
    "Priapulida":        "Priapulus caudatus",
    "Tardigrada":        "Hypsibius exemplaris",
    "Collembola":        "Folsomia candida",
    "Crustacea":         "Daphnia pulex",
    "Myriapoda":         "Strigamia maritima",
    "Chelicerata":       "Ixodes scapularis",
    "Insecta":           "Drosophila melanogaster",
    # Deuterostomia
    "Echinodermata":     "Strongylocentrotus purpuratus",
    "Hemichordata":      "Saccoglossus kowalevskii",
    "Chordata":          "Branchiostoma lanceolatum",
}


GAP_FILL_SUGGESTIONS: dict[str, list[tuple[str, str, str]]] = {
    "Bacteria↔Archaea (3500 Mya — no true fill possible, phylogenetic outgroup)": [],
    "Eukaryota↔Prokaryota (2100 Mya — add deep-branching excavates)": [
        ("giardia_intestinalis_gca000002435v2", "Giardia intestinalis",  "Metamonada"),
        ("naegleria_gruberi_gca000004985v1",    "Naegleria gruberi",     "Heterolobosea"),
    ],
    "Opisthokonta↔Amoebozoa (1500 Mya — add more excavates/SAR)": [
        ("guillardia_theta_gca000315265v1",     "Guillardia theta",      "Cryptista"),
        ("emiliania_huxleyi_gca000372725v1",    "Emiliania huxleyi",     "Haptophyta"),
    ],
    "Fungi↔Choanoflagellatea (900 Mya — choanoflagellates already added)": [],
    "Choanoflagellatea↔Porifera (800 Mya — ichthyosporeans fill ~750 Mya)": [
        ("sphaeroforma_arctica_gca000181245v2", "Sphaeroforma arctica",  "Ichthyosporea"),
        ("capsaspora_owczarzaki_gca000151315v2","Capsaspora owczarzaki", "Filasterea"),
    ],
    "Cnidaria↔Bilateria (650 Mya — myxozoa fill parasitic cnidarian branch)": [
        ("myxobolus_cerebralis_gca002909905v1", "Myxobolus cerebralis",  "Cnidaria"),
    ],
}

# ── TimeTree API ──────────────────────────────────────────────────────────────

def load_cache() -> dict:
    if CACHE_FILE.exists():
        return json.loads(CACHE_FILE.read_text())
    return {}


def save_cache(cache: dict):
    CACHE_FILE.write_text(json.dumps(cache, indent=2))


def query_timetree(sp_a: str, sp_b: str, cache: dict) -> float | None:
    """Query TimeTree API for divergence time between two species (Mya)."""
    key = f"{sp_a}|||{sp_b}"
    key_rev = f"{sp_b}|||{sp_a}"
    if key in cache:
        return cache[key]
    if key_rev in cache:
        return cache[key_rev]

    try:
        payload = {"inputA": sp_a, "inputB": sp_b}
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36",
            "Accept": "application/json"
        }
        
        r = requests.post(TIMETREE_API, json=payload, headers=headers, timeout=15)
        time.sleep(0.1)  # Minimal thread-level delay
        
        if r.status_code == 200:
            try:
                data = r.json()
            except ValueError:
                print(f"    [API non-JSON response] {sp_a} vs {sp_b} -> {r.text[:40]!r}...")
                return None

            mya = None
            if isinstance(data, dict):
                mya = data.get("median") or data.get("mean") or data.get("time")
            elif isinstance(data, list) and data:
                mya = data[0].get("median") or data[0].get("time")
            if mya is not None:
                mya = float(mya)
            
            return mya
        else:
            return None
    except Exception:
        return None


def backbone_time(phylum_a: str, phylum_b: str) -> float | None:
    return BACKBONE_TIMES.get((phylum_a, phylum_b)) or BACKBONE_TIMES.get((phylum_b, phylum_a))

# ── Parallel Processing Wrapper ───────────────────────────────────────────────

def fetch_pair_worker(pa, pb, sp_a, sp_b, cache):
    """Worker task to process a single pair query."""
    mya = None
    if sp_a and sp_b:
        mya = query_timetree(sp_a, sp_b, cache)
    
    if mya is None:
        mya = backbone_time(pa, pb)
        source = "backbone" if mya else "unknown"
    else:
        source = "timetree_api"
        
    return {"phylum_a": pa, "phylum_b": pb, "mya": mya, "source": source, "key": f"{sp_a}|||{sp_b}"}

# ── Node-time matrix (Parallel) ───────────────────────────────────────────────

def build_divergence_table(enr_df: pd.DataFrame, cache: dict) -> pd.DataFrame:
    present_phyla = [p for p in PHYLUM_ORDER if p in enr_df["phylum"].values]
    rows = []
    futures = []

    print(f"Querying TimeTree concurrently using {MAX_WORKERS} threads for {len(present_phyla)} phyla...")
    
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for i, pa in enumerate(present_phyla):
            for pb in present_phyla[i+1:]:
                sp_a = PHYLUM_REPRESENTATIVE.get(pa)
                sp_b = PHYLUM_REPRESENTATIVE.get(pb)
                
                # Check cache synchronously first to minimize thread overhead
                key = f"{sp_a}|||{sp_b}"
                key_rev = f"{sp_b}|||{sp_a}"
                if key in cache or key_rev in cache:
                    mya = cache.get(key) or cache.get(key_rev)
                    rows.append({
                        "phylum_a": pa, "phylum_b": pb, "mya": mya,
                        "source": "timetree_api" if mya else "backbone"
                    })
                    continue
                
                futures.append(executor.submit(fetch_pair_worker, pa, pb, sp_a, sp_b, cache))

        for future in as_completed(futures):
            res = future.result()
            # Update active runtime cache
            if res["source"] == "timetree_api":
                cache[res["key"]] = res["mya"]
                print(f"  {res['phylum_a']} vs {res['phylum_b']}: {res['mya']:.0f} Mya (API)")
            else:
                print(f"  {res['phylum_a']} vs {res['phylum_b']}: {f'{res['mya']:.0f} Mya' if res['mya'] else 'No result'} (Fallback)")
            
            rows.append({
                "phylum_a": res["phylum_a"],
                "phylum_b": res["phylum_b"],
                "mya": res["mya"],
                "source": res["source"]
            })

    return pd.DataFrame(rows)

# ── Species-level crown node times (Parallel) ─────────────────────────────────

def fetch_crown_worker(phylum, species, cache):
    sp_a, sp_b = species[0], species[-1]
    mya = query_timetree(sp_a, sp_b, cache)
    return {"phylum": phylum, "crown_mya": mya, "n_species": len(species), "key": f"{sp_a}|||{sp_b}"}


def species_crown_times(enr_df: pd.DataFrame, cache: dict) -> pd.DataFrame:
    rows = []
    futures = []
    
    print("\nQuerying Within-phylum crown-group ages concurrently...")
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for phylum, group in enr_df.groupby("phylum"):
            species = group["display_name"].tolist()
            if len(species) < 2:
                rows.append({"phylum": phylum, "crown_mya": None, "n_species": 1})
                continue
                
            futures.append(executor.submit(fetch_crown_worker, phylum, species, cache))

        for future in as_completed(futures):
            res = future.result()
            if res["crown_mya"]:
                cache[res["key"]] = res["crown_mya"]
                print(f"  Crown [{res['phylum']}]: {res['crown_mya']:.0f} Mya")
            else:
                print(f"  Crown [{res['phylum']}]: No API result")
                
            rows.append({
                "phylum": res["phylum"],
                "crown_mya": res["crown_mya"],
                "n_species": res["n_species"]
            })

    return pd.DataFrame(rows)

# ── Gap detection ─────────────────────────────────────────────────────────────

def detect_gaps(diverg_df: pd.DataFrame) -> pd.DataFrame:
    present_phyla = [p for p in PHYLUM_ORDER if p in set(diverg_df["phylum_a"]) | set(diverg_df["phylum_b"])]
    gaps = []
    for i in range(len(present_phyla) - 1):
        pa, pb = present_phyla[i], present_phyla[i + 1]
        row = diverg_df[((diverg_df["phylum_a"] == pa) & (diverg_df["phylum_b"] == pb)) |
                        ((diverg_df["phylum_a"] == pb) & (diverg_df["phylum_b"] == pa))]
        if row.empty:
            continue
        mya = row["mya"].values[0]
        if mya and mya >= GAP_THRESHOLD_MYA:
            suggestion_key = next((k for k in GAP_FILL_SUGGESTIONS if pa in k or pb in k), None)
            suggestion = GAP_FILL_SUGGESTIONS.get(suggestion_key, [])
            gaps.append({
                "phylum_a":   pa,
                "phylum_b":   pb,
                "gap_mya":    mya,
                "suggested_fill": "; ".join(f"{n} ({p})" for _, n, p in suggestion) if suggestion else "see notes",
            })
    return pd.DataFrame(gaps)

# ── Figure ────────────────────────────────────────────────────────────────────

def fig_timetree(enr_df: pd.DataFrame, diverg_df: pd.DataFrame, crown_df: pd.DataFrame):
    present_phyla = [p for p in PHYLUM_ORDER if p in enr_df["phylum"].values]
    crown_map = dict(zip(crown_df["phylum"], crown_df["crown_mya"]))

    fig, ax = plt.subplots(figsize=(12, len(present_phyla) * 0.45 + 2))
    y_pos = {p: i for i, p in enumerate(reversed(present_phyla))}

    for phylum in present_phyla:
        y = y_pos[phylum]
        colour = PHYLUM_COLOURS.get(phylum, "#aaaaaa")
        crown = crown_map.get(phylum) or 0
        ax.barh(y, crown, left=0, height=0.6, color=colour, edgecolor="black", linewidth=0.5, align="center")
        ax.text(-20, y, phylum, ha="right", va="center", fontsize=8)

    shown_nodes: set[float] = set()
    for _, row in diverg_df.iterrows():
        pa, pb = row["phylum_a"], row["phylum_b"]
        mya = row["mya"]
        if not mya or pa not in y_pos or pb not in y_pos:
            continue
        ya, yb = y_pos[pa], y_pos[pb]
        node_key = round(mya, 0)
        ax.vlines(mya, min(ya, yb) - 0.1, max(ya, yb) + 0.1, color="black", linewidth=0.8, alpha=0.5)
        if node_key not in shown_nodes:
            ax.text(mya + 5, max(ya, yb) + 0.2, f"{mya:.0f}", fontsize=6, color="grey", va="bottom")
            shown_nodes.add(node_key)

    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel("Divergence time (Mya, right = present)", fontsize=11)
    ax.set_xlim(4000, -200)
    ax.set_yticks([])
    ax.set_title("Dated phylogeny of species in LCR terminal enrichment study\n(bar = within-phylum crown age; node lines = inter-phylum divergence)", fontsize=11)

    sig_phyla = set(enr_df.groupby("phylum").apply(lambda g: g["significant"].any()).pipe(lambda s: s[s].index)) if "significant" in enr_df.columns else set()
    for phylum in sig_phyla:
        if phylum in y_pos:
            ax.text(2, y_pos[phylum], "*", ha="left", va="center", fontsize=10, color="red", fontweight="bold")

    legend_patches = [mpatches.Patch(facecolor=PHYLUM_COLOURS.get(p, "#aaaaaa"), edgecolor="black", label=p) for p in present_phyla if p in PHYLUM_COLOURS]
    ax.legend(handles=legend_patches, fontsize=7, loc="lower left", ncol=3, framealpha=0.7)

    plt.tight_layout()
    out = FIGURES_DIR / "fig_timetree.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nFigure saved: {out}")
    plt.close()

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    enr_path = RESULTS_DIR / "enrichment.tsv"
    if not enr_path.exists():
        print(f"ERROR: {enr_path} not found — run 03_analyse.py first.")
        return

    enr_df = pd.read_csv(enr_path, sep="\t")
    enr_df = enr_df[enr_df["n_lcr"].notna() & (enr_df["n_lcr"] > 0)].copy()

    cache = load_cache()

    diverg_df = build_divergence_table(enr_df, cache)
    save_cache(cache)

    crown_df = species_crown_times(enr_df, cache)
    save_cache(cache)

    print("\n=== Temporal coverage gaps ===")
    gap_df = detect_gaps(diverg_df)

    diverg_df.to_csv(OUT_DIVERG, sep="\t", index=False)
    print(f"Divergence table → {OUT_DIVERG}")

    gap_df.to_csv(OUT_GAPS, sep="\t", index=False)
    print(f"Gap report       → {OUT_GAPS}")

    if not gap_df.empty:
        print(f"\nGaps ≥{GAP_THRESHOLD_MYA} Mya between adjacent sampled phyla:")
        for _, row in gap_df.iterrows():
            print(f"  {row['phylum_a']} ↔ {row['phylum_b']}: {row['gap_mya']:.0f} Mya — fill: {row['suggested_fill']}")
    else:
        print("No gaps above threshold detected.")

    print("\n=== Generating figure ===")
    fig_timetree(enr_df, diverg_df, crown_df)

    n_api    = (diverg_df["source"] == "timetree_api").sum()
    n_back   = (diverg_df["source"] == "backbone").sum()
    print(f"\nDivergence estimates: {n_api} from TimeTree API, {n_back} from backbone literature values.")
    print(f"Species in analysis: {len(enr_df)}, phyla: {enr_df['phylum'].nunique()}")


if __name__ == "__main__":
    main()