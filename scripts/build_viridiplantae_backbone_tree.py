#!/usr/bin/env python3
"""
Build a time-calibrated backbone phylogeny for the Viridiplantae PGLS species set.

This produces results/viridiplantae_timetree.nwk so that 19_pgls_viridiplantae.py
can run a genuine PGLS instead of the OLS-on-tier fallback, WITHOUT requiring a
manual timetree.org upload.

Scope and caveat
----------------
The tree is an ultrametric backbone: the deep nodes (green-algae / streptophyte /
land-plant / angiosperm / monocot–eudicot splits) are fixed at well-established
divergence times (Ma), and the species within each evolutionary tier are placed
as a dated polytomy at that tier's approximate crown age. It is therefore a
coarse, conservative phylogeny — it captures the between-tier structure that the
PGLS tests and correctly stops treating the many eudicot/grass genera as fully
independent, but it does NOT resolve within-tier (family/genus) topology.

For final submission, a species-level dated tree from timetree.org ("Load a List
of Species" on results/pgls_species_list.txt) is preferable; 19_pgls_viridiplantae.py
consumes either tree identically (leaves are matched to species by genus).

Leaf names are the study species_key values (genus-resolvable), so the tree
matches the analysis species directly.

Divergence times (Ma) are rounded consensus values from the plant timescale
(e.g. Morris et al. 2018 PNAS; Kumar et al. 2017 TimeTree; APG IV topology).
The seed-plant crown (gymnosperms | angiosperms) is set to 350 Ma, the midpoint of the
330-370 Ma range in Morris et al. 2018; the euphyllophyte crown (ferns | seed plants) is
set to 400 Ma, between that and the 420-430 Ma lycophyte split.

Tiers 4 and 5 are narrowly sampled: both fern genera are Pteridaceae and all three
gymnosperms are conifers (no Ginkgo, cycads, or Gnetales — those genomes lack NCBI
protein annotation; see 01c_download_viridiplantae_expansion.py). Their within-tier
crown ages are set conservatively low, which increases shared covariance and therefore
down-weights them slightly rather than overstating their independence.

Output:
  results/viridiplantae_timetree.nwk
"""

import csv
from pathlib import Path

from config import RESULTS_DIR

PGLS_TABLE = RESULTS_DIR / "pgls_viridiplantae.tsv"
OUT_TREE   = RESULTS_DIR / "viridiplantae_timetree.nwk"

# Per-tier crown age (Ma) for the within-tier polytomy of extant sampled genera.
# Only used when a tier has >1 sampled species.
TIER_CROWN = {
    0: 750.0,   # Chlorophyta
    1: 400.0,   # Charophytes
    2: 430.0,   # Bryophytes
    3: 350.0,   # Lycophytes
    4: 150.0,   # Ferns — sampled genera are both Pteridaceae (Ceratopteris, Adiantum)
    5: 300.0,   # Gymnosperms — sampled genera are all conifers (Pinaceae | Cupressaceae/Taxaceae)
    6: 135.0,   # Basal angiosperms (ANA grade)
    7: 120.0,   # Eudicots
    8: 110.0,   # Non-grass monocots
    9: 65.0,    # Grasses (Poaceae)
}


def build_tier_clade(species: list[str], crown: float) -> tuple[str, float]:
    """
    Return (newick_substring, attach_age) for one tier.
    Single species  → a bare leaf reaching the present (attach_age 0).
    Multiple species → a polytomy at `crown`; each tip branch = crown.
    """
    if not species:
        return "", 0.0
    if len(species) == 1:
        return species[0], 0.0
    tips = ",".join(f"{s}:{crown:g}" for s in species)
    return f"({tips})", crown


def join(child_a: tuple[str, float], child_b: tuple[str, float],
         node_age: float) -> tuple[str, float]:
    """Join two subtrees at an internal node of the given age (ultrametric)."""
    (sa, aa), (sb, ab) = child_a, child_b
    return f"({sa}:{node_age - aa:g},{sb}:{node_age - ab:g})", node_age


def main():
    if not PGLS_TABLE.exists():
        raise SystemExit(f"ERROR: {PGLS_TABLE} missing — run 19_pgls_viridiplantae.py first.")

    tiers: dict[int, list[str]] = {}
    with open(PGLS_TABLE) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tiers.setdefault(int(row["tier"]), []).append(row["species_key"])

    def clade(t: int):
        return build_tier_clade(sorted(tiers.get(t, [])), TIER_CROWN[t])

    def join_opt(child_a, child_b, node_age: float):
        """Join two subtrees, collapsing the node if either tier has no sampled species."""
        if not child_a[0]:
            return child_b
        if not child_b[0]:
            return child_a
        return join(child_a, child_b, node_age)

    # Assemble the ladder from the youngest split outward, following land-plant
    # phylogeny: ferns sister to seed plants, gymnosperms sister to angiosperms,
    # monocots and grasses nested, ANA grade basal in angiosperms.
    monocots  = join_opt(clade(8), clade(9), 130.0)   # non-grass monocots | grasses
    mesangio  = join_opt(clade(7), monocots, 160.0)   # eudicots | monocots
    angio     = join_opt(clade(6), mesangio, 175.0)   # basal angiosperms | mesangiosperms
    spermato  = join_opt(clade(5), angio,    350.0)   # gymnosperms | angiosperms (seed-plant crown)
    euphyllo  = join_opt(clade(4), spermato, 400.0)   # ferns | seed plants (euphyllophyte crown)
    tracheo   = join_opt(clade(3), euphyllo, 430.0)   # lycophytes | euphyllophytes
    embryo    = join_opt(clade(2), tracheo, 480.0)    # bryophytes | vascular plants
    strepto   = join_opt(clade(1), embryo,  850.0)    # charophytes | land plants
    root      = join_opt(clade(0), strepto, 1100.0)   # chlorophytes | streptophytes

    newick = root[0] + ";\n"
    OUT_TREE.write_text(newick, encoding="utf-8")

    n_species = sum(len(v) for v in tiers.values())
    print(f"Wrote {OUT_TREE}")
    print(f"  {n_species} species across {len(tiers)} tiers; root age 1100 Ma (ultrametric).")
    print("  Tiers present:", ", ".join(f"{t}(n={len(tiers[t])})" for t in sorted(tiers)))
    print("\n  NOTE: coarse time-calibrated backbone (within-tier polytomies). For a "
          "species-level\n  tree, upload results/pgls_species_list.txt to timetree.org "
          "and overwrite this file.")


if __name__ == "__main__":
    main()
