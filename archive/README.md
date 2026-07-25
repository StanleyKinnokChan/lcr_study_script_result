# archive/ — deprecated deep-time (evolutionary) analysis

Contents here were retired on 2026-07-11 when the manuscript pivoted from an
evolutionary/deep-time framing (v6 and earlier) to a **mechanism-first**
architecture framing (`manuscript/manuscript_v7.md`). The LUCA/LECA origin
inference and Scenario A/B (ancestral vs. convergent) argument were dropped;
nothing downstream in the v7 analysis pipeline consumes these files.

Moved out of the live tree (not deleted — restore by moving back if the
evolutionary angle is ever revived):

| File | Was | Why retired |
|------|-----|-------------|
| `scripts/13_timetree_phylogeny.py` | Dated phylogeny: TimeTree API divergence-time queries, "Divergence time (Mya)" figure, temporal-gap fill-species analysis | Pure deep-time analysis; not referenced by `main.py`; its only figure (former Supp Fig 7) was removed from v7 |
| `results/timetree_divergences.tsv` | Inter-phylum divergence estimates (Mya) | Output of 13; not read by any other script |
| `results/timetree_gaps.tsv` | Temporal sampling gaps >200 Mya + suggested fill taxa | Output of 13; not read by any other script |
| `results/timetree_cache.json` | TimeTree API response cache for 13 | Output of 13 |
| `figures/fig_timetree.pdf` / `.png` | Dated phylogeny figure | Deep-time figure; dropped from v7 |

### Still live (NOT archived) — these belong to the kept PGLS analysis
- `scripts/build_viridiplantae_backbone_tree.py` — builds the conservative dated
  backbone (`results/viridiplantae_timetree.nwk`) that script 19 (PGLS) uses. Now
  wired into `main.py` Phase 5 before step 19.
- `scripts/build_timetree_from_api.py` — optional alternative tree builder
  (timetree.org species-level tree); cited in v7 as a future within-tier
  refinement. Kept as a standalone helper, not in `main.py`.
- `results/timetree_pairwise_cache.json` — cache for `build_timetree_from_api.py`.
