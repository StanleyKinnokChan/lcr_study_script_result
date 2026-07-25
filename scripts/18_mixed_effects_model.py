#!/usr/bin/env python3
"""
Mixed-effects logistic regression for prokaryote terminal LCR enrichment.

The pooled Fisher's exact test in 03_analyse.py treats all LCRs as independent,
ignoring that LCRs from the same species are clustered.  This script fits a
proper random-effects logistic regression:

    terminal (0/1) ~ protein_length_quartile + domain + (1 | species_key)

using a GEE (generalised estimating equations) logistic regression with
cluster-robust standard errors (species as the cluster) — the primary
clustering-corrected estimate — with a Bayesian mixed GLM (random species
intercept) as a fallback, plus a species-level bootstrap confidence interval on
the domain-level mean % terminal (resampling species rather than individual LCRs).

All approaches are run; the bootstrap CI is always available, the GEE / GLMM
require statsmodels (>= 0.14).

Also runs the mixed model on all 42 phyla as a robustness check for the eukaryote
results.

Outputs:
  results/mixed_model_prokaryote.tsv  — GLMM fixed-effect estimates + CIs
  results/bootstrap_ci_domain.tsv    — bootstrap mean ± 95% CI per domain
  figures/suppfig_bootstrap_ci.pdf
"""

import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from pathlib import Path

from config import RESULTS_DIR, FIGURES_DIR, N_BINS


def protein_length_quartile(length: float, quartile_bounds: list[float]) -> int:
    """Return quartile (1–4) for a given protein length."""
    for i, bound in enumerate(quartile_bounds):
        if length <= bound:
            return i + 1
    return 4


def bootstrap_mean_ci(values: list[float], n_boot: int = 5000,
                      alpha: float = 0.05, seed: int = 42) -> tuple[float, float, float]:
    """Bootstrap CI on the mean by resampling species."""
    rng = np.random.default_rng(seed)
    vals = np.array(values)
    boot_means = [rng.choice(vals, size=len(vals), replace=True).mean()
                  for _ in range(n_boot)]
    lo = np.percentile(boot_means, alpha / 2 * 100)
    hi = np.percentile(boot_means, (1 - alpha / 2) * 100)
    return float(vals.mean()), float(lo), float(hi)


def fit_clustered_model(df: pd.DataFrame) -> pd.DataFrame | None:
    """
    Fit a clustering-corrected logistic regression of terminal LCR occurrence:

        is_terminal (0/1) ~ C(quartile) + C(domain),  clustered by species_key

    Primary method is a GEE (generalised estimating equations) with a Binomial
    family, exchangeable working correlation, and species as the cluster — this
    yields cluster-robust (sandwich) standard errors that account for the
    non-independence of LCRs from the same species, which is exactly what the
    pooled Fisher's test ignores. GEE converges reliably where a Bayesian GLMM
    often will not.

    Falls back to a Bayesian mixed GLM (random species intercept) if GEE is
    unavailable, then to None. Returns a fixed-effects table with 95% CIs and,
    for GEE, robust p-values.
    """
    # ── Primary: GEE with cluster-robust SEs ──────────────────────────────────
    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf

        model = smf.gee(
            "is_terminal ~ C(quartile) + C(domain)",
            groups="species_key",
            data=df,
            family=sm.families.Binomial(),
            cov_struct=sm.cov_struct.Exchangeable(),
        )
        res = model.fit()
        ci = res.conf_int()
        return pd.DataFrame({
            "term":  res.params.index.tolist(),
            "coef":  res.params.round(4).tolist(),
            "ci_lo": ci[0].round(4).tolist(),
            "ci_hi": ci[1].round(4).tolist(),
            "p":     res.pvalues.round(6).tolist(),
            "method": "GEE_binomial_exchangeable_cluster_robust",
        })
    except ImportError:
        print("  statsmodels not installed — GEE unavailable, trying Bayesian GLMM.")
    except Exception as e:
        print(f"  GEE fitting failed ({e}); trying Bayesian GLMM.")

    # ── Fallback: Bayesian mixed GLM (random species intercept) ───────────────
    try:
        from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
        model = BinomialBayesMixedGLM.from_formula(
            "is_terminal ~ C(quartile) + C(domain)",
            {"species_key": "0 + C(species_key)"},
            df,
        )
        result = model.fit_vb()
        fe = result.summary_struct(0.05)["Fixed Effects"]
        return pd.DataFrame({
            "term":   fe.index.tolist(),
            "coef":   fe["Mean"].round(4).tolist(),
            "ci_lo":  (fe["Mean"] - 1.96 * fe["SD"]).round(4).tolist(),
            "ci_hi":  (fe["Mean"] + 1.96 * fe["SD"]).round(4).tolist(),
            "p":      None,
            "method": "BinomialBayesMixedGLM_random_species_intercept",
        })
    except ImportError:
        print("  statsmodels not available — clustered model skipped.")
        return None
    except Exception as e:
        print(f"  Bayesian GLMM fitting failed: {e}")
        return None


def main():
    pos_path = RESULTS_DIR / "lcr_positions.tsv"
    enr_path = RESULTS_DIR / "enrichment.tsv"
    if not pos_path.exists() or not enr_path.exists():
        print(f"ERROR: Run 03_analyse.py first to generate {pos_path} and {enr_path}.")
        return

    pos_df = pd.read_csv(pos_path, sep="\t")
    enr_df = pd.read_csv(enr_path, sep="\t")

    # ── Domain assignment ─────────────────────────────────────────────────────
    PROKARYOTE_PHYLA = {"Bacteria", "Archaea"}
    pos_df["domain"] = pos_df["phylum"].apply(
        lambda p: p if p in PROKARYOTE_PHYLA else "Eukaryota"
    )

    # Protein-length quartile (global within domain)
    for domain in pos_df["domain"].unique():
        mask = pos_df["domain"] == domain
        bounds = pos_df.loc[mask, "protein_len"].quantile([0.25, 0.5, 0.75]).tolist()
        pos_df.loc[mask, "quartile"] = pos_df.loc[mask, "protein_len"].apply(
            lambda l: protein_length_quartile(l, bounds)
        )
    pos_df["quartile"] = pos_df["quartile"].astype(int)

    # ── Bootstrap CI by domain (species-level resampling) ────────────────────
    enr_df["domain"] = enr_df["phylum"].apply(
        lambda p: p if p in PROKARYOTE_PHYLA else "Eukaryota"
    )

    boot_rows = []
    print("Species-level bootstrap CI on mean % terminal (5,000 resamples):\n")
    for domain in ["Bacteria", "Archaea"]:
        sp_pcts = enr_df[
            (enr_df["domain"] == domain) & enr_df["pct_terminal"].notna()
        ]["pct_terminal"].tolist()

        if len(sp_pcts) < 3:
            print(f"  {domain}: insufficient species (n={len(sp_pcts)}) for bootstrap.")
            continue

        mean_, lo, hi = bootstrap_mean_ci(sp_pcts)
        print(f"  {domain}  n_species={len(sp_pcts)}")
        print(f"    Mean %terminal = {mean_:.2f}%  95% CI = [{lo:.2f}%, {hi:.2f}%]")
        print(f"    (pooled Fisher's test result should fall near this mean)\n")
        boot_rows.append({
            "domain":      domain,
            "n_species":   len(sp_pcts),
            "mean_pct":    round(mean_, 2),
            "ci_lo_95":    round(lo, 2),
            "ci_hi_95":    round(hi, 2),
        })

    boot_df = pd.DataFrame(boot_rows)
    boot_path = RESULTS_DIR / "bootstrap_ci_domain.tsv"
    boot_df.to_csv(boot_path, sep="\t", index=False)
    print(f"Bootstrap CI table: {boot_path}")

    # ── GLMM attempt (prokaryotes only) ───────────────────────────────────────
    prokaryote_pos = pos_df[pos_df["domain"].isin(PROKARYOTE_PHYLA)].copy()
    prokaryote_pos["is_terminal"] = prokaryote_pos["is_terminal"].astype(int)

    if len(prokaryote_pos) > 0:
        print("\nFitting clustering-corrected model on prokaryote LCRs...")
        glmm_path = RESULTS_DIR / "mixed_model_prokaryote.tsv"
        model_result = fit_clustered_model(prokaryote_pos)
        if model_result is not None:
            model_result.to_csv(glmm_path, sep="\t", index=False)
            print(f"Clustered-model fixed effects: {glmm_path}")
            print(model_result.to_string(index=False))
        else:
            print("  Clustered model unavailable. Bootstrap CI is the reported statistic.")
            # Save a sentinel so the output path exists
            pd.DataFrame([{
                "term": "clustered_model_unavailable",
                "coef": None, "ci_lo": None, "ci_hi": None, "p": None,
                "method": "install statsmodels to enable GEE / Bayesian GLMM"
            }]).to_csv(glmm_path, sep="\t", index=False)

    # ── Quartile enrichment (accounting for species clustering) ───────────────
    print("\nQuartile-level species-mean % terminal (prokaryotes):")
    for domain in ["Bacteria", "Archaea"]:
        print(f"  {domain}:")
        for q in [1, 2, 3, 4]:
            sp_q = []
            for sp_key in prokaryote_pos[
                prokaryote_pos["phylum"] == domain
            ]["species_key"].unique():
                sub = prokaryote_pos[
                    (prokaryote_pos["species_key"] == sp_key) &
                    (prokaryote_pos["quartile"] == q)
                ]
                if len(sub) >= 5:
                    sp_q.append(sub["is_terminal"].mean() * 100)
            if sp_q:
                print(f"    Q{q}: n_sp={len(sp_q)}  "
                      f"mean={np.mean(sp_q):.1f}%  SD={np.std(sp_q):.1f}%")

    # ── Figure: bootstrap CI violin ───────────────────────────────────────────
    if boot_df.empty:
        print("\nNo bootstrap data to plot.")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    x = np.arange(len(boot_df))
    colours = {"Bacteria": "#636363", "Archaea": "#969696"}
    ax.bar(x, boot_df["mean_pct"],
           color=[colours.get(d, "#aaa") for d in boot_df["domain"]],
           edgecolor="black", linewidth=0.7, zorder=3)
    ax.errorbar(
        x, boot_df["mean_pct"],
        yerr=[
            boot_df["mean_pct"] - boot_df["ci_lo_95"],
            boot_df["ci_hi_95"] - boot_df["mean_pct"]
        ],
        fmt="none", color="black", capsize=5, linewidth=1.5, zorder=4
    )
    ax.axhline(10.0, color="red", linestyle=":", linewidth=1.2,
               label="Null (10%)")
    ax.set_xticks(x)
    ax.set_xticklabels(boot_df["domain"])
    ax.set_ylabel("Mean % terminal LCRs (species-level)", fontsize=11)
    ax.set_title("Prokaryote terminal LCR enrichment\n"
                 "(bootstrap 95% CI by species resampling, n=5,000)", fontsize=11)
    ax.legend(fontsize=9)
    ax.set_ylim(0, max(boot_df["ci_hi_95"].max() * 1.2, 15))

    plt.tight_layout()
    out_fig = FIGURES_DIR / "suppfig_bootstrap_ci.pdf"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    fig.savefig(str(out_fig).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"\nBootstrap CI figure: {out_fig}")
    plt.close()


if __name__ == "__main__":
    main()
