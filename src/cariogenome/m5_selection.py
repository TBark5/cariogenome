"""M5: selection analysis with pooled NG86 dN/dS and codon-bootstrap confidence intervals.

Primary analysis: all pairs of *S. mutans* strains (every gene is present in every strain,
so all genes share the same taxon set). Between-species comparisons are computed but are
flagged as saturated when pS exceeds the configured threshold (Jukes-Cantor becomes
unreliable). Within-species dN/dS measures polymorphism (pN/pS) and is known to be inflated
by slightly deleterious mutations that selection has not yet removed (Rocha et al. 2006),
so absolute values are compared between gene classes, not interpreted in isolation.
"""
from __future__ import annotations

import warnings
from itertools import combinations, product

import numpy as np
import pandas as pd

from .config import RESULTS, coding_genes, gene_class, housekeeping_genes, load_config, rng, virulence_genes
from .dnds import PooledNG86, bootstrap_totals, omega_from_totals
from .seqio import read_fasta, species_map
from .stats import benjamini_hochberg, compare_groups, compare_table

ALN_DIR = RESULTS / "alignments"
FOCAL = "S. mutans"


def codon_msa(gene: str) -> dict[str, str]:
    return read_fasta(ALN_DIR / f"{gene}.codon.fasta")


def _ci(tot_boot: np.ndarray) -> tuple[float, float]:
    om = omega_from_totals(tot_boot)
    om = om[np.isfinite(om)]
    if len(om) < 0.5 * len(tot_boot):
        return np.nan, np.nan
    lo, hi = np.percentile(om, [2.5, 97.5])
    return float(lo), float(hi)


def estimate(msa: dict, pairs: list, n_boot: int, gen: np.random.Generator,
             sat: float) -> tuple[dict, PooledNG86]:
    """Pooled NG86 estimate with bootstrap CI; omega is withheld if pS is saturated."""
    ng = PooledNG86(msa, pairs)
    est = ng.estimate()
    lo, hi = _ci(bootstrap_totals(ng.column_totals(), n_boot, gen))
    saturated = bool(est["pS"] > sat) if np.isfinite(est["pS"]) else False
    out = {**est, "omega_ci_low": lo, "omega_ci_high": hi, "saturated": saturated,
           "n_pairs": len(pairs)}
    if saturated:
        out.update({"omega": np.nan, "omega_ci_low": np.nan, "omega_ci_high": np.nan})
    return out, ng


def within_species(n_boot: int, sat: float) -> tuple[pd.DataFrame, dict]:
    """Pooled within-species omega for every gene and every species with >= 3 strains."""
    sp = species_map()
    rows, models = [], {}
    for i, gene in enumerate(coding_genes()):
        msa = codon_msa(gene)
        by_sp: dict[str, list[str]] = {}
        for k in msa:
            by_sp.setdefault(sp[k], []).append(k)
        for j, (species, labs) in enumerate(by_sp.items()):
            if len(labs) < 3:
                continue
            est, ng = estimate(msa, list(combinations(labs, 2)), n_boot, rng(300 + 10 * i + j), sat)
            rows.append({"gene": gene, "class": gene_class(gene), "species": species,
                         "n_sequences": len(labs), **est})
            models[(gene, species)] = ng
    return pd.DataFrame(rows), models


def between_species(n_boot: int, sat: float) -> pd.DataFrame:
    """S. mutans strains vs each other species (all cross pairs), flagged if saturated."""
    sp = species_map()
    rows = []
    for i, gene in enumerate(coding_genes()):
        msa = codon_msa(gene)
        smu = [k for k in msa if sp[k] == FOCAL]
        for j, other in enumerate(sorted({sp[k] for k in msa} - {FOCAL})):
            oth = [k for k in msa if sp[k] == other]
            est, _ = estimate(msa, list(product(smu, oth)), n_boot, rng(600 + 10 * i + j), sat)
            rows.append({"gene": gene, "class": gene_class(gene), "comparison": f"S. mutans vs {other}",
                         "n_sequences": len(smu) + len(oth), **est})
    return pd.DataFrame(rows)


def gene_vs_controls(models: dict, n_boot: int) -> pd.DataFrame:
    """Each virulence gene's omega relative to pooled S. mutans controls (bootstrap ratio)."""
    ctl_cols = np.hstack([models[(g, FOCAL)].column_totals() for g in housekeeping_genes()
                          if (g, FOCAL) in models])
    ctl_boot = bootstrap_totals(ctl_cols, n_boot, rng(900))
    ctl_om = omega_from_totals(ctl_cols.sum(1))
    rows = []
    for i, gene in enumerate(virulence_genes()):
        if (gene, FOCAL) not in models:
            continue
        cols = models[(gene, FOCAL)].column_totals()
        g_boot = bootstrap_totals(cols, n_boot, rng(910 + i))
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = omega_from_totals(g_boot) / omega_from_totals(ctl_boot)
        ratio = ratio[np.isfinite(ratio)]
        est = omega_from_totals(cols.sum(1)) / ctl_om
        p = 2 * min((ratio <= 1).mean(), (ratio >= 1).mean())
        rows.append({"gene": gene, "omega_ratio_to_controls": float(est),
                     "ratio_ci_low": float(np.percentile(ratio, 2.5)),
                     "ratio_ci_high": float(np.percentile(ratio, 97.5)),
                     "p_bootstrap": float(max(p, 1 / n_boot)), "pooled_control_omega": float(ctl_om)})
    out = pd.DataFrame(rows)
    out["q_bh"] = benjamini_hochberg(out["p_bootstrap"].values)
    return out


def sliding_windows(models: dict, n_boot: int = 200) -> pd.DataFrame:
    """Windowed pN, pS and omega within S. mutans with bootstrap CIs."""
    cfg = load_config()["selection"]
    win, step = cfg["window_codons"], cfg["step_codons"]
    rows = []
    for i, gene in enumerate(coding_genes()):
        cols = models[(gene, FOCAL)].column_totals()
        L = cols.shape[1]
        for start in range(0, max(1, L - win + 1), step):
            sub = cols[:, start:start + win]
            tot = sub.sum(1)
            lo, hi = _ci(bootstrap_totals(sub, n_boot, rng(1000 + i * 1000 + start)))
            rows.append({"gene": gene, "class": gene_class(gene), "start_codon": start + 1,
                         "end_codon": min(L, start + win), "mid_codon": start + win / 2,
                         "pN": tot[3] / tot[1], "pS": tot[2] / tot[0],
                         "omega": float(omega_from_totals(tot)), "omega_ci_low": lo, "omega_ci_high": hi})
    return pd.DataFrame(rows)


def replication(within: pd.DataFrame, n_boot: int) -> pd.DataFrame:
    """Virulence-homolog vs control omega within each species (does the pattern replicate?)."""
    rows = []
    for i, (species, sub) in enumerate(within.groupby("species", sort=False)):
        sub = sub.dropna(subset=["omega"])
        vir, ctl = sub[sub["class"] == "virulence"], sub[sub["class"] == "control"]
        if len(vir) < 2 or len(ctl) < 2:
            continue
        res = compare_groups(vir["omega"].values, ctl["omega"].values, rng(950 + i), n_boot)
        rows.append({"species": species, "virulence_genes_present": ";".join(vir["gene"]), **res})
    return pd.DataFrame(rows)


def run() -> dict[str, pd.DataFrame]:
    """Run all selection analyses, tests and figures."""
    from . import m5_figures
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    cfg = load_config()
    n_boot, sat = cfg["bootstrap"]["dnds_replicates"], cfg["selection"]["saturation_ps"]
    within, models = within_species(n_boot, sat)
    within.to_csv(RESULTS / "m5_dnds_within_species.csv", index=False)
    between = between_species(200, sat)
    between.to_csv(RESULTS / "m5_dnds_between_species.csv", index=False)
    smu = within[within["species"] == FOCAL].set_index("gene").reindex(coding_genes()).reset_index()
    smu.to_csv(RESULTS / "m5_dnds_Smutans.csv", index=False)
    tests = compare_table(smu, ["omega", "dN", "dS"], rng(7), n_boot=cfg["bootstrap"]["stats_replicates"])
    tests.to_csv(RESULTS / "m5_tests.csv", index=False)
    rep = replication(within, cfg["bootstrap"]["stats_replicates"])
    rep.to_csv(RESULTS / "m5_replication_by_species.csv", index=False)
    per_gene = gene_vs_controls(models, n_boot)
    per_gene.to_csv(RESULTS / "m5_gene_vs_controls.csv", index=False)
    windows = sliding_windows(models)
    windows.to_csv(RESULTS / "m5_sliding_windows.csv", index=False)
    m5_figures.plot_all(smu, within, between, tests, per_gene, windows)
    print(smu[["gene", "n_sequences", "pN", "pS", "omega", "omega_ci_low", "omega_ci_high"]].round(4).to_string(index=False))
    print(tests[["metric", "median_virulence", "median_control", "cliffs_delta", "delta_ci_low",
                 "delta_ci_high", "q_bh"]].round(4).to_string(index=False))
    print(per_gene.round(4).to_string(index=False))
    print(rep[["species", "virulence_genes_present", "median_virulence", "median_control",
               "cliffs_delta", "delta_ci_low", "delta_ci_high", "p_mannwhitney"]].round(3).to_string(index=False))
    print("between-species saturated:", int(between["saturated"].sum()), "of", len(between))
    return {"smutans": smu, "tests": tests, "per_gene": per_gene, "between": between}
