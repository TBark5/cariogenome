"""Statistics shared by all modules: group comparisons, effect sizes, bootstrap, FDR.

The unit of replication is the gene. Each comparison reports:
- the difference in medians (virulence minus control) with a bootstrap 95% CI,
- Cliff's delta (a rank effect size in [-1, 1]) with a bootstrap 95% CI,
- a two-sided Mann-Whitney U p-value,
and p-values from a family of tests are adjusted with Benjamini-Hochberg.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def benjamini_hochberg(pvals: list[float] | np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values), NaNs passed through."""
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    m = ok.sum()
    if m == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """P(a > b) - P(a < b) over all pairs; 0 = no difference, +/-1 = complete separation."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    diff = a[:, None] - b[None, :]
    return float((np.sign(diff)).mean())


def bootstrap_ci(
    stat: Callable[[np.ndarray, np.ndarray], float],
    a: np.ndarray,
    b: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 5000,
    level: float = 0.95,
) -> tuple[float, float]:
    """Percentile CI of ``stat(a, b)`` resampling each group independently."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    vals = np.empty(n_boot)
    for i in range(n_boot):
        vals[i] = stat(rng.choice(a, len(a)), rng.choice(b, len(b)))
    lo, hi = np.nanpercentile(vals, [(1 - level) / 2 * 100, (1 + level) / 2 * 100])
    return float(lo), float(hi)


def median_difference(x: np.ndarray, y: np.ndarray) -> float:
    """Median of x minus median of y."""
    return float(np.median(x) - np.median(y))


def compare_groups(
    vir: np.ndarray, ctl: np.ndarray, rng: np.random.Generator, n_boot: int = 5000
) -> dict[str, float]:
    """Full virulence-vs-control comparison for one metric."""
    vir = np.asarray(vir, float)[~np.isnan(np.asarray(vir, float))]
    ctl = np.asarray(ctl, float)[~np.isnan(np.asarray(ctl, float))]
    d_lo, d_hi = bootstrap_ci(median_difference, vir, ctl, rng, n_boot)
    c_lo, c_hi = bootstrap_ci(cliffs_delta, vir, ctl, rng, n_boot)
    p = mannwhitneyu(vir, ctl, alternative="two-sided").pvalue if len(vir) and len(ctl) else np.nan
    return {
        "n_virulence": len(vir),
        "n_control": len(ctl),
        "median_virulence": float(np.median(vir)),
        "median_control": float(np.median(ctl)),
        "median_difference": median_difference(vir, ctl),
        "diff_ci_low": d_lo,
        "diff_ci_high": d_hi,
        "cliffs_delta": cliffs_delta(vir, ctl),
        "delta_ci_low": c_lo,
        "delta_ci_high": c_hi,
        "p_mannwhitney": float(p),
    }


def compare_table(
    df: pd.DataFrame,
    metrics: list[str],
    rng: np.random.Generator,
    class_col: str = "class",
    n_boot: int = 5000,
) -> pd.DataFrame:
    """Run ``compare_groups`` for each metric column and add BH q-values across metrics."""
    rows = []
    for m in metrics:
        res = compare_groups(
            df.loc[df[class_col] == "virulence", m].values,
            df.loc[df[class_col] == "control", m].values,
            rng,
            n_boot,
        )
        rows.append({"metric": m, **res})
    out = pd.DataFrame(rows)
    out["q_bh"] = benjamini_hochberg(out["p_mannwhitney"].values)
    return out


def fmt_ci(est: float, lo: float, hi: float, digits: int = 3) -> str:
    """'est [lo, hi]' string for tables and figure annotations."""
    return f"{est:.{digits}f} [{lo:.{digits}f}, {hi:.{digits}f}]"
