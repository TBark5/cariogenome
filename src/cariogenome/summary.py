"""Cross-module synthesis: one table and one forest plot of every virulence-vs-control test."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .config import RESULTS
from .plotting import (FS_ANNOT, LW_THICK, LW_THIN, NOT_SIGNIFICANT, SIGNIFICANT, W_PAIR,
                       annotate_class_n, save)
from .seqio import mode_tag

TESTS = [
    ("M2 composition", "m2_tests.csv", {"gc": "GC content", "gc3": "GC3", "cai": "CAI",
                                        "gc_skew": "GC skew"}),
    ("M3 conservation", "m3_tests.csv", {"mean_pid_Smutans": "mean aa identity",
                                         "mean_entropy_Smutans": "mean per-site entropy"}),
    ("M4 phylogeny", "m4_tests.csv", {"nRF_within_Smutans": "RF distance to reference",
                                      "n_supported_conflicts": "supported conflicting splits"}),
    ("M5 selection", "m5_tests.csv", {"omega": "dN/dS (omega)", "dN": "dN", "dS": "dS"}),
]


def effect_table() -> pd.DataFrame:
    """All class comparisons in one table (Cliff's delta > 0 = higher in virulence genes)."""
    rows = []
    for module, fname, labels in TESTS:
        t = pd.read_csv(RESULTS / fname).set_index("metric")
        for key, label in labels.items():
            r = t.loc[key]
            rows.append({"module": module, "metric": key, "label": label,
                         "n_virulence": int(r.n_virulence), "n_control": int(r.n_control),
                         "median_virulence": r.median_virulence, "median_control": r.median_control,
                         "cliffs_delta": r.cliffs_delta, "ci_low": r.delta_ci_low,
                         "ci_high": r.delta_ci_high, "p": r.p_mannwhitney, "q_bh": r.q_bh})
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "summary_effect_sizes.csv", index=False)
    return out


def plot_forest(tab: pd.DataFrame) -> None:
    """Every class comparison on one Cliff's delta axis."""
    fig, ax = plt.subplots(figsize=(W_PAIR, 6.8))
    y = np.arange(len(tab))[::-1]
    for yi, r in zip(y, tab.itertuples(), strict=True):
        sig = r.q_bh < 0.05
        col = SIGNIFICANT if sig else NOT_SIGNIFICANT
        ax.plot([r.ci_low, r.ci_high], [yi, yi], color=col, lw=LW_THICK)
        ax.scatter(r.cliffs_delta, yi, color=col, s=45, zorder=3, marker="D" if sig else "o")
        ax.text(1.08, yi, f"q = {r.q_bh:.3f}" if np.isfinite(r.q_bh) else "q = n/a", va="center",
                fontsize=FS_ANNOT, color=col)
    ax.axvline(0, color="black", lw=LW_THIN)
    ax.set_yticks(y, [f"{r.module}: {r.label}" for r in tab.itertuples()])
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Cliff's delta, virulence-associated vs housekeeping (95% bootstrap CI)")
    ax.text(-1.0, len(tab) - 0.2, "lower in virulence-associated genes", fontsize=FS_ANNOT)
    ax.text(1.0, len(tab) - 0.2, "higher in virulence-associated genes", fontsize=FS_ANNOT,
            ha="right")
    ax.set_ylim(-0.7, len(tab) + 0.2)
    ax.scatter([], [], marker="D", color=SIGNIFICANT, label="BH q < 0.05 (within its module)")
    ax.scatter([], [], marker="o", color=NOT_SIGNIFICANT, label="not significant")
    ax.legend(loc="lower left", bbox_to_anchor=(0, -0.2), ncol=2)
    ax.set_title("All virulence-vs-housekeeping comparisons in S. mutans" + mode_tag())
    annotate_class_n(fig, y=-0.06)
    save(fig, "summary_effect_sizes")


def run() -> pd.DataFrame:
    """Build the effect-size table and the forest plot."""
    tab = effect_table()
    plot_forest(tab)
    cols = ["module", "label", "cliffs_delta", "ci_low", "ci_high", "q_bh"]
    print(tab[cols].round(3).to_string(index=False))
    return tab
