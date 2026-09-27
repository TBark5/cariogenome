"""Cross-module synthesis: one table and one forest plot of every virulence-vs-control test."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .config import RESULTS
from .plotting import CLASS_COLORS, OKABE_ITO, save
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
    fig, ax = plt.subplots(figsize=(10, 6.5))
    y = np.arange(len(tab))[::-1]
    for yi, r in zip(y, tab.itertuples()):
        sig = r.q_bh < 0.05
        col = OKABE_ITO["vermillion"] if sig else "#777777"
        ax.plot([r.ci_low, r.ci_high], [yi, yi], color=col, lw=2)
        ax.scatter(r.cliffs_delta, yi, color=col, s=45, zorder=3, marker="D" if sig else "o")
        ax.text(1.08, yi, f"q = {r.q_bh:.3f}" if np.isfinite(r.q_bh) else "q = n/a", va="center",
                fontsize=8.5, color=col)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y, [f"{r.module}: {r.label}" for r in tab.itertuples()], fontsize=9)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Cliff's delta (virulence vs housekeeping genes), 95% bootstrap CI")
    ax.text(-1.0, len(tab) - 0.2, "lower in virulence genes", fontsize=8, color=CLASS_COLORS["control"])
    ax.text(1.0, len(tab) - 0.2, "higher in virulence genes", fontsize=8, ha="right",
            color=CLASS_COLORS["virulence"])
    ax.set_ylim(-0.7, len(tab) + 0.2)
    n = f"{int(tab['n_virulence'].iloc[0])} virulence vs {int(tab['n_control'].iloc[0])} control genes"
    ax.set_title(f"All virulence-vs-control comparisons in S. mutans ({n})\n"
                 "red diamond = BH q < 0.05 within its module; grey = not significant" + mode_tag(),
                 fontsize=10.5)
    save(fig, "summary_effect_sizes")


def run() -> pd.DataFrame:
    tab = effect_table()
    plot_forest(tab)
    print(tab[["module", "label", "cliffs_delta", "ci_low", "ci_high", "q_bh"]].round(3).to_string(index=False))
    return tab
