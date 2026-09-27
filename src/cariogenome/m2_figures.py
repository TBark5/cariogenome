"""Figures for M2 (composition). All values plotted are read from the M2 result tables."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .config import coding_genes, genomes
from .plotting import CLASS_COLORS, SEQUENTIAL_CMAP, DIVERGING_CMAP, save
from .seqio import mode_tag

LABELS = {"gc": "GC content", "gc3": "GC at third codon position (GC3)",
          "gc_skew": "GC skew (G-C)/(G+C)", "cai": "Codon Adaptation Index (CAI)"}


def _class_panel(ax, gl: pd.DataFrame, metric: str, test: pd.Series) -> None:
    for i, cls in enumerate(["virulence", "control"]):
        sub = gl[gl["class"] == cls]
        vals = sub[metric].values
        ax.boxplot(vals, positions=[i], widths=0.5, showfliers=False,
                   medianprops={"color": "black"}, boxprops={"color": "#555555"},
                   whiskerprops={"color": "#555555"}, capprops={"color": "#555555"})
        x = i + np.linspace(-0.15, 0.15, len(vals))
        ax.scatter(x, vals, color=CLASS_COLORS[cls], s=30, zorder=3, edgecolor="white", lw=0.5)
        for xi, v, g in zip(x, vals, sub.index):
            ax.annotate(g, (xi, v), xytext=(4, 0), textcoords="offset points", fontsize=6.5,
                        va="center", color="#333333")
    ax.set_xticks([0, 1], [f"virulence\n(n={int((gl['class'] == 'virulence').sum())} genes)",
                           f"housekeeping\n(n={int((gl['class'] == 'control').sum())} genes)"])
    ax.set_xlim(-0.6, 1.8)
    ax.set_ylabel(LABELS[metric])
    ax.set_title(f"{LABELS[metric].split(' (')[0]}\nCliff's d = {test.cliffs_delta:.2f} "
                 f"[{test.delta_ci_low:.2f}, {test.delta_ci_high:.2f}], q = {test.q_bh:.3f}",
                 fontsize=9.5)


def plot_comparison(gl: pd.DataFrame, tests: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6))
    for ax, m in zip(axes, ["gc", "gc3", "cai", "gc_skew"]):
        _class_panel(ax, gl, m, tests.set_index("metric").loc[m])
    n_seq = int(gl["n_sequences"].min()), int(gl["n_sequences"].max())
    fig.suptitle("S. mutans: composition of virulence-associated vs housekeeping genes "
                 f"(each point = gene mean over {n_seq[0]}-{n_seq[1]} strains)" + mode_tag(),
                 fontweight="bold", y=1.04)
    save(fig, "m2_composition_comparison")


def plot_background(bg: pd.DataFrame, per_seq: pd.DataFrame) -> None:
    reps = [g for g in genomes() if g.label in set(bg["label"])]
    fig, axes = plt.subplots(1, len(reps), figsize=(4 * len(reps), 4), sharey=True)
    for ax, g in zip(np.atleast_1d(axes), reps):
        sub = bg[bg["label"] == g.label]
        ax.scatter(sub["gc3"], sub["cai"], s=3, color="#BBBBBB", alpha=0.5, lw=0,
                   label=f"all CDS (n={len(sub)})")
        mine = per_seq[(per_seq["label"] == g.label) & per_seq["gene"].isin(coding_genes())]
        for cls in ("virulence", "control"):
            m = mine[mine["class"] == cls]
            ax.scatter(m["gc3"], m["cai"], s=28, color=CLASS_COLORS[cls], edgecolor="black",
                       lw=0.4, label=f"{cls} (n={len(m)})", zorder=3)
            for r in m.itertuples():
                ax.annotate(r.gene, (r.gc3, r.cai), xytext=(3, 2), textcoords="offset points",
                            fontsize=6.5)
        ax.set_title(f"{g.species} {g.label.split('_', 1)[1]}", style="italic")
        ax.set_xlabel("GC3")
        ax.legend(loc="lower right", fontsize=7)
    np.atleast_1d(axes)[0].set_ylabel("CAI (reference: ribosomal-protein genes)")
    fig.suptitle("Panel genes against the genome-wide codon-usage background" + mode_tag(),
                 fontweight="bold", y=1.03)
    save(fig, "m2_genome_background")


def plot_rscu(rs: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 13))
    im = ax.imshow(np.log2(rs.values.astype(float) + 0.05), cmap=DIVERGING_CMAP, vmin=-2.5,
                   vmax=2.5, aspect="auto")
    ax.set_yticks(range(len(rs)), rs.index, fontsize=6.5, family="monospace")
    ax.set_xticks(range(rs.shape[1]), rs.columns, rotation=60, ha="right")
    for lab in ax.get_xticklabels():
        t = lab.get_text()
        lab.set_color(CLASS_COLORS["virulence"] if t in coding_genes()[:6] else
                      (CLASS_COLORS["control"] if t in coding_genes() else "black"))
    cb = fig.colorbar(im, ax=ax, shrink=0.4)
    cb.set_label("log2 RSCU (0 = no codon preference)")
    ax.set_title("Codon usage (RSCU) in S. mutans: virulence (red), housekeeping (blue),\n"
                 "ribosomal-protein reference set (black); pooled over strains" + mode_tag(),
                 fontsize=10)
    save(fig, "m2_codon_usage_rscu")


def plot_aa(aa: pd.DataFrame, aa_tests: pd.DataFrame) -> None:
    mat = aa.drop(columns="class")
    z = (mat - mat.mean()) / mat.std(ddof=0)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    im = ax.imshow(z.values, cmap=DIVERGING_CMAP, vmin=-2.5, vmax=2.5, aspect="auto")
    q = aa_tests.set_index("metric")["q_bh"]
    ax.set_xticks(range(mat.shape[1]), [f"{a}{'*' if q[a] < 0.05 else ''}" for a in mat.columns])
    ax.set_yticks(range(mat.shape[0]), mat.index)
    for lab in ax.get_yticklabels():
        lab.set_color(CLASS_COLORS[aa.loc[lab.get_text(), "class"]])
    ax.axhline((aa["class"] == "virulence").sum() - 0.5, color="black", lw=1)
    cb = fig.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label("z-score of amino-acid fraction (per column)")
    n_sig = int((q < 0.05).sum())
    ax.set_title("Amino-acid composition in S. mutans (gene means over strains)\n"
                 f"* = virulence vs control differs at BH q < 0.05 ({n_sig} of 20 amino acids)"
                 + mode_tag(), fontsize=10)
    save(fig, "m2_aa_composition")


def plot_species_heatmap(per_seq: pd.DataFrame) -> None:
    sub = per_seq[per_seq["gene"].isin(coding_genes())]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    for ax, m in zip(axes, ["gc3", "cai"]):
        piv = sub.pivot_table(index="gene", columns="species", values=m, aggfunc="mean")
        piv = piv.reindex(coding_genes())[[s for s in dict.fromkeys(g.species for g in genomes()) if s in piv]]
        im = ax.imshow(piv.values, cmap=SEQUENTIAL_CMAP, aspect="auto")
        ax.set_xticks(range(piv.shape[1]), piv.columns, rotation=40, ha="right", style="italic")
        ax.set_yticks(range(piv.shape[0]), piv.index)
        ax.axhline(5.5, color="white", lw=1.5)
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = piv.values[i, j]
                ax.text(j, i, "absent" if np.isnan(v) else f"{v:.2f}", ha="center", va="center",
                        fontsize=6.5, color="black" if np.isnan(v) or v > np.nanmean(piv.values) else "white")
        fig.colorbar(im, ax=ax, shrink=0.7).set_label(LABELS[m])
        ax.set_title(LABELS[m])
    fig.suptitle("Composition by species (mean over strains; rows above the line are "
                 "virulence-associated)" + mode_tag(), fontweight="bold")
    save(fig, "m2_species_composition")


def plot_all(per_seq, gl, tests, aa, aa_tests, rs, bg, pct) -> None:
    """Draw every M2 figure."""
    plot_comparison(gl, tests)
    plot_background(bg, per_seq)
    plot_rscu(rs)
    plot_aa(aa, aa_tests)
    plot_species_heatmap(per_seq)
