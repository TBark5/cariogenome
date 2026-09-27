"""Figures for M2 (composition). All values plotted are read from the M2 result tables."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import coding_genes, genomes, virulence_genes
from .plotting import (CLASS_COLORS, DARK_GREY, DIVERGING_CMAP, FS_SMALL, LW_THIN,
                       SEQUENTIAL_CMAP, W_FULL, annotate_class_n, delta_stats, save,
                       title_with_stats)
from .seqio import mode_tag

LABELS = {"gc": "GC content (fraction)", "gc3": "GC at third codon position (fraction)",
          "gc_skew": "GC skew, (G-C)/(G+C)", "cai": "Codon Adaptation Index (0-1)"}
TITLES = {"gc": "GC content", "gc3": "GC3", "gc_skew": "GC skew", "cai": "Codon Adaptation Index"}


def _color_gene_labels(labels) -> None:
    for lab in labels:
        g = lab.get_text()
        if g in coding_genes():
            lab.set_color(CLASS_COLORS["virulence" if g in virulence_genes() else "control"])


def _class_panel(ax, gl: pd.DataFrame, metric: str, test: pd.Series) -> None:
    for i, cls in enumerate(["virulence", "control"]):
        sub = gl[gl["class"] == cls]
        vals = sub[metric].to_numpy()
        ax.boxplot(vals, positions=[i], widths=0.5, showfliers=False,
                   medianprops={"color": "black"}, boxprops={"color": DARK_GREY},
                   whiskerprops={"color": DARK_GREY}, capprops={"color": DARK_GREY})
        x = i + np.linspace(-0.15, 0.15, len(vals))
        ax.scatter(x, vals, color=CLASS_COLORS[cls], s=30, zorder=3, edgecolor="white", lw=0.5)
        for xi, v, g in zip(x, vals, sub.index, strict=True):
            ax.annotate(g, (xi, v), xytext=(4, 0), textcoords="offset points", fontsize=FS_SMALL,
                        va="center", color=DARK_GREY)
    n_v, n_c = int((gl["class"] == "virulence").sum()), int((gl["class"] == "control").sum())
    ax.set_xticks([0, 1], [f"virulence\n(n={n_v} genes)", f"housekeeping\n(n={n_c} genes)"])
    ax.set_xlim(-0.6, 1.8)
    ax.set_xlabel("Gene class")
    ax.set_ylabel(LABELS[metric])
    title_with_stats(ax, TITLES[metric], delta_stats(test))


def plot_comparison(gl: pd.DataFrame, tests: pd.DataFrame) -> None:
    """Four composition metrics, virulence vs housekeeping genes."""
    fig, axes = plt.subplots(1, 4, figsize=(W_FULL, 5.0))
    for ax, m in zip(axes, ["gc", "gc3", "cai", "gc_skew"], strict=True):
        _class_panel(ax, gl, m, tests.set_index("metric").loc[m])
    lo, hi = int(gl["n_sequences"].min()), int(gl["n_sequences"].max())
    fig.suptitle("S. mutans: composition of virulence-associated vs housekeeping genes "
                 f"(each point = gene mean over {lo}-{hi} strains)" + mode_tag())
    fig.tight_layout()
    annotate_class_n(fig)
    save(fig, "m2_composition_comparison")


def plot_background(bg: pd.DataFrame, per_seq: pd.DataFrame) -> None:
    """GC3 vs CAI for every CDS of each species representative, panel genes overlaid."""
    reps = [g for g in genomes() if g.label in set(bg["label"])]
    fig, axes = plt.subplots(1, len(reps), figsize=(W_FULL, 4.2), sharey=True)
    for ax, g in zip(np.atleast_1d(axes), reps, strict=True):
        sub = bg[bg["label"] == g.label]
        ax.scatter(sub["gc3"], sub["cai"], s=3, color="#BBBBBB", alpha=0.5, lw=0,
                   label="all CDS of the genome")
        mine = per_seq[(per_seq["label"] == g.label) & per_seq["gene"].isin(coding_genes())]
        n = {}
        for cls in ("virulence", "control"):
            m = mine[mine["class"] == cls]
            n[cls] = len(m)
            name = "virulence-associated gene" if cls == "virulence" else "housekeeping gene"
            ax.scatter(m["gc3"], m["cai"], s=22, color=CLASS_COLORS[cls], edgecolor="black",
                       lw=0.4, label=name, zorder=3)
        title_with_stats(ax, f"{g.species} {g.label.split('_', 1)[1]}",
                         f"all CDS n={len(sub)}\n"
                         f"virulence n={n['virulence']}, housekeeping n={n['control']}")
        ax.title.set_style("italic")
        ax.set_xlabel("GC3 (fraction)")
    first = np.atleast_1d(axes)[0]
    first.set_ylabel("CAI (reference: ribosomal proteins)")
    fig.legend(*first.get_legend_handles_labels(), loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.04), markerscale=1.5)
    fig.suptitle("Panel genes against the genome-wide codon-usage background" + mode_tag())
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "m2_genome_background")


def plot_rscu(rs: pd.DataFrame) -> None:
    """Heatmap of log2 RSCU per gene, plus the ribosomal reference set."""
    fig, ax = plt.subplots(figsize=(W_FULL / 2, 13))
    im = ax.imshow(np.log2(rs.to_numpy(dtype=float) + 0.05), cmap=DIVERGING_CMAP, vmin=-2.5,
                   vmax=2.5, aspect="auto")
    ax.set_yticks(range(len(rs)), rs.index, fontsize=FS_SMALL, family="monospace")
    ax.set_xticks(range(rs.shape[1]), rs.columns, rotation=60, ha="right")
    _color_gene_labels(ax.get_xticklabels())
    ax.set_xlabel("Gene (vermillion = virulence-associated, blue = housekeeping)")
    ax.set_ylabel("Amino acid - codon")
    fig.colorbar(im, ax=ax, shrink=0.4).set_label("log2 RSCU (0 = no preference)")
    ax.set_title("Codon usage (RSCU) in S. mutans, pooled over strains\n"
                 "(white = amino acid absent from that gene)" + mode_tag())
    annotate_class_n(fig)
    save(fig, "m2_codon_usage_rscu")


def plot_aa(aa: pd.DataFrame, aa_tests: pd.DataFrame) -> None:
    """Heatmap of amino-acid composition (z-scored per amino acid)."""
    mat = aa.drop(columns="class")
    z = (mat - mat.mean()) / mat.std(ddof=0)
    fig, ax = plt.subplots(figsize=(W_FULL * 0.8, 6.0))
    im = ax.imshow(z.to_numpy(), cmap=DIVERGING_CMAP, vmin=-2.5, vmax=2.5, aspect="auto")
    q = aa_tests.set_index("metric")["q_bh"]
    ax.set_xticks(range(mat.shape[1]), [f"{a}{'*' if q[a] < 0.05 else ''}" for a in mat.columns])
    ax.set_yticks(range(mat.shape[0]), mat.index)
    _color_gene_labels(ax.get_yticklabels())
    ax.axhline((aa["class"] == "virulence").sum() - 0.5, color="black", lw=LW_THIN)
    ax.set_xlabel("Amino acid (* = BH q < 0.05, virulence vs housekeeping)")
    ax.set_ylabel("Gene")
    fig.colorbar(im, ax=ax, shrink=0.7).set_label("z-score of amino-acid fraction")
    n_sig = int((q < 0.05).sum())
    ax.set_title("Amino-acid composition in S. mutans (gene means over strains): "
                 f"{n_sig} of 20 amino acids differ" + mode_tag())
    annotate_class_n(fig)
    save(fig, "m2_aa_composition")


def plot_species_heatmap(per_seq: pd.DataFrame) -> None:
    """GC3 and CAI per gene and species."""
    sub = per_seq[per_seq["gene"].isin(coding_genes())]
    fig, axes = plt.subplots(1, 2, figsize=(W_FULL, 6.0))
    order = list(dict.fromkeys(g.species for g in genomes()))
    for ax, m in zip(axes, ["gc3", "cai"], strict=True):
        piv = sub.pivot_table(index="gene", columns="species", values=m, aggfunc="mean")
        piv = piv.reindex(coding_genes())[[s for s in order if s in piv]]
        vals = piv.to_numpy()
        im = ax.imshow(vals, cmap=SEQUENTIAL_CMAP, aspect="auto")
        ax.set_xticks(range(piv.shape[1]), piv.columns, rotation=40, ha="right", style="italic")
        ax.set_yticks(range(piv.shape[0]), piv.index)
        _color_gene_labels(ax.get_yticklabels())
        ax.axhline(len(virulence_genes()) - 0.5, color="white", lw=1.5)
        mid = np.nanmean(vals)
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = vals[i, j]
                ax.text(j, i, "absent" if np.isnan(v) else f"{v:.2f}", ha="center", va="center",
                        fontsize=FS_SMALL, color="black" if np.isnan(v) or v > mid else "white")
        ax.set_xlabel("Species")
        ax.set_ylabel("Gene")
        fig.colorbar(im, ax=ax, shrink=0.7).set_label(LABELS[m])
        ax.set_title(TITLES[m])
    fig.suptitle("Composition by species (mean over strains; genes above the white line are "
                 "virulence-associated)" + mode_tag())
    fig.tight_layout()
    annotate_class_n(fig)
    save(fig, "m2_species_composition")


def plot_all(per_seq, gl, tests, aa, aa_tests, rs, bg, pct) -> None:
    """Draw every M2 figure."""
    plot_comparison(gl, tests)
    plot_background(bg, per_seq)
    plot_rscu(rs)
    plot_aa(aa, aa_tests)
    plot_species_heatmap(per_seq)

