"""Figures for M4: bootstrapped trees, tanglegrams and discordance summaries."""

from __future__ import annotations

import copy
from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.Phylo.BaseTree import Tree
from matplotlib.axes import Axes
from matplotlib.lines import Line2D

from .config import housekeeping_genes, virulence_genes
from .m4_phylogeny import SUPPORTED, load_tree
from .plotting import (
    CLASS_COLORS,
    DARK_GREY,
    FS_ANNOT,
    FS_SMALL,
    FS_TICK,
    GREY,
    LW_THIN,
    W_FULL,
    W_XL,
    annotate_class_n,
    delta_stats,
    save,
    species_color,
    title_with_stats,
)
from .seqio import mode_tag, species_map


def _layout(tree: Tree) -> tuple[dict, dict]:
    """x = distance from the root, y = leaf order (internal nodes at the mean child y)."""
    x: dict = {tree.root: tree.root.branch_length or 0.0}
    for c in tree.find_clades(order="preorder"):
        for ch in c.clades:
            x[ch] = x[c] + (ch.branch_length or 0.0)
    y: dict = {t: float(i) for i, t in enumerate(tree.get_terminals())}
    for c in tree.find_clades(order="postorder"):
        if c.clades:
            y[c] = float(np.mean([y[ch] for ch in c.clades]))
    return x, y


def draw_tree(
    ax: Axes, tree: Tree, title: str, mirror: bool = False, fontsize: float = FS_SMALL
) -> dict:
    """Rectangular tree; support >= 70 in black, < 70 in grey italics ('unsupported')."""
    tree = copy.deepcopy(tree)
    tree.ladderize()  # consistent child order reduces arbitrary crossings in tanglegrams
    sp = species_map()
    x, y = _layout(tree)
    sign = -1 if mirror else 1
    xmax = max(x.values()) or 1.0
    for c in tree.find_clades():
        for ch in c.clades:
            ax.plot([sign * x[c], sign * x[c]], [y[c], y[ch]], color=DARK_GREY, lw=LW_THIN)
            ax.plot([sign * x[c], sign * x[ch]], [y[ch], y[ch]], color=DARK_GREY, lw=LW_THIN)
        if c.is_terminal():
            ax.text(
                sign * (x[c] + xmax * 0.02),
                y[c],
                c.name,
                va="center",
                ha="right" if mirror else "left",
                fontsize=fontsize,
                color=species_color(sp.get(c.name, "")),
            )
        elif c.confidence is not None and c != tree.root:
            ok = c.confidence >= SUPPORTED
            ax.text(
                sign * x[c] - sign * xmax * 0.01,
                y[c] - 0.25,
                f"{int(c.confidence)}",
                ha="left" if mirror else "right",
                va="bottom",
                fontsize=fontsize - 1,
                color="black" if ok else GREY,
                style="normal" if ok else "italic",
            )
    ax.set_ylim(len(tree.get_terminals()) - 0.5, -0.8)
    ax.set_xlim(*((-xmax * 1.6, xmax * 0.05) if mirror else (-xmax * 0.05, xmax * 1.6)))
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("Branch length (substitutions/site, K2P)", fontsize=FS_ANNOT)
    ax.tick_params(axis="x", labelsize=FS_SMALL)
    ax.set_title(title)
    return {t.name: y[t] for t in tree.get_terminals()}


def _support_note() -> str:
    return (
        f"Numbers = bootstrap support (%, 100 replicates); grey italic < {int(SUPPORTED)}% "
        "= unsupported branch. Tip color = species. Midpoint-rooted for display."
    )


def plot_reference() -> None:
    """NJ and UPGMA reference trees from the concatenated housekeeping genes."""
    fig, axes = plt.subplots(1, 2, figsize=(W_FULL, 7.5))
    for ax, m in zip(axes, ("nj", "upgma"), strict=True):
        t = load_tree(f"reference_housekeeping_{m}")
        draw_tree(ax, t, f"{m.upper()} (n={len(t.get_terminals())} genomes)", fontsize=FS_ANNOT)
    fig.suptitle(
        "Reference tree: concatenated housekeeping genes "
        f"({', '.join(housekeeping_genes())})" + mode_tag()
    )
    fig.text(0.5, 0.0, _support_note(), ha="center", fontsize=FS_ANNOT)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    save(fig, "m4_reference_tree")


def plot_gene_trees(genes: list[str], name: str, ncols: int = 3) -> None:
    """Grid of NJ gene trees."""
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(W_XL, 5.6 * nrows))
    for ax, gene in zip(axes.flat, genes, strict=False):  # the grid may have spare panels
        t = load_tree(f"{gene}_nj")
        draw_tree(ax, t, f"{gene} NJ (n={len(t.get_terminals())})")
        cls = "virulence" if gene in virulence_genes() else "control"
        ax.title.set_color(CLASS_COLORS[cls])
    for ax in list(axes.flat)[len(genes) :]:
        ax.axis("off")
    label = "virulence-associated" if name.endswith("virulence") else "housekeeping (and 16S)"
    fig.suptitle(f"Neighbor-joining gene trees: {label} genes" + mode_tag())
    fig.text(0.5, 0.0, _support_note(), ha="center", fontsize=FS_ANNOT)
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    save(fig, name)


def plot_tanglegram(gene: str, ax_pair: Sequence[Axes]) -> None:
    """Reference tree (pruned to the gene's taxa) facing the gene tree."""
    ref = load_tree("reference_housekeeping_nj")
    gt = load_tree(f"{gene}_nj")
    keep = {t.name for t in gt.get_terminals()}
    for t in list(ref.get_terminals()):
        if t.name not in keep:
            ref.prune(t)
    sp = species_map()
    y1 = draw_tree(ax_pair[0], ref, f"reference (housekeeping), n={len(keep)}")
    y2 = draw_tree(ax_pair[1], gt, f"{gene} gene tree, n={len(keep)}", mirror=True)
    fig = ax_pair[0].figure
    to_fig = fig.transFigure.inverted()
    for name in sorted(keep):
        p1 = to_fig.transform(ax_pair[0].transData.transform((ax_pair[0].get_xlim()[1], y1[name])))
        p2 = to_fig.transform(ax_pair[1].transData.transform((ax_pair[1].get_xlim()[0], y2[name])))
        fig.add_artist(
            Line2D(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                transform=fig.transFigure,
                color=species_color(sp[name]),
                lw=LW_THIN,
                alpha=0.8,
            )
        )


def plot_tanglegrams(disc: pd.DataFrame) -> None:
    """Most discordant multi-species virulence gene and control gene vs the reference."""
    multi = disc[disc["n_species"] >= 3]
    vir = multi[multi["class"] == "virulence"].sort_values("nRF_vs_reference").iloc[-1]["gene"]
    ctl = multi[multi["class"] == "control"].sort_values("nRF_vs_reference").iloc[-1]["gene"]
    fig, axes = plt.subplots(1, 4, figsize=(W_XL, 7.5), gridspec_kw={"wspace": 0.35})
    plot_tanglegram(vir, axes[:2])
    plot_tanglegram(ctl, axes[2:])
    fig.suptitle(
        f"Tanglegrams: most discordant multi-species virulence gene ({vir}) and "
        f"control gene ({ctl}) vs the reference tree" + mode_tag()
    )
    fig.text(
        0.5,
        0.0,
        "Crossing lines = taxa placed differently. " + _support_note(),
        ha="center",
        fontsize=FS_ANNOT,
    )
    save(fig, "m4_tanglegrams")


def plot_discordance(disc: pd.DataFrame, tests: pd.DataFrame) -> None:
    """Per-gene RF distance and supported conflicting splits, colored by gene class."""
    coding = disc[disc["class"] != "rRNA control"].reset_index(drop=True)
    t = tests.set_index("metric")
    fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 5.2))
    colors = [CLASS_COLORS[c] for c in coding["class"]]
    x = np.arange(len(coding))
    axes[0].bar(x, coding["nRF_vs_reference"], color=colors)
    axes[0].set_ylabel("Normalized RF distance (0-1)")
    title_with_stats(axes[0], "All taxa of each gene", "taxon sets differ: not tested")
    axes[1].bar(x, coding["nRF_within_Smutans"], color=colors)
    axes[1].set_ylabel("Normalized RF distance (0-1)")
    r = t.loc["nRF_within_Smutans"]
    title_with_stats(axes[1], "Within S. mutans only", delta_stats(r))
    within = coding["n_supported_conflicts"] - coding["n_supported_conflicts_between_species"]
    axes[2].bar(x, within, color=colors, alpha=0.5, label="regroups strains of one species")
    axes[2].bar(
        x,
        coding["n_supported_conflicts_between_species"],
        bottom=within,
        color=colors,
        hatch="///",
        edgecolor="white",
        label="clade mixes species",
    )
    r = t.loc["n_supported_conflicts"]
    title_with_stats(axes[2], f"Supported (>= {int(SUPPORTED)}%) conflicts", delta_stats(r))
    axes[2].set_ylabel("Conflicting splits (count)")
    axes[2].set_ylim(0, coding["n_supported_conflicts"].max() + 2)  # headroom for the legend
    axes[2].legend(loc="upper right")
    for ax in axes:
        ax.set_xticks(
            x,
            [f"{g} (n={n})" for g, n in zip(coding["gene"], coding["n_taxa"], strict=True)],
            rotation=90,
            fontsize=FS_TICK,
        )
        ax.set_xlabel("Gene (n = taxa in the gene tree)")
        for lab, c in zip(ax.get_xticklabels(), coding["class"], strict=True):
            lab.set_color(CLASS_COLORS[c])
    fig.suptitle(
        "Gene-tree discordance with the housekeeping reference tree (suggestive of "
        "transfer or recombination, not proof)" + mode_tag()
    )
    fig.tight_layout()
    annotate_class_n(fig)
    save(fig, "m4_discordance")


def plot_all(disc: pd.DataFrame, tests: pd.DataFrame) -> None:
    """Draw every M4 figure."""
    plot_reference()
    plot_gene_trees(virulence_genes(), "m4_gene_trees_virulence")
    plot_gene_trees([*housekeeping_genes(), "16S"], "m4_gene_trees_controls")
    plot_discordance(disc, tests)
    plot_tanglegrams(disc)
