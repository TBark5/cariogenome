"""Figures for M4: bootstrapped trees, a tanglegram and discordance summaries."""
from __future__ import annotations

import copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .config import housekeeping_genes, virulence_genes
from .m4_phylogeny import SUPPORTED, load_tree
from .plotting import CLASS_COLORS, OKABE_ITO, SPECIES_COLORS, save
from .seqio import mode_tag, species_map


def _layout(tree) -> tuple[dict, dict]:
    """x = distance from root, y = leaf order (internal nodes at mean child y)."""
    x, y = {}, {}
    leaves = tree.get_terminals()
    for i, t in enumerate(leaves):
        y[t] = i
    for c in tree.find_clades(order="preorder"):
        parent_x = x.get(getattr(c, "_parent", None), 0.0)
        x[c] = parent_x + (c.branch_length or 0.0)
        for ch in c.clades:
            ch._parent = c
    for c in tree.find_clades(order="postorder"):
        if c.clades:
            y[c] = float(np.mean([y[ch] for ch in c.clades]))
    return x, y


def draw_tree(ax, tree, title: str, mirror: bool = False, fontsize: float = 7) -> dict:
    """Rectangular tree; support >= 70 in black, < 70 in grey italics ('unsupported')."""
    tree = copy.deepcopy(tree)
    tree.ladderize()  # consistent child order reduces arbitrary crossings in tanglegrams
    sp = species_map()
    x, y = _layout(tree)
    sign = -1 if mirror else 1
    xmax = max(x.values()) or 1.0
    for c in tree.find_clades():
        for ch in c.clades:
            ax.plot([sign * x[c], sign * x[c]], [y[c], y[ch]], color="#444444", lw=0.8)
            ax.plot([sign * x[c], sign * x[ch]], [y[ch], y[ch]], color="#444444", lw=0.8)
        if c.is_terminal():
            ax.text(sign * (x[c] + xmax * 0.02), y[c], c.name, va="center",
                    ha="right" if mirror else "left", fontsize=fontsize,
                    color=SPECIES_COLORS.get(sp.get(c.name, ""), "black"))
        elif c.confidence is not None and c != tree.root:
            ok = c.confidence >= SUPPORTED
            ax.text(sign * x[c] - sign * xmax * 0.01, y[c] - 0.25, f"{int(c.confidence)}",
                    ha="left" if mirror else "right", va="bottom", fontsize=fontsize - 1.5,
                    color="black" if ok else "#999999", style="normal" if ok else "italic")
    ax.set_ylim(len(tree.get_terminals()) - 0.5, -0.8)
    lo, hi = (-xmax * 1.6, xmax * 0.05) if mirror else (-xmax * 0.05, xmax * 1.6)
    ax.set_xlim(lo, hi)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("substitutions/site (K2P)", fontsize=7)
    ax.tick_params(axis="x", labelsize=6)
    ax.set_title(title, fontsize=9.5)
    return {t.name: y[t] for t in tree.get_terminals()}


def _support_note() -> str:
    return (f"Numbers = bootstrap support (%, 100 replicates); grey italic < {int(SUPPORTED)}% "
            "= unsupported branch. Midpoint-rooted for display.")


def plot_reference() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, m in zip(axes, ("nj", "upgma")):
        t = load_tree(f"reference_housekeeping_{m}")
        draw_tree(ax, t, f"{m.upper()} (n={len(t.get_terminals())} genomes)", fontsize=8)
    fig.suptitle(f"Reference tree: concatenated housekeeping genes ({', '.join(housekeeping_genes())})"
                 + mode_tag(), fontweight="bold")
    fig.text(0.5, 0.005, _support_note(), ha="center", fontsize=8)
    save(fig, "m4_reference_tree")


def plot_gene_trees(genes: list[str], name: str, ncols: int = 3) -> None:
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 5.2 * nrows))
    for ax, gene in zip(axes.flat, genes):
        t = load_tree(f"{gene}_nj")
        cls = "virulence" if gene in virulence_genes() else "control"
        draw_tree(ax, t, f"{gene} NJ (n={len(t.get_terminals())})", fontsize=6.5)
        ax.title.set_color(CLASS_COLORS[cls])
    for ax in list(axes.flat)[len(genes):]:
        ax.axis("off")
    label = "virulence-associated" if name.endswith("virulence") else "control"
    fig.suptitle(f"Neighbor-joining gene trees: {label} genes" + mode_tag(), fontweight="bold")
    fig.text(0.5, 0.0, _support_note(), ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    save(fig, name)


def plot_tanglegram(gene: str, ax_pair) -> None:
    ref = load_tree("reference_housekeeping_nj")
    gt = load_tree(f"{gene}_nj")
    keep = {t.name for t in gt.get_terminals()}
    for t in list(ref.get_terminals()):
        if t.name not in keep:
            ref.prune(t)
    sp = species_map()
    y1 = draw_tree(ax_pair[0], ref, f"reference (housekeeping), n={len(keep)}", fontsize=6.5)
    y2 = draw_tree(ax_pair[1], gt, f"{gene} gene tree, n={len(keep)}", mirror=True, fontsize=6.5)
    fig = ax_pair[0].figure
    for name in keep:
        p1 = fig.transFigure.inverted().transform(ax_pair[0].transData.transform((ax_pair[0].get_xlim()[1], y1[name])))
        p2 = fig.transFigure.inverted().transform(ax_pair[1].transData.transform((ax_pair[1].get_xlim()[0], y2[name])))
        fig.add_artist(plt.Line2D([p1[0], p2[0]], [p1[1], p2[1]], transform=fig.transFigure,
                                  color=SPECIES_COLORS.get(sp[name], "grey"), lw=0.8, alpha=0.8))


def plot_discordance(disc: pd.DataFrame, tests: pd.DataFrame) -> None:
    coding = disc[disc["class"] != "rRNA control"].reset_index(drop=True)
    t = tests.set_index("metric")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    colors = [CLASS_COLORS[c] for c in coding["class"]]
    x = np.arange(len(coding))
    axes[0].bar(x, coding["nRF_vs_reference"], color=colors)
    axes[0].set_ylabel("normalised RF distance")
    axes[0].set_title("Gene tree vs reference, all taxa of the gene
(taxon sets differ between genes, so not tested)", fontsize=9.5)
    axes[1].bar(x, coding["nRF_within_Smutans"], color=colors)
    r = t.loc["nRF_within_Smutans"]
    axes[1].set_title(f"Within S. mutans only: Cliff's d = {r.cliffs_delta:.2f} "
                      f"[{r.delta_ci_low:.2f}, {r.delta_ci_high:.2f}], q = {r.q_bh:.3f}", fontsize=9.5)
    within = coding["n_supported_conflicts"] - coding["n_supported_conflicts_between_species"]
    axes[2].bar(x, within, color=colors, alpha=0.5, label="regroups strains of one species")
    axes[2].bar(x, coding["n_supported_conflicts_between_species"], bottom=within, color=colors,
                hatch="///", edgecolor="white", label="clade mixes species")
    r = t.loc["n_supported_conflicts"]
    axes[2].set_title(f"Supported (>= {int(SUPPORTED)}%) conflicting splits: Cliff's d = "
                      f"{r.cliffs_delta:.2f} [{r.delta_ci_low:.2f}, {r.delta_ci_high:.2f}]", fontsize=9.5)
    axes[2].set_ylabel("number of splits")
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.set_xticks(x, [f"{g}\n(n={n})" for g, n in zip(coding["gene"], coding["n_taxa"])],
                      rotation=90, fontsize=7.5)
        for lab, c in zip(ax.get_xticklabels(), coding["class"]):
            lab.set_color(CLASS_COLORS[c])
    fig.suptitle("Gene-tree discordance with the housekeeping reference tree (red = virulence, "
                 "blue = control; discordance is suggestive of HGT/recombination, not proof)"
                 + mode_tag(), fontweight="bold", fontsize=11)
    fig.tight_layout()
    save(fig, "m4_discordance")


def plot_all(disc: pd.DataFrame, tests: pd.DataFrame) -> None:
    """Draw every M4 figure."""
    plot_reference()
    plot_gene_trees(virulence_genes(), "m4_gene_trees_virulence")
    plot_gene_trees(housekeeping_genes() + ["16S"], "m4_gene_trees_controls")
    plot_discordance(disc, tests)
    multi = disc[(disc["n_species"] >= 3)]
    vir = multi[multi["class"] == "virulence"].sort_values("nRF_vs_reference").iloc[-1]["gene"]
    ctl = multi[multi["class"] == "control"].sort_values("nRF_vs_reference").iloc[-1]["gene"]
    fig, axes = plt.subplots(1, 4, figsize=(20, 7.5), gridspec_kw={"wspace": 0.35})
    plot_tanglegram(vir, axes[:2])
    plot_tanglegram(ctl, axes[2:])
    fig.suptitle(f"Tanglegrams: most discordant multi-species virulence gene ({vir}) and "
                 f"control gene ({ctl}) vs the reference tree" + mode_tag(), fontweight="bold")
    fig.text(0.5, 0.0, "Crossing lines = taxa placed differently. " + _support_note(),
             ha="center", fontsize=8)
    save(fig, "m4_tanglegrams")
