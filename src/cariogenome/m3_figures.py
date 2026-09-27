"""Figures for M3 (alignment and conservation)."""

from __future__ import annotations

from collections.abc import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from .config import RESULTS, all_genes, genomes, virulence_genes
from .plotting import (
    CLASS_COLORS,
    DARK_GREY,
    FS_ANNOT,
    FS_SMALL,
    LIGHT_GREY,
    LW,
    LW_HAIR,
    OKABE_ITO,
    RESIDUE_COLORS,
    RESIDUE_GROUPS,
    SEQUENTIAL_CMAP,
    W_FULL,
    W_PAIR,
    W_XL,
    annotate_class_n,
    delta_stats,
    save,
    species_color,
    title_with_stats,
)
from .seqio import mode_tag, read_fasta, species_map

MOST, LEAST = OKABE_ITO["green"], OKABE_ITO["orange"]


def _order(labels: list[str]) -> list[str]:
    rank = {g.label: i for i, g in enumerate(genomes())}
    return sorted(labels, key=lambda k: rank.get(k, 999))


def _gene_color(gene: str) -> str:
    return CLASS_COLORS["virulence" if gene in virulence_genes() else "control"]


def _color_species_labels(labels: Iterable[Text], sp: dict[str, str]) -> None:
    for t in labels:
        t.set_color(species_color(sp[t.get_text()]))


def _species_blocks(order: list[str], sp: dict[str, str]) -> list[tuple[str, int, int]]:
    """(species, first index, last index) for consecutive runs of one species."""
    blocks: list[tuple[str, int, int]] = []
    for i, lab in enumerate(order):
        if blocks and blocks[-1][0] == sp[lab]:
            blocks[-1] = (blocks[-1][0], blocks[-1][1], i)
        else:
            blocks.append((sp[lab], i, i))
    return blocks


def plot_identity(summary: pd.DataFrame) -> None:
    """Grid of pairwise percent-identity heatmaps, one per gene, labelled by species block."""
    fig, axes = plt.subplots(3, 5, figsize=(W_XL, 10.5))
    fig.subplots_adjust(wspace=0.55, hspace=0.35, left=0.08, right=0.88)
    sp = species_map()
    im = None
    for ax, gene in zip(axes.flat, all_genes(), strict=True):
        m = pd.read_csv(RESULTS / "identity" / f"{gene}.csv", index_col=0)
        order = _order(list(m.index))
        im = ax.imshow(m.loc[order, order].to_numpy(), cmap=SEQUENTIAL_CMAP, vmin=50, vmax=100)
        blocks = _species_blocks(order, sp)
        for _, _, last in blocks[:-1]:
            ax.axhline(last + 0.5, color="white", lw=LW_HAIR)
            ax.axvline(last + 0.5, color="white", lw=LW_HAIR)
        centers = [(first + last) / 2 for _, first, last in blocks]
        names = [f"{spc} ({last - first + 1})" for spc, first, last in blocks]
        ax.set_yticks(centers, names, fontsize=FS_SMALL)
        for t, (spc, _, _) in zip(ax.get_yticklabels(), blocks, strict=True):
            t.set_color(species_color(spc))
        ax.set_xticks([])
        kind = "nt" if gene == "16S" else "aa"
        ax.set_title(f"{gene} ({kind}, n={len(order)})", color=_gene_color(gene))
        ax.set_xlabel("same strains as rows", fontsize=FS_SMALL)
        ax.spines[:].set_visible(False)
    cax = fig.add_axes((0.91, 0.3, 0.012, 0.4))
    assert im is not None  # the loop above always draws at least one heatmap
    fig.colorbar(im, cax=cax).set_label("Percent identity over aligned positions (%)")
    fig.suptitle(
        "Pairwise percent identity per gene (strains ordered by species; "
        "vermillion titles = virulence-associated, blue = housekeeping)" + mode_tag(),
        y=0.98,
    )
    annotate_class_n(fig, y=0.05)
    save(fig, "m3_identity_heatmaps")


def plot_profiles(summary: pd.DataFrame, cons_all: dict, regions: pd.DataFrame) -> None:
    """Conservation along each virulence protein and two controls."""
    genes = [*virulence_genes(), "recA", "rpoB"]
    fig, axes = plt.subplots(len(genes), 1, figsize=(W_FULL, 2.0 * len(genes)))
    s = summary.set_index("gene")
    for ax, gene in zip(axes, genes, strict=True):
        c = cons_all[gene].dropna(subset=["ref_position"])
        win = int(regions.loc[regions["gene"] == gene, "window"].iloc[0])
        ax.plot(
            c["ref_position"],
            c["conservation"],
            color="#CCCCCC",
            lw=LW_HAIR,
            label="per-site conservation",
        )
        sm = (
            c.set_index("ref_position")["conservation"]
            .rolling(win, center=True, min_periods=win // 2)
            .mean()
        )
        ax.plot(
            sm.index,
            sm.to_numpy(),
            color=_gene_color(gene),
            lw=LW,
            label=f"{win}-residue running mean",
        )
        for r in regions[regions["gene"] == gene].itertuples():
            ax.axvspan(
                r.start,
                r.end,
                color=MOST if r.type == "most conserved" else LEAST,
                alpha=0.25,
                lw=0,
            )
        ax.set_ylim(-0.02, 1.05)
        ax.set_xlim(0, c["ref_position"].max())
        ax.set_ylabel("conservation\n(0-1)")
        ax.set_title(
            f"{gene}: n={s.at[gene, 'n_sequences']} sequences from "
            f"{s.at[gene, 'n_species']} species",
            loc="left",
            color=_gene_color(gene),
        )
    axes[-1].set_xlabel("Residue position in S. mutans UA159 (amino acids)")
    handles = [
        *axes[0].get_legend_handles_labels()[0],
        Rectangle((0, 0), 1, 1, color=MOST, alpha=0.4),
        Rectangle((0, 0), 1, 1, color=LEAST, alpha=0.4),
    ]
    labels = [
        "per-site conservation",
        "running mean (class color)",
        "3 most conserved windows",
        "3 least conserved windows",
    ]
    fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(
        "Per-site conservation (1 - H / log2 20, Henikoff-weighted) along each "
        "virulence-associated protein and two housekeeping controls" + mode_tag()
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.975))
    save(fig, "m3_conservation_profiles")


def plot_alignment_residues(gene: str = "luxS", block: int = 80) -> None:
    """Residue-level alignment view colored by amino-acid class."""
    msa = read_fasta(RESULTS / "alignments" / f"{gene}.aa.fasta")
    labels = _order(list(msa))
    n_cols = len(msa[labels[0]])
    n_blocks = int(np.ceil(n_cols / block))
    fig, axes = plt.subplots(n_blocks, 1, figsize=(W_XL, 0.19 * len(labels) * n_blocks + 1.6))
    sp = species_map()
    for b, ax in enumerate(np.atleast_1d(axes)):
        lo, hi = b * block, min(n_cols, (b + 1) * block)
        for i, lab in enumerate(labels):
            for j in range(lo, hi):
                ch = msa[lab][j]
                if ch != "-":
                    ax.add_patch(
                        Rectangle(
                            (j - 0.5, i - 0.5), 1, 1, lw=0, color=RESIDUE_COLORS.get(ch, LIGHT_GREY)
                        )
                    )
                ax.text(
                    j, i, ch, ha="center", va="center", fontsize=FS_SMALL - 1.3, family="monospace"
                )
        ax.set_xlim(lo - 0.5, lo + block - 0.5)
        ax.set_ylim(len(labels) - 0.5, -0.5)
        ax.set_yticks(range(len(labels)), labels, fontsize=FS_SMALL)
        _color_species_labels(ax.get_yticklabels(), sp)
        ticks = list(range(lo, hi, 10))
        ax.set_xticks(ticks, [str(t + 1) for t in ticks])
        ax.set_xlabel("Alignment column")
        ax.spines[:].set_visible(False)
    fig.legend(
        [Rectangle((0, 0), 1, 1, color=c) for _, c in RESIDUE_GROUPS.values()],
        list(RESIDUE_GROUPS),
        loc="lower center",
        ncol=7,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(
        f"{gene} protein alignment (center-star progressive; n={len(labels)} sequences, "
        f"{len({sp[k] for k in labels})} species; label color = species)" + mode_tag()
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    save(fig, f"m3_alignment_{gene}")


def plot_alignment_overview(gene: str = "gtfD") -> None:
    """Whole-alignment map: same as UA159, different, or gap."""
    msa = read_fasta(RESULTS / "alignments" / f"{gene}.aa.fasta")
    labels = _order(list(msa))
    ref = msa["Smu_UA159"]
    mat = np.array(
        [
            [0 if ch == "-" else (2 if ch == ref[j] else 1) for j, ch in enumerate(msa[k])]
            for k in labels
        ]
    )
    fig, ax = plt.subplots(figsize=(W_FULL, 0.25 * len(labels) + 1.8))
    colors = ["white", DARK_GREY, LIGHT_GREY]
    ax.imshow(
        mat, cmap=ListedColormap(colors), aspect="auto", vmin=0, vmax=2, interpolation="nearest"
    )
    ax.set_yticks(range(len(labels)), labels, fontsize=FS_ANNOT)
    _color_species_labels(ax.get_yticklabels(), species_map())
    ax.set_xlabel("Alignment column")
    ax.set_ylabel("Strain (label color = species)")
    ax.legend(
        [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=DARK_GREY) for c in colors[::-1]],
        ["same residue as UA159", "different residue", "gap"],
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
    )
    ax.set_title(
        f"{gene} alignment overview (n={len(labels)} sequences) relative to "
        "S. mutans UA159" + mode_tag()
    )
    save(fig, f"m3_alignment_overview_{gene}")


def plot_summary(summary: pd.DataFrame, tests: pd.DataFrame) -> None:
    """Class comparison of within-S. mutans identity and entropy."""
    coding = summary[summary["class"] != "rRNA control"]
    fig, axes = plt.subplots(1, 2, figsize=(W_PAIR, 4.8))
    t = tests.set_index("metric")
    metrics = [
        ("mean_pid_Smutans", "Mean pairwise aa identity (%)"),
        ("mean_entropy_Smutans", "Mean per-site entropy (bits)"),
    ]
    for ax, (m, lab) in zip(axes, metrics, strict=True):
        for i, cls in enumerate(["virulence", "control"]):
            sub = coding[coding["class"] == cls]
            x = i + np.linspace(-0.15, 0.15, len(sub))
            ax.scatter(x, sub[m], color=CLASS_COLORS[cls], s=32, zorder=3)
            for xi, v, g in zip(x, sub[m], sub["gene"], strict=True):
                ax.annotate(
                    g, (xi, v), xytext=(4, 0), textcoords="offset points", fontsize=FS_SMALL
                )
            ax.hlines(
                np.median(sub[m]), i - 0.3, i + 0.3, color="black", label="median" if i else None
            )
        r = t.loc[m]
        ax.set_xticks(
            [0, 1],
            [
                f"virulence\n(n={int(r.n_virulence)} genes)",
                f"housekeeping\n(n={int(r.n_control)} genes)",
            ],
        )
        ax.set_xlim(-0.6, 1.8)
        ax.set_xlabel("Gene class")
        ax.set_ylabel(lab)
        title_with_stats(ax, lab.split(" (")[0], delta_stats(r))
        ax.legend(loc="best")
    n = coding["n_Smutans"]
    fig.suptitle(
        f"Conservation within S. mutans ({n.min()}-{n.max()} strains per gene)" + mode_tag()
    )
    fig.tight_layout()
    annotate_class_n(fig)
    save(fig, "m3_conservation_summary")


def plot_all(
    summary: pd.DataFrame,
    cons_all: dict[str, pd.DataFrame],
    regions: pd.DataFrame,
    tests: pd.DataFrame,
) -> None:
    """Draw every M3 figure."""
    plot_identity(summary)
    plot_profiles(summary, cons_all, regions)
    plot_alignment_residues("luxS")
    plot_alignment_overview("gtfD")
    plot_summary(summary, tests)
