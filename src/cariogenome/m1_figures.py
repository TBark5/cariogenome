"""Figures for M1: ortholog presence/absence and record-length QC."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle

from .config import all_genes, genomes, load_config, virulence_genes
from .plotting import (
    CLASS_COLORS,
    DARK_GREY,
    FS_ANNOT,
    FS_SMALL,
    LIGHT_GREY,
    LW,
    OKABE_ITO,
    W_FULL,
    save,
    species_color,
)
from .seqio import mode_tag

PRESENT, EXCLUDED, ABSENT = DARK_GREY, OKABE_ITO["yellow"], "#F2F2F2"


def _species_of(label: str) -> str:
    return next(g.species for g in genomes() if g.label == label)


def plot_presence(mat: pd.DataFrame) -> None:
    """Heatmap of ortholog presence (included / excluded / absent) across genomes."""
    fig, ax = plt.subplots(figsize=(W_FULL, 6.0))
    cmap = ListedColormap([ABSENT, EXCLUDED, PRESENT])
    ax.imshow(mat.values, cmap=cmap, vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=60, ha="right")
    ax.set_yticks(range(mat.shape[0]), mat.index)
    n_vir = len(virulence_genes())
    ax.axhline(n_vir - 0.5, color="black", lw=LW)
    ax.set_xticks(np.arange(-0.5, mat.shape[1]), minor=True)
    ax.set_yticks(np.arange(-0.5, mat.shape[0]), minor=True)
    ax.grid(which="minor", color="white", lw=1)  # cell separators, not a data grid
    ax.tick_params(which="minor", length=0)
    for lab in ax.get_xticklabels():
        lab.set_color(species_color(_species_of(lab.get_text())))
    for lab in ax.get_yticklabels():
        cls = "virulence" if lab.get_text() in virulence_genes() else "control"
        lab.set_color(CLASS_COLORS[cls])
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=LIGHT_GREY)
        for c in (PRESENT, EXCLUDED, ABSENT)
    ]
    handles += [Rectangle((0, 0), 1, 1, color=CLASS_COLORS[c]) for c in ("virulence", "control")]
    ax.legend(
        handles,
        [
            "ortholog, passed QC",
            "ortholog, excluded by QC",
            "no ortholog (RBH)",
            "virulence-associated gene (label)",
            "housekeeping gene (label)",
        ],
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
    )
    ax.set_xlabel("Genome (label color = species)")
    ax.set_ylabel("Gene")
    ax.set_title(
        "Reciprocal-best-hit orthologs of the UA159 gene panel across 22 genomes" + mode_tag()
    )
    save(fig, "m1_presence_absence")


def plot_lengths(cat: pd.DataFrame) -> None:
    """Length of every record relative to UA159, with the QC acceptance band."""
    genes = [g for g in all_genes() if g in set(cat["gene"])]
    fig, ax = plt.subplots(figsize=(W_FULL, 5.0))
    cfg = load_config()["qc"]
    ax.axhspan(
        cfg["min_length_fraction"],
        cfg["max_length_fraction"],
        color="#EEEEEE",
        zorder=0,
        label="expected range (80-120% of UA159)",
    )
    jitter = np.random.default_rng(0)  # fixed seed: jitter is cosmetic but reproducible
    for i, gene in enumerate(genes):
        sub = cat[cat["gene"] == gene]
        x = i + jitter.uniform(-0.18, 0.18, len(sub))
        colors = [species_color(s) for s in sub["species"]]
        ax.scatter(x, sub["length_ratio_to_UA159"], c=colors, s=22, edgecolor="none", alpha=0.9)
        bad = ~sub["included"].to_numpy()
        ax.scatter(
            x[bad], sub.loc[bad, "length_ratio_to_UA159"], marker="x", color="black", s=40, lw=1.2
        )
        ax.text(
            i, 0.77, f"n={int(sub['included'].sum())}", ha="center", va="bottom", fontsize=FS_SMALL
        )
    ax.set_ylim(0.75, 1.25)
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right")
    for lab in ax.get_xticklabels():
        g = lab.get_text()
        lab.set_color(
            CLASS_COLORS["virulence"] if g in virulence_genes() else CLASS_COLORS["control"]
        )
    ax.set_xlabel("Gene (n = sequences passing QC)")
    ax.set_ylabel("Record length / UA159 length (ratio)")
    ax.set_title("Record length QC (x = excluded record)" + mode_tag())
    for sp in dict.fromkeys(g.species for g in genomes()):
        ax.scatter([], [], color=species_color(sp), label=sp)
    ax.scatter([], [], marker="x", color="black", label="excluded by QC")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=FS_ANNOT)
    save(fig, "m1_length_qc")
