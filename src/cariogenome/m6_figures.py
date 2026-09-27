"""Figures for M6: family conservation map with catalytic residues, and sequence logos."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

from .m6_motifs import information_content
from .plotting import (
    CATALYTIC,
    FS_ANNOT,
    FS_SMALL,
    LIGHT_GREY,
    LW,
    LW_HAIR,
    LW_THIN,
    OKABE_ITO,
    RESIDUE_COLORS,
    RESIDUE_GROUPS,
    W_FULL,
    save,
)
from .seqio import mode_tag

FONT = FontProperties(family="DejaVu Sans Mono", weight="bold")
MOTIF = OKABE_ITO["green"]


def draw_letter(
    ax: Axes, letter: str, x: float, y: float, width: float, height: float, color: str
) -> None:
    """Draw one letter scaled to fill the box (x, y, width, height)."""
    if height <= 0:
        return
    tp = TextPath((0, 0), letter, size=1, prop=FONT)
    bb = tp.get_extents()
    t = (
        Affine2D()
        .translate(-bb.x0, -bb.y0)
        .scale(width / bb.width * 0.95, height / bb.height)
        .translate(x - width / 2 * 0.95, y)
    )
    ax.add_patch(PathPatch(t.transform_path(tp), color=color, lw=0))


def draw_logo(ax: Axes, freq: pd.DataFrame, title: str, highlight: int | None = None) -> None:
    """Information-content sequence logo (letter height = frequency x IC)."""
    ic = information_content(freq)
    for i, (col, row) in enumerate(freq.iterrows()):
        y = 0.0
        for aa, f in sorted(row.items(), key=lambda kv: kv[1]):
            h = f * ic[col]
            if h > 0.01:
                draw_letter(ax, str(aa), i, y, 0.9, h, RESIDUE_COLORS.get(str(aa), LIGHT_GREY))
                y += h
    if highlight is not None:
        ax.axvspan(
            highlight - 0.5, highlight + 0.5, color=OKABE_ITO["yellow"], alpha=0.35, zorder=0
        )
    ax.set_xlim(-0.6, len(freq) - 0.4)
    ax.set_ylim(0, np.log2(20))
    ax.set_xticks(range(len(freq)), [str(c) for c in freq.index], fontsize=FS_SMALL, rotation=90)
    ax.set_xlabel("Alignment column", fontsize=FS_ANNOT)
    ax.set_ylabel("Information (bits)", fontsize=FS_ANNOT)
    ax.set_title(title)


def plot_map(
    cons: pd.DataFrame, ranks: pd.DataFrame, motifs: pd.DataFrame, n: int, n_species: int
) -> None:
    """Family conservation along GtfC with motifs and catalytic residues."""
    c = cons.dropna(subset=["ref_position"])
    fig, ax = plt.subplots(figsize=(W_FULL, 4.8))
    ax.plot(
        c["ref_position"],
        c["conservation"],
        color="#BBBBBB",
        lw=LW_HAIR,
        label="per-site conservation",
    )
    sm = c.set_index("ref_position")["conservation"].rolling(15, center=True, min_periods=5).mean()
    ax.plot(sm.index, sm.to_numpy(), color="black", lw=LW, label="15-residue running mean")
    ax.axvspan(
        250,
        1050,
        color=OKABE_ITO["sky"],
        alpha=0.1,
        label="catalytic region (UniProt P13470, approx.)",
    )
    for k, m in enumerate(motifs.sort_values("start").itertuples()):
        ax.axvspan(m.start, m.end, color=MOTIF, alpha=0.35, lw=0)
        ax.text(
            (m.start + m.end) / 2,
            1.05 + 0.06 * (k % 2),
            m.motif,
            ha="center",
            fontsize=FS_SMALL,
            color=MOTIF,
        )
    for k, r in enumerate(ranks.itertuples()):
        ax.axvline(r.position, color=CATALYTIC, lw=LW_THIN, ls="--", ymax=0.9)
        ax.text(
            r.position,
            0.02 + 0.08 * k,
            f" {r.residue}{r.position}",
            ha="left",
            fontsize=FS_ANNOT,
            color=CATALYTIC,
            fontweight="bold",
        )
    ax.set_ylim(0, 1.16)
    ax.set_xlim(0, c["ref_position"].max())
    ax.set_xlabel("Residue position in S. mutans UA159 GtfC (amino acids)")
    ax.set_ylabel("Conservation, 1 - H / log2 20 (0-1)")
    ax.plot(
        [],
        [],
        color=CATALYTIC,
        ls="--",
        lw=LW_THIN,
        label="catalytic residues (canonical GH70 motifs)",
    )
    ax.add_patch(
        Rectangle((0, 0), 0, 0, color=MOTIF, alpha=0.35, label="top-6 conserved 10-residue motifs")
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
    ax.set_title(
        f"Glucansucrase (GH70) family conservation mapped on GtfC (n={n} sequences, "
        f"{n_species} species)" + mode_tag()
    )
    save(fig, "m6_gtf_family_conservation")


def plot_logos(pwms: dict[str, pd.DataFrame], n: int) -> None:
    """Grid of sequence logos: top motifs, then the three catalytic regions."""
    names = list(pwms)
    ncols = 3
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(W_FULL, 2.8 * nrows))
    for ax, name in zip(axes.flat, names, strict=False):  # the grid may have spare panels
        freq = pwms[name]
        draw_logo(ax, freq, name, len(freq) // 2 if "around" in name else None)
    for ax in list(axes.flat)[len(names) :]:
        ax.axis("off")
    handles = [Rectangle((0, 0), 1, 1, color=c) for _, c in RESIDUE_GROUPS.values()]
    handles.append(Rectangle((0, 0), 1, 1, color=OKABE_ITO["yellow"], alpha=0.35))
    fig.legend(
        handles,
        [*RESIDUE_GROUPS, "catalytic residue"],
        loc="lower center",
        ncol=8,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(
        f"Sequence logos of the GH70 family (n={n} sequences, Henikoff-weighted): "
        "most conserved motifs (top) and catalytic motifs (bottom)" + mode_tag()
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, "m6_motif_logos")


def plot_all(
    cons: pd.DataFrame,
    ranks: pd.DataFrame,
    motifs: pd.DataFrame,
    pwms: dict[str, pd.DataFrame],
    n: int,
    n_species: int,
) -> None:
    """Draw every M6 figure."""
    plot_map(cons, ranks, motifs, n, n_species)
    plot_logos(pwms, n)
