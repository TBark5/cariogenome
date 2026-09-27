"""Figures for M6: family conservation map with catalytic residues, and sequence logos."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

from .m3_figures import AA_GROUPS
from .m6_motifs import CATALYTIC, information_content
from .plotting import OKABE_ITO, save
from .seqio import mode_tag

FONT = FontProperties(family="DejaVu Sans Mono", weight="bold")


def draw_letter(ax, letter: str, x: float, y: float, width: float, height: float, color: str) -> None:
    """Draw one letter scaled to fill the box (x, y, width, height)."""
    if height <= 0:
        return
    tp = TextPath((0, 0), letter, size=1, prop=FONT)
    bb = tp.get_extents()
    t = (Affine2D().translate(-bb.x0, -bb.y0).scale(width / bb.width * 0.95, height / bb.height)
         .translate(x - width / 2 * 0.95, y))
    ax.add_patch(PathPatch(t.transform_path(tp), color=color, lw=0))


def draw_logo(ax, freq: pd.DataFrame, title: str, highlight: int | None = None) -> None:
    """Information-content sequence logo (letter height = frequency x IC)."""
    ic = information_content(freq)
    for i, (col, row) in enumerate(freq.iterrows()):
        y = 0.0
        for aa, f in sorted(row.items(), key=lambda kv: kv[1]):
            h = f * ic[col]
            if h > 0.01:
                draw_letter(ax, aa, i, y, 0.9, h, AA_GROUPS.get(aa, "grey"))
                y += h
    if highlight is not None:
        ax.axvspan(highlight - 0.5, highlight + 0.5, color=OKABE_ITO["yellow"], alpha=0.35, zorder=0)
    ax.set_xlim(-0.6, len(freq) - 0.4)
    ax.set_ylim(0, np.log2(20))
    ax.set_xticks(range(len(freq)), [str(c) for c in freq.index], fontsize=6, rotation=90)
    ax.set_ylabel("bits", fontsize=8)
    ax.set_title(title, fontsize=9)


def plot_map(cons: pd.DataFrame, ranks: pd.DataFrame, motifs: pd.DataFrame, n: int,
             n_species: int) -> None:
    c = cons.dropna(subset=["ref_position"])
    fig, ax = plt.subplots(figsize=(14, 4.2))
    ax.plot(c["ref_position"], c["conservation"], color="#BBBBBB", lw=0.5)
    sm = c.set_index("ref_position")["conservation"].rolling(15, center=True, min_periods=5).mean()
    ax.plot(sm.index, sm.values, color=OKABE_ITO["blue"], lw=1.3, label="15-residue running mean")
    ax.axvspan(250, 1050, color=OKABE_ITO["sky"], alpha=0.08, label="catalytic region (UniProt P13470, approx.)")
    for k, m in enumerate(motifs.sort_values("start").itertuples()):
        ax.axvspan(m.start, m.end, color=OKABE_ITO["green"], alpha=0.35, lw=0)
        ax.text((m.start + m.end) / 2, 1.05 + 0.06 * (k % 2), m.motif, ha="center", fontsize=7,
                color=OKABE_ITO["green"])
    for k, r in enumerate(ranks.itertuples()):
        ax.axvline(r.position, color=OKABE_ITO["vermillion"], lw=1, ymax=0.9)
        ax.text(r.position, 0.02 + 0.08 * k, f" {r.residue}{r.position}", ha="left", fontsize=7.5,
                color=OKABE_ITO["vermillion"], fontweight="bold")
    ax.set_ylim(0, 1.16)
    ax.set_xlim(0, c["ref_position"].max())
    ax.set_xlabel("Residue position in S. mutans UA159 GtfC")
    ax.set_ylabel("conservation (1 - H / log2 20)")
    ax.plot([], [], color=OKABE_ITO["vermillion"], label="catalytic residues (Ito et al. 2011)")
    ax.fill_between([], [], color=OKABE_ITO["green"], alpha=0.35, label="top-6 conserved 10-residue motifs")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=7.5, ncol=4)
    ax.set_title(f"Glucansucrase (GH70) family conservation mapped on GtfC (n={n} sequences, "
                 f"{n_species} species)" + mode_tag(), fontsize=10.5)
    save(fig, "m6_gtf_family_conservation")


def plot_logos(pwms: dict[str, pd.DataFrame], n: int) -> None:
    names = list(pwms)
    ncols = 3
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 2.6 * nrows))
    for ax, name in zip(axes.flat, names):
        freq = pwms[name]
        hl = None
        if "around" in name:
            hl = len(freq) // 2
        draw_logo(ax, freq, name, hl)
    for ax in list(axes.flat)[len(names):]:
        ax.axis("off")
    fig.suptitle(f"Sequence logos of the GH70 family (n={n} sequences, Henikoff-weighted). Top: "
                 "most conserved motifs; bottom: catalytic motifs (yellow = catalytic residue).\n"
                 "x-axis = alignment column" + mode_tag(), fontweight="bold", fontsize=10.5)
    fig.tight_layout()
    save(fig, "m6_motif_logos")


def plot_all(cons, ranks, motifs, pwms, n, n_species) -> None:
    """Draw every M6 figure."""
    plot_map(cons, ranks, motifs, n, n_species)
    plot_logos(pwms, n)
