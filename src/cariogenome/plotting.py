"""The one shared figure style. Every figure is created and saved through this module.

- Fonts, line widths, DPI, margins and backgrounds are set once here (``apply_style``);
  figure modules use the named constants below instead of literal sizes.
- Color semantics are fixed: a gene class, a species or a significance level has the
  same color in every figure.
- ``save`` writes ``figures/NN_<name>.png`` where NN is the figure's position in the
  pipeline (registry in ``captions.py``), always at 300 dpi on an opaque white background
  so figures stay readable on GitHub's dark theme.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (the backend must be chosen before pyplot loads)
from matplotlib.figure import Figure  # noqa: E402

from .captions import figure_filename  # noqa: E402
from .config import FIGURES, housekeeping_genes, virulence_genes  # noqa: E402

# Okabe-Ito palette (Okabe & Ito 2008), safe for the common forms of color blindness.
OKABE_ITO = {
    "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
    "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000",
}
GREY = "#7F7F7F"
LIGHT_GREY = "#D9D9D9"
DARK_GREY = "#404040"

# Fixed color semantics used by every figure.
CLASS_COLORS = {"virulence": OKABE_ITO["vermillion"], "control": OKABE_ITO["blue"]}
SPECIES_COLORS = {  # never vermillion or blue, which are reserved for the gene classes
    "S. mutans": OKABE_ITO["black"],
    "S. sanguinis": OKABE_ITO["sky"],
    "S. gordonii": OKABE_ITO["green"],
    "S. mitis": OKABE_ITO["purple"],
    "S. salivarius": OKABE_ITO["orange"],
}
SIGNIFICANT, NOT_SIGNIFICANT = OKABE_ITO["black"], GREY
CATALYTIC = OKABE_ITO["black"]
# Amino-acid classes in alignments and logos: Paul Tol's colorblind-safe "light" scheme,
# kept separate from the saturated gene-class colors above.
RESIDUE_GROUPS = {"hydrophobic AILMV": ("AILMV", "#77AADD"), "aromatic FWY": ("FWY", "#99DDFF"),
                  "basic KRH": ("KRH", "#EE8866"), "acidic DE": ("DE", "#FFAABB"),
                  "polar STNQ": ("STNQ", "#44BB99"), "G/P": ("GP", "#EEDD88"),
                  "C": ("C", "#BBCC33")}
RESIDUE_COLORS = {aa: col for letters, col in RESIDUE_GROUPS.values() for aa in letters}
SEQUENTIAL_CMAP = "viridis"
DIVERGING_CMAP = "PuOr_r"  # purple-orange: avoids the red/blue of the gene classes
DPI = 300

# Font sizes (points): the only sizes any figure may use.
FS_SUPTITLE = 12.0
FS_TITLE = 10.5
FS_LABEL = 10.0
FS_TICK = 8.5
FS_ANNOT = 8.0
FS_SMALL = 6.5

# Line widths (points).
LW_HAIR = 0.5
LW_THIN = 0.8
LW = 1.4
LW_THICK = 2.0

# Figure widths (inches) per figure class; heights follow the content.
W_SINGLE = 7.0   # one panel
W_PAIR = 10.5    # two side-by-side panels (class comparisons, maps)
W_FULL = 14.0    # a row of panels or one wide panel
W_XL = 18.0      # dense grids (identity heatmaps, tree grids, tanglegrams)


def apply_style() -> None:
    """Set the matplotlib defaults shared by every figure."""
    plt.rcParams.update({
        "figure.dpi": 100, "savefig.dpi": DPI,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "savefig.transparent": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.15,
        "font.family": "DejaVu Sans", "font.size": FS_LABEL,
        "figure.titlesize": FS_SUPTITLE, "figure.titleweight": "bold",
        "axes.titlesize": FS_TITLE, "axes.titleweight": "bold", "axes.labelsize": FS_LABEL,
        "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK, "legend.fontsize": FS_ANNOT,
        "legend.frameon": False, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": LW_THIN, "lines.linewidth": LW, "patch.linewidth": LW_THIN,
        "axes.grid": False,
    })


def class_n_text() -> str:
    """Standard sample-size note for figures that compare gene classes."""
    return (f"n = {len(virulence_genes())} virulence-associated, "
            f"n = {len(housekeeping_genes())} housekeeping genes")


def annotate_class_n(fig: Figure, y: float = -0.01) -> None:
    """Write the class sample sizes under the figure, in the same place everywhere."""
    fig.text(0.5, y, class_n_text(), ha="center", va="top", fontsize=FS_ANNOT, color=DARK_GREY)


def title_with_stats(ax, title: str, stats: str) -> None:
    """Bold panel title with a smaller, non-bold statistics line underneath."""
    ax.set_title(title, pad=16 + 11 * stats.count(chr(10)))  # chr(10): one more line of stats
    ax.text(0.5, 1.012, stats, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=FS_ANNOT)


def delta_stats(r) -> str:
    """'Cliff's δ = d [lo, hi], q = q' from a stats row (Cliff's delta with CI and BH q)."""
    return (f"Cliff's δ = {r.cliffs_delta:.2f} [{r.delta_ci_low:.2f}, {r.delta_ci_high:.2f}], "
            f"q = {r.q_bh:.3f}")


def species_color(species: str) -> str:
    return SPECIES_COLORS.get(species, DARK_GREY)


def figure_path(name: str, ext: str = "png") -> Path:
    """Path of a registered figure (``figures/NN_<name>.<ext>``)."""
    return FIGURES / f"{figure_filename(name)}.{ext}"


def save(fig: Figure, name: str) -> Path:
    """Save ``fig`` as ``figures/NN_<name>.png`` at 300 dpi on white, then close it."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = figure_path(name)
    fig.savefig(path, dpi=DPI, facecolor="white", transparent=False)
    plt.close(fig)
    return path


apply_style()
