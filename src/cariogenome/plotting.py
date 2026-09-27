"""Shared figure style: one colorblind-safe palette and one save function for all modules."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .config import FIGURES  # noqa: E402

# Okabe-Ito palette (Okabe & Ito 2008), safe for the common forms of color blindness.
OKABE_ITO = {
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#7F7F7F",
}
CLASS_COLORS = {"virulence": OKABE_ITO["vermillion"], "control": OKABE_ITO["blue"]}
SPECIES_COLORS = {
    "S. mutans": OKABE_ITO["vermillion"],
    "S. sanguinis": OKABE_ITO["blue"],
    "S. gordonii": OKABE_ITO["green"],
    "S. mitis": OKABE_ITO["purple"],
    "S. salivarius": OKABE_ITO["orange"],
}
SEQUENTIAL_CMAP = "viridis"
DIVERGING_CMAP = "RdBu_r"
DPI = 300


def apply_style() -> None:
    """Set matplotlib defaults used by every figure in the project."""
    plt.rcParams.update(
        {
            "figure.dpi": 100,
            "savefig.dpi": DPI,
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "legend.frameon": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "font.family": "DejaVu Sans",
        }
    )


def save(fig: plt.Figure, name: str) -> Path:
    """Save ``fig`` as ``figures/<name>.png`` at 300 dpi and close it."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / f"{name}.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


apply_style()
