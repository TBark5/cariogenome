"""Figures for M7: conservation on the structure (static + interactive) and Ramachandran."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import py3Dmol
from matplotlib import colormaps

from .config import FIGURES
from .m6_motifs import CATALYTIC
from .plotting import OKABE_ITO, SEQUENTIAL_CMAP, save
from .seqio import mode_tag


def _viridis_hex(n: int = 9) -> list[str]:
    cmap = colormaps[SEQUENTIAL_CMAP]
    return ["#%02x%02x%02x" % tuple(int(255 * c) for c in cmap(i / (n - 1))[:3]) for i in range(n)]


def build_view(pdb_text: str, width: int = 900, height: int = 600) -> py3Dmol.view:
    """Cartoon colored by conservation (B-factor, viridis 0-100), unscored residues grey,
    catalytic residues as magenta sticks."""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "linear", "min": 0,
                                                "max": 100, "colors": _viridis_hex()}}})
    unscored = sorted({int(line[22:26]) for line in pdb_text.splitlines()
                       if line.startswith("ATOM") and float(line[60:66]) < 0})
    if unscored:
        view.setStyle({"resi": unscored}, {"cartoon": {"color": "#BBBBBB"}})
    view.addStyle({"resi": list(CATALYTIC)}, {"stick": {"color": OKABE_ITO["purple"], "radius": 0.3}})
    for pos, (aa, role, _) in CATALYTIC.items():
        view.addLabel(f"{aa}{pos} {role}", {"fontSize": 11, "backgroundColor": "white",
                                            "fontColor": "black", "backgroundOpacity": 0.7},
                      {"resi": pos, "atom": "CA"})
    view.zoomTo({"resi": list(CATALYTIC)})
    view.zoom(0.35)
    return view


def write_html(pdb_text: str, pid: str) -> None:
    """Stand-alone interactive HTML page (open in any browser)."""
    html = build_view(pdb_text)._make_html()
    legend = (f"<p style='font-family:sans-serif'><b>{pid} chain A (S. mutans GtfC)</b>: cartoon "
              "colored by GH70-family conservation (viridis: purple = variable, yellow = "
              "invariant; grey = not scored). Magenta sticks = catalytic residues D477, E515, "
              "D588. Drag to rotate, scroll to zoom.</p>")
    (FIGURES / "m7_structure_conservation.html").write_text(legend + html, encoding="utf-8")


def plot_structure_map(table: pd.DataFrame, summary: dict, pid: str) -> None:
    """2D principal-component projection of C-alpha atoms + conservation vs distance."""
    xyz = table[["x", "y", "z"]].values
    centered = xyz - xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    proj = centered @ vt[:2].T
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.8))
    a1.plot(proj[:, 0], proj[:, 1], color="#DDDDDD", lw=0.5, zorder=1)
    ok = table["conservation"].notna().values
    sc = a1.scatter(proj[ok, 0], proj[ok, 1], c=table.loc[ok, "conservation"], cmap=SEQUENTIAL_CMAP,
                    s=9, vmin=0, vmax=1, zorder=2)
    a1.scatter(proj[~ok, 0], proj[~ok, 1], color="#BBBBBB", s=6, zorder=2, label="not scored")
    cat = table["pdb_resnum"].isin(list(CATALYTIC)).values
    a1.scatter(proj[cat, 0], proj[cat, 1], marker="*", s=220, color=OKABE_ITO["purple"],
               edgecolor="black", zorder=3, label="catalytic D477 / E515 / D588")
    fig.colorbar(sc, ax=a1, shrink=0.8).set_label("GH70 family conservation")
    a1.set_aspect("equal")
    a1.set_xlabel("PC1 of C-alpha coordinates (A)")
    a1.set_ylabel("PC2 (A)")
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), fontsize=8, ncol=2)
    a1.set_title(f"{pid} chain A, C-alpha atoms in 2D (n={len(table)})", fontsize=10)
    s = table.dropna(subset=["conservation"])
    a2.scatter(s["distance_to_active_site"], s["conservation"], s=6, alpha=0.4, color=OKABE_ITO["blue"])
    bins = np.arange(0, s["distance_to_active_site"].max() + 5, 5)
    mids = (bins[:-1] + bins[1:]) / 2
    means = [s.loc[(s["distance_to_active_site"] >= a) & (s["distance_to_active_site"] < b), "conservation"].mean()
             for a, b in zip(bins[:-1], bins[1:])]
    a2.plot(mids, means, color=OKABE_ITO["vermillion"], lw=2, marker="o", ms=4, label="mean per 5 A bin")
    a2.set_xlabel("C-alpha distance to catalytic center (A)")
    a2.set_ylabel("conservation")
    a2.legend(fontsize=8)
    a2.set_title(f"Spearman rho = {summary['spearman_rho_conservation_vs_distance']:.2f} "
                 f"[{summary['rho_ci_low']:.2f}, {summary['rho_ci_high']:.2f}] "
                 f"(n={summary['n_residues_scored']} residues)", fontsize=10)
    fig.suptitle("Conservation in 3D: are conserved residues near the active site?" + mode_tag(),
                 fontweight="bold")
    save(fig, "m7_structure_conservation_map")


def plot_ramachandran(rama: pd.DataFrame, pid: str) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    styles = {"general": OKABE_ITO["blue"], "Gly": OKABE_ITO["orange"], "Pro": OKABE_ITO["green"]}
    for ax, (t, col) in zip(axes, styles.items()):
        sub = rama[rama["type"] == t]
        ax.hist2d(rama["phi"], rama["psi"], bins=72, range=[[-180, 180], [-180, 180]], cmap="Greys",
                  cmin=1, alpha=0.35)
        ax.scatter(sub["phi"], sub["psi"], s=4, color=col, alpha=0.7)
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(0, color="grey", lw=0.5)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xticks(range(-180, 181, 90))
        ax.set_yticks(range(-180, 181, 90))
        ax.set_xlabel("phi (degrees)")
        ax.set_aspect("equal")
        ax.set_title(f"{t} residues (n={len(sub)})", fontsize=10)
    axes[0].set_ylabel("psi (degrees)")
    frac = rama["region"].value_counts(normalize=True) * 100
    txt = "; ".join(f"{k}: {v:.1f}%" for k, v in frac.items())
    fig.suptitle(f"Ramachandran plot, {pid} chain A (grey = all residues)" + mode_tag(),
                 fontweight="bold")
    fig.text(0.5, -0.02, f"Broad regions (simple phi/psi boxes, not a validation score): {txt}",
             ha="center", fontsize=8.5)
    save(fig, "m7_ramachandran")


def plot_all(table, rama, summary, pdb_text, pid) -> None:
    """Draw every M7 figure and the interactive HTML view."""
    plot_structure_map(table, summary, pid)
    plot_ramachandran(rama, pid)
    write_html(pdb_text, pid)
