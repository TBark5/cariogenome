"""Figures for M7: conservation on the structure (static + interactive) and Ramachandran.

The interactive view is a single self-contained HTML file: the 3Dmol.js library (v2.5.5,
BSD-3-Clause, vendored in data/vendor/) is inlined, so it works without internet. A static
300 dpi PNG of the same model and coloring is written next to it.
"""

from __future__ import annotations

from itertools import pairwise

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import py3Dmol
from matplotlib import colormaps
from matplotlib.artist import Artist
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from .config import ROOT
from .m6_motifs import CATALYTIC as CATALYTIC_SITES
from .plotting import (
    CATALYTIC,
    DARK_GREY,
    FS_ANNOT,
    GREY,
    LIGHT_GREY,
    LW_HAIR,
    LW_THICK,
    OKABE_ITO,
    SEQUENTIAL_CMAP,
    W_FULL,
    W_PAIR,
    figure_path,
    save,
)
from .seqio import mode_tag

VENDOR_JS = ROOT / "data" / "vendor" / "3Dmol-min.js"
CDN_LOADER = "loadScriptAsync('https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js')"


def _viridis_hex(n: int = 9) -> list[str]:
    cmap = colormaps[SEQUENTIAL_CMAP]
    return [
        "#{:02x}{:02x}{:02x}".format(*(int(255 * c) for c in cmap(i / (n - 1))[:3]))
        for i in range(n)
    ]


def build_view(pdb_text: str, width: int = 900, height: int = 600) -> py3Dmol.view:
    """Cartoon colored by conservation (B-factor, viridis 0-100), unscored residues grey,
    catalytic residues as black sticks."""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")
    view.setStyle(
        {
            "cartoon": {
                "colorscheme": {
                    "prop": "b",
                    "gradient": "linear",
                    "min": 0,
                    "max": 100,
                    "colors": _viridis_hex(),
                }
            }
        }
    )
    unscored = sorted(
        {
            int(line[22:26])
            for line in pdb_text.splitlines()
            if line.startswith("ATOM") and float(line[60:66]) < 0
        }
    )
    if unscored:
        view.setStyle({"resi": unscored}, {"cartoon": {"color": LIGHT_GREY}})
    view.addStyle({"resi": list(CATALYTIC_SITES)}, {"stick": {"color": CATALYTIC, "radius": 0.3}})
    for pos, (aa, role, _) in CATALYTIC_SITES.items():
        view.addLabel(
            f"{aa}{pos} {role}",
            {
                "fontSize": 11,
                "backgroundColor": "white",
                "fontColor": "black",
                "backgroundOpacity": 0.7,
            },
            {"resi": pos, "atom": "CA"},
        )
    view.zoomTo({"resi": list(CATALYTIC_SITES)})
    view.zoom(0.35)
    return view


def write_html(pdb_text: str, pid: str) -> None:
    """Self-contained interactive HTML page (no network needed to open it)."""
    view = build_view(pdb_text)
    html = view._make_html().replace(str(view.uniqueid), "cariogenome3d")  # deterministic file
    if CDN_LOADER not in html:
        raise RuntimeError("py3Dmol HTML template changed; cannot inline 3Dmol.js")
    html = html.replace(CDN_LOADER, "Promise.resolve()")  # library is inlined below instead
    library = VENDOR_JS.read_text(encoding="utf-8")
    png = figure_path("m7_structure_3d").name
    fallback = (
        "<script>setTimeout(function(){if(!window.viewer_cariogenome3d){"
        "var d=document.getElementById('3dmolviewer_cariogenome3d');"
        'd.innerHTML=\'<p style="font-family:sans-serif;background:#ffe9cc;padding:8px">'
        "The 3D viewer could not start in this browser (WebGL may be disabled). "
        f"See the static image {png} in the same folder.</p>';}}}},3000);</script>"
    )
    legend = (
        f"<p style='font-family:sans-serif'><b>{pid} chain A (S. mutans GtfC)</b>: "
        "cartoon colored by GH70-family conservation (viridis: purple = variable, "
        "yellow = invariant; grey = not scored). Black sticks = catalytic residues D477, "
        "E515, D588. Drag to rotate, scroll to zoom. Works offline (3Dmol.js 2.5.5, "
        "BSD-3-Clause, is embedded in this file).</p>"
    )
    page = (
        f"<!DOCTYPE html>\n<html><head><meta charset='utf-8'><title>{pid} conservation"
        f"</title>\n<script>{library}</script>\n</head><body>\n{legend}\n{html}\n{fallback}"
        "\n</body></html>\n"
    )
    figure_path("m7_structure_3d", "html").write_text(page, encoding="utf-8")


def _principal_axes(xyz: np.ndarray) -> np.ndarray:
    """Coordinates rotated onto their principal axes (largest spread first)."""
    centered = xyz - xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return centered @ vt.T


def plot_structure_3d(table: pd.DataFrame, pid: str) -> None:
    """Static 3D rendering of the C-alpha trace colored by conservation."""
    xyz = _principal_axes(table[["x", "y", "z"]].to_numpy(dtype=float))
    resnum = table["pdb_resnum"].to_numpy()
    cons = table["conservation"].to_numpy(dtype=float)
    link = (np.diff(resnum) == 1) & (np.linalg.norm(np.diff(xyz, axis=0), axis=1) < 4.3)
    segs = np.stack([xyz[:-1], xyz[1:]], axis=1)[link]
    a, b = cons[:-1][link], cons[1:][link]
    # mean of the scored endpoints; NaN (drawn grey) only if neither endpoint is scored
    seg_cons = np.where(np.isnan(a), b, np.where(np.isnan(b), a, (a + b) / 2))
    cmap = colormaps[SEQUENTIAL_CMAP]
    colors = [LIGHT_GREY if np.isnan(v) else cmap(v) for v in seg_cons]
    fig = plt.figure(figsize=(W_PAIR, 8.0))
    ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Line3DCollection(segs, colors=colors, linewidths=LW_THICK))
    cat = np.isin(resnum, list(CATALYTIC_SITES))
    ax.scatter(
        *xyz[cat].T,
        marker="*",
        s=420,
        color=CATALYTIC,
        edgecolor="white",
        depthshade=False,
        label="catalytic D477 / E515 / D588",
    )
    lo, hi = xyz.min(axis=0) - 2, xyz.max(axis=0) + 2
    ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), zlim=(lo[2], hi[2]))
    ax.set_box_aspect(tuple(hi - lo), zoom=1.25)  # true proportions, filling the panel
    ax.view_init(elev=22, azim=-58)
    ax.set_xlabel("PC1 (A)")
    ax.set_ylabel("PC2 (A)")
    ax.set_zlabel("PC3 (A)")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    fig.colorbar(sm, ax=ax, orientation="horizontal", shrink=0.5, pad=0.02).set_label(
        "GH70 family conservation (0-1; grey = not scored)"
    )
    ax.legend(loc="upper left")
    ax.set_title(
        f"{pid} chain A (S. mutans GtfC), C-alpha trace colored by conservation "
        f"(n={len(table)} residues)" + mode_tag()
    )
    save(fig, "m7_structure_3d")


def plot_structure_map(table: pd.DataFrame, summary: dict, pid: str) -> None:
    """2D principal-component projection of C-alpha atoms + conservation vs distance."""
    proj = _principal_axes(table[["x", "y", "z"]].to_numpy(dtype=float))[:, :2]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W_FULL, 6.0))
    a1.plot(proj[:, 0], proj[:, 1], color="#DDDDDD", lw=LW_HAIR, zorder=1)
    ok = table["conservation"].notna().to_numpy()
    sc = a1.scatter(
        proj[ok, 0],
        proj[ok, 1],
        c=table.loc[ok, "conservation"],
        cmap=SEQUENTIAL_CMAP,
        s=9,
        vmin=0,
        vmax=1,
        zorder=2,
    )
    a1.scatter(proj[~ok, 0], proj[~ok, 1], color=LIGHT_GREY, s=6, zorder=2, label="not scored")
    cat = table["pdb_resnum"].isin(list(CATALYTIC_SITES)).to_numpy()
    a1.scatter(
        proj[cat, 0],
        proj[cat, 1],
        marker="*",
        s=220,
        color=CATALYTIC,
        edgecolor="white",
        zorder=3,
        label="catalytic D477 / E515 / D588",
    )
    fig.colorbar(sc, ax=a1, shrink=0.8).set_label("GH70 family conservation (0-1)")
    a1.set_aspect("equal")
    a1.set_xlabel("PC1 of C-alpha coordinates (A)")
    a1.set_ylabel("PC2 of C-alpha coordinates (A)")
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
    a1.set_title(f"{pid} chain A, C-alpha atoms in 2D (n={len(table)})")
    s = table.dropna(subset=["conservation"])
    d = s["distance_to_active_site"]
    a2.scatter(d, s["conservation"], s=6, alpha=0.4, color=GREY, label="residue")
    bins = np.arange(0, d.max() + 5, 5)
    means = [s.loc[(d >= lo) & (d < hi), "conservation"].mean() for lo, hi in pairwise(bins)]
    a2.plot(
        (bins[:-1] + bins[1:]) / 2,
        means,
        color="black",
        lw=LW_THICK,
        marker="o",
        ms=4,
        label="mean per 5 A bin",
    )
    a2.set_xlabel("C-alpha distance to catalytic center (A)")
    a2.set_ylabel("GH70 family conservation (0-1)")
    a2.legend(loc="lower left")
    a2.set_title(
        f"Spearman rho = {summary['spearman_rho_conservation_vs_distance']:.2f} "
        f"[{summary['rho_ci_low']:.2f}, {summary['rho_ci_high']:.2f}] "
        f"(n={summary['n_residues_scored']} residues)"
    )
    fig.suptitle("Conservation in 3D: are conserved residues near the active site?" + mode_tag())
    fig.tight_layout()
    save(fig, "m7_structure_conservation_map")


def plot_ramachandran(rama: pd.DataFrame, pid: str) -> None:
    """phi/psi scatter for general, glycine and proline residues."""
    fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 5.2), sharey=True)
    styles = {"general": DARK_GREY, "Gly": OKABE_ITO["orange"], "Pro": OKABE_ITO["green"]}
    for ax, (t, col) in zip(axes, styles.items(), strict=True):
        sub = rama[rama["type"] == t]
        ax.hist2d(
            rama["phi"],
            rama["psi"],
            bins=72,
            range=[[-180, 180], [-180, 180]],
            cmap="Greys",
            cmin=1,
            alpha=0.35,
        )
        ax.scatter(sub["phi"], sub["psi"], s=4, color=col, alpha=0.7, label=f"{t} residues")
        ax.axhline(0, color=GREY, lw=LW_HAIR)
        ax.axvline(0, color=GREY, lw=LW_HAIR)
        ax.set(
            xlim=(-180, 180),
            ylim=(-180, 180),
            xticks=range(-180, 181, 90),
            yticks=range(-180, 181, 90),
        )
        ax.set_xlabel("phi (degrees)")
        ax.set_aspect("equal")
        ax.set_title(f"{t} residues (n={len(sub)})")
    axes[0].set_ylabel("psi (degrees)")
    frac = rama["region"].value_counts(normalize=True) * 100
    txt = "; ".join(f"{k}: {v:.1f}%" for k, v in frac.items())
    fig.suptitle(f"Ramachandran plot, {pid} chain A (grey background = all residues)" + mode_tag())
    handles: list[Artist] = [plt.Rectangle((0, 0), 1, 1, color="#BBBBBB")]
    handles += [plt.Line2D([], [], marker="o", ls="", color=c) for c in styles.values()]
    fig.legend(
        handles,
        ["all residues (density)", *[f"{t} residues" for t in styles]],
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 0.0),
    )
    fig.text(
        0.5,
        -0.07,
        f"Broad regions (simple phi/psi boxes, not a validation score): {txt}",
        ha="center",
        fontsize=FS_ANNOT,
    )
    fig.tight_layout()
    save(fig, "m7_ramachandran")


def plot_all(
    table: pd.DataFrame, rama: pd.DataFrame, summary: dict, pdb_text: str, pid: str
) -> None:
    """Draw every M7 figure and the interactive HTML view."""
    plot_structure_map(table, summary, pid)
    plot_structure_3d(table, pid)
    plot_ramachandran(rama, pid)
    write_html(pdb_text, pid)
