"""Figures for M5 (selection)."""

from __future__ import annotations

from collections.abc import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from .config import coding_genes, genomes, load_config, virulence_genes
from .plotting import (
    CLASS_COLORS,
    DARK_GREY,
    FS_ANNOT,
    FS_TICK,
    GREY,
    LW,
    LW_HAIR,
    LW_THIN,
    OKABE_ITO,
    W_FULL,
    annotate_class_n,
    delta_stats,
    save,
    species_color,
    title_with_stats,
)
from .seqio import mode_tag


def _color_gene_labels(labels: Iterable[Text]) -> None:
    for lab in labels:
        g = lab.get_text().split(" ")[0].split("\n")[0]
        lab.set_color(CLASS_COLORS["virulence" if g in virulence_genes() else "control"])


def plot_dnds_bar(smu: pd.DataFrame, tests: pd.DataFrame, per_gene: pd.DataFrame) -> None:
    """Per-gene omega with CIs, plus the class comparison."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(W_FULL, 5.4), gridspec_kw={"width_ratios": [3, 1]})
    x = np.arange(len(smu))
    om = smu["omega"].to_numpy()
    lo = np.clip(om - smu["omega_ci_low"].to_numpy(), 0, None)
    hi = np.clip(smu["omega_ci_high"].to_numpy() - om, 0, None)
    ax.bar(
        x,
        om,
        color=[CLASS_COLORS[c] for c in smu["class"]],
        yerr=[lo, hi],
        capsize=3,
        error_kw={"lw": LW_THIN, "ecolor": DARK_GREY},
    )
    neutral = ax.axhline(1, color="black", ls="--", lw=LW_THIN)
    ctl_med = float(np.median(smu.loc[smu["class"] == "control", "omega"]))
    median = ax.axhline(ctl_med, color=CLASS_COLORS["control"], ls=":", lw=LW_THIN)
    q = per_gene.set_index("gene")["q_bh"]
    for xi, g, v, h in zip(x, smu["gene"], om, hi, strict=True):
        if g in q.index and q[g] < 0.05:
            ax.text(xi, v + h + 0.01, "*", ha="center", fontsize=FS_ANNOT * 1.6)
    ax.set_xticks(
        x,
        [f"{g}\n(n={n})" for g, n in zip(smu["gene"], smu["n_sequences"], strict=True)],
        fontsize=FS_TICK,
    )
    _color_gene_labels(ax.get_xticklabels())
    ax.set_xlabel("Gene (n = strains)")
    ax.set_ylabel("dN/dS (omega), pooled over strain pairs")
    ax.set_title(
        "dN/dS within S. mutans with 95% codon-bootstrap CI (1000 replicates)\n"
        "* = differs from pooled housekeeping genes, BH q < 0.05"
    )
    ax.set_ylim(0, max(1.1, float(np.nanmax(om + hi)) * 1.1))
    ax.legend(
        [
            *[Rectangle((0, 0), 1, 1, color=CLASS_COLORS[c]) for c in ("virulence", "control")],
            neutral,
            median,
        ],
        [
            "virulence-associated",
            "housekeeping",
            "omega = 1 (neutral)",
            f"housekeeping median ({ctl_med:.3f})",
        ],
        loc="upper left",
        bbox_to_anchor=(0.02, 0.9),
    )
    t = tests.set_index("metric").loc["omega"]
    for i, cls in enumerate(["virulence", "control"]):
        v = smu.loc[smu["class"] == cls, "omega"]
        ax2.scatter(
            i + np.linspace(-0.15, 0.15, len(v)), v, color=CLASS_COLORS[cls], s=30, zorder=3
        )
        ax2.hlines(np.median(v), i - 0.3, i + 0.3, color="black", label="median" if i else None)
    ax2.set_xticks(
        [0, 1],
        [
            f"virulence\n(n={int(t.n_virulence)} genes)",
            f"housekeeping\n(n={int(t.n_control)} genes)",
        ],
    )
    ax2.set_xlim(-0.6, 1.6)
    ax2.set_xlabel("Gene class")
    ax2.set_ylabel("dN/dS (omega)")
    ax2.legend(loc="upper right")
    title_with_stats(ax2, "Class comparison", delta_stats(t).replace(" [", "\n[", 1))
    fig.suptitle(
        "Selection pressure: virulence-associated vs housekeeping genes "
        "(omega < 1 = purifying selection)" + mode_tag()
    )
    fig.tight_layout()
    annotate_class_n(fig)
    save(fig, "m5_dnds_comparison")


def plot_species(within: pd.DataFrame) -> None:
    """Within-species omega of every gene in every species."""
    species = list(dict.fromkeys(g.species for g in genomes()))
    fig, ax = plt.subplots(figsize=(W_FULL, 5.4))
    genes = coding_genes()
    width = 0.8 / len(species)
    for j, sp in enumerate(species):
        sub = within[within["species"] == sp].set_index("gene").reindex(genes)
        x = np.arange(len(genes)) + (j - len(species) / 2 + 0.5) * width
        om = sub["omega"].to_numpy()
        err = [
            np.clip(om - sub["omega_ci_low"].to_numpy(), 0, None),
            np.clip(sub["omega_ci_high"].to_numpy() - om, 0, None),
        ]
        n = int(sub["n_sequences"].max()) if sub["n_sequences"].notna().any() else 0
        ax.errorbar(
            x,
            om,
            yerr=err,
            fmt="o",
            ms=4,
            color=species_color(sp),
            lw=LW_THIN,
            capsize=1.5,
            label=f"{sp} (n={n} strains)",
        )
    ax.axhline(1, color="black", ls="--", lw=LW_THIN)
    ax.axvline(len(virulence_genes()) - 0.5, color=GREY, lw=LW_THIN)
    ax.text(
        len(virulence_genes()) - 0.6,
        -0.035,
        "virulence-associated genes",
        ha="right",
        va="center",
        color=CLASS_COLORS["virulence"],
        fontsize=FS_ANNOT,
    )
    ax.text(
        len(virulence_genes()) - 0.4,
        -0.035,
        "housekeeping genes",
        ha="left",
        va="center",
        color=CLASS_COLORS["control"],
        fontsize=FS_ANNOT,
    )
    ax.set_xticks(range(len(genes)), genes)
    _color_gene_labels(ax.get_xticklabels())
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Within-species dN/dS (omega, symlog scale)")
    ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    ax.set_title(
        "Within-species dN/dS in every species (missing point = gene absent, fewer than "
        "3 strains, or no synonymous differences)" + mode_tag()
    )
    annotate_class_n(fig, y=-0.12)
    save(fig, "m5_dnds_by_species")


def plot_windows(windows: pd.DataFrame) -> None:
    """Sliding-window pN, pS and omega within S. mutans."""
    genes = [*virulence_genes(), "recA", "rpoB"]
    fig, axes = plt.subplots(len(genes), 1, figsize=(W_FULL, 2.0 * len(genes)))
    cfg = load_config()["selection"]
    for ax, gene in zip(axes, genes, strict=True):
        w = windows[windows["gene"] == gene]
        cls = "virulence" if gene in virulence_genes() else "control"
        ax.plot(w["mid_codon"], w["pS"], color=GREY, lw=LW_THIN, label="pS (synonymous)")
        ax.plot(
            w["mid_codon"],
            w["pN"],
            color=CLASS_COLORS[cls],
            lw=LW,
            label="pN (nonsynonymous, class color)",
        )
        ax.set_ylabel("p-distance\n(per site)")
        ax2 = ax.twinx()
        ok = w["omega"].notna()
        ax2.scatter(
            w.loc[ok, "mid_codon"],
            w.loc[ok, "omega"],
            s=10,
            color="black",
            zorder=3,
            label="omega (where pS > 0)",
        )
        ax2.axhline(1, color="black", ls="--", lw=LW_HAIR)
        ax2.set_ylabel("omega", fontsize=FS_ANNOT)
        ax2.spines["right"].set_visible(True)
        pos = w[(w["omega_ci_low"] > 1)]
        for r in pos.itertuples():
            ax.axvspan(r.start_codon, r.end_codon, color=OKABE_ITO["orange"], alpha=0.3)
        ax.set_title(
            f"{gene}: {len(w)} windows of {cfg['window_codons']} codons, "
            f"{int(ok.sum())} with defined omega, {len(pos)} with CI > 1",
            loc="left",
            color=CLASS_COLORS[cls],
        )
    h1, l1 = axes[0].get_legend_handles_labels()
    h2 = [plt.Line2D([], [], marker="o", ls="", color="black")]
    fig.legend(
        [*h1, *h2],
        [*l1, "omega (right axis, where pS > 0)"],
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.01),
    )
    axes[-1].set_xlabel("Codon position in the alignment (window midpoint)")
    fig.suptitle(
        "Sliding-window polymorphism within S. mutans (orange shading would mark a "
        "window with omega CI entirely above 1)" + mode_tag()
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.975))
    save(fig, "m5_sliding_window")


def plot_saturation(between: pd.DataFrame) -> None:
    """Between-species synonymous distance per gene and comparison."""
    sat = load_config()["selection"]["saturation_ps"]
    fig, ax = plt.subplots(figsize=(W_FULL, 5.0))
    comps = sorted(between["comparison"].unique())
    width = 0.8 / len(comps)
    genes = coding_genes()
    for j, c in enumerate(comps):
        sub = between[between["comparison"] == c].set_index("gene").reindex(genes)
        x = np.arange(len(genes)) + (j - len(comps) / 2 + 0.5) * width
        ax.bar(x, sub["pS"], width, color=species_color(c.replace("S. mutans vs ", "")), label=c)
    ax.axhline(sat, color="black", ls="--", lw=LW_THIN, label=f"saturation threshold pS = {sat}")
    ax.axhline(0.75, color=DARK_GREY, ls=":", lw=LW, label="Jukes-Cantor limit pS = 0.75")
    ax.set_xticks(range(len(genes)), genes)
    _color_gene_labels(ax.get_xticklabels())
    ax.set_xlabel("Gene (no bar = no ortholog in that species)")
    ax.set_ylabel("Synonymous p-distance, pS (per site)")
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    ax.set_title(
        f"Why dN/dS is estimated within species: {int(between['saturated'].sum())} of "
        f"{len(between)} between-species comparisons have saturated synonymous sites" + mode_tag()
    )
    save(fig, "m5_synonymous_saturation")


def plot_all(
    smu: pd.DataFrame,
    within: pd.DataFrame,
    between: pd.DataFrame,
    tests: pd.DataFrame,
    per_gene: pd.DataFrame,
    windows: pd.DataFrame,
) -> None:
    """Draw every M5 figure."""
    plot_dnds_bar(smu, tests, per_gene)
    plot_species(within)
    plot_windows(windows)
    plot_saturation(between)
