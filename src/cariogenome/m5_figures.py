"""Figures for M5 (selection)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .config import coding_genes, genomes, load_config, virulence_genes
from .plotting import CLASS_COLORS, OKABE_ITO, SPECIES_COLORS, save
from .seqio import mode_tag


def plot_dnds_bar(smu: pd.DataFrame, tests: pd.DataFrame, per_gene: pd.DataFrame) -> None:
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 5.2), gridspec_kw={"width_ratios": [3, 1]})
    x = np.arange(len(smu))
    cols = [CLASS_COLORS[c] for c in smu["class"]]
    om = smu["omega"].values
    lo = np.clip(om - smu["omega_ci_low"].values, 0, None)
    hi = np.clip(smu["omega_ci_high"].values - om, 0, None)
    ax.bar(x, om, color=cols, yerr=[lo, hi], capsize=3, error_kw={"lw": 1, "ecolor": "#333333"})
    ax.axhline(1, color="black", ls="--", lw=0.8)
    ax.text(len(smu) - 0.5, 1.02, "omega = 1 (neutral)", ha="right", va="bottom", fontsize=8)
    ctl_med = np.median(smu.loc[smu["class"] == "control", "omega"])
    ax.axhline(ctl_med, color=CLASS_COLORS["control"], ls=":", lw=1)
    ax.text(-0.4, ctl_med, f"control median {ctl_med:.3f}", color=CLASS_COLORS["control"],
            fontsize=8, va="bottom")
    q = per_gene.set_index("gene")["q_bh"]
    for xi, g, v, h in zip(x, smu["gene"], om, hi):
        if g in q.index and q[g] < 0.05:
            ax.text(xi, v + h + 0.01, "*", ha="center", fontsize=14)
    ax.set_xticks(x, [f"{g}\n(n={n})" for g, n in zip(smu["gene"], smu["n_sequences"])], fontsize=8)
    for lab, c in zip(ax.get_xticklabels(), smu["class"]):
        lab.set_color(CLASS_COLORS[c])
    ax.set_ylabel("dN/dS (omega), pooled over strain pairs")
    ax.set_title("dN/dS within S. mutans, 95% codon-bootstrap CI (1000 replicates);\n"
                 "* = differs from pooled controls, BH q < 0.05", fontsize=10)
    ax.set_ylim(0, max(1.1, np.nanmax(om + hi) * 1.1))
    t = tests.set_index("metric").loc["omega"]
    for i, cls in enumerate(["virulence", "control"]):
        v = smu.loc[smu["class"] == cls, "omega"]
        xs = i + np.linspace(-0.15, 0.15, len(v))
        ax2.scatter(xs, v, color=CLASS_COLORS[cls], s=30, zorder=3)
        ax2.hlines(np.median(v), i - 0.3, i + 0.3, color="black")
    ax2.set_xticks([0, 1], [f"virulence\n(n={int(t.n_virulence)})", f"control\n(n={int(t.n_control)})"])
    ax2.set_xlim(-0.6, 1.6)
    ax2.set_ylabel("omega")
    ax2.set_title(f"Cliff's d = {t.cliffs_delta:.2f}\n[{t.delta_ci_low:.2f}, {t.delta_ci_high:.2f}], "
                  f"q = {t.q_bh:.3f}", fontsize=10)
    fig.suptitle("Selection pressure: virulence-associated vs housekeeping genes (omega < 1 = "
                 "purifying selection)" + mode_tag(), fontweight="bold")
    fig.tight_layout()
    save(fig, "m5_dnds_comparison")


def plot_species(within: pd.DataFrame) -> None:
    species = list(dict.fromkeys(g.species for g in genomes()))
    fig, ax = plt.subplots(figsize=(13, 5))
    genes = coding_genes()
    width = 0.8 / len(species)
    for j, sp in enumerate(species):
        sub = within[within["species"] == sp].set_index("gene").reindex(genes)
        x = np.arange(len(genes)) + (j - len(species) / 2 + 0.5) * width
        om = sub["omega"].values
        err = [np.clip(om - sub["omega_ci_low"].values, 0, None), np.clip(sub["omega_ci_high"].values - om, 0, None)]
        ax.errorbar(x, om, yerr=err, fmt="o", ms=4, color=SPECIES_COLORS[sp], lw=0.8,
                    capsize=1.5, label=f"{sp} (n={int(sub['n_sequences'].max()) if sub['n_sequences'].notna().any() else 0} strains)")
    ax.axhline(1, color="black", ls="--", lw=0.8)
    ax.axvline(len(virulence_genes()) - 0.5, color="grey", lw=0.8)
    ax.set_xticks(range(len(genes)), genes)
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_ylabel("within-species omega (symlog scale)")
    ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=8)
    ax.set_title("Within-species dN/dS in every species (missing point = gene absent, fewer than 3 "
                 "strains, or no synonymous differences)" + mode_tag(), fontsize=10)
    save(fig, "m5_dnds_by_species")


def plot_windows(windows: pd.DataFrame) -> None:
    genes = virulence_genes() + ["recA", "rpoB"]
    fig, axes = plt.subplots(len(genes), 1, figsize=(12, 1.9 * len(genes)))
    cfg = load_config()["selection"]
    for ax, gene in zip(axes, genes):
        w = windows[windows["gene"] == gene]
        cls = "virulence" if gene in virulence_genes() else "control"
        ax.plot(w["mid_codon"], w["pS"], color=OKABE_ITO["grey"], lw=1, label="pS (synonymous)")
        ax.plot(w["mid_codon"], w["pN"], color=CLASS_COLORS[cls], lw=1.4, label="pN (nonsynonymous)")
        ax.set_ylabel("p-distance")
        ax2 = ax.twinx()
        ok = w["omega"].notna()
        ax2.scatter(w.loc[ok, "mid_codon"], w.loc[ok, "omega"], s=10, color="black", zorder=3,
                    label="omega (where pS > 0)")
        ax2.axhline(1, color="black", ls="--", lw=0.6)
        ax2.set_ylabel("omega", fontsize=8)
        pos = w[(w["omega_ci_low"] > 1)]
        for r in pos.itertuples():
            ax.axvspan(r.start_codon, r.end_codon, color=OKABE_ITO["orange"], alpha=0.3)
        ax.set_title(f"{gene}: {len(w)} windows of {cfg['window_codons']} codons, "
                     f"{int(ok.sum())} with defined omega, {len(pos)} with CI > 1",
                     loc="left", fontsize=9, color=CLASS_COLORS[cls])
    axes[0].legend(loc="upper left", fontsize=7, ncol=2)
    axes[-1].set_xlabel("codon position (alignment)")
    fig.suptitle("Sliding-window polymorphism within S. mutans (orange = window with omega CI "
                 "entirely above 1)" + mode_tag(), fontweight="bold")
    fig.tight_layout()
    save(fig, "m5_sliding_window")


def plot_saturation(between: pd.DataFrame) -> None:
    sat = load_config()["selection"]["saturation_ps"]
    fig, ax = plt.subplots(figsize=(12, 4.5))
    comps = sorted(between["comparison"].unique())
    width = 0.8 / len(comps)
    genes = coding_genes()
    for j, c in enumerate(comps):
        sub = between[between["comparison"] == c].set_index("gene").reindex(genes)
        x = np.arange(len(genes)) + (j - len(comps) / 2 + 0.5) * width
        sp = c.replace("S. mutans vs ", "")
        ax.bar(x, sub["pS"], width, color=SPECIES_COLORS.get(sp, "grey"), label=c)
    ax.axhline(sat, color="black", ls="--", lw=0.8)
    ax.axhline(0.75, color="red", ls=":", lw=0.8)
    ax.text(len(genes) - 0.5, sat, f"saturation threshold pS = {sat}", ha="right", va="bottom", fontsize=8)
    ax.text(len(genes) - 0.5, 0.75, "JC limit 0.75", ha="right", va="bottom", fontsize=8, color="red")
    ax.set_xticks(range(len(genes)), genes)
    ax.set_ylabel("synonymous p-distance (pS)")
    ax.legend(fontsize=8, ncol=2)
    ax.set_title(f"Why dN/dS is estimated within species: {int(between['saturated'].sum())} of "
                 f"{len(between)} between-species comparisons have saturated synonymous sites"
                 + mode_tag(), fontsize=10)
    save(fig, "m5_synonymous_saturation")


def plot_all(smu, within, between, tests, per_gene, windows) -> None:
    """Draw every M5 figure."""
    plot_dnds_bar(smu, tests, per_gene)
    plot_species(within)
    plot_windows(windows)
    plot_saturation(between)
