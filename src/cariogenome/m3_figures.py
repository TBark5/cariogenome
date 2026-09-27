"""Figures for M3 (alignment and conservation)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from .config import RESULTS, all_genes, genomes, virulence_genes
from .plotting import CLASS_COLORS, OKABE_ITO, SEQUENTIAL_CMAP, SPECIES_COLORS, save
from .seqio import mode_tag, read_fasta, species_map

AA_GROUPS = {**dict.fromkeys("AILMV", OKABE_ITO["sky"]), **dict.fromkeys("FWY", OKABE_ITO["blue"]),
             **dict.fromkeys("KRH", OKABE_ITO["vermillion"]), **dict.fromkeys("DE", OKABE_ITO["purple"]),
             **dict.fromkeys("STNQ", OKABE_ITO["green"]), **dict.fromkeys("GP", OKABE_ITO["orange"]),
             "C": OKABE_ITO["yellow"]}


def _order(labels: list[str]) -> list[str]:
    rank = {g.label: i for i, g in enumerate(genomes())}
    return sorted(labels, key=lambda k: rank.get(k, 999))


def plot_identity(summary: pd.DataFrame) -> None:
    genes = all_genes()
    fig, axes = plt.subplots(3, 5, figsize=(22, 13))
    fig.subplots_adjust(wspace=0.55, hspace=0.25)
    sp = species_map()
    for ax, gene in zip(axes.flat, genes):
        m = pd.read_csv(RESULTS / "identity" / f"{gene}.csv", index_col=0)
        order = _order(list(m.index))
        m = m.loc[order, order]
        im = ax.imshow(m.values, cmap=SEQUENTIAL_CMAP, vmin=50, vmax=100)
        ax.set_xticks([])
        ax.set_yticks(range(len(order)), order, fontsize=5.5)
        for lab in ax.get_yticklabels():
            lab.set_color(SPECIES_COLORS.get(sp[lab.get_text()], "black"))
        kind = "nt" if gene == "16S" else "aa"
        cls = "virulence" if gene in virulence_genes() else "control"
        ax.set_title(f"{gene} ({kind}, n={len(order)})", color=CLASS_COLORS[cls], fontsize=10)
        ax.spines[:].set_visible(False)
    cb = fig.colorbar(im, ax=axes, shrink=0.5, location="right")
    cb.set_label("Percent identity (aligned positions)")
    fig.suptitle("Pairwise percent identity per gene (rows ordered by species; red titles = "
                 "virulence-associated, blue = controls)" + mode_tag(), fontweight="bold", y=0.93)
    save(fig, "m3_identity_heatmaps")


def plot_profiles(summary: pd.DataFrame, cons_all: dict, regions: pd.DataFrame) -> None:
    genes = virulence_genes() + ["recA", "rpoB"]
    fig, axes = plt.subplots(len(genes), 1, figsize=(12, 2.0 * len(genes)))
    s = summary.set_index("gene")
    for ax, gene in zip(axes, genes):
        c = cons_all[gene].dropna(subset=["ref_position"])
        win = int(regions.loc[regions["gene"] == gene, "window"].iloc[0])
        ax.plot(c["ref_position"], c["conservation"], color="#CCCCCC", lw=0.5)
        sm = c.set_index("ref_position")["conservation"].rolling(win, center=True, min_periods=win // 2).mean()
        cls = "virulence" if gene in virulence_genes() else "control"
        ax.plot(sm.index, sm.values, color=CLASS_COLORS[cls], lw=1.5)
        for r in regions[regions["gene"] == gene].itertuples():
            col = OKABE_ITO["green"] if r.type == "most conserved" else OKABE_ITO["orange"]
            ax.axvspan(r.start, r.end, color=col, alpha=0.25, lw=0)
        ax.set_ylim(-0.02, 1.05)
        ax.set_xlim(0, c["ref_position"].max())
        ax.set_ylabel("conservation")
        ax.set_title(f"{gene}: n={s.at[gene, 'n_sequences']} sequences from "
                     f"{s.at[gene, 'n_species']} species", loc="left", fontsize=9.5,
                     color=CLASS_COLORS[cls])
    axes[-1].set_xlabel("Residue position in S. mutans UA159")
    handles = [plt.Rectangle((0, 0), 1, 1, color=OKABE_ITO["green"], alpha=0.4),
               plt.Rectangle((0, 0), 1, 1, color=OKABE_ITO["orange"], alpha=0.4)]
    axes[0].legend(handles, ["3 most conserved windows", "3 least conserved windows"],
                   loc="lower right", ncol=2, fontsize=8)
    fig.suptitle("Per-site conservation (1 - Shannon entropy / log2 20; Henikoff-weighted), "
                 "running mean over the region window (30 aa; 16 aa for luxS)" + mode_tag(), fontweight="bold", y=1.0)
    fig.tight_layout()
    save(fig, "m3_conservation_profiles")


def plot_alignment_residues(gene: str = "luxS", block: int = 80) -> None:
    msa = read_fasta(RESULTS / "alignments" / f"{gene}.aa.fasta")
    labels = _order(list(msa))
    L = len(msa[labels[0]])
    n_blocks = int(np.ceil(L / block))
    fig, axes = plt.subplots(n_blocks, 1, figsize=(16, 0.19 * len(labels) * n_blocks + 1.2))
    sp = species_map()
    for b, ax in enumerate(np.atleast_1d(axes)):
        lo, hi = b * block, min(L, (b + 1) * block)
        for i, lab in enumerate(labels):
            for j in range(lo, hi):
                ch = msa[lab][j]
                if ch != "-":
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                               color=AA_GROUPS.get(ch, "#DDDDDD"), alpha=0.75, lw=0))
                ax.text(j, i, ch, ha="center", va="center", fontsize=5.2, family="monospace")
        ax.set_xlim(lo - 0.5, lo + block - 0.5)
        ax.set_ylim(len(labels) - 0.5, -0.5)
        ax.set_yticks(range(len(labels)), labels, fontsize=6)
        for t in ax.get_yticklabels():
            t.set_color(SPECIES_COLORS.get(sp[t.get_text()], "black"))
        ticks = list(range(lo, hi, 10))
        ax.set_xticks(ticks, [str(t + 1) for t in ticks], fontsize=7)
        ax.spines[:].set_visible(False)
    names = ["hydrophobic AILMV", "aromatic FWY", "basic KRH", "acidic DE", "polar STNQ", "G/P", "C"]
    cols = [OKABE_ITO[k] for k in ("sky", "blue", "vermillion", "purple", "green", "orange", "yellow")]
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.75) for c in cols], names,
               loc="lower center", ncol=7, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"{gene} protein alignment (center-star progressive; n={len(labels)} sequences, "
                 f"{len({sp[k] for k in labels})} species)" + mode_tag(), fontweight="bold")
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    save(fig, f"m3_alignment_{gene}")


def plot_alignment_overview(gene: str = "gtfD") -> None:
    msa = read_fasta(RESULTS / "alignments" / f"{gene}.aa.fasta")
    labels = _order(list(msa))
    ref = msa["Smu_UA159"]
    mat = np.array([[0 if ch == "-" else (2 if ch == ref[j] else 1) for j, ch in enumerate(msa[k])]
                    for k in labels])
    fig, ax = plt.subplots(figsize=(14, 0.25 * len(labels) + 1.5))
    cmap = ListedColormap(["white", OKABE_ITO["vermillion"], "#9ECAE1"])
    ax.imshow(mat, cmap=cmap, aspect="auto", vmin=0, vmax=2, interpolation="nearest")
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    sp = species_map()
    for t in ax.get_yticklabels():
        t.set_color(SPECIES_COLORS.get(sp[t.get_text()], "black"))
    ax.set_xlabel("Alignment column")
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=c, ec="grey") for c in cmap.colors[::-1]],
              ["same residue as UA159", "different residue", "gap"], loc="upper left",
              bbox_to_anchor=(1.01, 1))
    ax.set_title(f"{gene} alignment overview (n={len(labels)} sequences) relative to "
                 "S. mutans UA159" + mode_tag())
    save(fig, f"m3_alignment_overview_{gene}")


def plot_summary(summary: pd.DataFrame, tests: pd.DataFrame) -> None:
    coding = summary[summary["class"] != "rRNA control"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    t = tests.set_index("metric")
    for ax, (m, lab) in zip(axes, [("mean_pid_Smutans", "Mean pairwise aa identity (%)"),
                                   ("mean_entropy_Smutans", "Mean per-site entropy (bits)")]):
        for i, cls in enumerate(["virulence", "control"]):
            sub = coding[coding["class"] == cls]
            x = i + np.linspace(-0.15, 0.15, len(sub))
            ax.scatter(x, sub[m], color=CLASS_COLORS[cls], s=32, zorder=3)
            for xi, v, g in zip(x, sub[m], sub["gene"]):
                ax.annotate(g, (xi, v), xytext=(4, 0), textcoords="offset points", fontsize=7)
            ax.hlines(np.median(sub[m]), i - 0.3, i + 0.3, color="black")
        r = t.loc[m]
        ax.set_xticks([0, 1], [f"virulence (n={int(r.n_virulence)})", f"housekeeping (n={int(r.n_control)})"])
        ax.set_xlim(-0.6, 1.8)
        ax.set_ylabel(lab)
        ax.set_title(f"Cliff's d = {r.cliffs_delta:.2f} [{r.delta_ci_low:.2f}, {r.delta_ci_high:.2f}]"
                     f", q = {r.q_bh:.3f}", fontsize=9.5)
    n = coding["n_Smutans"]
    fig.suptitle(f"Conservation within S. mutans ({n.min()}-{n.max()} strains per gene); "
                 "black line = median" + mode_tag(), fontweight="bold", y=1.02)
    save(fig, "m3_conservation_summary")


def plot_all(summary, cons_all, regions, tests) -> None:
    """Draw every M3 figure."""
    plot_identity(summary)
    plot_profiles(summary, cons_all, regions)
    plot_alignment_residues("luxS")
    plot_alignment_overview("gtfD")
    plot_summary(summary, tests)
