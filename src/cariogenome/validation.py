"""Synthetic recovery test: does the pipeline recover a known tree and known omega values?

Codon sequences are simulated along the known tree in ``synthetic.py`` (no indels, so the
simulated sequences are already aligned). We then run the same M4 (K2P + NJ) and M5
(pooled NG86) code used on the real data and compare with the truth.
"""
from __future__ import annotations

import io
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from scipy.stats import spearmanr

from . import phylo
from .config import FIGURES, RESULTS, genomes, load_config
from .dnds import PooledNG86
from .plotting import OKABE_ITO, save
from .synthetic import simulate_codons, true_tree_newick

TEST_OMEGAS = [0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0]


def true_splits() -> set[frozenset]:
    return phylo.splits(Phylo.read(io.StringIO(true_tree_newick()), "newick"))


def recover_topology(n_genes: int, n_codons: int, gen: np.random.Generator) -> dict:
    """Simulate housekeeping-like genes, concatenate, build an NJ tree, compare with truth."""
    newick = true_tree_newick()
    genes = [simulate_codons(newick, n_codons, 0.1, gen) for _ in range(n_genes)]
    concat = {k: "".join(g[k] for g in genes) for k in genes[0]}
    tree, sup = phylo.run_bootstrap_tree(concat, "nj", 50, gen)
    rf, nrf = phylo.robinson_foulds(set(sup), true_splits())
    return {"n_genes": n_genes, "n_codons_per_gene": n_codons, "rf": rf, "nRF": nrf,
            "mean_support": float(np.mean(list(sup.values())))}


def recover_omega(n_codons: int, gen: np.random.Generator, reps: int = 3) -> pd.DataFrame:
    """Estimate omega for genes simulated with known omega (pairs within the S. mutans clade)."""
    smu = [g.label for g in genomes() if g.species == "S. mutans"]
    rows = []
    for om in TEST_OMEGAS:
        for r in range(reps):
            seqs = simulate_codons(true_tree_newick(), n_codons, om, gen, scale=3.0)
            est = PooledNG86(seqs, list(combinations(smu, 2))).estimate()
            rows.append({"true_omega": om, "replicate": r, "estimated_omega": est["omega"],
                         "pS": est["pS"]})
    return pd.DataFrame(rows)


def run(save_outputs: bool = True) -> dict:
    """Run both recovery checks and write results/validation_synthetic_*.csv + a figure."""
    gen = np.random.default_rng(load_config()["seed"] + 42)
    topo = pd.DataFrame([recover_topology(8, 400, gen)])
    omega = recover_omega(500, gen)
    rho = spearmanr(omega["true_omega"], omega["estimated_omega"]).statistic
    rel_err = float(np.median(np.abs(omega["estimated_omega"] / omega["true_omega"] - 1)))
    summary = {"topology_nRF": float(topo["nRF"].iloc[0]), "topology_mean_support": float(topo["mean_support"].iloc[0]),
               "omega_spearman": float(rho), "omega_median_relative_error": rel_err}
    if save_outputs:
        topo.to_csv(RESULTS / "validation_synthetic_topology.csv", index=False)
        omega.to_csv(RESULTS / "validation_synthetic_omega.csv", index=False)
        pd.DataFrame([summary]).to_csv(RESULTS / "validation_synthetic_summary.csv", index=False)
        fig, ax = plt.subplots(figsize=(5.5, 5))
        ax.plot([0, 1.1], [0, 1.1], color="grey", ls="--", lw=0.8, label="perfect recovery")
        ax.scatter(omega["true_omega"], omega["estimated_omega"], color=OKABE_ITO["blue"], s=28,
                   label=f"estimate ({len(omega)} simulated genes)")
        ax.set_xlabel("true omega (simulation)")
        ax.set_ylabel("estimated omega (pooled NG86)")
        ax.set_title(f"SYNTHETIC validation: omega recovery (Spearman {rho:.2f})\n"
                     f"NJ topology vs true tree: normalised RF = {summary['topology_nRF']:.2f}",
                     fontsize=10)
        ax.legend(fontsize=8)
        save(fig, "validation_synthetic_recovery")
    print("  synthetic validation:", {k: round(v, 3) for k, v in summary.items()})
    return summary
