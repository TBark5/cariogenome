"""M4: per-gene NJ and UPGMA trees with bootstrap support, a concatenated housekeeping
reference tree, and gene-tree vs reference-tree discordance.

Discordance is summarised by the normalised Robinson-Foulds distance and by the number of
gene-tree splits with bootstrap >= 70% that are incompatible with the reference tree.
Discordance can reflect horizontal transfer or recombination, but also incomplete lineage
sorting, paralogy, alignment error or simply too little signal. It is suggestive, not proof.
"""
from __future__ import annotations

from io import StringIO

import numpy as np
import pandas as pd
from Bio import Phylo

from . import phylo
from .config import RESULTS, all_genes, gene_class, housekeeping_genes, load_config, rng
from .seqio import read_fasta, species_map
from .stats import compare_table

TREE_DIR = RESULTS / "trees"
ALN_DIR = RESULTS / "alignments"
SUPPORTED = 70.0
FOCAL = "S. mutans"


def gene_alignment(gene: str) -> dict[str, str]:
    """Codon alignment for coding genes, nucleotide MSA for 16S."""
    name = f"{gene}.nt.fasta" if gene == "16S" else f"{gene}.codon.fasta"
    return read_fasta(ALN_DIR / name)


def concatenated_housekeeping() -> dict[str, str]:
    """Concatenate the housekeeping codon alignments; missing genes are filled with gaps."""
    parts = {g: gene_alignment(g) for g in housekeeping_genes()}
    taxa = list(dict.fromkeys(t for aln in parts.values() for t in aln))
    out = {t: "" for t in taxa}
    for aln in parts.values():
        L = len(next(iter(aln.values())))
        for t in taxa:
            out[t] += aln.get(t, "-" * L)
    (ALN_DIR / "housekeeping_concatenated.fasta").write_text(
        "".join(f">{k}\n{v}\n" for k, v in out.items()))
    return out


def save_tree(tree, name: str) -> None:
    TREE_DIR.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(TREE_DIR / f"{name}.nwk"), "newick")


def load_tree(name: str):
    return Phylo.read(str(TREE_DIR / f"{name}.nwk"), "newick")


def _smaller_side(split: frozenset, taxa: frozenset) -> frozenset:
    other = taxa - split
    return split if len(split) <= len(other) else other


def discordance(gene: str, gene_sup: dict, ref_sup: dict, sp: dict[str, str]) -> tuple[dict, list]:
    """Compare one gene tree (split -> support) with the reference tree."""
    taxa = frozenset(gene_alignment(gene))
    ref_all = frozenset().union(*ref_sup.keys())
    ref_r = phylo.restrict(ref_sup, taxa, ref_all)
    _, nrf = phylo.robinson_foulds(set(gene_sup), set(ref_r))
    conflicts = []
    for s, sup in gene_sup.items():
        bad = [r for r in ref_r if not phylo.compatible(s, r, taxa)]
        if bad and sup >= SUPPORTED:
            side = _smaller_side(s, taxa)
            # A conflicting clade that mixes species is the pattern relevant to HGT between
            # species; one that only regroups strains of one species is within-species.
            kind = "within-species" if len({sp[t] for t in side}) == 1 else "between-species"
            conflicts.append({"gene": gene, "gene_support": sup, "category": kind,
                              "max_conflicting_reference_support": max(ref_r[r] for r in bad),
                              "taxa_on_smaller_side": ";".join(sorted(side))})
    focal = frozenset(t for t in taxa if sp[t] == FOCAL)
    g_f = phylo.restrict(gene_sup, focal, taxa)
    r_f = phylo.restrict(ref_sup, focal, ref_all)
    _, nrf_focal = phylo.robinson_foulds(set(g_f), set(r_f))
    species_sets = {}
    for t in taxa:
        species_sets.setdefault(sp[t], set()).add(t)
    mono = [phylo.canonical(frozenset(v), taxa) in gene_sup or len(v) < 2 or len(taxa - v) < 2
            for v in species_sets.values()]
    sups = np.array(list(gene_sup.values())) if gene_sup else np.array([np.nan])
    row = {"gene": gene, "class": "rRNA control" if gene == "16S" else gene_class(gene),
           "n_taxa": len(taxa), "n_species": len(species_sets), "n_splits": len(gene_sup),
           "mean_support": float(np.nanmean(sups)),
           "n_unsupported_splits": int((sups < SUPPORTED).sum()),
           "nRF_vs_reference": nrf, "nRF_within_Smutans": nrf_focal,
           "n_supported_conflicts": len(conflicts),
           "n_supported_conflicts_between_species": sum(c["category"] == "between-species" for c in conflicts),
           "n_strong_conflicts": sum(c["max_conflicting_reference_support"] >= SUPPORTED for c in conflicts),
           "all_species_monophyletic": all(mono)}
    return row, conflicts


def run() -> dict[str, pd.DataFrame]:
    """Build every tree, compute discordance, and test virulence vs control."""
    from . import m4_figures
    cfg = load_config()["bootstrap"]
    n_rep = cfg["tree_replicates"]
    sp = species_map()
    supports: dict[str, dict] = {}
    rows_meth = []
    ref_msa = concatenated_housekeeping()
    for method in ("nj", "upgma"):
        tree, sup = phylo.run_bootstrap_tree(ref_msa, method, n_rep, rng(100))
        save_tree(tree, f"reference_housekeeping_{method}")
        supports[f"reference_{method}"] = sup
    for i, gene in enumerate(all_genes()):
        msa = gene_alignment(gene)
        for method in ("nj", "upgma"):
            tree, sup = phylo.run_bootstrap_tree(msa, method, n_rep, rng(200 + i))
            save_tree(tree, f"{gene}_{method}")
            supports[f"{gene}_{method}"] = sup
        _, nrf_m = phylo.robinson_foulds(set(supports[f"{gene}_nj"]), set(supports[f"{gene}_upgma"]))
        rows_meth.append({"gene": gene, "nRF_NJ_vs_UPGMA": nrf_m})
    rows, conflicts = [], []
    for gene in all_genes():
        r, c = discordance(gene, supports[f"{gene}_nj"], supports["reference_nj"], sp)
        rows.append(r)
        conflicts += c
    disc = pd.DataFrame(rows).merge(pd.DataFrame(rows_meth), on="gene")
    disc.to_csv(RESULTS / "m4_discordance.csv", index=False)
    pd.DataFrame(conflicts, columns=["gene", "gene_support", "category",
                                     "max_conflicting_reference_support",
                                     "taxa_on_smaller_side"]).sort_values(
        ["gene", "taxa_on_smaller_side"]).to_csv(RESULTS / "m4_conflicts.csv", index=False)
    ref_sup = pd.DataFrame([{"method": m, "split": ";".join(sorted(s)), "support": v}
                            for m in ("nj", "upgma") for s, v in supports[f"reference_{m}"].items()])
    ref_sup.sort_values(["method", "split"]).to_csv(RESULTS / "m4_reference_supports.csv", index=False)
    coding = disc[disc["class"] != "rRNA control"]
    tests = compare_table(coding, ["nRF_within_Smutans", "n_supported_conflicts",
                                   "n_supported_conflicts_between_species", "mean_support"],
                          rng(6), n_boot=cfg["stats_replicates"])
    tests.to_csv(RESULTS / "m4_tests.csv", index=False)
    m4_figures.plot_all(disc, tests)
    print(disc[["gene", "n_taxa", "mean_support", "n_unsupported_splits", "nRF_vs_reference",
                "nRF_within_Smutans", "n_supported_conflicts", "n_supported_conflicts_between_species",
                "all_species_monophyletic", "nRF_NJ_vs_UPGMA"]].round(3).to_string(index=False))
    print(tests[["metric", "median_virulence", "median_control", "cliffs_delta", "delta_ci_low",
                 "delta_ci_high", "q_bh"]].round(4).to_string(index=False))
    return {"discordance": disc, "tests": tests}


def tree_from_string(newick: str):
    """Parse a Newick string (used by tests)."""
    return Phylo.read(StringIO(newick), "newick")
