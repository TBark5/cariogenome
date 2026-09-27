"""SYNTHETIC fallback data: codon sequences simulated along a known tree with known omega.

Used when NCBI is unreachable (so the pipeline never blocks) and by the test suite to
check that the pipeline recovers a known topology and known selection pressures.

Simulation model: mutations are proposed at rate ``t`` per nucleotide along each branch
(uniform over the three alternative bases). Mutations creating stop codons are rejected,
synonymous mutations are accepted, and nonsynonymous mutations are accepted with
probability omega. The expected dN/dS of the output is therefore close to omega.
"""
from __future__ import annotations

import gzip
import io
import json

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Seq import Seq

from .config import DATA, GENES_DIR, GENOMES_DIR, coding_genes, genomes, load_config, virulence_genes

CODE = {a + b + c: str(Seq(a + b + c).translate(table=11))
        for a in "ACGT" for b in "ACGT" for c in "ACGT"}
SENSE = [c for c, aa in CODE.items() if aa != "*"]
TRUE_OMEGA = {"gtfB": 0.45, "gtfC": 0.40, "gtfD": 0.35, "spaP": 0.90, "ftf": 0.60,
              "luxS": 0.15, "recA": 0.06, "rpoB": 0.04, "gyrB": 0.05, "gyrA": 0.07,
              "sodA": 0.08, "pheS": 0.10, "atpD": 0.05, "tuf": 0.03}
TRUTH_FILE = DATA / "synthetic_truth.json"


def _clade(labels: list[str], bl: float) -> str:
    """Balanced, fully resolved Newick subtree over ``labels``."""
    if len(labels) == 1:
        return f"{labels[0]}:{bl}"
    mid = len(labels) // 2
    return f"({_clade(labels[:mid], bl)},{_clade(labels[mid:], bl)}):{bl}"


def true_tree_newick() -> str:
    """Known species/strain tree over the configured genome labels."""
    by_sp: dict[str, list[str]] = {}
    for g in genomes():
        by_sp.setdefault(g.species, []).append(g.label)
    c = {sp: _clade(lbls, 0.02) for sp, lbls in by_sp.items()}
    return (f"(({c['S. mutans']}:0.25,{c['S. salivarius']}:0.25):0.1,"
            f"(({c['S. sanguinis']}:0.12,{c['S. gordonii']}:0.12):0.1,{c['S. mitis']}:0.2):0.1);")


def simulate_codons(newick: str, n_codons: int, omega: float, rng: np.random.Generator,
                    scale: float = 1.0) -> dict[str, str]:
    """Evolve a random coding sequence (ATG ... TAA) down the tree; return leaf sequences."""
    tree = Phylo.read(io.StringIO(newick), "newick")
    root = [SENSE[i] for i in rng.integers(0, len(SENSE), n_codons)]
    out: dict[str, str] = {}

    def evolve(clade, seq: list[str]) -> None:
        seq = list(seq)
        t = (clade.branch_length or 0.0) * scale
        for _ in range(rng.poisson(t * 3 * n_codons)):
            ci, pos = int(rng.integers(n_codons)), int(rng.integers(3))
            old = seq[ci]
            new = old[:pos] + rng.choice([b for b in "ACGT" if b != old[pos]]) + old[pos + 1:]
            if CODE[new] == "*":
                continue
            if CODE[new] == CODE[old] or rng.random() < omega:
                seq[ci] = new
        if clade.is_terminal():
            out[clade.name] = "ATG" + "".join(seq) + "TAA"
        for child in clade.clades:
            evolve(child, seq)

    evolve(tree.root, root)
    return out


def simulate_nucleotides(newick: str, length: int, rng: np.random.Generator,
                         scale: float = 0.3) -> dict[str, str]:
    """Jukes-Cantor evolution of a non-coding sequence (16S stand-in)."""
    tree = Phylo.read(io.StringIO(newick), "newick")
    root = "".join(rng.choice(list("ACGT"), length))
    out: dict[str, str] = {}

    def evolve(clade, seq: str) -> None:
        s = list(seq)
        for _ in range(rng.poisson((clade.branch_length or 0.0) * scale * length)):
            i = int(rng.integers(length))
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
        if clade.is_terminal():
            out[clade.name] = "".join(s)
        for child in clade.clades:
            evolve(child, "".join(s))

    evolve(tree.root, root)
    return out


def generate_dataset(seed_offset: int = 7) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Write a complete SYNTHETIC gene panel to ``data/`` and return (hits, gtf family)."""
    rng = np.random.default_rng(load_config()["seed"] + seed_offset)
    newick = true_tree_newick()
    gl = genomes()
    GENES_DIR.mkdir(parents=True, exist_ok=True)
    GENOMES_DIR.mkdir(parents=True, exist_ok=True)
    rows, protein = [], {}
    panel = {g: simulate_codons(newick, int(rng.integers(300, 600)), TRUE_OMEGA[g], rng)
             for g in coding_genes()}
    panel["16S"] = simulate_nucleotides(newick, 1500, rng)
    for gene, seqs in panel.items():
        with open(GENES_DIR / f"{gene}.fna", "w") as fna, open(GENES_DIR / f"{gene}.faa", "w") as faa:
            for g in gl:
                nt = seqs[g.label]
                head = f">{g.label} gene={gene} SYNTHETIC"
                fna.write(f"{head}\n{nt}\n")
                if gene != "16S":
                    aa = str(Seq(nt).translate(table=11)).rstrip("*")
                    faa.write(f"{head}\n{aa}\n")
                    protein[(gene, g.label)] = aa
                rows.append({"gene": gene, "label": g.label, "species": g.species,
                             "accession": "SYNTHETIC", "row_index": -1,
                             "locus_tag": f"SYN_{gene}", "old_locus_tag": "", "protein_id": "",
                             "annotation_gene": gene, "product": "SYNTHETIC", "start": 1,
                             "end": len(nt), "strand": "+", "pseudo": False,
                             "identity_to_query": 1.0, "query_coverage": 1.0, "hit_coverage": 1.0,
                             "reciprocal_best_hit": True, "ortholog_call": True})
    for g in gl:
        ribo = [simulate_codons("(a:0.01,b:0.01);", 150, 0.05, rng)["a"] for _ in range(20)]
        (GENOMES_DIR / f"{g.label}_ribosomal.fna").write_text(
            "".join(f">SYN_ribo_{i} SYNTHETIC ribosomal protein\n{s}\n" for i, s in enumerate(ribo)))
        gz = gzip.GzipFile(GENOMES_DIR / f"{g.label}_cds.fna.gz", "wb", mtime=0)  # no timestamp
        with io.TextIOWrapper(gz, newline="\n") as fh:
            for i in range(200):
                s = simulate_codons("(a:0.01,b:0.01);", int(rng.integers(100, 500)), 0.2, rng)["a"]
                fh.write(f">SYN_cds_{i} SYNTHETIC\n{s}\n")
    fam_rows = []
    with open(GENES_DIR / "gtf_family.faa", "w") as fh:
        for gene in ("gtfB", "gtfC", "gtfD"):
            for g in gl[::3]:
                sid = f"{g.label}|SYN_{gene}"
                fh.write(f">{sid} SYNTHETIC\n{protein[(gene, g.label)]}\n")
                fam_rows.append({"id": sid, "label": g.label, "species": g.species,
                                 "locus_tag": f"SYN_{gene}", "old_locus_tag": "", "protein_id": "",
                                 "product": "SYNTHETIC", "length": len(protein[(gene, g.label)]),
                                 "best_identity_to_Smu_gtf": 1.0})
    TRUTH_FILE.write_text(json.dumps({"newick": newick, "omega": TRUE_OMEGA,
                                      "virulence": virulence_genes()}, indent=2))
    return pd.DataFrame(rows), pd.DataFrame(fam_rows)
