"""M2: compositional analysis (GC, GC3, GC skew, CAI, codon usage, amino-acid composition).

Primary test: 6 virulence vs 8 housekeeping genes in *S. mutans* (gene-level means over
strains). The CAI reference set is the ribosomal-protein genes of the same genome.
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
from Bio import SeqIO

from .codon import CODE, SENSE_CODONS, SYNONYMS, cai, cai_weights, codon_counts, rscu
from .config import GENOMES_DIR, RESULTS, coding_genes, all_genes, gene_class, genomes, load_config, rng
from .seqio import load_gene, species_map
from .stats import benjamini_hochberg, compare_table

METRICS = ["gc", "gc3", "gc_skew", "cai"]
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
FOCAL = "S. mutans"


def gc_skew(seq: str) -> float:
    """(G - C) / (G + C) on the given strand."""
    g, c = seq.count("G"), seq.count("C")
    return (g - c) / (g + c) if g + c else float("nan")


def sequence_metrics(nt: str, weights: dict[str, float] | None) -> dict[str, float]:
    """GC, GC at codon positions 1-3, GC skew and CAI of one sequence."""
    gc = lambda s: (s.count("G") + s.count("C")) / len(s) if s else float("nan")  # noqa: E731
    out = {"length": len(nt), "gc": gc(nt), "gc_skew": gc_skew(nt)}
    if weights is None:  # non-coding (16S)
        out.update({"gc1": np.nan, "gc2": np.nan, "gc3": np.nan, "cai": np.nan})
        return out
    body = nt[: len(nt) - len(nt) % 3]
    out.update({"gc1": gc(body[0::3]), "gc2": gc(body[1::3]), "gc3": gc(body[2::3]),
                "cai": cai(nt, weights)})
    return out


def _weights(label: str) -> dict[str, float]:
    ribo = [str(r.seq) for r in SeqIO.parse(str(GENOMES_DIR / f"{label}_ribosomal.fna"), "fasta")]
    return cai_weights(ribo)


def per_sequence_table() -> pd.DataFrame:
    """Metrics for every QC-passing record of every gene."""
    sp = species_map()
    wcache: dict[str, dict] = {}
    rows = []
    for gene in all_genes():
        for label, nt in load_gene(gene, "nt").items():
            w = None if gene == "16S" else wcache.setdefault(label, _weights(label))
            rows.append({"gene": gene, "label": label, "species": sp[label],
                         "class": "rRNA control" if gene == "16S" else gene_class(gene),
                         **sequence_metrics(nt, w)})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "m2_per_sequence.csv", index=False)
    return df


def gene_level(df: pd.DataFrame, species: str = FOCAL) -> pd.DataFrame:
    """Mean of each metric per gene within one species (coding genes only)."""
    sub = df[(df["species"] == species) & (df["gene"].isin(coding_genes()))]
    out = sub.groupby("gene")[METRICS + ["gc1", "gc2", "length"]].mean()
    out["n_sequences"] = sub.groupby("gene").size()
    out["class"] = [gene_class(g) for g in out.index]
    return out.reindex([g for g in coding_genes() if g in out.index])


def aa_composition(species: str = FOCAL) -> pd.DataFrame:
    """Mean amino-acid fractions per gene over the strains of one species."""
    sp = species_map()
    rows = {}
    for gene in coding_genes():
        fracs = []
        for label, aa in load_gene(gene, "aa").items():
            if sp[label] == species:
                fracs.append([aa.count(a) / len(aa) for a in AMINO_ACIDS])
        if fracs:
            rows[gene] = np.mean(fracs, axis=0)
    return pd.DataFrame.from_dict(rows, orient="index", columns=list(AMINO_ACIDS))


def rscu_table(species: str = FOCAL) -> pd.DataFrame:
    """RSCU per gene (pooled over strains) plus the ribosomal-protein reference set."""
    sp = species_map()
    cols = {}
    for gene in coding_genes():
        seqs = [s for lab, s in load_gene(gene, "nt").items() if sp[lab] == species]
        if seqs:
            cols[gene] = rscu(codon_counts(seqs))
    rep = next(g.label for g in genomes() if g.species == species)
    ribo = [str(r.seq) for r in SeqIO.parse(str(GENOMES_DIR / f"{rep}_ribosomal.fna"), "fasta")]
    cols["ribosomal proteins"] = rscu(codon_counts(ribo))
    tab = pd.DataFrame(cols).loc[[c for c in SENSE_CODONS if len(SYNONYMS[CODE[c]]) > 1]]
    tab.index = [f"{CODE[c]}-{c}" for c in tab.index]
    return tab.sort_index()


def genome_background() -> pd.DataFrame:
    """GC3 and CAI of every CDS in each species representative; panel genes marked."""
    rows = []
    for g in genomes():
        path = GENOMES_DIR / f"{g.label}_cds.fna.gz"
        if not path.exists():
            continue
        w = _weights(g.label)
        with gzip.open(path, "rt") as fh:
            for r in SeqIO.parse(fh, "fasta"):
                s = str(r.seq).upper()
                if len(s) < 300 or len(s) % 3:
                    continue
                m = sequence_metrics(s, w)
                rows.append({"label": g.label, "species": g.species, "id": r.id,
                             "gc3": m["gc3"], "cai": m["cai"], "gc": m["gc"]})
    bg = pd.DataFrame(rows)
    bg.to_csv(RESULTS / "m2_genome_background.csv", index=False)
    return bg


def percentiles(per_seq: pd.DataFrame, bg: pd.DataFrame) -> pd.DataFrame:
    """Percentile of each panel gene's CAI and GC3 within its own genome's distribution."""
    rows = []
    for label, sub in bg.groupby("label"):
        mine = per_seq[(per_seq["label"] == label) & per_seq["gene"].isin(coding_genes())]
        for r in mine.to_dict("records"):
            rows.append({"label": label, "species": r["species"], "gene": r["gene"],
                         "class": r["class"], "cai": r["cai"],
                         "cai_percentile": (sub["cai"] < r["cai"]).mean() * 100,
                         "gc3_percentile": (sub["gc3"] < r["gc3"]).mean() * 100})
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "m2_genome_percentiles.csv", index=False)
    return out


def run() -> dict[str, pd.DataFrame]:
    """Compute all M2 tables, statistics and figures."""
    from . import m2_figures
    n_boot = load_config()["bootstrap"]["stats_replicates"]
    per_seq = per_sequence_table()
    gl = gene_level(per_seq)
    gl.to_csv(RESULTS / "m2_gene_level_Smutans.csv")
    tests = compare_table(gl.reset_index(), METRICS, rng(2), n_boot=n_boot)
    tests.insert(0, "species", FOCAL)
    tests.to_csv(RESULTS / "m2_tests.csv", index=False)
    by_sp = []
    for sp in sorted(set(per_seq["species"]) - {FOCAL}):
        g2 = gene_level(per_seq, sp)
        if (g2["class"] == "virulence").sum() >= 2:
            t = compare_table(g2.reset_index(), METRICS, rng(3), n_boot=n_boot)
            t.insert(0, "species", sp)
            by_sp.append(t)
    if by_sp:
        pd.concat(by_sp).to_csv(RESULTS / "m2_tests_other_species.csv", index=False)
    aa = aa_composition()
    aa["class"] = [gene_class(g) for g in aa.index]
    aa_tests = compare_table(aa.reset_index(), list(AMINO_ACIDS), rng(4), n_boot=1000)
    aa_tests["q_bh"] = benjamini_hochberg(aa_tests["p_mannwhitney"])
    aa_tests.to_csv(RESULTS / "m2_aa_tests.csv", index=False)
    aa.to_csv(RESULTS / "m2_aa_composition_Smutans.csv")
    rs = rscu_table()
    rs.to_csv(RESULTS / "m2_rscu_Smutans.csv")
    bg = genome_background()
    pct = percentiles(per_seq, bg)
    m2_figures.plot_all(per_seq, gl, tests, aa, aa_tests, rs, bg, pct)
    print(tests[["metric", "median_virulence", "median_control", "cliffs_delta",
                 "delta_ci_low", "delta_ci_high", "p_mannwhitney", "q_bh"]].round(4).to_string(index=False))
    return {"per_sequence": per_seq, "gene_level": gl, "tests": tests, "aa_tests": aa_tests}
