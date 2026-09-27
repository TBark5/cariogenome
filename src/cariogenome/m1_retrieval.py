"""M1 (part 1): download genomes, parse annotations and extract the gene panel.

Orthologs of each UA159 query gene are called by reciprocal best hit, not by annotation
name. The first 16S rRNA copy of each genome is taken. Every extracted record is cached
in ``data/genes`` (FASTA) and ``data/genbank`` (GenBank slice of the source record).
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

from . import entrez_client
from .config import (GENBANK_DIR, GENES_DIR, GENOMES_DIR, RAW, Genome, genomes,
                     load_config)
from .homology import ProteinIndex, all_hits, local_stats, reciprocal_best_hit

PARSED = RAW / "parsed"
RIBO = re.compile(r"^(30S|50S) ribosomal protein")


def parse_genome(genome: Genome) -> pd.DataFrame:
    """One row per CDS / rRNA feature with sequences; cached as a pickle after first parse."""
    cache = PARSED / f"{genome.accession}.pkl"
    if cache.exists():
        return pd.read_pickle(cache)
    record = SeqIO.read(entrez_client.fetch_genome(genome.accession), "genbank")
    rows = []
    for f in record.features:
        if f.type not in ("CDS", "rRNA"):
            continue
        q = f.qualifiers
        rows.append({
            "type": f.type,
            "locus_tag": q.get("locus_tag", [""])[0],
            "old_locus_tag": q.get("old_locus_tag", [""])[0],
            "protein_id": q.get("protein_id", [""])[0],
            "gene_name": q.get("gene", [""])[0],
            "product": q.get("product", [""])[0],
            "start": int(f.location.start),
            "end": int(f.location.end),
            "strand": int(f.location.strand or 1),
            "pseudo": "pseudo" in q or "pseudogene" in q,
            "nt": str(f.extract(record.seq)).upper(),
            "aa": q.get("translation", [""])[0],
        })
    df = pd.DataFrame(rows)
    PARSED.mkdir(parents=True, exist_ok=True)
    df.to_pickle(cache)
    meta = {"title": record.description, "length": len(record.seq),
            "gc": round(gc_fraction(record.seq), 4)}
    (PARSED / f"{genome.accession}.json").write_text(json.dumps(meta))
    return df


def genome_meta(genome: Genome) -> dict:
    """Title, length and GC of a parsed genome."""
    path = PARSED / f"{genome.accession}.json"
    if not path.exists():
        parse_genome(genome)
    return json.loads(path.read_text())


def _write_genbank_slice(record: SeqRecord, genome: Genome, row: pd.Series, gene: str) -> Path:
    """Save the source-record region covering one gene as a small GenBank file."""
    sub: SeqRecord = record[row.start:row.end]
    sub.id = f"{genome.accession}:{row.start + 1}-{row.end}"
    sub.name = f"{genome.label}_{gene}"[:16]
    sub.description = f"{genome.label} {gene} {row.locus_tag} region of {genome.accession}"
    sub.annotations = {"molecule_type": "DNA", "organism": genome.species,
                       "source": record.annotations.get("source", genome.species)}
    out = GENBANK_DIR / gene / f"{genome.label}.gb"
    out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(sub, out, "genbank")
    return out


def _query_proteins(ref: pd.DataFrame) -> dict[str, int]:
    """Map each coding query gene to its row index in the UA159 table."""
    genes_cfg = load_config()["genes"]
    tags = {**genes_cfg["virulence"], **genes_cfg["housekeeping"]}
    out = {}
    for gene, tag in tags.items():
        hit = ref.index[(ref["type"] == "CDS") & (ref["old_locus_tag"] == tag)]
        if len(hit) != 1:
            raise ValueError(f"query {gene} ({tag}) not found once in the reference")
        out[gene] = int(hit[0])
    return out


def extract_panel() -> pd.DataFrame:
    """Find every panel gene in every genome and cache sequences. Returns raw hits."""
    cfg = load_config()["homology"]
    gl = genomes()
    tables = {g.label: parse_genome(g) for g in gl}
    ref_label = gl[0].label
    ref = tables[ref_label]
    ref_cds = ref[(ref["type"] == "CDS") & (ref["aa"] != "")]
    ref_index = ProteinIndex(list(ref_cds["aa"]), cfg["kmer"])
    ref_pos = {i: p for p, i in enumerate(ref_cds.index)}
    queries = _query_proteins(ref)

    rows = []
    for g in gl:
        df = tables[g.label]
        cds = df[(df["type"] == "CDS") & (df["aa"] != "")]
        index = ProteinIndex(list(cds["aa"]), cfg["kmer"])
        for gene, qi in queries.items():
            qseq = ref.at[qi, "aa"]
            hit, rbh = reciprocal_best_hit(qseq, ref_pos[qi], index, ref_index,
                                           cfg["prefilter_top"])
            if hit is None:
                continue
            row = cds.iloc[hit.index]
            ok = rbh and hit.identity >= cfg["min_identity"] and min(hit.qcov, hit.tcov) >= cfg["min_coverage"]
            rows.append(_hit_row(g, gene, row, hit.identity, hit.qcov, hit.tcov, rbh, ok))
        rrna = df[(df["type"] == "rRNA") & df["product"].str.contains("16S")].sort_values("start")
        if len(rrna):
            rows.append(_hit_row(g, "16S", rrna.iloc[0], 1.0, 1.0, 1.0, True, True))
        _save_genome_sets(g, df, is_representative=_is_rep(g, gl))
    hits = pd.DataFrame(rows)
    _write_gene_fastas(hits, {g.label: g for g in gl}, tables)
    return hits


def _is_rep(genome: Genome, gl: list[Genome]) -> bool:
    return next(x for x in gl if x.species == genome.species).label == genome.label


def _hit_row(g: Genome, gene: str, row: pd.Series, ident: float, qcov: float, tcov: float,
             rbh: bool, ortholog: bool) -> dict:
    return {
        "gene": gene, "label": g.label, "species": g.species, "accession": g.accession,
        "row_index": int(row.name), "locus_tag": row.locus_tag,
        "old_locus_tag": row.old_locus_tag, "protein_id": row.protein_id,
        "annotation_gene": row.gene_name, "product": row["product"],
        "start": int(row.start) + 1, "end": int(row.end), "strand": "+" if row.strand == 1 else "-",
        "pseudo": bool(row.pseudo), "identity_to_query": round(float(ident), 4),
        "query_coverage": round(float(qcov), 4), "hit_coverage": round(float(tcov), 4),
        "reciprocal_best_hit": bool(rbh), "ortholog_call": bool(ortholog),
    }


def _write_gene_fastas(hits: pd.DataFrame, gmap: dict[str, Genome],
                       tables: dict[str, pd.DataFrame]) -> None:
    """Write nucleotide and protein FASTA per gene plus one GenBank slice per record."""
    GENES_DIR.mkdir(parents=True, exist_ok=True)
    calls = hits[hits["ortholog_call"]]
    for gene, sub in calls.groupby("gene", sort=False):
        with open(GENES_DIR / f"{gene}.fna", "w") as fna, open(GENES_DIR / f"{gene}.faa", "w") as faa:
            for r in sub.itertuples():
                src = tables[r.label].loc[r.row_index]
                head = f">{r.label} gene={gene} locus_tag={r.locus_tag} {r.accession}:{r.start}-{r.end}({r.strand})"
                fna.write(f"{head}\n{src.nt}\n")
                if src.aa:
                    faa.write(f"{head} protein_id={r.protein_id}\n{src.aa}\n")
    for label, sub in calls.groupby("label", sort=False):
        genome = gmap[label]
        if all((GENBANK_DIR / r.gene / f"{label}.gb").exists() for r in sub.itertuples()):
            continue  # slices already cached; skip re-reading the 5 MB genome
        record = SeqIO.read(entrez_client.fetch_genome(genome.accession), "genbank")
        for r in sub.itertuples():
            _write_genbank_slice(record, genome, tables[label].loc[r.row_index], r.gene)


def _save_genome_sets(g: Genome, df: pd.DataFrame, is_representative: bool) -> None:
    """Ribosomal-protein CDS (CAI reference set) for every genome; all CDS for representatives."""
    GENOMES_DIR.mkdir(parents=True, exist_ok=True)
    cds = df[(df["type"] == "CDS") & (~df["pseudo"]) & (df["aa"] != "")]
    ribo = cds[cds["product"].str.match(RIBO)]
    with open(GENOMES_DIR / f"{g.label}_ribosomal.fna", "w") as fh:
        for r in ribo.itertuples():
            fh.write(f">{r.locus_tag} {r.product}\n{r.nt}\n")
    if is_representative:
        with gzip.open(GENOMES_DIR / f"{g.label}_cds.fna.gz", "wt") as fh:
            for r in cds.itertuples():
                fh.write(f">{r.locus_tag} {r.product}\n{r.nt}\n")


def extract_gtf_family() -> pd.DataFrame:
    """All glucansucrase (GH70) homologs of GtfB/C/D in each species representative."""
    cfg = load_config()["homology"]
    gl = genomes()
    ref = parse_genome(gl[0])
    queries = _query_proteins(ref)
    gtf_seqs = [ref.at[queries[g], "aa"] for g in ("gtfB", "gtfC", "gtfD")]
    rows, seen = [], set()
    for g in (x for x in gl if _is_rep(x, gl)):
        df = parse_genome(g)
        cds = df[(df["type"] == "CDS") & (df["aa"] != "")]
        index = ProteinIndex(list(cds["aa"]), cfg["kmer"])
        for q in gtf_seqs:
            for h in all_hits(q, index, 10, cfg["gtf_family_min_identity"], 0.5):
                r = cds.iloc[h.index]
                if r.aa in seen:
                    continue
                seen.add(r.aa)
                best = max((local_stats(x, r.aa) for x in gtf_seqs), key=lambda s: s.score)
                rows.append({"id": f"{g.label}|{r.locus_tag}", "label": g.label,
                             "species": g.species, "locus_tag": r.locus_tag,
                             "old_locus_tag": r.old_locus_tag, "protein_id": r.protein_id,
                             "product": r["product"], "length": len(r.aa),
                             "best_identity_to_Smu_gtf": round(float(best.identity), 4),
                             "aa": r.aa})
    fam = pd.DataFrame(rows)
    with open(GENES_DIR / "gtf_family.faa", "w") as fh:
        for r in fam.itertuples():
            fh.write(f">{r.id} {r.product} protein_id={r.protein_id}\n{r.aa}\n")
    return fam.drop(columns="aa")
