"""M1 (part 2): quality control, catalog, presence/absence matrix and ACCESSIONS.md.

Exclusion rule (applied to every record, documented in the catalog ``flags`` column):
a record is excluded if it has any non-ACGT base; if it is a protein-coding gene whose
length is not a multiple of 3, lacks a start (ATG/GTG/TTG) or stop codon, contains an
internal stop codon, is annotated as a pseudogene, or is shorter than 80% of the UA159
query ("truncated"); or if it is a 16S rRNA gene shorter than 1400 nt.
Records longer than 120% of the query are flagged "extended" as a warning only: a complete
CDS with valid start and stop codons that is longer than the query reflects biology (extra
domains), not a broken record. See DECISIONS.md.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from .config import GENES_DIR, ROOT, RESULTS, all_genes, gene_class, genomes, load_config
from .entrez_client import download_log
from .plotting import CLASS_COLORS, OKABE_ITO, SPECIES_COLORS, save
from .seqio import CATALOG, mode_tag, read_fasta

STARTS = {"ATG", "GTG", "TTG"}
STOPS = {"TAA", "TAG", "TGA"}
WARNING_ONLY = {"extended"}


def qc_flags(nt: str, gene: str, ref_len: int, pseudo: bool = False) -> list[str]:
    """Return the list of QC problems of one record (empty list = passes)."""
    cfg = load_config()["qc"]
    flags = []
    if sum(c not in "ACGT" for c in nt) > cfg["max_ambiguous"]:
        flags.append("ambiguous_bases")
    if gene == "16S":
        if len(nt) < cfg["min_16s_length"]:
            flags.append("short_16S")
        return flags
    if pseudo:
        flags.append("annotated_pseudogene")
    if len(nt) % 3:
        flags.append("length_not_multiple_of_3")
    if nt[:3] not in STARTS:
        flags.append("no_start_codon")
    if nt[-3:] not in STOPS:
        flags.append("no_stop_codon")
    if "*" in str(Seq(nt[: len(nt) - len(nt) % 3 - 3]).translate(table=11)):
        flags.append("internal_stop")
    ratio = len(nt) / ref_len
    if ratio < cfg["min_length_fraction"]:
        flags.append("truncated")
    elif ratio > cfg["max_length_fraction"]:
        flags.append("extended")
    return flags


def build_catalog(hits: pd.DataFrame) -> pd.DataFrame:
    """Apply QC to every ortholog call and write ``results/m1_catalog.csv``."""
    cat = hits[hits["ortholog_call"]].copy()
    ref_label = genomes()[0].label
    nt_len, flags = [], []
    for gene, sub in cat.groupby("gene", sort=False):
        seqs = read_fasta(GENES_DIR / f"{gene}.fna")
        ref_len = len(seqs[ref_label])
        for idx, r in sub.iterrows():
            nt = seqs[r.label]
            nt_len.append((idx, len(nt), round(len(nt) / ref_len, 4)))
            flags.append((idx, ";".join(qc_flags(nt, gene, ref_len, bool(r.pseudo)))))
    cat.loc[[i for i, *_ in nt_len], "nt_length"] = [n for _, n, _ in nt_len]
    cat.loc[[i for i, *_ in nt_len], "length_ratio_to_UA159"] = [x for *_, x in nt_len]
    cat.loc[[i for i, _ in flags], "flags"] = [f for _, f in flags]
    cat["flags"] = cat["flags"].fillna("")
    cat["included"] = [not (set(f.split(";")) - WARNING_ONLY - {""}) for f in cat["flags"]]
    cat["class"] = [("rRNA control" if g == "16S" else gene_class(g)) for g in cat["gene"]]
    cat["nt_length"] = cat["nt_length"].astype(int)
    RESULTS.mkdir(parents=True, exist_ok=True)
    cat.to_csv(CATALOG, index=False)
    hits[~hits["ortholog_call"]].to_csv(RESULTS / "m1_rejected_hits.csv", index=False)
    return cat


def presence_matrix(cat: pd.DataFrame) -> pd.DataFrame:
    """Gene x genome matrix: 2 = included, 1 = found but excluded by QC, 0 = no ortholog."""
    labels = [g.label for g in genomes()]
    mat = pd.DataFrame(0, index=all_genes(), columns=labels)
    for r in cat.itertuples():
        mat.at[r.gene, r.label] = 2 if r.included else 1
    mat.to_csv(RESULTS / "m1_presence_matrix.csv")
    return mat


def qc_summary(cat: pd.DataFrame) -> pd.DataFrame:
    """Per-gene counts and length statistics."""
    g = cat.groupby("gene", sort=False)
    out = pd.DataFrame({
        "class": g["class"].first(),
        "n_found": g.size(),
        "n_included": g["included"].sum(),
        "n_species": g.apply(lambda d: d.loc[d.included, "species"].nunique()),
        "median_nt_length": g["nt_length"].median(),
        "min_nt_length": g["nt_length"].min(),
        "max_nt_length": g["nt_length"].max(),
    }).reindex(all_genes())
    out.index.name = "gene"
    out.to_csv(RESULTS / "m1_qc_summary.csv")
    return out


def plot_presence(mat: pd.DataFrame) -> None:
    """Heatmap of ortholog presence (included / excluded / absent) across genomes."""
    fig, ax = plt.subplots(figsize=(11, 5.2))
    cmap = ListedColormap(["#F2F2F2", OKABE_ITO["yellow"], OKABE_ITO["blue"]])
    ax.imshow(mat.values, cmap=cmap, vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(mat.shape[0]), mat.index)
    n_vir = len(load_config()["genes"]["virulence"])
    ax.axhline(n_vir - 0.5, color="black", lw=1)
    ax.set_xticks(np.arange(-0.5, mat.shape[1]), minor=True)
    ax.set_yticks(np.arange(-0.5, mat.shape[0]), minor=True)
    ax.grid(which="minor", color="white", lw=1)
    ax.tick_params(which="minor", length=0)
    for lab in ax.get_xticklabels():
        sp = next(g.species for g in genomes() if g.label == lab.get_text())
        lab.set_color(SPECIES_COLORS.get(sp, "black"))
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in cmap.colors[::-1]]
    ax.legend(handles, ["ortholog, passed QC", "ortholog, excluded by QC", "no ortholog (RBH)"],
              loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title("Reciprocal-best-hit orthologs of the UA159 gene panel\n(red labels = virulence-associated, blue = housekeeping controls)" + mode_tag())
    for lab in ax.get_yticklabels():
        cls = "virulence" if lab.get_text() in load_config()["genes"]["virulence"] else "control"
        lab.set_color(CLASS_COLORS[cls])
    save(fig, "m1_presence_absence")


def plot_lengths(cat: pd.DataFrame) -> None:
    """Length of every record relative to UA159, with the QC acceptance band."""
    genes = [g for g in all_genes() if g in set(cat["gene"])]
    fig, ax = plt.subplots(figsize=(11, 4.5))
    cfg = load_config()["qc"]
    ax.axhspan(cfg["min_length_fraction"], cfg["max_length_fraction"], color="#EEEEEE", zorder=0)
    rng = np.random.default_rng(0)
    for i, gene in enumerate(genes):
        sub = cat[cat["gene"] == gene]
        x = i + rng.uniform(-0.18, 0.18, len(sub))
        colors = [SPECIES_COLORS.get(s, "grey") for s in sub["species"]]
        ax.scatter(x, sub["length_ratio_to_UA159"], c=colors, s=22,
                   marker="o", edgecolor="none", alpha=0.9)
        bad = sub[~sub["included"]]
        ax.scatter(x[~sub["included"].values], bad["length_ratio_to_UA159"], marker="x",
                   color="black", s=40, lw=1.2)
        ax.text(i, ax.get_ylim()[0], f"n={int(sub['included'].sum())}", ha="center",
                va="bottom", fontsize=7)
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right")
    ax.set_ylabel("Length / UA159 length")
    ax.set_title("Record length QC (grey band = expected range; x = excluded)" + mode_tag())
    for sp, c in SPECIES_COLORS.items():
        ax.scatter([], [], color=c, label=sp)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1))
    save(fig, "m1_length_qc")


def write_accessions(cat: pd.DataFrame, fam: pd.DataFrame | None = None) -> None:
    """Write ACCESSIONS.md: every genome, gene record and structure used, with access dates."""
    from .m1_retrieval import genome_meta
    log = download_log()
    lines = ["# Accessions", "",
             "Every sequence and structure used in this project, with the date it was "
             "downloaded. Coordinates are 1-based and inclusive.", "",
             "## Genomes (NCBI nuccore, RefSeq)", "",
             "| Label | Species | Accession | Length (bp) | GC | Title | Accessed |",
             "|---|---|---|---|---|---|---|"]
    for g in genomes():
        meta = genome_meta(g)
        lines.append(f"| {g.label} | *{g.species}* | {g.accession} | {meta['length']:,} | "
                     f"{meta['gc']:.3f} | {meta['title']} | {log.get(g.accession, {}).get('accessed', '')} |")
    lines += ["", "## Gene records", "",
              "| Gene | Label | Locus tag | Protein ID | Location | Identity to UA159 | Included | QC flags |",
              "|---|---|---|---|---|---|---|---|"]
    for r in cat.itertuples():
        pid = r.protein_id if isinstance(r.protein_id, str) else ""
        lines.append(f"| {r.gene} | {r.label} | {r.locus_tag} | {pid} | {r.accession}:{r.start}-{r.end}"
                     f"({r.strand}) | {r.identity_to_query:.3f} | {'yes' if r.included else 'no'} | {r.flags} |")
    if fam is not None and len(fam):
        lines += ["", "## Glucansucrase (GH70) family set used in M6", "",
                  "| ID | Species | Protein ID | Product | Length (aa) |", "|---|---|---|---|---|"]
        for r in fam.itertuples():
            lines.append(f"| {r.id} | *{r.species}* | {r.protein_id} | {r.product} | {r.length} |")
    other = {k: v for k, v in log.items() if not k.startswith(("NC_", "NZ_"))}
    if other:
        lines += ["", "## Other downloads", "",
                  "Structures of *S. mutans* GtfC. 3AIE (2.1 A, highest resolution) chain A is "
                  "the one analysed in M7; the others were inspected when choosing it.", "",
                  "| File | Source | URL | Accessed |", "|---|---|---|---|"]
        for k, v in sorted(other.items()):
            lines.append(f"| {k} | {v.get('source', '')} | {v.get('url', '')} | {v.get('accessed', '')} |")
    (ROOT / "ACCESSIONS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(force_synthetic: bool = False) -> pd.DataFrame:
    """Run M1 end to end. Falls back to simulated data if NCBI cannot be used."""
    from . import m1_retrieval, synthetic
    from .entrez_client import ncbi_reachable
    from .seqio import MODE_FILE
    RESULTS.mkdir(parents=True, exist_ok=True)
    mode, fam = "REAL", None
    cached = all((m1_retrieval.PARSED / f"{g.accession}.pkl").exists() for g in genomes())
    if force_synthetic or not (cached or ncbi_reachable()):
        mode = "SYNTHETIC"
    else:
        try:
            hits = m1_retrieval.extract_panel()
            fam = m1_retrieval.extract_gtf_family()
        except Exception as exc:  # any network or parsing failure triggers the fallback
            print(f"  NCBI retrieval failed ({exc}); switching to SYNTHETIC data")
            mode = "SYNTHETIC"
    if mode == "SYNTHETIC":
        hits, fam = synthetic.generate_dataset()
    MODE_FILE.write_text(mode + "\n")
    fam.to_csv(RESULTS / "m1_gtf_family.csv", index=False)
    cat = build_catalog(hits)
    mat = presence_matrix(cat)
    summary = qc_summary(cat)
    plot_presence(mat)
    plot_lengths(cat)
    if mode == "REAL":
        write_accessions(cat, fam)
    print(f"  M1 [{mode}]: {int(cat['included'].sum())}/{len(cat)} records pass QC")
    print(summary[["n_found", "n_included", "n_species"]].to_string())
    return cat
