"""M3: pairwise and progressive alignments, identity matrices and per-site conservation.

Conservation of an alignment column is 1 - H / log2(K), where H is the Shannon entropy
(bits) of the Henikoff-weighted residue frequencies among non-gap rows and K is the
alphabet size (20 amino acids, 4 nucleotides). 1 = invariant, 0 = all residues equally
frequent. Columns with more than 50% gaps are not scored.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .alignment import all_pairs, center_star, codon_alignment, henikoff_weights
from .config import GENES_DIR, RESULTS, all_genes, gene_class, load_config, rng
from .seqio import load_gene, read_fasta, species_map
from .stats import compare_table

ALN_DIR = RESULTS / "alignments"
CONS_DIR = RESULTS / "conservation"
ID_DIR = RESULTS / "identity"
REF = "Smu_UA159"
FOCAL = "S. mutans"
WINDOW = 30
GTFC_ID = "Smu_UA159|SMU_RS04625"  # UA159 GtfC, numbering reference for M6/M7


def write_fasta(path, seqs: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))


def column_conservation(msa: dict[str, str], kind: str = "aa", ref: str | None = REF,
                        max_gap: float = 0.5) -> pd.DataFrame:
    """Per-column entropy, conservation, gap fraction and reference-residue numbering."""
    labels = list(msa)
    rows = np.array([list(msa[k]) for k in labels])
    w = henikoff_weights(msa) if len(labels) > 1 else np.ones(1)
    K = 20 if kind == "aa" else 4
    ref_row = msa.get(ref) if ref else None
    out, ref_pos = [], 0
    for j in range(rows.shape[1]):
        col = rows[:, j]
        mask = col != "-"
        gap = 1 - mask.mean()
        ent = cons = np.nan
        if mask.any() and gap <= max_gap:
            ww = w[mask] / w[mask].sum()
            _, inv = np.unique(col[mask], return_inverse=True)
            freqs = np.bincount(inv, weights=ww)
            ent = float(-(freqs * np.log2(freqs)).sum()) + 0.0
            cons = 1 - ent / np.log2(K)
        rp, rres = np.nan, ""
        if ref_row is not None and ref_row[j] != "-":
            ref_pos += 1
            rp, rres = ref_pos, ref_row[j]
        vals, cnt = np.unique(col[mask], return_counts=True) if mask.any() else (["-"], [1])
        out.append({"column": j + 1, "ref_position": rp, "ref_residue": rres,
                    "consensus": vals[int(np.argmax(cnt))], "gap_fraction": round(gap, 4),
                    "entropy_bits": ent, "conservation": cons})
    return pd.DataFrame(out)


def region_window(length: int) -> int:
    """Window size: ~10% of the protein, between 10 and 30 residues."""
    return int(min(WINDOW, max(10, length // 10)))


def extreme_regions(cons: pd.DataFrame, ref_seq: str, n: int = 3) -> list[dict]:
    """Top-n most and least conserved non-overlapping windows in reference coordinates.

    Windows of the two kinds never overlap each other either.
    """
    window = region_window(len(ref_seq))
    s = cons.dropna(subset=["ref_position"]).set_index("ref_position")["conservation"]
    if len(s) < window:
        return []
    sm = s.rolling(window, min_periods=int(window * 0.7)).mean().dropna()
    found, used = [], []
    for kind, order in (("most conserved", False), ("least conserved", True)):
        taken: list[int] = []
        for pos, val in sm.sort_values(ascending=order).items():
            end = int(pos)
            if any(abs(end - t) < window for t in used):
                continue
            taken.append(end)
            used.append(end)
            start = end - window + 1
            found.append({"type": kind, "rank": len(taken), "window": window,
                          "start": start, "end": end,
                          "mean_conservation": round(float(val), 4),
                          "ref_sequence": ref_seq[start - 1:end]})
            if len(taken) == n:
                break
    return found


def analyse_gene(gene: str) -> tuple[dict, pd.DataFrame, list[dict]]:
    """Align one gene, save outputs, and return (summary, conservation table, regions)."""
    kind = "nt" if gene == "16S" else "aa"
    seqs = load_gene(gene, kind)
    sp = species_map()
    pid, scores = all_pairs(seqs, kind)
    labels = list(seqs)
    pd.DataFrame(pid, index=labels, columns=labels).round(3).to_csv(ID_DIR / f"{gene}.csv")
    msa = center_star(seqs, kind, scores)
    write_fasta(ALN_DIR / f"{gene}.{kind}.fasta", msa)
    if kind == "aa":
        write_fasta(ALN_DIR / f"{gene}.codon.fasta", codon_alignment(msa, load_gene(gene, "nt")))
    cons = column_conservation(msa, kind)
    cons.to_csv(CONS_DIR / f"{gene}.csv", index=False)
    focal = [i for i, k in enumerate(labels) if sp[k] == FOCAL]
    sub_pid = pid[np.ix_(focal, focal)][np.triu_indices(len(focal), 1)]
    focal_msa = {k: msa[k] for k in labels if sp[k] == FOCAL}
    focal_cons = column_conservation(focal_msa, kind)
    summary = {
        "gene": gene, "class": "rRNA control" if gene == "16S" else gene_class(gene),
        "n_sequences": len(seqs), "n_species": len({sp[k] for k in labels}),
        "n_Smutans": len(focal), "alignment_columns": len(next(iter(msa.values()))),
        "mean_pid_all": float(pid[np.triu_indices(len(labels), 1)].mean()),
        "min_pid_all": float(pid.min()),
        "mean_pid_Smutans": float(sub_pid.mean()),
        "mean_conservation_all": float(cons["conservation"].mean()),
        "mean_entropy_Smutans": float(focal_cons["entropy_bits"].mean()),
        "variable_sites_Smutans": int((focal_cons["entropy_bits"] > 0).sum()),
    }
    regions = [{"gene": gene, **r} for r in extreme_regions(cons, seqs[REF])] if REF in seqs else []
    return summary, cons, regions


def gtf_family() -> pd.DataFrame:
    """Align the GH70 glucansucrase family (M6 input) and score its conservation."""
    fam = read_fasta(GENES_DIR / "gtf_family.faa")
    _, scores = all_pairs(fam, "aa")
    msa = center_star(fam, "aa", scores)
    write_fasta(ALN_DIR / "gtf_family.aa.fasta", msa)
    ref = GTFC_ID if GTFC_ID in fam else next(iter(fam))
    cons = column_conservation(msa, "aa", ref=ref)
    cons.to_csv(CONS_DIR / "gtf_family.csv", index=False)
    return cons


def run() -> dict[str, pd.DataFrame]:
    """Run M3 for every gene and the GH70 family; test virulence vs control conservation."""
    from . import m3_figures
    for d in (ALN_DIR, CONS_DIR, ID_DIR):
        d.mkdir(parents=True, exist_ok=True)
    summaries, regions, cons_all = [], [], {}
    for gene in all_genes():
        s, cons, reg = analyse_gene(gene)
        summaries.append(s)
        regions += reg
        cons_all[gene] = cons
    summary = pd.DataFrame(summaries)
    summary.to_csv(RESULTS / "m3_summary.csv", index=False)
    reg_df = pd.DataFrame(regions)
    reg_df.to_csv(RESULTS / "m3_extreme_regions.csv", index=False)
    coding = summary[summary["class"] != "rRNA control"]
    n_boot = load_config()["bootstrap"]["stats_replicates"]
    tests = compare_table(coding, ["mean_pid_Smutans", "mean_entropy_Smutans"], rng(5), n_boot=n_boot)
    tests.to_csv(RESULTS / "m3_tests.csv", index=False)
    fam_cons = gtf_family()
    m3_figures.plot_all(summary, cons_all, reg_df, tests)
    print(summary[["gene", "n_sequences", "n_species", "mean_pid_all", "mean_pid_Smutans",
                   "mean_entropy_Smutans"]].round(3).to_string(index=False))
    print(tests[["metric", "median_virulence", "median_control", "cliffs_delta", "delta_ci_low",
                 "delta_ci_high", "q_bh"]].round(4).to_string(index=False))
    return {"summary": summary, "regions": reg_df, "tests": tests, "gtf_family": fam_cons}
