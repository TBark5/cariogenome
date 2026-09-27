"""Pairwise alignment, center-star progressive MSA and codon back-translation.

No external aligner is used. The multiple alignment is the center-star approximation
(Gusfield 1993): the sequence with the highest summed pairwise score is the center, every
other sequence is aligned to it, and the pairwise alignments are merged so that any gap
inserted in the center is propagated to all rows ("once a gap, always a gap"). This is
faster and simpler than MUSCLE/MAFFT but can place gaps less well among non-center
sequences; it is an approximation and is reported as such.
"""
from __future__ import annotations

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices


def _aligner(kind: str) -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    if kind == "aa":
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score, a.extend_gap_score = -10, -0.5
    else:
        a.match_score, a.mismatch_score = 2, -3
        a.open_gap_score, a.extend_gap_score = -5, -2
    # Terminal gaps are free so length differences at the ends are not over-penalised.
    a.end_gap_score = 0
    return a


ALIGNERS = {"aa": _aligner("aa"), "nt": _aligner("nt")}


def pairwise(a: str, b: str, kind: str = "aa") -> tuple[str, str, float]:
    """Global alignment of two sequences; returns gapped a, gapped b and the score."""
    aln = ALIGNERS[kind].align(a, b)[0]
    return aln[0], aln[1], float(aln.score)


def identity(ga: str, gb: str) -> float:
    """Percent identity = identical pairs / aligned (non-gap) pairs x 100."""
    pairs = [(x, y) for x, y in zip(ga, gb, strict=True) if x != "-" and y != "-"]
    return 100.0 * sum(x == y for x, y in pairs) / len(pairs) if pairs else float("nan")


def all_pairs(seqs: dict[str, str], kind: str) -> tuple[np.ndarray, np.ndarray]:
    """Percent-identity and score matrices from all pairwise global alignments.

    Identical sequences are aligned only once (strains often share alleles).
    """
    labels = list(seqs)
    uniq = list(dict.fromkeys(seqs.values()))
    ui = {s: i for i, s in enumerate(uniq)}
    n_u = len(uniq)
    pid_u, score_u = np.full((n_u, n_u), 100.0), np.zeros((n_u, n_u))
    for i in range(n_u):
        score_u[i, i] = ALIGNERS[kind].score(uniq[i], uniq[i])
        for j in range(i + 1, n_u):
            ga, gb, sc = pairwise(uniq[i], uniq[j], kind)
            pid_u[i, j] = pid_u[j, i] = identity(ga, gb)
            score_u[i, j] = score_u[j, i] = sc
    idx = [ui[seqs[k]] for k in labels]
    return pid_u[np.ix_(idx, idx)], score_u[np.ix_(idx, idx)]


def center_star(seqs: dict[str, str], kind: str = "aa",
                scores: np.ndarray | None = None) -> dict[str, str]:
    """Center-star multiple alignment. Returns gapped sequences in input order."""
    labels = list(seqs)
    if len(labels) == 1:
        return dict(seqs)
    if scores is None:
        _, scores = all_pairs(seqs, kind)
    c = labels[int(np.argmax(scores.sum(axis=1)))]
    center = seqs[c]
    L = len(center)
    ins: dict[str, list[str]] = {}
    match: dict[str, list[str]] = {}
    for lab in labels:
        gc, gs = (center, center) if lab == c else pairwise(center, seqs[lab], kind)[:2]
        # ins[i] = residues of `lab` inserted before center position i (i == L: after end)
        slots, aligned, pos = [""] * (L + 1), ["-"] * L, 0
        for x, y in zip(gc, gs, strict=True):
            if x == "-":
                slots[pos] += y
            else:
                aligned[pos] = y
                pos += 1
        ins[lab], match[lab] = slots, aligned
    width = [max(len(ins[lab][i]) for lab in labels) for i in range(L + 1)]
    out = {}
    for lab in labels:
        parts = []
        for i in range(L + 1):
            parts.append(ins[lab][i].ljust(width[i], "-"))
            if i < L:
                parts.append(match[lab][i])
        out[lab] = "".join(parts)
    return out


def codon_alignment(protein_msa: dict[str, str], nts: dict[str, str]) -> dict[str, str]:
    """Thread each CDS onto its aligned protein (like PAL2NAL); gaps become '---'.

    Each protein must be the translation of its CDS (checked; the start codon is skipped
    because alternative starts GTG/TTG are translated as Met).
    """
    from .codon import CODE
    out = {}
    for lab, row in protein_msa.items():
        nt = nts[lab]
        cods = [nt[i:i + 3] for i in range(0, len(nt) - 2, 3)]
        n_res = sum(ch != "-" for ch in row)
        cods = cods[:n_res]
        k, parts = 0, []
        for ch in row:
            if ch == "-":
                parts.append("---")
            else:
                if k > 0 and CODE.get(cods[k], "X") != ch:
                    raise ValueError(f"{lab}: codon {k} {cods[k]} does not encode {ch}")
                parts.append(cods[k])
                k += 1
        out[lab] = "".join(parts)
    return out


def henikoff_weights(msa: dict[str, str]) -> np.ndarray:
    """Position-based sequence weights (Henikoff & Henikoff 1994), normalised to sum 1.

    Down-weights near-duplicate sequences so that many similar strains do not dominate
    the per-column frequencies.
    """
    rows = np.array([list(s) for s in msa.values()])
    n, L = rows.shape
    w = np.zeros(n)
    for j in range(L):
        col = rows[:, j]
        vals, counts = np.unique(col, return_counts=True)
        k = len(vals)
        lookup = dict(zip(vals, counts, strict=True))
        w += np.array([1.0 / (k * lookup[x]) for x in col])
    return w / w.sum()
