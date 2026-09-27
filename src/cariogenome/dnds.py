"""Nei-Gojobori (1986) dN/dS, pooled over sequence pairs, with a codon bootstrap.

Sites: for each codon, each of the 9 single-nucleotide neighbours is counted as a
synonymous or nonsynonymous change (a change to a stop codon counts as nonsynonymous,
as in Biopython); S = synonymous count / 3 and N = 3 - S. Sites of a codon pair are the
average of the two codons.

Differences: codon pairs differing at one position are classified directly. For two or
three differences, all mutational pathways are enumerated and the synonymous and
nonsynonymous counts are averaged over pathways that avoid stop codons.

Pooled estimate over pairs: pS = sum Sd / sum S, pN = sum Nd / sum N, each Jukes-Cantor
corrected, d = -3/4 ln(1 - 4p/3); omega = dN / dS. Codons with a gap, an ambiguous base or
a stop in either sequence are skipped for that pair (pairwise deletion).
"""

from __future__ import annotations

from functools import cache
from itertools import combinations, permutations

import numpy as np

from .codon import CODE, STOP_CODONS

BASES = "ACGT"


@cache
def codon_sites(codon: str) -> tuple[float, float]:
    """(synonymous sites, nonsynonymous sites) of one sense codon."""
    syn = 0.0
    for pos in range(3):
        for b in BASES:
            if b == codon[pos]:
                continue
            mut = codon[:pos] + b + codon[pos + 1 :]
            if mut not in STOP_CODONS and CODE[mut] == CODE[codon]:
                syn += 1
    s = syn / 3.0
    return s, 3.0 - s


@cache
def codon_diffs(c1: str, c2: str) -> tuple[float, float]:
    """(synonymous, nonsynonymous) differences between two sense codons (NG86 pathways)."""
    pos = [i for i in range(3) if c1[i] != c2[i]]
    if not pos:
        return 0.0, 0.0
    tot_s = tot_n = 0.0
    n_paths = 0
    for order in permutations(pos):
        cur, s, n, ok = c1, 0.0, 0.0, True
        for p in order:
            nxt = cur[:p] + c2[p] + cur[p + 1 :]
            if nxt in STOP_CODONS:
                ok = False
                break
            if CODE[nxt] == CODE[cur]:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            tot_s, tot_n, n_paths = tot_s + s, tot_n + n, n_paths + 1
    if n_paths == 0:  # every pathway passes a stop codon: treat all as nonsynonymous
        return 0.0, float(len(pos))
    return tot_s / n_paths, tot_n / n_paths


def jc(p: np.ndarray | float) -> np.ndarray:
    """Jukes-Cantor correction; NaN where p >= 0.75 (saturated)."""
    p = np.asarray(p, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        # "+ 0.0" turns -0.0 (from p = 0) into 0.0 so tables never print "-0.000"
        return np.where(p < 0.75, -0.75 * np.log(1 - 4 * p / 3), np.nan) + 0.0


def _split(seq: str) -> list[str]:
    return [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]


def _valid(c: str) -> bool:
    return c in CODE and c not in STOP_CODONS


class PooledNG86:
    """Per-pair, per-codon NG86 counts for one codon alignment and a set of pairs."""

    def __init__(self, msa: dict[str, str], pairs: list[tuple[str, str]]) -> None:
        cods = {k: _split(v) for k, v in msa.items()}
        self.pairs = pairs
        L = len(next(iter(cods.values())))
        shape = (len(pairs), L)
        self.S, self.N = np.zeros(shape), np.zeros(shape)
        self.Sd, self.Nd = np.zeros(shape), np.zeros(shape)
        for i, (a, b) in enumerate(pairs):
            for j, (x, y) in enumerate(zip(cods[a], cods[b], strict=True)):
                if _valid(x) and _valid(y):
                    s1, n1 = codon_sites(x)
                    s2, n2 = codon_sites(y)
                    self.S[i, j], self.N[i, j] = (s1 + s2) / 2, (n1 + n2) / 2
                    self.Sd[i, j], self.Nd[i, j] = codon_diffs(x, y)
        self.n_codons = L

    def estimate(self, weights: np.ndarray | None = None, cols: slice | None = None) -> dict:
        """Pooled pN, pS, dN, dS and omega (optionally codon-weighted or windowed)."""
        S, N, Sd, Nd = self.S, self.N, self.Sd, self.Nd
        if cols is not None:
            S, N, Sd, Nd = S[:, cols], N[:, cols], Sd[:, cols], Nd[:, cols]
        w = np.ones(S.shape[1]) if weights is None else weights
        tS, tN, tSd, tNd = S.sum(0) @ w, N.sum(0) @ w, Sd.sum(0) @ w, Nd.sum(0) @ w
        pS = tSd / tS if tS else np.nan
        pN = tNd / tN if tN else np.nan
        dS, dN = float(jc(pS)) + 0.0, float(jc(pN)) + 0.0  # +0.0 avoids printing -0.0
        omega = dN / dS if dS and dS > 0 else np.nan
        return {
            "pN": pN,
            "pS": pS,
            "dN": dN,
            "dS": dS,
            "omega": omega,
            "Nd": tNd,
            "Sd": tSd,
            "N_sites": tN,
            "S_sites": tS,
        }

    def column_totals(self) -> np.ndarray:
        """(4, L) array of per-codon totals summed over pairs: S, N, Sd, Nd."""
        return np.vstack([self.S.sum(0), self.N.sum(0), self.Sd.sum(0), self.Nd.sum(0)])


def omega_from_totals(tot: np.ndarray) -> np.ndarray:
    """Vectorised omega from (..., 4) totals [S, N, Sd, Nd] (NaN if undefined)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        dS = jc(tot[..., 2] / tot[..., 0])
        dN = jc(tot[..., 3] / tot[..., 1])
        return np.where(dS > 0, dN / dS, np.nan)


def bootstrap_totals(col_tot: np.ndarray, n_rep: int, rng: np.random.Generator) -> np.ndarray:
    """(n_rep, 4) totals from resampling codon columns with replacement."""
    L = col_tot.shape[1]
    counts = np.stack([np.bincount(rng.integers(0, L, L), minlength=L) for _ in range(n_rep)])
    return counts @ col_tot.T


def all_pairs(labels: list[str]) -> list[tuple[str, str]]:
    """Every unordered pair of labels."""
    return list(combinations(labels, 2))
