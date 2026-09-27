"""Distance-based phylogenetics: JC69/K2P distances, NJ/UPGMA, splits and bootstrap support.

Bootstrap support is computed on unrooted bipartitions ("splits"), which is correct for
NJ trees whose root position is arbitrary.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from numpy.typing import ArrayLike

ENC = {"A": 0, "C": 1, "G": 2, "T": 3}
MAX_DIST = 3.0  # cap for saturated pairs (log argument <= 0)


def encode(msa: dict[str, str]) -> np.ndarray:
    """Alignment to an (n_seq, n_col) uint8 array; gaps and ambiguity codes -> 4."""
    return np.array([[ENC.get(c, 4) for c in s] for s in msa.values()], dtype=np.uint8)


def _is_purine(x: np.ndarray) -> np.ndarray:
    """True for encoded A (0) and G (2)."""
    return (x == 0) | (x == 2)


class PairCounts:
    """Per-pair, per-column indicators so distances for any column weighting are fast."""

    def __init__(self, arr: np.ndarray) -> None:
        n = arr.shape[0]
        self.n = n
        self.iu = np.triu_indices(n, 1)
        a, b = arr[self.iu[0]], arr[self.iu[1]]
        valid = (a < 4) & (b < 4)
        diff = valid & (a != b)
        transition = diff & (_is_purine(a) == _is_purine(b))
        self.valid = valid.astype(np.float64)
        self.diff = diff.astype(np.float64)
        self.ts = transition.astype(np.float64)
        self.tv = (diff & ~transition).astype(np.float64)

    def distances(self, weights: np.ndarray | None = None, model: str = "K2P") -> np.ndarray:
        """Full symmetric distance matrix under JC69 or K2P (pairwise deletion of gaps)."""
        w = np.ones(self.valid.shape[1]) if weights is None else weights
        sites = self.valid @ w
        sites[sites == 0] = np.nan
        if model == "JC":
            d = jukes_cantor(self.diff @ w / sites)
        else:
            d = kimura_2p(self.ts @ w / sites, self.tv @ w / sites)
        d = np.nan_to_num(d, nan=MAX_DIST)
        mat = np.zeros((self.n, self.n))
        mat[self.iu] = d
        return mat + mat.T


def jukes_cantor(p: ArrayLike) -> np.ndarray:
    """JC69: d = -3/4 ln(1 - 4p/3)."""
    arg = 1 - 4 * np.asarray(p, float) / 3
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(arg > 0, -0.75 * np.log(np.clip(arg, 1e-12, None)), MAX_DIST)


def kimura_2p(P: ArrayLike, Q: ArrayLike) -> np.ndarray:
    """K2P: d = -1/2 ln(1 - 2P - Q) - 1/4 ln(1 - 2Q); P = transitions, Q = transversions."""
    a1 = 1 - 2 * np.asarray(P, float) - np.asarray(Q, float)
    a2 = 1 - 2 * np.asarray(Q, float)
    ok = (a1 > 0) & (a2 > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = -0.5 * np.log(np.clip(a1, 1e-12, None)) - 0.25 * np.log(np.clip(a2, 1e-12, None))
    return np.where(ok, d, MAX_DIST)


def build_tree(dist: np.ndarray, labels: list[str], method: str = "nj") -> Tree:
    """NJ or UPGMA tree from a distance matrix (negative NJ branch lengths set to 0)."""
    lower = [[float(dist[i, j]) for j in range(i + 1)] for i in range(len(labels))]
    dm = DistanceMatrix(list(labels), lower)
    tree = getattr(DistanceTreeConstructor(), method)(dm)
    for c in tree.find_clades():
        if c.branch_length is not None and c.branch_length < 0:
            c.branch_length = 0.0
    return tree


def canonical(side: frozenset, taxa: frozenset) -> frozenset:
    """Represent a bipartition by the side that does not contain the first taxon."""
    anchor = min(taxa)
    return taxa - side if anchor in side else side


def clade_split(clade: Clade, taxa: frozenset) -> frozenset | None:
    """Canonical split induced by the edge above ``clade`` (None if trivial)."""
    side = frozenset(t.name for t in clade.get_terminals())
    if len(side) < 2 or len(taxa - side) < 2:
        return None
    return canonical(side, taxa)


def splits(tree: Tree) -> set[frozenset]:
    """All non-trivial bipartitions of a tree."""
    taxa = frozenset(t.name for t in tree.get_terminals())
    out = set()
    for c in tree.find_clades():
        s = clade_split(c, taxa)
        if s is not None:
            out.add(s)
    return out


def restrict(split_set: set[frozenset] | dict, subset: frozenset, full: frozenset) -> dict:
    """Project splits onto a subset of taxa; values keep the best support if given."""
    items = split_set.items() if isinstance(split_set, dict) else ((s, 1.0) for s in split_set)
    out: dict[frozenset, float] = {}
    for s, sup in items:
        side = s & subset
        if len(side) >= 2 and len(subset - side) >= 2:
            key = canonical(side, subset)
            out[key] = max(out.get(key, 0.0), sup)
    return out


def compatible(a: frozenset, b: frozenset, taxa: frozenset) -> bool:
    """Two splits can co-exist in one tree iff one of the four intersections is empty."""
    ac, bc = taxa - a, taxa - b
    return not (a & b) or not (a & bc) or not (ac & b) or not (ac & bc)


def robinson_foulds(s1: set, s2: set) -> tuple[int, float]:
    """RF distance and its normalised form |A xor B| / (|A| + |B|)."""
    rf = len(s1 ^ s2)
    denom = len(s1) + len(s2)
    return rf, (rf / denom if denom else 0.0)


def bootstrap_support(
    arr: np.ndarray,
    labels: list[str],
    method: str,
    n_rep: int,
    rng: np.random.Generator,
    model: str = "K2P",
) -> tuple[Tree, dict]:
    """Reference tree plus split support (%) from ``n_rep`` column-resampling replicates."""
    pc = PairCounts(arr)
    ref = build_tree(pc.distances(model=model), labels, method)
    ref_splits = splits(ref)
    counts = dict.fromkeys(ref_splits, 0)
    L = arr.shape[1]
    for _ in range(n_rep):
        w = np.bincount(rng.integers(0, L, L), minlength=L).astype(float)
        rep = splits(build_tree(pc.distances(w, model), labels, method))
        for s in ref_splits & rep:
            counts[s] += 1
    return ref, {s: 100.0 * c / n_rep for s, c in counts.items()}


def annotate(tree: Tree, support: dict, midpoint: bool = True) -> Tree:
    """Midpoint-root (for display) and store split support as clade.confidence."""
    if midpoint:
        tree.root_at_midpoint()
    taxa = frozenset(t.name for t in tree.get_terminals())
    for c in tree.find_clades():
        c.confidence = None
        if not c.is_terminal():
            c.name = None
            s = clade_split(c, taxa)
            if s is not None and s in support:
                c.confidence = round(support[s])
    return tree


def run_bootstrap_tree(
    msa: dict[str, str], method: str, n_rep: int, rng: np.random.Generator, model: str = "K2P"
) -> tuple[Tree, dict]:
    """Convenience wrapper: alignment dict -> annotated tree and split supports."""
    labels = list(msa)
    tree, sup = bootstrap_support(encode(msa), labels, method, n_rep, rng, model)
    return annotate(tree, sup), sup


SupportFn = Callable[[frozenset], float]
