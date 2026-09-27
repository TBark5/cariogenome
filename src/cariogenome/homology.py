"""Pure-Python homology search used in place of BLAST.

1. A k-mer prefilter ranks all proteins of a genome by the number of amino-acid k-mers
   they share with the query (an inverted index keeps this fast).
2. The top candidates are scored with Smith-Waterman local alignment (BLOSUM62,
   gap open -11, gap extend -1, the BLASTP defaults).
3. Orthologs are called by reciprocal best hit (RBH): the best hit in genome X must
   itself hit the original query best when searched back against the reference genome.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices


@dataclass(frozen=True)
class Hit:
    """Summary of one local alignment between a query and a target protein."""

    index: int          # position of the target in the searched index
    score: float
    identity: float     # identical positions / alignment columns
    qcov: float         # fraction of the query covered by the local alignment
    tcov: float         # fraction of the target covered


def _local_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


LOCAL = _local_aligner()


def kmers(seq: str, k: int = 3) -> set[str]:
    """Set of all overlapping k-mers of ``seq``."""
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


class ProteinIndex:
    """Inverted k-mer index over the proteins of one genome."""

    def __init__(self, seqs: list[str], k: int = 3) -> None:
        self.seqs = seqs
        self.k = k
        inv: dict[str, list[int]] = defaultdict(list)
        for i, s in enumerate(seqs):
            for km in kmers(s, k):
                inv[km].append(i)
        self.inv = {km: np.array(ix, dtype=np.int32) for km, ix in inv.items()}

    def candidates(self, query: str, top: int) -> list[int]:
        """Indices of the ``top`` proteins sharing the most k-mers with ``query``."""
        counts = np.zeros(len(self.seqs), dtype=np.int32)
        for km in kmers(query, self.k):
            ix = self.inv.get(km)
            if ix is not None:
                counts[ix] += 1
        order = np.argsort(-counts, kind="stable")[:top]
        return [int(i) for i in order if counts[i] > 0]


def local_stats(query: str, target: str, index: int = -1) -> Hit:
    """Align two proteins locally and return identity and coverage."""
    aln = LOCAL.align(query, target)[0]
    q_segs, t_segs = aln.aligned
    ident = aligned = gaps = 0
    for (qs, qe), (ts, te) in zip(q_segs, t_segs):
        aligned += qe - qs
        ident += sum(a == b for a, b in zip(query[qs:qe], target[ts:te]))
    for j in range(1, len(q_segs)):
        gaps += (q_segs[j][0] - q_segs[j - 1][1]) + (t_segs[j][0] - t_segs[j - 1][1])
    length = aligned + gaps
    if length == 0:
        return Hit(index, 0.0, 0.0, 0.0, 0.0)
    qcov = (q_segs[-1][1] - q_segs[0][0]) / len(query)
    tcov = (t_segs[-1][1] - t_segs[0][0]) / len(target)
    return Hit(index, float(aln.score), ident / length, qcov, tcov)


def best_hit(query: str, index: ProteinIndex, top: int = 5) -> Hit | None:
    """Best Smith-Waterman hit among the k-mer candidates, or None if no candidate."""
    cands = index.candidates(query, top)
    if not cands:
        return None
    scores = [LOCAL.score(query, index.seqs[i]) for i in cands]
    best = cands[int(np.argmax(scores))]
    return local_stats(query, index.seqs[best], best)


def all_hits(query: str, index: ProteinIndex, top: int, min_identity: float,
             min_coverage: float) -> list[Hit]:
    """Every candidate passing identity and (query and target) coverage thresholds."""
    hits = [local_stats(query, index.seqs[i], i) for i in index.candidates(query, top)]
    return [h for h in hits if h.identity >= min_identity
            and h.qcov >= min_coverage and h.tcov >= min_coverage]


def reciprocal_best_hit(query: str, query_index_in_ref: int, target_index: ProteinIndex,
                        ref_index: ProteinIndex, top: int = 5) -> tuple[Hit | None, bool]:
    """Forward best hit of ``query`` in the target genome, and whether it is reciprocal."""
    fwd = best_hit(query, target_index, top)
    if fwd is None:
        return None, False
    back = best_hit(target_index.seqs[fwd.index], ref_index, top)
    return fwd, back is not None and back.index == query_index_in_ref
