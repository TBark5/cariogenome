"""Genetic code (NCBI table 11), codon usage, RSCU and the Codon Adaptation Index."""
from __future__ import annotations

import math
from collections import Counter

from Bio.Seq import Seq

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODE: dict[str, str] = {c: str(Seq(c).translate(table=11)) for c in CODONS}
STOP_CODONS = {c for c, aa in CODE.items() if aa == "*"}
SENSE_CODONS = [c for c in CODONS if c not in STOP_CODONS]
SYNONYMS: dict[str, list[str]] = {}
for _c in SENSE_CODONS:
    SYNONYMS.setdefault(CODE[_c], []).append(_c)
# Met and Trp have a single codon and carry no information about codon choice.
INFORMATIVE_AA = sorted(aa for aa, cs in SYNONYMS.items() if len(cs) > 1)


def codons_of(seq: str) -> list[str]:
    """Split a CDS into codons, dropping a trailing stop and any incomplete/ambiguous codon."""
    cods = [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    if cods and cods[-1] in STOP_CODONS:
        cods = cods[:-1]
    return [c for c in cods if c in CODE]


def codon_counts(seqs: list[str]) -> Counter:
    """Pooled codon counts over several CDS (start codon included, stop excluded)."""
    counts: Counter = Counter()
    for s in seqs:
        counts.update(codons_of(s))
    return counts


def rscu(counts: Counter, pseudocount: float = 0.0) -> dict[str, float]:
    """Relative synonymous codon usage: observed / expected under equal synonymous use."""
    out = {}
    for aa, cods in SYNONYMS.items():
        vals = [counts.get(c, 0) + pseudocount for c in cods]
        total = sum(vals)
        for c, v in zip(cods, vals):
            out[c] = v * len(cods) / total if total > 0 else float("nan")
    return out


def cai_weights(reference: list[str], pseudocount: float = 0.5) -> dict[str, float]:
    """Relative adaptiveness w = RSCU / max RSCU of its amino acid (Sharp & Li 1987).

    A pseudocount of 0.5 keeps codons absent from the reference set from receiving w = 0,
    which would force CAI to 0.
    """
    r = rscu(codon_counts(reference), pseudocount)
    w = {}
    for aa in INFORMATIVE_AA:
        mx = max(r[c] for c in SYNONYMS[aa])
        for c in SYNONYMS[aa]:
            w[c] = r[c] / mx
    return w


def cai(seq: str, weights: dict[str, float]) -> float:
    """Geometric mean of codon weights over informative codons (skips the start codon)."""
    logs = [math.log(weights[c]) for c in codons_of(seq)[1:] if c in weights]
    return math.exp(sum(logs) / len(logs)) if logs else float("nan")
