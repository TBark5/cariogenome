"""Unit tests for codon usage, CAI, NG86 dN/dS, alignment, conservation and QC."""
import math
import warnings

import numpy as np
import pytest

from cariogenome.alignment import center_star, codon_alignment, henikoff_weights, identity, pairwise
from cariogenome.codon import CODE, cai, cai_weights, codon_counts, codons_of, rscu
from cariogenome.dnds import PooledNG86, codon_diffs, codon_sites, jc
from cariogenome.homology import ProteinIndex, best_hit
from cariogenome.m1_qc import qc_flags
from cariogenome.m3_conservation import column_conservation


def test_genetic_code_table11():
    assert CODE["ATG"] == "M" and CODE["TAA"] == "*" and CODE["TGG"] == "W"


def test_codons_of_drops_stop_and_partial():
    assert codons_of("ATGAAATAA") == ["ATG", "AAA"]
    assert codons_of("ATGAAAT") == ["ATG", "AAA"]


def test_rscu_equal_use_is_one():
    counts = codon_counts(["TTTTTC"])  # Phe: TTT and TTC once each
    r = rscu(counts)
    assert r["TTT"] == pytest.approx(1.0) and r["TTC"] == pytest.approx(1.0)


def test_cai_is_one_for_preferred_codons_only():
    ref = ["ATG" + "GCT" * 20 + "AAA" * 20 + "TAA"]
    w = cai_weights(ref)
    assert w["GCT"] == pytest.approx(1.0)
    assert cai("ATG" + "GCTAAA" * 10 + "TAA", w) == pytest.approx(1.0)
    assert cai("ATG" + "GCCAAG" * 10 + "TAA", w) < 0.2


def test_ng86_sites_known_codons():
    assert codon_sites("ATG") == (0.0, 3.0)           # Met: no synonymous changes
    assert codon_sites("TTT")[0] == pytest.approx(1 / 3)  # Phe: only TTC is synonymous
    assert codon_sites("GGG")[0] == pytest.approx(1.0)    # Gly: 4-fold degenerate 3rd position


def test_ng86_differences():
    assert codon_diffs("TTT", "TTC") == (1.0, 0.0)
    assert codon_diffs("TTT", "TTA") == (0.0, 1.0)
    s, n = codon_diffs("TTT", "CTC")  # two differences, pathways averaged
    assert s + n == pytest.approx(2.0)


def test_jukes_cantor_formula():
    assert float(jc(0.0)) == 0.0
    assert float(jc(0.1)) == pytest.approx(-0.75 * math.log(1 - 4 * 0.1 / 3))
    assert np.isnan(float(jc(0.8)))


def test_pooled_ng86_matches_biopython():
    """Our NG86 must agree with Bio.codonalign on a gap-free pair."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
    rng = np.random.default_rng(1)
    sense = [c for c, a in CODE.items() if a != "*"]
    a = "".join(rng.choice(sense, 300))
    b = list(a)
    for i in rng.choice(len(a), 25, replace=False):
        cand = b[:]
        cand[i] = "ACGT"[(("ACGT".index(b[i])) + 1) % 4]
        cod = "".join(cand[i - i % 3: i - i % 3 + 3])
        if CODE[cod] != "*":
            b = cand
    b = "".join(b)
    ours = PooledNG86({"a": a, "b": b}, [("a", "b")]).estimate()
    dn, ds = cal_dn_ds(CodonSeq(a), CodonSeq(b), method="NG86")
    assert ours["dN"] == pytest.approx(dn, rel=0.02)
    assert ours["dS"] == pytest.approx(ds, rel=0.02)


def test_pairwise_identity_and_center_star():
    seqs = {"a": "MKVLAAGIVG", "b": "MKVLAGIVG", "c": "MKILAAGIVGK"}
    msa = center_star(seqs, "aa")
    assert len({len(v) for v in msa.values()}) == 1
    assert all(msa[k].replace("-", "") == seqs[k] for k in seqs)
    ga, gb, _ = pairwise("MKVLAAGIVG", "MKVLAAGIVG")
    assert identity(ga, gb) == 100.0


def test_codon_alignment_threads_codons_and_checks_translation():
    msa = {"x": "MK-A", "y": "MKLA"}
    nts = {"x": "ATGAAAGCTTAA", "y": "ATGAAACTTGCT"}
    out = codon_alignment(msa, nts)
    assert out["x"] == "ATGAAA---GCT"
    with pytest.raises(ValueError):
        codon_alignment({"x": "MW"}, {"x": "ATGAAA"})


def test_conservation_invariant_and_uniform_columns():
    aa = "ACDEFGHIKLMNPQRSTVWY"
    msa = {f"s{i}": "M" + aa[i] for i in range(20)}
    cons = column_conservation(msa, "aa", ref=None)
    assert cons.loc[0, "conservation"] == pytest.approx(1.0)
    assert cons.loc[1, "conservation"] == pytest.approx(0.0, abs=1e-9)


def test_henikoff_weights_downweight_duplicates():
    w = henikoff_weights({"a": "AAAA", "b": "AAAA", "c": "CCCC"})
    assert w.sum() == pytest.approx(1.0)
    assert w[2] > w[0]


def test_qc_flags_rules():
    good = "ATG" + "GCT" * 30 + "TAA"
    assert qc_flags(good, "recA", len(good)) == []
    assert "internal_stop" in qc_flags("ATG" + "TAA" + "GCT" * 29 + "TAA", "recA", len(good))
    assert "truncated" in qc_flags("ATG" + "GCT" * 10 + "TAA", "recA", len(good))
    assert "ambiguous_bases" in qc_flags(good[:-6] + "NNNTAA", "recA", len(good))


def test_homology_finds_itself():
    prots = ["MKVLAAGIVGLLLAASSQA", "MSTNPKPQRKTKRNTNRRPQDVKFPGG", "MAHHHHHHSSGLVPRGSHM"]
    hit = best_hit(prots[1], ProteinIndex(prots), top=3)
    assert hit.index == 1 and hit.identity == pytest.approx(1.0)
