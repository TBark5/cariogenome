"""Unit tests for distances, trees, splits, bootstrap support and statistics."""

import math

import numpy as np
import pytest

from cariogenome import phylo
from cariogenome.m4_phylogeny import tree_from_string
from cariogenome.stats import benjamini_hochberg, cliffs_delta, compare_groups


def test_kimura_reduces_to_known_values():
    assert float(phylo.kimura_2p(0.0, 0.0)) == 0.0
    P, Q = 0.1, 0.05
    expected = -0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q)
    assert float(phylo.kimura_2p(P, Q)) == pytest.approx(expected)


def test_paircounts_transitions_and_transversions():
    arr = phylo.encode({"a": "AAAA", "b": "GAAT"})  # 1 transition (A-G), 1 transversion (A-T)
    pc = phylo.PairCounts(arr)
    assert pc.ts.sum() == 1 and pc.tv.sum() == 1 and pc.valid.sum() == 4


def test_splits_and_rf():
    t1 = tree_from_string("((A,B),(C,D),E);")
    t2 = tree_from_string("((A,C),(B,D),E);")
    assert phylo.robinson_foulds(phylo.splits(t1), phylo.splits(t1)) == (0, 0.0)
    assert phylo.robinson_foulds(phylo.splits(t1), phylo.splits(t2))[1] == 1.0


def test_compatibility():
    taxa = frozenset("ABCDE")
    assert phylo.compatible(frozenset("AB"), frozenset("ABC"), taxa)
    assert not phylo.compatible(frozenset("AB"), frozenset("AC"), taxa)


def test_nj_recovers_additive_tree():
    # Distances from the tree ((A:1,B:2):1,(C:1,D:3):1); are additive, NJ must recover AB|CD.
    labels = ["A", "B", "C", "D"]
    d = np.array([[0, 3, 4, 6], [3, 0, 5, 7], [4, 5, 0, 4], [6, 7, 4, 0]], float)
    tree = phylo.build_tree(d, labels, "nj")
    assert frozenset({"C", "D"}) in phylo.splits(tree) or frozenset({"A", "B"}) in phylo.splits(
        tree
    )


def test_bootstrap_support_is_high_for_clear_signal():
    rng = np.random.default_rng(0)
    base = rng.choice(list("ACGT"), 600)
    seqs = {}
    for name, flips in [("A", 0), ("B", 5), ("C", 90), ("D", 95)]:
        s = base.copy()
        idx = (
            rng.choice(600, flips, replace=False)
            if name in "AB"
            else np.r_[np.arange(60), rng.choice(np.arange(60, 600), flips - 60, replace=False)]
        )
        s[idx] = [{"A": "C", "C": "G", "G": "T", "T": "A"}[x] for x in s[idx]]
        seqs[name] = "".join(s)
    _, sup = phylo.run_bootstrap_tree(seqs, "nj", 50, np.random.default_rng(1))
    assert max(sup.values()) >= 95


def test_benjamini_hochberg_known_values():
    q = benjamini_hochberg([0.01, 0.04, 0.03, 0.5])
    assert q == pytest.approx([0.04, 0.0533333, 0.0533333, 0.5], rel=1e-4)


def test_cliffs_delta_extremes():
    assert cliffs_delta(np.array([5, 6]), np.array([1, 2])) == 1.0
    assert cliffs_delta(np.array([1, 2]), np.array([1, 2])) == 0.0


def test_compare_groups_reports_effect_and_ci():
    res = compare_groups(
        np.array([5, 6, 7.0]), np.array([1, 2, 3.0]), np.random.default_rng(0), 500
    )
    assert res["cliffs_delta"] == 1.0 and res["delta_ci_low"] <= 1.0 and res["n_control"] == 3
