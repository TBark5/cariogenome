"""The pipeline must recover a known tree and known selection pressures from SYNTHETIC data."""

import numpy as np
from scipy.stats import spearmanr

from cariogenome import validation
from cariogenome.synthetic import CODE, simulate_codons, true_tree_newick


def test_simulated_sequences_are_valid_cds():
    seqs = simulate_codons(true_tree_newick(), 120, 0.3, np.random.default_rng(3))
    assert len(seqs) == 22
    for s in seqs.values():
        assert s.startswith("ATG") and s.endswith("TAA") and len(s) % 3 == 0
        assert all(CODE[s[i : i + 3]] != "*" for i in range(0, len(s) - 3, 3))


def test_recovers_true_topology():
    res = validation.recover_topology(8, 300, np.random.default_rng(11))
    assert res["nRF"] == 0.0


def test_recovers_omega_ranking_and_scale():
    om = validation.recover_omega(400, np.random.default_rng(12), reps=2)
    rho = spearmanr(om["true_omega"], om["estimated_omega"]).statistic
    rel = np.median(np.abs(om["estimated_omega"] / om["true_omega"] - 1))
    assert rho > 0.9
    assert rel < 0.25
