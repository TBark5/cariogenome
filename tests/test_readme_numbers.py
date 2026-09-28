"""The numbers quoted in README.md must match the tables in results/."""

import re
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
pytestmark = pytest.mark.skipif(
    not (RES / "summary_effect_sizes.csv").exists(), reason="run run_all.py first"
)


def _fmt(x: float, d: int) -> str:
    # README uses the typographic minus sign, so the test must too.
    return f"{x:.{d}f}".replace("-", "−")


def _readme() -> str:
    return (ROOT / "README.md").read_text(encoding="utf-8")


def test_effect_size_table_matches_results():
    text = _readme()
    eff = pd.read_csv(RES / "summary_effect_sizes.csv")
    for r in eff.itertuples():
        cell = f"{_fmt(r.cliffs_delta, 2)} [{_fmt(r.ci_low, 2)}, {_fmt(r.ci_high, 2)}]"
        assert cell in text, f"README missing effect size for {r.metric}: {cell}"
        q = f"{r.q_bh:.4f}" if r.q_bh < 0.01 else f"{r.q_bh:.3f}"
        assert q in text, f"README missing q for {r.metric}: {q}"


def test_dnds_table_matches_results():
    text = _readme()
    smu = pd.read_csv(RES / "m5_dnds_Smutans.csv")
    for r in smu.rename(columns={"class": "gene_class"}).itertuples():
        lo, hi = max(r.omega_ci_low, 0.0) + 0.0, max(r.omega_ci_high, 0.0) + 0.0  # no "-0.000"
        cls = "virulence" if r.gene_class == "virulence" else "housekeeping"
        row = f"| {r.gene} | {cls} | {r.omega:.3f} [{lo:.3f}, {hi:.3f}] |"
        assert row in text, f"README dN/dS row mismatch: expected {row!r}"


def test_headline_numbers_match_results():
    text = _readme()
    val = pd.read_csv(RES / "validation_synthetic_summary.csv").iloc[0]
    m7 = pd.read_csv(RES / "m7_summary.csv").iloc[0]
    m6 = pd.read_csv(RES / "m6_summary.csv").iloc[0]
    cat = pd.read_csv(RES / "m1_catalog.csv")
    assert f"ρ = {val.omega_spearman:.3f}" in text
    assert f"{100 * val.omega_median_relative_error:.1f}%" in text
    assert (
        f"ρ = {_fmt(m7.spearman_rho_conservation_vs_distance, 2)} "
        f"[{_fmt(m7.rho_ci_low, 2)}, {_fmt(m7.rho_ci_high, 2)}]"
    ) in text
    assert f"{100 * m6.fraction_invariant_columns:.1f}%" in text
    assert f"{int(cat['included'].sum())} of {len(cat)} records pass" in text
    genomes = re.search(r"(\d+) complete genomes", text)
    assert genomes is not None and genomes.group(1) == str(cat["label"].nunique())
