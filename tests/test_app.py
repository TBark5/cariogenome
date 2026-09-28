"""Headless dashboard test: every page renders, and the gene pages work for every gene type."""

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(
    not (ROOT / "results" / "m1_catalog.csv").exists(),
    reason="run run_all.py first to create results/",
)

PAGES = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / "app_pages").glob("*.py"))
GENE_PAGES = [p for p in PAGES if p.split("/")[-1][:2] in {"m1", "m2", "m3", "m4", "m5"}]


def test_dashboard_all_pages_and_genes():
    at = AppTest.from_file(str(ROOT / "app.py"), default_timeout=120).run()
    assert not at.exception, [e.value for e in at.exception]
    assert len(PAGES) == 8
    for page in PAGES:
        at.switch_page(page).run()
        assert not at.exception, (page, [e.value for e in at.exception])
    for page in GENE_PAGES:
        at.switch_page(page).run()
        for gene in ["gtfB", "spaP", "recA", "16S"]:  # S. mutans-only, multi-species, control, rRNA
            at.selectbox(key="gene").select(gene).run()
            assert not at.exception, (page, gene, [e.value for e in at.exception])
    at.switch_page("app_pages/m4_phylogeny.py").run()
    at.segmented_control(key="method").set_value("upgma").run()
    assert not at.exception
