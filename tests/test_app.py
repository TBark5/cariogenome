"""Headless dashboard test: every tab renders for every gene without an exception."""
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest


ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(not (ROOT / "results" / "m1_catalog.csv").exists(),
                                reason="run run_all.py first to create results/")


def test_dashboard_all_genes_and_tabs():
    at = AppTest.from_file(str(ROOT / "app.py"), default_timeout=120).run()
    assert not at.exception, [e.value for e in at.exception]
    assert len(at.tabs) == 8
    for gene in ["gtfB", "spaP", "recA", "16S"]:  # S. mutans-only, multi-species, control, rRNA
        at.selectbox(key="gene").select(gene).run()
        assert not at.exception, (gene, [e.value for e in at.exception])
    at.segmented_control(key="method").set_value("upgma").run()
    assert not at.exception
