"""CARIOGENOME dashboard: an overview page plus one page per analysis module.

Run with:  .venv\\Scripts\\streamlit run app.py
Hosted on Streamlit Community Cloud with app.py as the entry point.
All numbers shown are read from results/ (produced by run_all.py); nothing is recomputed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent / "src"))

from cariogenome.dashboard import REPO_URL, asset, is_dark, results_ready

st.set_page_config(
    page_title="CARIOGENOME · oral streptococcus comparative genomics",
    page_icon=":material/genetics:",
    layout="wide",
    menu_items={
        "Get help": REPO_URL,
        "Report a bug": f"{REPO_URL}/issues",
        "About": "Comparative genomics of cariogenic and commensal oral streptococci. "
        f"Code and data: {REPO_URL}",
    },
)

mode = "dark" if is_dark() else "light"
st.logo(
    str(asset(f"logo_{mode}.svg")),
    icon_image=str(asset(f"icon_{mode}.svg")),
    link=REPO_URL,
    size="large",
)

if not results_ready():
    st.error("No results found. Run `python run_all.py` first.", icon=":material/error:")
    st.stop()

page = st.navigation(
    {
        "": [
            st.Page(
                "app_pages/overview.py",
                title="Overview",
                icon=":material/space_dashboard:",
                default=True,
            ),
        ],
        "Sequences": [
            st.Page(
                "app_pages/m1_data_qc.py", title="M1 · Data & QC", icon=":material/fact_check:"
            ),
            st.Page(
                "app_pages/m2_composition.py", title="M2 · Composition", icon=":material/pie_chart:"
            ),
            st.Page(
                "app_pages/m3_conservation.py",
                title="M3 · Conservation",
                icon=":material/align_horizontal_left:",
            ),
        ],
        "Evolution": [
            st.Page(
                "app_pages/m4_phylogeny.py", title="M4 · Phylogeny", icon=":material/account_tree:"
            ),
            st.Page(
                "app_pages/m5_selection.py", title="M5 · Selection", icon=":material/trending_up:"
            ),
        ],
        "Protein": [
            st.Page("app_pages/m6_motifs.py", title="M6 · Motifs", icon=":material/pattern:"),
            st.Page(
                "app_pages/m7_structure.py", title="M7 · Structure", icon=":material/view_in_ar:"
            ),
        ],
    },
    position="top",
)
page.run()
