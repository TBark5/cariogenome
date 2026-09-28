"""M7: sequence conservation mapped onto the GtfC crystal structure."""

import streamlit as st

from cariogenome.captions import CAPTIONS
from cariogenome.dashboard import csv, figure_gallery, page_header, section, viewer_html
from cariogenome.plotting import figure_path

page_header(
    "Structure",
    ":material/view_in_ar:",
    "Do the most conserved residues of GtfC cluster around its active site in 3D?",
    ["P6"],
)

html = figure_path("m7_structure_3d", "html")
if not html.exists():
    st.info("M7 needs the real GtfC structure and is skipped for synthetic data.")
    st.stop()

s = csv("m7_summary.csv").iloc[0]
with st.container(horizontal=True):
    st.metric("PDB entry", f"{s.pdb_id} · chain {s.chain}", border=True)
    st.metric("Residues scored", f"{int(s.n_residues_scored):,}", border=True)
    st.metric(
        "Spearman ρ, conservation vs distance",
        f"{s.spearman_rho_conservation_vs_distance:.2f}",
        f"CI {s.rho_ci_low:.2f} – {s.rho_ci_high:.2f}",
        delta_color="green",
        delta_arrow="off",
        border=True,
        help="Negative: conservation falls with distance from the catalytic residues.",
    )
    st.metric(
        "Mean conservation within 12 Å",
        f"{s.mean_conservation_within_12A:.2f}",
        f"vs {s.mean_conservation_beyond_12A:.2f} beyond",
        delta_color="green",
        delta_arrow="off",
        border=True,
    )

with st.container(border=True):
    section(
        "Interactive GtfC model",
        ":material/3d_rotation:",
        "Drag to rotate, scroll to zoom. " + CAPTIONS["m7_structure_3d"],
    )
    st.iframe(viewer_html(str(html)), height=660)

with (
    st.expander("Ramachandran statistics", icon=":material/query_stats:"),
    st.container(horizontal=True),
):
    st.metric("α region", f"{s.rama_fraction_alpha:.1%}", border=True)
    st.metric("β region", f"{s.rama_fraction_beta:.1%}", border=True)
    st.metric("Left-handed", f"{s['rama_fraction_left-handed']:.1%}", border=True)
    st.metric("Other", f"{s.rama_fraction_other:.2%}", border=True)

figure_gallery(
    [
        ("m7_structure_conservation_map", "Conservation vs distance"),
        ("m7_structure_3d", "Structure (static)"),
        ("m7_ramachandran", "Ramachandran"),
    ],
    key="m7_figures",
)
