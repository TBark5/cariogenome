"""M6: conserved motifs and catalytic residues across the GH70 glucansucrase family."""

import pandas as pd
import streamlit as st

from cariogenome.config import RESULTS
from cariogenome.dashboard import csv, figure_gallery, page_header, section

page_header(
    "Motifs",
    ":material/pattern:",
    "Are the catalytic residues of the GtfC glucansucrase among the most conserved columns "
    "of the GH70 family alignment?",
    ["P5"],
)

if not (RESULTS / "m6_catalytic_residues.csv").exists():
    st.info("M6 needs the real GtfC sequence and is skipped for synthetic data.")
    st.stop()

res = csv("m6_catalytic_residues.csv")
summ = csv("m6_summary.csv").iloc[0] if (RESULTS / "m6_summary.csv").exists() else None
with st.container(horizontal=True):
    st.metric("Catalytic residues", len(res), border=True)
    st.metric(
        "Invariant across the family",
        f"{int((res['conservation'] >= 1.0).sum())} / {len(res)}",
        border=True,
    )
    st.metric(
        "Mid-rank percentile",
        f"{res['top_percent_midrank'].iloc[0]:.1f}%",
        "threshold was the top 10%",
        delta_color="red",
        delta_arrow="off",
        border=True,
        help="With ties shared: invariant columns are tied, so each gets the mid rank.",
    )
    st.metric(
        "Columns that are invariant",
        f"{res['top_percent_strict'].iloc[0]:.1f}%",
        border=True,
        help="The share of all alignment columns with conservation 1.",
    )

st.warning(
    "**P5 failed as written.** D477, E515 and D588 are invariant, but so is about a fifth of "
    "the alignment, so with ties shared they rank just outside the pre-registered top 10%. "
    "The rule was fixed before the run and was not changed afterwards.",
    icon=":material/gavel:",
)

left, right = st.columns([1, 1], gap="medium")
with left.container(border=True, height="stretch"):
    section("Catalytic residues of GtfC", ":material/target:")
    st.dataframe(
        res[["position", "residue", "role", "alignment_column", "conservation", "gap_fraction"]],
        hide_index=True,
        column_config={
            "position": "GtfC position",
            "residue": "Residue",
            "role": "Role",
            "alignment_column": "Column",
            "conservation": st.column_config.ProgressColumn(
                "Conservation", min_value=0.0, max_value=1.0, format="%.2f"
            ),
            "gap_fraction": st.column_config.NumberColumn("Gaps", format="percent"),
        },
    )
with right.container(border=True, height="stretch"):
    section("Checked against the published motifs", ":material/verified:")
    st.dataframe(
        csv("m6_catalytic_verification.csv"),
        hide_index=True,
        column_config={
            "position": "Position",
            "expected": "Expected",
            "observed": "Observed",
            "role": "Role",
            "motif": "Motif",
            "in_motif": st.column_config.CheckboxColumn("In motif"),
            "context": st.column_config.TextColumn("Sequence context"),
        },
    )

with st.container(border=True):
    section(
        "Most conserved motifs in the family",
        ":material/format_list_numbered:",
        "Windows of 10 alignment columns with the highest mean conservation.",
    )
    st.dataframe(
        csv("m6_motifs.csv").assign(
            catalytic_residues_inside=lambda d: d["catalytic_residues_inside"].map(
                lambda v: "" if pd.isna(v) else str(v).removesuffix(".0")
            )
        ),
        hide_index=True,
        column_config={
            "motif": "Motif",
            "start": "GtfC start",
            "end": "GtfC end",
            "first_column": "First column",
            "last_column": "Last column",
            "mean_conservation": st.column_config.ProgressColumn(
                "Mean conservation", min_value=0.0, max_value=1.0, format="%.3f"
            ),
            "GtfC_sequence": st.column_config.TextColumn("GtfC"),
            "consensus": st.column_config.TextColumn("Consensus"),
            "catalytic_residues_inside": "Catalytic residues",
        },
    )

figure_gallery(
    [
        ("m6_gtf_family_conservation", "Family conservation"),
        ("m6_motif_logos", "Motif logos"),
    ],
    key="m6_figures",
)
