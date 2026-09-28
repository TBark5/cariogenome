"""M3: alignments, per-column conservation and pairwise identity."""

import streamlit as st

from cariogenome.dashboard import (
    conservation_chart,
    csv,
    figure_gallery,
    identity_heatmap,
    page_header,
    section,
    tests_table,
)

gene = page_header(
    "Conservation",
    ":material/align_horizontal_left:",
    "Are virulence-associated proteins less conserved among S. mutans strains than "
    "housekeeping proteins?",
    ["P2"],
    gene_select=True,
)
assert gene is not None

r = csv("m3_summary.csv").set_index("gene").loc[gene]
with st.container(horizontal=True):
    st.metric("Sequences aligned", int(r.n_sequences), border=True)
    st.metric("Species", int(r.n_species), border=True)
    st.metric("Mean identity, all (%)", f"{r.mean_pid_all:.2f}", border=True)
    st.metric("Mean identity, S. mutans (%)", f"{r.mean_pid_Smutans:.2f}", border=True)
    st.metric("Variable sites, S. mutans", int(r.variable_sites_Smutans), border=True)

with st.container(border=True):
    section(
        f"{gene}: conservation along the alignment",
        ":material/show_chart:",
        "Shannon-entropy conservation of each alignment column (1 = invariant), with a "
        "15-column rolling mean. Hover to read the UA159 residue.",
    )
    cons = csv(f"conservation/{gene}.csv").dropna(subset=["conservation"])
    st.altair_chart(conservation_chart(cons, gene))

left, right = st.columns([1.1, 1], gap="medium")
with left.container(border=True, height="stretch"):
    section("Pairwise percent identity", ":material/grid_view:")
    st.altair_chart(identity_heatmap(csv(f"identity/{gene}.csv", index_col=0)))
with right.container(border=True, height="stretch"):
    section(
        "Most and least conserved 30-residue windows",
        ":material/straighten:",
    )
    reg = csv("m3_extreme_regions.csv")
    st.dataframe(
        reg[reg["gene"] == gene].drop(columns=["gene", "window"]),
        hide_index=True,
        column_config={
            "type": "Region",
            "rank": st.column_config.NumberColumn("Rank", width="small"),
            "start": "Start",
            "end": "End",
            "mean_conservation": st.column_config.ProgressColumn(
                "Conservation", min_value=0.0, max_value=1.0, format="%.3f"
            ),
            "ref_sequence": st.column_config.TextColumn("UA159 sequence"),
        },
    )

with st.container(border=True):
    section(
        "Virulence vs housekeeping in S. mutans",
        ":material/balance:",
        "Gene-level values compared with a Mann-Whitney U test; BH-adjusted.",
    )
    tests_table(csv("m3_tests.csv"))

figure_gallery(
    [
        ("m3_conservation_summary", "Summary"),
        ("m3_conservation_profiles", "Profiles"),
        ("m3_identity_heatmaps", "Identity heatmaps"),
        ("m3_alignment_luxS", "luxS alignment"),
        ("m3_alignment_overview_gtfD", "gtfD overview"),
    ],
    key="m3_figures",
)
