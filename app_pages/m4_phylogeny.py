"""M4: gene trees, bootstrap support and discordance with the reference tree."""

import streamlit as st

from cariogenome.dashboard import (
    csv,
    figure_gallery,
    page_header,
    section,
    tests_table,
    tree_png,
)

gene = page_header(
    "Phylogeny",
    ":material/account_tree:",
    "Do virulence-gene trees disagree with the species tree more often than housekeeping "
    "trees, as horizontal transfer would predict?",
    ["P4b"],
    gene_select=True,
)
assert gene is not None

d = csv("m4_discordance.csv").set_index("gene").loc[gene]
with st.container(horizontal=True):
    st.metric("Taxa", int(d.n_taxa), border=True)
    st.metric("Mean bootstrap support", f"{d.mean_support:.1f}%", border=True)
    st.metric("Unsupported splits (<70%)", int(d.n_unsupported_splits), border=True)
    st.metric("nRF vs reference", f"{d.nRF_vs_reference:.2f}", border=True)
    st.metric(
        "Species-mixing conflicts",
        int(d.n_supported_conflicts_between_species),
        border=True,
        help="Well-supported splits that group strains of different species together.",
    )

left, right = st.columns([1, 1.15], gap="medium")
with left.container(border=True, height="stretch"):
    section(f"{gene} gene tree", ":material/park:")
    method = st.segmented_control(
        "Tree method",
        ["nj", "upgma"],
        default="nj",
        required=True,
        key="method",
        format_func=lambda m: {"nj": "Neighbor joining", "upgma": "UPGMA"}[m],
        label_visibility="collapsed",
    )
    st.image(tree_png(gene, method or "nj"), width="stretch")
    st.caption("Bootstrap support (%) on internal branches; tips colored by species.")
with right.container(border=True, height="stretch"):
    section(
        "Splits that conflict with the reference tree",
        ":material/call_split:",
        "Gene-tree splits with ≥70% bootstrap support that are incompatible with at least "
        "one split of the reference tree.",
    )
    conf = csv("m4_conflicts.csv")
    conf = conf[conf["gene"] == gene].drop(columns=["gene"])
    if conf.empty:
        st.success("No well-supported conflicts for this gene.", icon=":material/check:")
    else:
        st.dataframe(
            conf,
            hide_index=True,
            column_config={
                "gene_support": st.column_config.NumberColumn("Support (%)", format="%.0f"),
                "category": "Conflict",
                "max_conflicting_reference_support": st.column_config.NumberColumn(
                    "Reference support (%)", format="%.0f"
                ),
                "taxa_on_smaller_side": st.column_config.TextColumn("Taxa on smaller side"),
            },
        )

with st.container(border=True):
    section(
        "Virulence vs housekeeping",
        ":material/balance:",
        "Normalized Robinson-Foulds distance saturates near 1 within S. mutans, so these "
        "comparisons have little power.",
    )
    tests_table(csv("m4_tests.csv"))

figure_gallery(
    [
        ("m4_reference_tree", "Reference tree"),
        ("m4_gene_trees_virulence", "Virulence trees"),
        ("m4_gene_trees_controls", "Control trees"),
        ("m4_discordance", "Discordance"),
        ("m4_tanglegrams", "Tanglegrams"),
    ],
    key="m4_figures",
)
