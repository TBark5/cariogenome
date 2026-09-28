"""M1: ortholog retrieval, quality control and the presence/absence matrix."""

import streamlit as st

from cariogenome.dashboard import (
    csv,
    figure_gallery,
    page_header,
    presence_heatmap,
    section,
)

gene = page_header(
    "Data & quality control",
    ":material/fact_check:",
    "Which genomes carry an ortholog of each gene, and which records pass quality control? "
    "Orthologs are reciprocal best hits of the S. mutans UA159 query.",
    ["P4a"],
    gene_select=True,
)
assert gene is not None

cat = csv("m1_catalog.csv")
sub = cat[cat["gene"] == gene]
inc = sub[sub["included"]]
with st.container(horizontal=True):
    st.metric("Records found", len(sub), border=True)
    st.metric("Pass QC", f"{len(inc)} / {len(sub)}", border=True)
    st.metric("Species with an ortholog", sub["species"].nunique(), border=True)
    st.metric(
        "Median length (nt)",
        f"{int(inc['nt_length'].median()):,}" if len(inc) else "–",
        border=True,
    )

with st.container(border=True):
    section(
        "Ortholog presence across 22 genomes",
        ":material/grid_on:",
        f"Rows are loci, columns are genomes grouped by species; {gene} is highlighted. "
        "Hover for details.",
    )
    st.altair_chart(presence_heatmap(csv("m1_presence_matrix.csv", index_col=0), gene))

with st.container(border=True):
    section(f"{gene} records", ":material/list_alt:")
    st.dataframe(
        sub[
            [
                "label",
                "species",
                "accession",
                "locus_tag",
                "product",
                "nt_length",
                "identity_to_query",
                "reciprocal_best_hit",
                "flags",
                "included",
            ]
        ],
        hide_index=True,
        column_config={
            "label": "Strain",
            "species": "Species",
            "accession": st.column_config.TextColumn("Accession"),
            "locus_tag": "Locus tag",
            "product": "Product",
            "nt_length": st.column_config.NumberColumn("Length (nt)", format="%d"),
            "identity_to_query": st.column_config.ProgressColumn(
                "Identity to UA159", min_value=0.0, max_value=1.0, format="percent"
            ),
            "reciprocal_best_hit": st.column_config.CheckboxColumn("RBH"),
            "flags": "QC flags",
            "included": st.column_config.CheckboxColumn("Passes QC"),
        },
    )

with st.expander("QC summary for every locus", icon=":material/summarize:"):
    st.dataframe(
        csv("m1_qc_summary.csv"),
        hide_index=True,
        column_config={
            "gene": "Gene",
            "class": "Class",
            "n_found": "Found",
            "n_included": "Pass QC",
            "n_species": "Species",
            "median_nt_length": st.column_config.NumberColumn("Median length", format="%d"),
            "min_nt_length": "Min length",
            "max_nt_length": "Max length",
        },
    )

figure_gallery(
    [
        ("m1_presence_absence", "Presence / absence"),
        ("m1_length_qc", "Length QC"),
    ],
    key="m1_figures",
)
