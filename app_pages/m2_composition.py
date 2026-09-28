"""M2: nucleotide composition, codon usage and amino-acid composition."""

import pandas as pd
import streamlit as st

from cariogenome.dashboard import (
    METRIC_LABELS,
    csv,
    figure_gallery,
    page_header,
    section,
    strain_strip_chart,
    tests_table,
)

gene = page_header(
    "Composition",
    ":material/pie_chart:",
    "Do virulence-associated genes carry a compositional signature (GC, GC3, codon "
    "adaptation, GC skew) that would hint at recent horizontal transfer?",
    ["P3"],
    gene_select=True,
)
assert gene is not None

per = csv("m2_per_sequence.csv")
g = per[(per["gene"] == gene) & (per["species"] == "S. mutans")]


def mean(col: str, fmt: str) -> str:
    v = g[col].mean()
    return "–" if pd.isna(v) else format(v, fmt)


with st.container(horizontal=True):
    st.metric("GC content", mean("gc", ".1%"), border=True)
    st.metric("GC3", mean("gc3", ".1%"), border=True)
    st.metric("Codon adaptation index", mean("cai", ".3f"), border=True)
    st.metric("GC skew", mean("gc_skew", "+.3f"), border=True)
st.caption(
    f"Means over {len(g)} S. mutans strains."
    + (" Codon metrics are undefined for the non-coding 16S rRNA." if gene == "16S" else "")
)

with st.container(border=True):
    section(
        "Every S. mutans strain, every gene",
        ":material/scatter_plot:",
        "Dots are strains, ticks are gene medians. The two classes overlap on every metric.",
    )
    metric = st.segmented_control(
        "Metric",
        ["gc", "gc3", "cai", "gc_skew"],
        default="gc3",
        required=True,
        format_func=lambda m: METRIC_LABELS[m],
        key="m2_metric",
        label_visibility="collapsed",
    )
    st.altair_chart(strain_strip_chart(per, metric or "gc3", gene))

with st.container(border=True):
    section(
        "Virulence vs housekeeping in S. mutans",
        ":material/balance:",
        "Gene-level medians compared with a Mann-Whitney U test; BH-adjusted.",
    )
    tests_table(csv("m2_tests.csv"))

with st.expander("Same comparison in the commensal species", icon=":material/public:"):
    other = csv("m2_tests_other_species.csv")
    for species, df in other.groupby("species", sort=False):
        st.markdown(f"**{species}**")
        tests_table(df)

with st.expander(f"{gene}: composition of every sequence", icon=":material/table:"):
    st.dataframe(
        per[per["gene"] == gene].drop(columns=["gene", "class"]),
        hide_index=True,
        column_config={
            "label": "Strain",
            "species": "Species",
            "length": st.column_config.NumberColumn("Length (nt)", format="%d"),
            "gc": st.column_config.NumberColumn("GC", format="percent"),
            "gc_skew": st.column_config.NumberColumn("GC skew", format="%.3f"),
            "gc1": st.column_config.NumberColumn("GC1", format="percent"),
            "gc2": st.column_config.NumberColumn("GC2", format="percent"),
            "gc3": st.column_config.NumberColumn("GC3", format="percent"),
            "cai": st.column_config.NumberColumn("CAI", format="%.3f"),
        },
    )

figure_gallery(
    [
        ("m2_composition_comparison", "Class comparison"),
        ("m2_genome_background", "Genome background"),
        ("m2_species_composition", "Across species"),
        ("m2_aa_composition", "Amino acids"),
        ("m2_codon_usage_rscu", "Codon usage"),
    ],
    key="m2_figures",
)
