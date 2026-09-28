"""M5: dN/dS within species, per gene and along each gene."""

import streamlit as st

from cariogenome.dashboard import (
    csv,
    figure_gallery,
    gene_interval_chart,
    page_header,
    section,
    tests_table,
    window_chart,
)

gene = page_header(
    "Selection",
    ":material/trending_up:",
    "Are virulence-associated genes under different selection pressure? dN/dS (ω) is the "
    "ratio of non-synonymous to synonymous divergence among S. mutans strains; ω < 1 means "
    "purifying selection.",
    ["P1"],
    gene_select=True,
)
assert gene is not None

smu = csv("m5_dnds_Smutans.csv").set_index("gene")
if gene in smu.index:
    g = smu.loc[gene]
    vs = csv("m5_gene_vs_controls.csv").set_index("gene")
    with st.container(horizontal=True):
        st.metric("dN/dS within S. mutans", f"{g.omega:.3f}", border=True)
        st.metric("95% CI", f"{g.omega_ci_low:.3f} – {g.omega_ci_high:.3f}", border=True)
        st.metric("pN / pS", f"{g.pN:.4f} / {g.pS:.4f}", border=True)
        if gene in vs.index:
            v = vs.loc[gene]
            st.metric(
                "Ratio to pooled controls",
                f"{v.omega_ratio_to_controls:.1f}×",
                f"CI {v.ratio_ci_low:.1f} – {v.ratio_ci_high:.1f}, q = {v.q_bh:.4f}",
                delta_color="orange",
                delta_arrow="off",
                border=True,
            )
        st.metric("Strains", int(g.n_sequences), border=True)
else:
    st.info(
        "dN/dS is defined for protein-coding genes only; 16S is a non-coding control. "
        "The comparison across genes is shown below.",
        icon=":material/info:",
    )

with st.container(border=True):
    section(
        "dN/dS of every gene within S. mutans",
        ":material/bar_chart:",
        "Pooled NG86 ω with 95% codon-bootstrap CIs. Every gene is well below 1.",
    )
    st.altair_chart(
        gene_interval_chart(
            smu, "omega", "omega_ci_low", "omega_ci_high", gene, "dN/dS (ω)", fmt=".3f"
        )
    )

if gene in smu.index:
    with st.container(border=True):
        section(
            f"{gene}: sliding window along the gene",
            ":material/timeline:",
            "pN and pS in 60-codon windows (step 20) within S. mutans.",
        )
        win = csv("m5_sliding_windows.csv")
        st.altair_chart(window_chart(win[win["gene"] == gene], gene))

with st.container(border=True):
    section(
        "Virulence vs housekeeping in S. mutans",
        ":material/balance:",
        "The difference comes from dN, not dS: relaxed constraint rather than faster mutation.",
    )
    tests_table(csv("m5_tests.csv"))

with st.expander("dN/dS in every species where the gene occurs", icon=":material/public:"):
    within = csv("m5_dnds_within_species.csv")
    st.dataframe(
        within[within["gene"] == gene][
            ["species", "n_sequences", "pN", "pS", "omega", "omega_ci_low", "omega_ci_high"]
        ],
        hide_index=True,
        column_config={
            "species": "Species",
            "n_sequences": "Strains",
            "pN": st.column_config.NumberColumn(format="%.4f"),
            "pS": st.column_config.NumberColumn(format="%.4f"),
            "omega": st.column_config.NumberColumn("ω", format="%.3f"),
            "omega_ci_low": st.column_config.NumberColumn("CI low", format="%.3f"),
            "omega_ci_high": st.column_config.NumberColumn("CI high", format="%.3f"),
        },
    )

with st.expander("Each virulence gene vs the pooled controls", icon=":material/compare:"):
    st.dataframe(
        csv("m5_gene_vs_controls.csv"),
        hide_index=True,
        column_config={
            "gene": "Gene",
            "omega_ratio_to_controls": st.column_config.NumberColumn("ω ratio", format="%.2f"),
            "ratio_ci_low": st.column_config.NumberColumn("CI low", format="%.2f"),
            "ratio_ci_high": st.column_config.NumberColumn("CI high", format="%.2f"),
            "p_bootstrap": st.column_config.NumberColumn("p (bootstrap)", format="%.3f"),
            "pooled_control_omega": st.column_config.NumberColumn(
                "Pooled control ω", format="%.4f"
            ),
            "q_bh": st.column_config.NumberColumn("q (BH)", format="%.4f"),
        },
    )

figure_gallery(
    [
        ("m5_dnds_comparison", "dN/dS comparison"),
        ("m5_dnds_by_species", "By species"),
        ("m5_sliding_window", "Sliding windows"),
        ("m5_synonymous_saturation", "Saturation"),
        ("validation_synthetic_recovery", "Method validation"),
    ],
    key="m5_figures",
)
