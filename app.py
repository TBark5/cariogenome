"""CARIOGENOME dashboard: one tab per analysis module, with a gene selector.

Run with:  .venv\\Scripts\\streamlit run app.py
All numbers shown are read from results/ (produced by run_all.py).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent / "src"))

from cariogenome.captions import CAPTIONS
from cariogenome.config import RESULTS, all_genes, gene_class
from cariogenome.m4_figures import draw_tree
from cariogenome.m4_phylogeny import load_tree
from cariogenome.plotting import W_SINGLE, figure_path

st.set_page_config(page_title="CARIOGENOME", page_icon=":material/biotech:", layout="wide")


@st.cache_data(max_entries=200)
def csv(name: str, index_col: int | None = None) -> pd.DataFrame:
    """Read one results table (cached)."""
    return pd.read_csv(RESULTS / name, index_col=index_col)


def figure(name: str) -> None:
    """Show a registered figure (figures/NN_<name>.png) with its caption."""
    path = figure_path(name)
    if path.exists():
        st.image(str(path), caption=CAPTIONS[name], width="stretch")
    else:
        st.info(f"Figure {path.name} not found; run `python run_all.py`.")


if not (RESULTS / "m1_catalog.csv").exists():
    st.error("No results found. Run `python run_all.py` first.")
    st.stop()

mode = (RESULTS / "data_mode.txt").read_text().strip()
st.title("CARIOGENOME")
st.caption("Comparative genomics of cariogenic and commensal oral streptococci: are genes "
           "associated with cariogenicity different from housekeeping genes in conservation, "
           "selection or composition?")
if mode != "REAL":
    st.warning("These results were produced from SYNTHETIC simulated sequences, not NCBI data.",
               icon=":material/warning:")

with st.sidebar:
    st.header("Controls")
    gene = st.selectbox("Gene", all_genes(), index=0, key="gene",
                        format_func=lambda g: f"{g} ({'rRNA control' if g == '16S' else gene_class(g)})")
    st.caption("The gene selector drives the M1-M5 tabs. M6 and M7 analyse the glucansucrase "
               "(GH70) family and the GtfC structure.")
    st.caption(f"Data mode: **{mode}**")

tabs = st.tabs(["Overview", "M1 Data & QC", "M2 Composition", "M3 Conservation",
                "M4 Phylogeny", "M5 Selection", "M6 Motifs", "M7 Structure"])

with tabs[0]:
    eff = csv("summary_effect_sizes.csv")
    cat = csv("m1_catalog.csv")
    with st.container(horizontal=True):
        st.metric("Genomes", cat["label"].nunique(), border=True)
        st.metric("Species", cat["species"].nunique(), border=True)
        st.metric("Records passing QC", f"{int(cat['included'].sum())}/{len(cat)}", border=True)
        om = eff.set_index("metric").loc["omega"]
        st.metric("Median dN/dS, virulence vs control",
                  f"{om.median_virulence:.3f} vs {om.median_control:.3f}", border=True)
    figure("summary_effect_sizes")
    st.subheader("All virulence-vs-control tests")
    st.dataframe(eff.drop(columns=["metric"]).round(4), hide_index=True)

with tabs[1]:
    cat = csv("m1_catalog.csv")
    sub = cat[cat["gene"] == gene]
    st.subheader(f"{gene}: {int(sub['included'].sum())} of {len(sub)} records pass QC")
    st.dataframe(sub[["label", "species", "accession", "locus_tag", "product", "nt_length",
                      "identity_to_query", "reciprocal_best_hit", "flags", "included"]],
                 hide_index=True)
    figure("m1_presence_absence")
    figure("m1_length_qc")
    st.dataframe(csv("m1_qc_summary.csv"), hide_index=True)

with tabs[2]:
    per = csv("m2_per_sequence.csv")
    st.subheader(f"{gene}: composition per strain")
    st.dataframe(per[per["gene"] == gene].drop(columns=["gene", "class"]).round(4), hide_index=True)
    st.subheader("Virulence vs housekeeping (S. mutans, gene level)")
    st.dataframe(csv("m2_tests.csv").round(4), hide_index=True)
    for name in ["m2_composition_comparison", "m2_genome_background", "m2_species_composition",
                 "m2_aa_composition", "m2_codon_usage_rscu"]:
        figure(name)

with tabs[3]:
    summ = csv("m3_summary.csv").set_index("gene")
    r = summ.loc[gene]
    with st.container(horizontal=True):
        st.metric("Sequences", int(r.n_sequences), border=True)
        st.metric("Species", int(r.n_species), border=True)
        st.metric("Mean identity, all (%)", f"{r.mean_pid_all:.1f}", border=True)
        st.metric("Mean identity, S. mutans (%)", f"{r.mean_pid_Smutans:.2f}", border=True)
    cons = csv(f"conservation/{gene}.csv").dropna(subset=["conservation"])
    st.subheader("Per-column conservation (1 - H / log2 K)")
    st.line_chart(cons.set_index("column")["conservation"], x_label="alignment column",
                  y_label="conservation")
    reg = csv("m3_extreme_regions.csv")
    st.dataframe(reg[reg["gene"] == gene], hide_index=True)
    st.subheader("Pairwise percent identity")
    st.dataframe(csv(f"identity/{gene}.csv", index_col=0).style.background_gradient(cmap="viridis"))
    st.dataframe(csv("m3_tests.csv").round(4), hide_index=True)
    for name in ["m3_conservation_summary", "m3_conservation_profiles", "m3_identity_heatmaps",
                 "m3_alignment_luxS", "m3_alignment_overview_gtfD"]:
        figure(name)

with tabs[4]:
    disc = csv("m4_discordance.csv").set_index("gene")
    d = disc.loc[gene]
    with st.container(horizontal=True):
        st.metric("Taxa", int(d.n_taxa), border=True)
        st.metric("Mean bootstrap support (%)", f"{d.mean_support:.1f}", border=True)
        st.metric("Unsupported splits (<70%)", int(d.n_unsupported_splits), border=True)
        st.metric("nRF vs reference", f"{d.nRF_vs_reference:.2f}", border=True)
        st.metric("Supported conflicting splits", int(d.n_supported_conflicts), border=True)
    method = st.segmented_control("Tree method", ["nj", "upgma"], default="nj", key="method",
                                  format_func=str.upper)
    fig, ax = plt.subplots(figsize=(W_SINGLE, 0.32 * d.n_taxa + 1.5))
    draw_tree(ax, load_tree(f"{gene}_{method or 'nj'}"), f"{gene} {(method or 'nj').upper()} tree")
    st.pyplot(fig, width="content")
    plt.close(fig)
    conf = csv("m4_conflicts.csv")
    st.dataframe(conf[conf["gene"] == gene], hide_index=True)
    st.dataframe(csv("m4_tests.csv").round(4), hide_index=True)
    for name in ["m4_reference_tree", "m4_gene_trees_virulence", "m4_gene_trees_controls",
                 "m4_discordance", "m4_tanglegrams"]:
        figure(name)

with tabs[5]:
    smu = csv("m5_dnds_Smutans.csv").set_index("gene")
    if gene in smu.index:
        g = smu.loc[gene]
        with st.container(horizontal=True):
            st.metric("dN/dS within S. mutans", f"{g.omega:.3f}", border=True)
            st.metric("95% CI", f"[{g.omega_ci_low:.3f}, {g.omega_ci_high:.3f}]", border=True)
            st.metric("pN / pS", f"{g.pN:.4f} / {g.pS:.4f}", border=True)
            st.metric("Strains", int(g.n_sequences), border=True)
        win = csv("m5_sliding_windows.csv")
        w = win[win["gene"] == gene].set_index("mid_codon")[["pN", "pS"]]
        st.subheader("Sliding window (60 codons) within S. mutans")
        st.line_chart(w, x_label="codon position", y_label="p-distance")
    else:
        st.info("dN/dS is defined for protein-coding genes only (16S is a non-coding control).")
    within = csv("m5_dnds_within_species.csv")
    st.dataframe(within[within["gene"] == gene].round(4), hide_index=True)
    st.dataframe(csv("m5_tests.csv").round(4), hide_index=True)
    st.dataframe(csv("m5_gene_vs_controls.csv").round(4), hide_index=True)
    for name in ["m5_dnds_comparison", "m5_dnds_by_species", "m5_sliding_window",
                 "m5_synonymous_saturation", "validation_synthetic_recovery"]:
        figure(name)

with tabs[6]:
    if (RESULTS / "m6_catalytic_residues.csv").exists():
        st.subheader("Catalytic residues of GtfC in the GH70 family alignment")
        st.dataframe(csv("m6_catalytic_residues.csv").round(3), hide_index=True)
        st.dataframe(csv("m6_catalytic_verification.csv"), hide_index=True)
        st.subheader("Most conserved motifs")
        st.dataframe(csv("m6_motifs.csv").round(3), hide_index=True)
        figure("m6_gtf_family_conservation")
        figure("m6_motif_logos")
    else:
        st.info("M6 needs the real GtfC sequence and is skipped for synthetic data.")

with tabs[7]:
    html = figure_path("m7_structure_3d", "html")
    if html.exists():
        summ7 = csv("m7_summary.csv").iloc[0]
        with st.container(horizontal=True):
            st.metric("PDB entry", f"{summ7.pdb_id} chain {summ7.chain}", border=True)
            st.metric("Residues scored", int(summ7.n_residues_scored), border=True)
            st.metric("Spearman rho (conservation vs distance)",
                      f"{summ7.spearman_rho_conservation_vs_distance:.2f}", border=True)
        st.caption("Interactive view (self-contained; works offline). "
                   + CAPTIONS["m7_structure_3d"])
        st.iframe(html, height=720)
        figure("m7_structure_3d")
        figure("m7_structure_conservation_map")
        figure("m7_ramachandran")
    else:
        st.info("M7 needs the real GtfC structure and is skipped for synthetic data.")
