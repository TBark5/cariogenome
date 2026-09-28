"""Overview: the headline result, the effect-size summary and the prediction scorecard."""

import streamlit as st

from cariogenome.dashboard import (
    ALPHA,
    PREDICTIONS,
    REPO_URL,
    csv,
    data_mode,
    figure,
    forest_chart,
    tests_table,
)

eff = csv("summary_effect_sizes.csv")
tests = eff.set_index("metric")
cat = csv("m1_catalog.csv")
smu = csv("m5_dnds_Smutans.csv")
cat_res = csv("m6_catalytic_residues.csv")
m7 = csv("m7_summary.csv").iloc[0]
om, pid = tests.loc["omega"], tests.loc["mean_pid_Smutans"]
n_genomes, n_species = cat["label"].nunique(), cat["species"].nunique()
n_genes = cat["gene"].nunique()

# ------------------------------------------------------------------------------ hero

st.caption("CARIOGENOME · comparative genomics of cariogenic and commensal oral streptococci")
st.title("Do the genes behind tooth decay evolve differently?")
st.markdown(
    "*Streptococcus mutans* is the bacterium most strongly linked to dental caries; its close "
    "relatives in the mouth are mostly harmless. This project compares "
    f"**6 virulence-associated genes** with **8 housekeeping controls** and 16S rRNA across "
    f"**{n_genomes} complete genomes** from **{n_species} species**, testing seven predictions "
    "that were written down before any data were analysed."
)
with st.container(horizontal=True, gap="small"):
    if data_mode() == "REAL":
        st.badge("Real NCBI RefSeq data", icon=":material/verified:", color="green")
    else:
        st.badge("Synthetic data", icon=":material/warning:", color="red")
    st.badge(f"{n_genomes} genomes", icon=":material/genetics:", color="gray")
    st.badge(f"{n_species} species", icon=":material/bubble_chart:", color="gray")
    st.badge(f"{n_genes} loci", icon=":material/format_list_bulleted:", color="gray")
    st.badge("Pre-registered", icon=":material/lock_clock:", color="gray")
    st.badge("Pure Python, no external binaries", icon=":material/code:", color="gray")

with st.container(horizontal=True, gap="small"):
    st.link_button("Source on GitHub", REPO_URL, icon=":material/code:", type="primary")
    st.link_button("Methods", f"{REPO_URL}/blob/master/METHODS.md", icon=":material/menu_book:")
    st.link_button(
        "Full results",
        f"{REPO_URL}/blob/master/RESULTS_DISCUSSION.md",
        icon=":material/article:",
    )

st.space("small")
st.info(
    f"**Headline.** Virulence-associated genes in *S. mutans* are under **weaker purifying "
    f"selection** than housekeeping genes (median dN/dS {om.median_virulence:.3f} vs "
    f"{om.median_control:.3f}). That is relaxed constraint, **not positive selection**: every "
    f"gene has dN/dS below 1. The same elevation appears in the commensal homologs, so the "
    "signal most likely reflects the biology of secreted and surface proteins, not "
    "cariogenicity itself.",
    icon=":material/lightbulb:",
)

with st.container(horizontal=True):
    st.metric(
        "Median dN/dS, virulence genes",
        f"{om.median_virulence:.3f}",
        f"{om.median_virulence / om.median_control:.1f}× the housekeeping median",
        delta_color="orange",
        delta_arrow="up",
        border=True,
        help="Within S. mutans; Cliff's δ and BH q in the table below.",
    )
    st.metric(
        "Mean protein identity, virulence",
        f"{pid.median_virulence:.2f}%",
        f"vs {pid.median_control:.2f}% for housekeeping",
        delta_color="blue",
        delta_arrow="down",
        border=True,
        help="Median over genes of the mean pairwise amino-acid identity among S. mutans strains.",
    )
    st.metric(
        "Highest dN/dS of any gene",
        f"{smu['omega'].max():.3f}",
        "every gene below 1 → no positive selection",
        delta_color="gray",
        delta_arrow="off",
        border=True,
    )
    st.metric(
        "Records passing QC",
        f"{int(cat['included'].sum())} / {len(cat)}",
        f"{n_genomes} genomes × {n_genes} loci",
        delta_color="gray",
        delta_arrow="off",
        border=True,
    )

# ----------------------------------------------------------------- results summary

comp_q = tests.loc[["gc", "gc3", "cai", "gc_skew"], "q_bh"].min()
disc = tests.loc["n_supported_conflicts"]
p5 = cat_res["top_percent_midrank"].iloc[0]
key_numbers = {
    "P1": f"median ω {om.median_virulence:.3f} vs {om.median_control:.3f}; "
    f"δ = {om.cliffs_delta:.2f} [{om.ci_low:.2f}, {om.ci_high:.2f}], q = {om.q_bh:.3f}",
    "P2": f"identity {pid.median_virulence:.2f}% vs {pid.median_control:.2f}%; "
    f"δ = {pid.cliffs_delta:.2f}, q = {pid.q_bh:.4f}",
    "P3": f"GC, GC3, CAI and GC skew: all q ≥ {comp_q:.2f}, every CI includes 0",
    "P4a": "gtfB and gtfC have one-to-one orthologs only in S. mutans; "
    "gtfD, spaP and ftf occur in 2 to 4 species",
    "P4b": f"δ = {disc.cliffs_delta:.2f} [{disc.ci_low:.2f}, {disc.ci_high:.2f}], "
    f"q = {disc.q_bh:.3f}",
    "P5": f"all {len(cat_res)} catalytic residues invariant, but mid-rank percentile "
    f"{p5:.1f}% (too many invariant columns)",
    "P6": f"Spearman ρ = {m7.spearman_rho_conservation_vs_distance:.2f} "
    f"[{m7.rho_ci_low:.2f}, {m7.rho_ci_high:.2f}] between conservation and distance",
}

left, right = st.columns([1.25, 1], gap="medium")
with left.container(border=True, height="stretch"):
    st.subheader("Effect sizes across every test", icon=":material/insights:")
    st.caption(
        f"Cliff's δ with 95% bootstrap CI, virulence vs housekeeping. Filled color: "
        f"Benjamini-Hochberg q < {ALPHA}. Only conservation and selection differ."
    )
    st.altair_chart(forest_chart(eff), width="stretch")
with right.container(border=True, height="stretch"):
    st.subheader("Prediction scorecard", icon=":material/checklist:")
    st.caption("Written in HYPOTHESIS.md before the analysis; thresholds were never changed.")
    for code, number in key_numbers.items():
        p = PREDICTIONS[code]
        with st.container(gap=None):
            st.markdown(f":{p.color}-badge[{p.icon} {p.verdict}] **{p.code}** · {p.claim}")
            st.caption(number)

with st.expander("All 11 virulence-vs-housekeeping tests", icon=":material/table:"):
    tests_table(eff)
    st.space("small")
    figure("summary_effect_sizes")

# ------------------------------------------------------------------------- caveats

st.subheader("What this project does not show", icon=":material/do_not_disturb_on:")
caveats = [
    (
        ":material/biotech:",
        "A cariogenicity-specific signature",
        "Commensal homologs of the virulence-associated genes also have elevated dN/dS, so "
        "the pattern is better explained by the kind of protein (secreted, surface-exposed) "
        "than by a role in caries.",
    ),
    (
        ":material/account_tree:",
        "Tree discordance within S. mutans",
        "Strain relationships are so poorly resolved that the Robinson-Foulds distance to the "
        "reference is about 1 for most genes. The metric is saturated.",
    ),
    (
        ":material/compare_arrows:",
        "dN/dS between species",
        "Synonymous sites between S. mutans and every commensal are saturated, so selection "
        "is measured only within species.",
    ),
    (
        ":material/medical_services:",
        "Causation or clinical relevance",
        "No gene or variant is shown to cause disease, and no diagnostic or therapeutic "
        "claim is made.",
    ),
]
for row in (caveats[:2], caveats[2:]):
    for col, (icon, head, body) in zip(st.columns(2, border=True), row, strict=True):
        col.markdown(f"{icon} **{head}**")
        col.caption(body)

# --------------------------------------------------------------------- module index

st.subheader("Explore the analysis", icon=":material/explore:")
modules = [
    ("app_pages/m1_data_qc.py", "M1 · Data & QC", "Ortholog calls, QC and the locus panel"),
    ("app_pages/m2_composition.py", "M2 · Composition", "GC, codon usage and amino acids"),
    ("app_pages/m3_conservation.py", "M3 · Conservation", "Alignments, identity and entropy"),
    ("app_pages/m4_phylogeny.py", "M4 · Phylogeny", "Gene trees, support and discordance"),
    ("app_pages/m5_selection.py", "M5 · Selection", "dN/dS within species and along genes"),
    ("app_pages/m6_motifs.py", "M6 · Motifs", "The GH70 glucansucrase active site"),
    ("app_pages/m7_structure.py", "M7 · Structure", "Conservation on the GtfC crystal structure"),
    (REPO_URL, "Code & methods", "The pipeline, tests and write-up on GitHub"),
]
for row in (modules[:4], modules[4:]):
    for col, (path, title, blurb) in zip(st.columns(4, border=True), row, strict=True):
        col.page_link(path, label=f"**{title}**", icon=":material/arrow_forward:")
        col.caption(blurb)

st.space("medium")
st.caption(
    "Code under the MIT License. Sequences are NCBI RefSeq records and the structure is RCSB "
    "PDB entry 3AIE, each under its provider's terms "
    f"([LICENSE-DATA]({REPO_URL}/blob/master/LICENSE-DATA)). "
    f"[Repository]({REPO_URL}) · [Accessions]({REPO_URL}/blob/master/ACCESSIONS.md) · "
    f"[Decisions log]({REPO_URL}/blob/master/DECISIONS.md)"
)
