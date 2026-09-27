"""Figure registry: pipeline order, file names and one-sentence captions.

The position of a figure in ``CAPTIONS`` defines its file name prefix, so
``figure_filename("m1_presence_absence") == "01_m1_presence_absence"``. The dashboard,
the README, figures/CAPTIONS.md and the tests all use this registry.
"""

from __future__ import annotations

from .config import FIGURES

CAPTIONS: dict[str, str] = {
    "m1_presence_absence": (
        "Reciprocal-best-hit orthologs of the 15-locus panel across 22 genomes, showing that "
        "gtfB and gtfC have one-to-one orthologs only in S. mutans while luxS and all "
        "housekeeping genes are found in every species."
    ),
    "m1_length_qc": (
        "Length of every extracted record relative to its UA159 query, colored by species, "
        "with QC-excluded records marked by a cross and the number of retained sequences per "
        "gene."
    ),
    "m2_composition_comparison": (
        "GC content, GC3, codon adaptation index and GC skew of the virulence-associated and "
        "housekeeping genes in S. mutans (gene means over strains), with Cliff's delta, its "
        "95% CI and the BH q-value."
    ),
    "m2_genome_background": (
        "GC3 and CAI of every CDS in each species representative (grey) with the panel genes "
        "overlaid, placing virulence and housekeeping genes within their genome's codon-usage"
        " distribution."
    ),
    "m2_codon_usage_rscu": (
        "Relative synonymous codon usage of every panel gene in S. mutans (log2 scale) "
        "alongside the ribosomal-protein reference set used for CAI."
    ),
    "m2_aa_composition": (
        "Amino-acid composition of each S. mutans panel protein (z-scored per amino acid), "
        "with asterisks marking amino acids that differ between gene classes after BH "
        "correction."
    ),
    "m2_species_composition": (
        "Mean GC3 and CAI of every panel gene in each species, with absent orthologs labelled."
    ),
    "m3_identity_heatmaps": (
        "Pairwise percent identity (from global pairwise alignments) for every gene, with "
        "rows ordered by species and the number of sequences in each panel title."
    ),
    "m3_conservation_profiles": (
        "Henikoff-weighted Shannon-entropy conservation along each virulence protein and two "
        "housekeeping controls, with the three most and least conserved windows shaded."
    ),
    "m3_alignment_luxS": (
        "Residue-level view of the center-star protein alignment of LuxS across 22 strains of"
        " five species, colored by amino-acid class."
    ),
    "m3_alignment_overview_gtfD": (
        "Overview of the GtfD alignment showing, for every sequence and column, whether the "
        "residue matches S. mutans UA159, differs, or is a gap."
    ),
    "m3_conservation_summary": (
        "Mean pairwise amino-acid identity and mean per-site entropy within S. mutans for "
        "virulence-associated vs housekeeping genes, showing lower conservation of the "
        "virulence-associated proteins."
    ),
    "m4_reference_tree": (
        "Neighbor-joining and UPGMA trees of the concatenated housekeeping genes (the "
        "reference species tree) with bootstrap support from 100 replicates; values below 70%"
        " are unsupported."
    ),
    "m4_gene_trees_virulence": (
        "Neighbor-joining gene trees of the six virulence-associated genes with bootstrap "
        "support (100 replicates; grey italic = unsupported)."
    ),
    "m4_gene_trees_controls": (
        "Neighbor-joining gene trees of the eight housekeeping genes and 16S rRNA with "
        "bootstrap support (100 replicates; grey italic = unsupported)."
    ),
    "m4_discordance": (
        "Normalized Robinson-Foulds distance of each gene tree to the reference tree and the "
        "number of well-supported conflicting splits, split into within-species and species-"
        "mixing conflicts."
    ),
    "m4_tanglegrams": (
        "Tanglegrams of the most discordant multi-species virulence gene and control gene "
        "against the reference tree, where crossing lines indicate taxa placed differently."
    ),
    "m5_dnds_comparison": (
        "Pooled NG86 dN/dS within S. mutans for each gene with 95% codon-bootstrap CIs, and "
        "the virulence-vs-control comparison, showing purifying selection in all genes but "
        "weaker constraint on virulence-associated genes."
    ),
    "m5_dnds_by_species": (
        "Within-species dN/dS for every gene in every species with at least three strains, "
        "showing that virulence-gene homologs in commensal species also have elevated omega."
    ),
    "m5_sliding_window": (
        "Sliding-window (60 codons) nonsynonymous and synonymous diversity and omega within "
        "S. mutans, with no window showing omega significantly above 1."
    ),
    "m5_synonymous_saturation": (
        "Synonymous p-distance between S. mutans and each commensal species per gene, showing"
        " the saturation that makes between-species dN/dS unreliable."
    ),
    "validation_synthetic_recovery": (
        "SYNTHETIC validation: omega estimated by the pipeline vs the true omega used to "
        "simulate codon evolution, plus the Robinson-Foulds distance between the recovered "
        "and true trees."
    ),
    "m6_gtf_family_conservation": (
        "Conservation of the 11-member glucansucrase (GH70) family mapped onto GtfC "
        "numbering, with the six most conserved motifs and the three catalytic residues "
        "marked."
    ),
    "m6_motif_logos": (
        "Information-content sequence logos of the six most conserved GH70 motifs and of the "
        "catalytic regions II, III and IV centered on D477, E515 and D588."
    ),
    "m7_structure_conservation_map": (
        "GtfC (PDB 3AIE chain A) C-alpha atoms projected to 2D and colored by GH70 "
        "conservation, and conservation vs distance to the catalytic center (Spearman rho "
        "with 95% CI)."
    ),
    "m7_structure_3d": (
        "GtfC (PDB 3AIE chain A) C-alpha trace in 3D colored by GH70 conservation with the "
        "catalytic residues marked; the .png is a static rendering and the .html is the "
        "interactive, self-contained py3Dmol view."
    ),
    "m7_ramachandran": (
        "Ramachandran plot of GtfC (PDB 3AIE chain A) for general, glycine and proline residues."
    ),
    "summary_effect_sizes": (
        "Forest plot of every virulence-vs-housekeeping comparison in S. mutans (Cliff's "
        "delta with 95% bootstrap CI; black diamonds = BH q < 0.05), showing that only "
        "conservation and dN/dS differ between the gene classes."
    ),
}
HTML_FIGURES = {"m7_structure_3d"}  # registered figures that also have an interactive .html


def figure_filename(name: str) -> str:
    """``NN_<name>`` with NN = 1-based position of ``name`` in the registry."""
    if name not in CAPTIONS:
        raise KeyError(f"figure {name!r} is not registered in captions.CAPTIONS")
    return f"{list(CAPTIONS).index(name) + 1:02d}_{name}"


def expected_files() -> list[str]:
    """Every figure file the pipeline produces (PNG for all, plus HTML where registered)."""
    out = []
    for name in CAPTIONS:
        out.append(f"{figure_filename(name)}.png")
        if name in HTML_FIGURES:
            out.append(f"{figure_filename(name)}.html")
    return out


def write_captions() -> None:
    """Write figures/CAPTIONS.md: exactly one caption per figure file, in pipeline order."""
    lines = [
        "# Figure captions",
        "",
        "All PNG figures are 300 dpi on a white background. Colors follow one scheme "
        "throughout (Okabe-Ito): vermillion = virulence-associated genes, blue = "
        "housekeeping controls; species are black (S. mutans), sky blue (S. sanguinis), "
        "green (S. gordonii), purple (S. mitis) and orange (S. salivarius).",
        "",
    ]
    for name, cap in CAPTIONS.items():
        stem = figure_filename(name)
        if name in HTML_FIGURES:
            lines.append(f"- **{stem}.png** / **{stem}.html**: {cap}")
        else:
            lines.append(f"- **{stem}.png**: {cap}")
    (FIGURES / "CAPTIONS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
