"""One-sentence captions for every figure (used by figures/CAPTIONS.md and the dashboard)."""
from __future__ import annotations

from .config import FIGURES

CAPTIONS: dict[str, str] = {
    "summary_effect_sizes": "Forest plot of every virulence-vs-housekeeping comparison in S. mutans (Cliff's delta with 95% bootstrap CI; red = BH q < 0.05), showing that only conservation and dN/dS differ between the gene classes.",
    "m1_presence_absence": "Reciprocal-best-hit orthologs of the 15-locus panel across 22 genomes, showing that gtfB and gtfC have one-to-one orthologs only in S. mutans while luxS and all housekeeping genes are found in every species.",
    "m1_length_qc": "Length of every extracted record relative to its UA159 query, colored by species, with QC-excluded records marked by a cross and the number of retained sequences per gene.",
    "m2_composition_comparison": "GC content, GC3, codon adaptation index and GC skew of 6 virulence-associated vs 8 housekeeping genes in S. mutans (gene means over strains), with Cliff's delta, its 95% CI and BH q-value.",
    "m2_genome_background": "GC3 and CAI of every CDS in each species representative (grey) with the panel genes overlaid, placing virulence and housekeeping genes within their genome's codon-usage distribution.",
    "m2_codon_usage_rscu": "Relative synonymous codon usage per gene in S. mutans (log2 scale) alongside the ribosomal-protein reference set used for CAI.",
    "m2_aa_composition": "Amino-acid composition of each S. mutans panel protein (z-scored per amino acid), with asterisks marking amino acids that differ between gene classes after BH correction.",
    "m2_species_composition": "Mean GC3 and CAI of every panel gene in each species, with absent orthologs labelled.",
    "m3_identity_heatmaps": "Pairwise percent identity (from global pairwise alignments) for every gene, with rows ordered by species and the number of sequences in each panel title.",
    "m3_conservation_profiles": "Henikoff-weighted Shannon-entropy conservation along each virulence protein and two housekeeping controls, with the three most and least conserved windows shaded.",
    "m3_alignment_luxS": "Residue-level view of the center-star protein alignment of LuxS across 22 strains of five species, colored by amino-acid class.",
    "m3_alignment_overview_gtfD": "Overview of the GtfD alignment showing, for every sequence and column, whether the residue matches S. mutans UA159, differs, or is a gap.",
    "m3_conservation_summary": "Mean pairwise amino-acid identity and mean per-site entropy within S. mutans for virulence vs housekeeping genes, showing lower conservation of the virulence-associated proteins.",
    "m4_reference_tree": "Neighbor-joining and UPGMA trees of the concatenated housekeeping genes (the reference species tree) with bootstrap support from 100 replicates; values below 70% are unsupported.",
    "m4_gene_trees_virulence": "Neighbor-joining gene trees of the six virulence-associated genes with bootstrap support (100 replicates; grey italic = unsupported).",
    "m4_gene_trees_controls": "Neighbor-joining gene trees of the eight housekeeping genes and 16S rRNA with bootstrap support (100 replicates; grey italic = unsupported).",
    "m4_discordance": "Normalised Robinson-Foulds distance of each gene tree to the reference tree and the number of well-supported conflicting splits, split into within-species and species-mixing conflicts.",
    "m4_tanglegrams": "Tanglegrams of the most discordant multi-species virulence gene and control gene against the reference tree, where crossing lines indicate taxa placed differently.",
    "m5_dnds_comparison": "Pooled NG86 dN/dS within S. mutans for each gene with 95% codon-bootstrap CIs, and the virulence-vs-control comparison, showing purifying selection in all genes but weaker constraint on virulence genes.",
    "m5_dnds_by_species": "Within-species dN/dS for every gene in every species with at least three strains, showing that virulence-gene homologs in commensal species also have elevated omega.",
    "m5_sliding_window": "Sliding-window (60 codons) nonsynonymous and synonymous diversity and omega within S. mutans, with no window showing omega significantly above 1.",
    "m5_synonymous_saturation": "Synonymous p-distance between S. mutans and each commensal species per gene, showing saturation that makes between-species dN/dS unreliable.",
    "m6_gtf_family_conservation": "Conservation of the 11-member glucansucrase (GH70) family mapped onto GtfC numbering, with the six most conserved motifs and the three catalytic residues marked.",
    "m6_motif_logos": "Information-content sequence logos of the six most conserved GH70 motifs and of the catalytic regions II, III and IV centered on D477, E515 and D588.",
    "m7_structure_conservation_map": "GtfC (PDB 3AIE chain A) C-alpha atoms projected to 2D and colored by GH70 conservation, and conservation vs distance to the catalytic center (Spearman rho with 95% CI).",
    "m7_ramachandran": "Ramachandran plot of GtfC (PDB 3AIE chain A) for general, glycine and proline residues.",
    "m7_structure_conservation": "Interactive py3Dmol view of GtfC (PDB 3AIE chain A) colored by GH70 conservation with catalytic residues shown as sticks (open the HTML file in a browser).",
    "validation_synthetic_recovery": "SYNTHETIC validation: omega estimated by the pipeline vs the true omega used to simulate codon evolution, plus the Robinson-Foulds distance between the recovered and true trees.",
}


def write_captions() -> None:
    """Write figures/CAPTIONS.md listing every figure that exists, in pipeline order."""
    lines = ["# Figure captions", "", "All PNG figures are 300 dpi. Colors use the Okabe-Ito "
             "colorblind-safe palette (red/orange = virulence-associated genes, blue = "
             "housekeeping controls).", ""]
    for name, cap in CAPTIONS.items():
        ext = "html" if name == "m7_structure_conservation" else "png"
        if (FIGURES / f"{name}.{ext}").exists():
            lines.append(f"- **{name}.{ext}**: {cap}")
    (FIGURES / "CAPTIONS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
