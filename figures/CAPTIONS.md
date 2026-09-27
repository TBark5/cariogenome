# Figure captions

All PNG figures are 300 dpi on a white background. Colors follow one scheme throughout (Okabe-Ito): vermillion = virulence-associated genes, blue = housekeeping controls; species are black (S. mutans), sky blue (S. sanguinis), green (S. gordonii), purple (S. mitis) and orange (S. salivarius).

- **01_m1_presence_absence.png**: Reciprocal-best-hit orthologs of the 15-locus panel across 22 genomes, showing that gtfB and gtfC have one-to-one orthologs only in S. mutans while luxS and all housekeeping genes are found in every species.
- **02_m1_length_qc.png**: Length of every extracted record relative to its UA159 query, colored by species, with QC-excluded records marked by a cross and the number of retained sequences per gene.
- **03_m2_composition_comparison.png**: GC content, GC3, codon adaptation index and GC skew of the virulence-associated and housekeeping genes in S. mutans (gene means over strains), with Cliff's delta, its 95% CI and the BH q-value.
- **04_m2_genome_background.png**: GC3 and CAI of every CDS in each species representative (grey) with the panel genes overlaid, placing virulence and housekeeping genes within their genome's codon-usage distribution.
- **05_m2_codon_usage_rscu.png**: Relative synonymous codon usage of every panel gene in S. mutans (log2 scale) alongside the ribosomal-protein reference set used for CAI.
- **06_m2_aa_composition.png**: Amino-acid composition of each S. mutans panel protein (z-scored per amino acid), with asterisks marking amino acids that differ between gene classes after BH correction.
- **07_m2_species_composition.png**: Mean GC3 and CAI of every panel gene in each species, with absent orthologs labelled.
- **08_m3_identity_heatmaps.png**: Pairwise percent identity (from global pairwise alignments) for every gene, with rows ordered by species and the number of sequences in each panel title.
- **09_m3_conservation_profiles.png**: Henikoff-weighted Shannon-entropy conservation along each virulence protein and two housekeeping controls, with the three most and least conserved windows shaded.
- **10_m3_alignment_luxS.png**: Residue-level view of the center-star protein alignment of LuxS across 22 strains of five species, colored by amino-acid class.
- **11_m3_alignment_overview_gtfD.png**: Overview of the GtfD alignment showing, for every sequence and column, whether the residue matches S. mutans UA159, differs, or is a gap.
- **12_m3_conservation_summary.png**: Mean pairwise amino-acid identity and mean per-site entropy within S. mutans for virulence-associated vs housekeeping genes, showing lower conservation of the virulence-associated proteins.
- **13_m4_reference_tree.png**: Neighbor-joining and UPGMA trees of the concatenated housekeeping genes (the reference species tree) with bootstrap support from 100 replicates; values below 70% are unsupported.
- **14_m4_gene_trees_virulence.png**: Neighbor-joining gene trees of the six virulence-associated genes with bootstrap support (100 replicates; grey italic = unsupported).
- **15_m4_gene_trees_controls.png**: Neighbor-joining gene trees of the eight housekeeping genes and 16S rRNA with bootstrap support (100 replicates; grey italic = unsupported).
- **16_m4_discordance.png**: Normalized Robinson-Foulds distance of each gene tree to the reference tree and the number of well-supported conflicting splits, split into within-species and species-mixing conflicts.
- **17_m4_tanglegrams.png**: Tanglegrams of the most discordant multi-species virulence gene and control gene against the reference tree, where crossing lines indicate taxa placed differently.
- **18_m5_dnds_comparison.png**: Pooled NG86 dN/dS within S. mutans for each gene with 95% codon-bootstrap CIs, and the virulence-vs-control comparison, showing purifying selection in all genes but weaker constraint on virulence-associated genes.
- **19_m5_dnds_by_species.png**: Within-species dN/dS for every gene in every species with at least three strains, showing that virulence-gene homologs in commensal species also have elevated omega.
- **20_m5_sliding_window.png**: Sliding-window (60 codons) nonsynonymous and synonymous diversity and omega within S. mutans, with no window showing omega significantly above 1.
- **21_m5_synonymous_saturation.png**: Synonymous p-distance between S. mutans and each commensal species per gene, showing the saturation that makes between-species dN/dS unreliable.
- **22_validation_synthetic_recovery.png**: SYNTHETIC validation: omega estimated by the pipeline vs the true omega used to simulate codon evolution, plus the Robinson-Foulds distance between the recovered and true trees.
- **23_m6_gtf_family_conservation.png**: Conservation of the 11-member glucansucrase (GH70) family mapped onto GtfC numbering, with the six most conserved motifs and the three catalytic residues marked.
- **24_m6_motif_logos.png**: Information-content sequence logos of the six most conserved GH70 motifs and of the catalytic regions II, III and IV centered on D477, E515 and D588.
- **25_m7_structure_conservation_map.png**: GtfC (PDB 3AIE chain A) C-alpha atoms projected to 2D and colored by GH70 conservation, and conservation vs distance to the catalytic center (Spearman rho with 95% CI).
- **26_m7_structure_3d.png** / **26_m7_structure_3d.html**: GtfC (PDB 3AIE chain A) C-alpha trace in 3D colored by GH70 conservation with the catalytic residues marked; the .png is a static rendering and the .html is the interactive, self-contained py3Dmol view.
- **27_m7_ramachandran.png**: Ramachandran plot of GtfC (PDB 3AIE chain A) for general, glycine and proline residues.
- **28_summary_effect_sizes.png**: Forest plot of every virulence-vs-housekeeping comparison in S. mutans (Cliff's delta with 95% bootstrap CI; black diamonds = BH q < 0.05), showing that only conservation and dN/dS differ between the gene classes.
