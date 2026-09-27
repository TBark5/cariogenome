# Hypotheses (written before any analysis was run)

## Biological question

*Streptococcus mutans* is strongly associated with dental caries. Its close relatives in
the mouth (*S. sanguinis*, *S. gordonii*, *S. mitis*, *S. salivarius*) are mostly
commensal. Several *S. mutans* genes are associated with cariogenicity: the
glucosyltransferases that build sticky extracellular glucan (**gtfB, gtfC, gtfD**), the
surface adhesin antigen I/II (**spaP**), fructosyltransferase (**ftf**), and the
quorum-sensing enzyme **luxS**.

**Do these cariogenicity-associated genes differ from housekeeping genes in sequence
conservation, evolutionary selection pressure, or compositional signature?**

Every comparison is made against a housekeeping control set: **recA, rpoB, gyrB, gyrA,
sodA, pheS, atpD, tuf** (protein-coding) and the **16S rRNA** gene (used for
composition and phylogeny only, because dN/dS needs codons).

The unit of replication is the **gene** (6 virulence vs 8 coding controls). Sequences from
different strains of the same gene are not independent, so gene-level summaries are
compared. With n = 6 vs 8, statistical power is low, and effect sizes with confidence
intervals are the primary result, not p-values.

## Predictions and what would falsify them

**P1. Selection.** Virulence genes will show a higher dN/dS (omega) than housekeeping
controls. Housekeeping genes encode core cellular machinery and should be under strong
purifying selection (omega well below 1). Surface-exposed and secreted virulence proteins
interact with the host and other microbes and may show relaxed constraint or diversifying
selection.
*Falsified if* the median omega of virulence genes is not higher than the control median,
or if the bootstrap 95% CI of the effect size includes zero (then "not supported").

**P2. Conservation.** Virulence proteins will be less conserved than housekeeping proteins
among *S. mutans* strains (lower mean pairwise amino-acid identity, higher mean per-site
Shannon entropy).
*Falsified if* virulence proteins are as conserved as, or more conserved than, the controls.

**P3. Composition.** Virulence genes will have a different compositional signature from
housekeeping genes: lower Codon Adaptation Index (codon usage less like highly expressed
genes) and different GC / GC3 content. Such a signature would be consistent with lower
expression or with acquisition by horizontal gene transfer.
*Falsified if* none of the compositional metrics differ after multiple-testing
correction and the effect-size CIs include zero.

**P4. Distribution and phylogeny.** The glucan- and fructan-related genes (gtfB, gtfC,
gtfD, ftf) and spaP will be restricted to *S. mutans* or highly divergent in commensals,
while luxS and all housekeeping genes will be present in every species. Virulence gene
trees will conflict with the housekeeping species tree more often than control gene trees
do (discordance may suggest horizontal gene transfer or recombination).
*Falsified if* virulence genes are present across species as often as controls, or if
virulence gene trees have no more strongly supported (bootstrap >= 70%) conflicting
branches than control gene trees.

**P5. Functional sites.** In an alignment of the glucansucrase (GH70) family across the
five species, the three catalytic residues of GtfC (nucleophile, acid/base,
transition-state stabilizer) will fall within the 10% most conserved alignment positions.
*Falsified if* any of the three catalytic residues falls outside the top 10%.

**P6. Structural context.** In the GtfC crystal structure, more conserved residues will lie
closer to the catalytic site than less conserved residues (negative Spearman correlation
between conservation score and distance to the active site).
*Falsified if* the correlation is zero or positive, or its 95% CI includes zero.

## Framing

All results are associations between sequence patterns and gene categories. A dN/dS
value or a conservation score is evidence about evolutionary history. It does not
establish a gene's role in disease.
