# Results and discussion

Every number below is read from a file in `results/` (named in brackets). The
predictions were written in `HYPOTHESIS.md` before any analysis was run; the git history
shows the order. "Supported" means the effect was in the predicted direction with a 95%
CI excluding zero and a Benjamini-Hochberg q < 0.05. Sequences were **real NCBI data**
(`results/data_mode.txt` = REAL).

## Summary

| Prediction | Verdict | Key number |
|---|---|---|
| P1 Virulence genes have higher dN/dS | **Supported** (relaxed constraint, not positive selection) | median ω 0.156 vs 0.019; Cliff's δ = 0.88 [0.50, 1.00], q = 0.012 |
| P2 Virulence proteins are less conserved | **Supported** | identity 99.06% vs 99.90%; δ = −0.96 [−1.00, −0.75], q = 0.0055 |
| P3 Compositional signature differs | **Not supported** | all 4 metrics q ≥ 0.33, all CIs include 0 |
| P4a gtf/ftf/spaP restricted to *S. mutans*; luxS and controls universal | **Partly supported** | gtfB, gtfC one-to-one orthologs only in *S. mutans*; gtfD, spaP, ftf present in 2–4 species |
| P4b Virulence gene trees more discordant | **Not supported** | δ = 0.42 [−0.12, 0.88], q = 0.366; no coding gene has a supported species-mixing conflict |
| P5 Catalytic residues in the top 10% most conserved columns | **Narrowly falsified as stated** | all 3 invariant, but mid-rank percentile 10.2% (20.4% of columns are invariant) |
| P6 Conserved residues cluster near the active site | **Supported** | Spearman ρ = −0.30 [−0.36, −0.23] |

## P1. Selection (M5)

[`m5_tests.csv`, `m5_dnds_Smutans.csv`, `m5_gene_vs_controls.csv`]

Within *S. mutans* (10 strains, 45 strain pairs per gene; 9 strains and 36 pairs for
gtfB), every gene has ω < 1, so **all genes are under purifying selection**. Virulence
genes have a higher ω than housekeeping genes:

| Gene | Class | ω | 95% CI |
|---|---|---|---|
| gtfB | virulence | 0.095 | 0.067–0.134 |
| gtfC | virulence | 0.081 | 0.056–0.117 |
| gtfD | virulence | 0.085 | 0.050–0.128 |
| spaP | virulence | 0.218 | 0.148–0.303 |
| ftf | virulence | 0.247 | 0.127–0.446 |
| luxS | virulence | 0.235 | 0.041–1.011 |
| recA | control | 0.000 | 0.000–0.000 |
| rpoB | control | 0.007 | 0.000–0.020 |
| gyrB | control | 0.009 | 0.000–0.026 |
| gyrA | control | 0.065 | 0.017–0.153 |
| sodA | control | 0.171 | 0.037–0.703 |
| pheS | control | 0.029 | 0.011–0.057 |
| atpD | control | 0.047 | 0.000–0.229 |
| tuf | control | 0.000 | 0.000–0.000 |

- Class comparison: median 0.156 vs 0.019, Cliff's δ = 0.88 [0.50, 1.00], q = 0.012.
- The difference is **driven by dN** (δ = 0.92 [0.62, 1.00], q = 0.012), **not dS**
  (δ = 0.42 [−0.21, 0.92], q = 0.228). Virulence genes accumulate more amino-acid
  polymorphism; they do not have unusual synonymous diversity.
- Each virulence gene individually has a higher ω than the pooled controls (pooled
  control ω = 0.024). Ratios range from 3.4 (gtfC) to 10.4 (ftf), all BH q ≤ 0.018.
- **No evidence of positive selection.** No gene has ω > 1, and none of 579 sliding
  windows has a CI entirely above 1 [`m5_sliding_windows.csv`]. 11 windows have a point
  estimate above 1, but their CIs include 1.

**Replication in commensals** [`m5_replication_by_species.csv`]. The same direction
appears in every other species with at least two virulence homologs:
- *S. gordonii*: median 0.138 vs 0.011
- *S. salivarius*: 0.173 vs 0.009
- *S. sanguinis*: 0.065 vs 0.023

None of these is significant (n = 2–3 genes; p = 0.13–0.89). The consistent direction
suggests that elevated ω is a property of these secreted and surface-protein families in
general, **not a signature specific to the cariogenic species**.

**Interpretation.** P1 as written (higher ω) is supported. The pattern is best described
as *weaker purifying selection* (relaxed constraint) on the virulence-associated genes.
Diversifying selection would require ω > 1, which is not observed. Several caveats apply:
- Within-species ω measures polymorphism. It is inflated by slightly deleterious variants
  that selection has not yet removed (Rocha et al. 2006), so the absolute values should
  not be read as long-term selection coefficients. Both classes are measured the same way,
  so the comparison is fair.
- recA and tuf have zero nonsynonymous polymorphisms among the 10 strains, so their CIs
  are degenerate. Their true ω is small but not exactly zero.
- luxS is classed as virulence-associated, but it is a core metabolic enzyme in every
  species. Its ω CI is very wide (5 variable amino-acid sites).
- Between-species dN/dS could not be used. Synonymous sites are saturated in 37 of 42
  *S. mutans*-vs-commensal comparisons (pS 0.34–0.81) [`m5_dnds_between_species.csv`].

## P2. Conservation (M3)

[`m3_tests.csv`, `m3_summary.csv`]

Within *S. mutans*, virulence proteins are less conserved:
- Mean pairwise amino-acid identity: 99.06% vs 99.90% (δ = −0.96 [−1.00, −0.75],
  q = 0.0055).
- Mean per-site entropy: 0.0193 vs 0.0021 bits (δ = 0.92 [0.62, 1.00], q = 0.0055).

**Supported.** This agrees with P1 (more amino-acid polymorphism). The absolute
differences are small: all *S. mutans* proteins are more than 98% identical between
strains.

Across species, where orthologs exist, the least conserved regions of GtfD are the
N-terminal residues 53–147 (conservation 0.56–0.67) [`m3_extreme_regions.csv`]. Two of
the three most conserved GtfD windows contain the catalytic motifs GVRVDAVDNV (452–481)
and ILEAW (498–527). This matches the expectation that catalytic cores are conserved and
N-terminal variable regions are not.

## P3. Composition (M2)

[`m2_tests.csv`, `m2_aa_tests.csv`]

| Metric | Virulence median | Control median | Cliff's δ [95% CI] | q |
|---|---|---|---|---|
| GC | 0.397 | 0.404 | −0.46 [−1.00, 0.17] | 0.362 |
| GC3 | 0.299 | 0.267 | 0.38 [−0.25, 0.92] | 0.377 |
| CAI | 0.558 | 0.578 | −0.29 [−0.83, 0.33] | 0.414 |
| GC skew | 0.055 | 0.124 | −0.58 [−1.00, 0.00] | 0.325 |

**Not supported.** No compositional metric differs after correction, and every CI
includes zero. With 6 vs 8 genes, only very large effects are detectable, so this is
"no evidence of a difference", not evidence of no difference. The genome-background plot
shows that the virulence genes sit inside the normal range of their genome's codon usage.
There is no strong atypical signature of recent horizontal acquisition.

Two amino acids differ after BH correction across 20 tests. Glutamate is lower in
virulence proteins (4.7% vs 8.5%, q = 0.027) and threonine is higher (8.2% vs 6.0%,
q = 0.027). Thr/Ser enrichment is typical of secreted and cell-surface proteins. This is
a description of protein class, not evidence about cariogenicity.

GC skew within genes mostly reflects leading- vs lagging-strand position, so it is not
interpreted further.

## P4. Distribution and phylogeny (M1, M4)

[`m1_presence_matrix.csv`, `m4_discordance.csv`, `m4_tests.csv`, `m4_conflicts.csv`]

**Distribution (partly supported).**
- gtfB and gtfC have reciprocal-best-hit orthologs only in *S. mutans*, as predicted.
- gtfD has orthologs in *S. sanguinis* (2 of 3 strains), *S. gordonii* and
  *S. salivarius*. These are the single-copy GH70 enzymes of those species.
- spaP has orthologs in *S. sanguinis* and *S. gordonii* (the antigen I/II family).
- ftf has an ortholog in *S. salivarius*.
- luxS and all 8 housekeeping genes are found in every species.

So the prediction holds for gtfB/C but not for gtfD, spaP or ftf. "Absent" means "no
one-to-one ortholog": *S. sanguinis* and *S. gordonii* do carry a related GH70 enzyme,
which RBH assigns to gtfD.

**Discordance (not supported).**
- The reference tree (concatenated housekeeping genes) recovers every species as a clade
  with 100% bootstrap support [`m4_reference_supports.csv`].
- Virulence gene trees have a median of 4 supported conflicting splits vs 3 for controls
  (δ = 0.42 [−0.12, 0.88], q = 0.366).
- **No protein-coding gene tree has a supported conflict that mixes species.** The only
  such conflict is in 16S rRNA (*S. mitis* grouping with *S. sanguinis*, 94% support),
  where 16S is known to resolve the mitis group poorly.
- In the pheS tree, *S. sanguinis* and *S. gordonii* are not separately monophyletic, but
  the branches involved are unsupported (< 70%).

Within *S. mutans*, gene trees almost never match the reference topology (normalized RF
≈ 1 for 12 of 14 coding genes). Strain relationships are poorly resolved and differ
between genes, as expected with frequent homologous recombination and very low
divergence. This saturated metric cannot separate the gene classes. It is a limitation of
the design, not a finding about virulence genes.

In short, with this sampling there is **no phylogenetic signal of horizontal transfer of
the virulence genes between these species**. Discordance could suggest transfer or
recombination, but it would not prove it, and here there is little to explain.

## P5. Catalytic residues in the GH70 family (M6)

[`m6_catalytic_residues.csv`, `m6_summary.csv`, `m6_catalytic_verification.csv`]

The family alignment has 11 glucansucrases from 4 species (*S. mitis* B6 has none) and
1545 scored columns.
- The three catalytic residues of GtfC (Ito et al. 2011) are D477 (nucleophile,
  region II, `SIRVDAVDNV`), E515 (acid/base, region III, `ILEAW`) and D588
  (transition-state stabilizer, region IV, `FIRAHD`). All three are verified to sit in
  their canonical GH70 motifs.
- **All three are invariant** across the 11 enzymes (conservation = 1.0).
- If 20.4% of columns are invariant, the chance that three random columns are all
  invariant is 0.204³ = 0.0085.
- However, because so many columns are invariant, the mid-rank percentile of each
  residue is 10.2% of the most conserved columns. That is just outside the pre-registered
  "top 10%" threshold.

**Verdict: narrowly falsified as operationalized, although the biological expectation is
met.** The residues reach the maximum possible conservation, but a percentile threshold
was the wrong test when a fifth of the columns are tied at the maximum. The tie rule
(mid-rank) was fixed in code before the first run. We report the failure rather than
change the rule after seeing the result.

The region II motif containing D477 is also one of the six most conserved 10-residue
windows (motif M5, residues 476–485) [`m6_motifs.csv`].

## P6. Structural context (M7)

[`m7_summary.csv`, `m7_residue_conservation.csv`, `m7_ramachandran.csv`]

On the GtfC crystal structure (PDB 3AIE chain A, residues 244–1087, 839 scored residues),
conservation decreases with distance from the catalytic center: Spearman ρ = −0.30
[−0.36, −0.23]. The 34 residues within 12 Å of the center have mean conservation 0.954,
compared with 0.767 elsewhere.

**Supported.** Caveat: neighboring residues are not independent, so the bootstrap CI is
optimistic. The direction and size of ρ are the result, not the p-value.

The Ramachandran plot (842 residues) places 51.0% in the broad α region, 42.3% in the
β/polyproline region, and 6.5% with φ > 0, mostly glycines. These are simple descriptive
boxes, not a validation score.

## Method validation (synthetic data)

[`validation_synthetic_summary.csv`]

On sequences simulated along a known 22-taxon tree with known ω:
- The K2P + NJ pipeline recovers the true topology exactly (normalized RF = 0.0, mean
  support 100%).
- Pooled NG86 recovers ω from 0.02 to 1.0 with Spearman ρ = 0.991 and a median relative
  error of 8.8%.

The unit tests also check our NG86 against Biopython's `cal_dn_ds` (agreement within 2%).

## What this means

The cariogenicity-associated genes of *S. mutans* differ from housekeeping genes mainly
in **how strongly they are constrained**. They carry more amino-acid polymorphism
(higher dN, higher ω, lower identity), but they are still under purifying selection.
They do **not** differ detectably in codon usage or GC, and their gene trees show no
supported evidence of inter-species transfer. The same elevated ω appears in commensal
homologs, so it is probably a general feature of secreted, surface-exposed proteins
rather than a hallmark of cariogenicity. The catalytic core of the glucansucrases is
invariant and spatially clustered with the most conserved residues, as expected for an
enzyme active site.

These are associations between sequence patterns and gene categories. They do not show
that any gene causes caries or that any difference in selection is functionally
important.

## Limitations

- Small gene panel (6 vs 8 genes), so only large class differences are detectable.
- Distance-based trees (NJ/UPGMA with K2P), not maximum likelihood or Bayesian.
- Approximate center-star alignment; no MUSCLE/MAFFT.
- Limited strain sampling (10 *S. mutans*, 3 per commensal); within-species ω
  reflects polymorphism.
- The virulence set is heterogeneous: luxS is a metabolic enzyme found in all species.
- Sequence signatures do not establish function or a role in disease.
