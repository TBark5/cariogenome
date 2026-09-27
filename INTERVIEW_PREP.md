# Interview preparation

Twenty questions you are likely to be asked, with short answers you can defend. The
numbers come from `results/`.

---

**1. Describe the project in 30 seconds.**
I compared 6 genes linked to tooth decay in *Streptococcus mutans* with 8 housekeeping
genes, across 22 genomes from 5 oral streptococcal species. I measured composition,
conservation, phylogeny and selection (dN/dS). The virulence genes are under *weaker*
purifying selection (median dN/dS 0.156 vs 0.019), but they show no compositional or
phylogenetic sign of horizontal transfer. Their commensal homologs show the same
elevated dN/dS, so it looks like a property of these protein families rather than of
cariogenicity itself.

**2. What does dN/dS measure, and why is it a ratio?**
dN is the rate of nonsynonymous (amino-acid-changing) substitutions per nonsynonymous
site. dS is the rate of synonymous (silent) substitutions per synonymous site. Silent
changes are roughly invisible to selection, so dS estimates the neutral rate: mutation
rate × time. Dividing by it cancels mutation rate and time, leaving how much faster or
slower amino acids change than neutral. ω < 1 means purifying selection removes
protein changes, ω ≈ 1 neutral, ω > 1 positive selection favors change.

**3. Why are housekeeping genes the right control?**
They are in the same genomes and strains, so they share the same population history,
mutation rate, recombination and sampling. They have well-understood, essential
functions under strong constraint. Comparing virulence genes to them asks whether the
virulence genes are *unusual for this genome*. Without a control, "ω = 0.2" has no
reference point, especially for within-species data, where every gene's ω is inflated.

**4. What does a bootstrap value actually mean?**
The alignment columns are resampled with replacement, and the tree is rebuilt 100 times.
A bootstrap value of 90 means 90% of those trees contain that same split. It measures how
consistently *the data* support the branch under this method. It is not the probability
that the branch is true, and it cannot fix a wrong model. I call a branch unsupported
below 70%.

**5. What is the difference between a gene tree and a species tree?**
A gene tree is the history of one locus. A species tree is the history of the organisms.
They can differ because of horizontal gene transfer, recombination, incomplete lineage
sorting, gene duplication and loss (paralogy), or estimation error. I approximated the
species tree by concatenating 8 housekeeping genes. It recovered every species as a clade
with 100% support.

**6. Why would discordance suggest horizontal gene transfer?**
If a gene moved from species A into species B, B's copy will group with A in the gene
tree, even though B groups elsewhere in the species tree. It is only *suggestive*: the
same conflict can come from recombination, paralogs, incomplete lineage sorting or just
poor signal. That is why I only counted conflicts with ≥ 70% bootstrap support. I found
no supported species-mixing conflict in any protein-coding gene; the only one was in 16S
rRNA.

**7. What are the weaknesses of distance-based phylogenetics?**
It compresses the alignment into one number per pair, losing site-by-site information.
It uses simple substitution models (K2P here), with no rate variation among sites. It
cannot compare models statistically. It is sensitive to saturation (UPGMA also assumes
a molecular clock). Maximum likelihood (e.g., IQ-TREE with model selection) or Bayesian
methods use all the data and a better model. The advantage of NJ is speed: 100 bootstraps
for 32 trees take about 20 seconds in pure Python.

**8. What would you do next with more resources?**
- Sample hundreds of *S. mutans* genomes and more commensal strains.
- Build the pangenome (Roary/Panaroo) instead of a fixed panel.
- Align with MAFFT and build ML trees with IQ-TREE.
- Use codon models (PAML codeml or HyPhy) for branch-site and site-level tests of
  positive selection.
- Add a larger random set of genes as controls.
- Detect recombination (ClonalFrameML, Gubbins), since *S. mutans* recombines a lot.
- Link the sequence variation to phenotype data (biofilm formation, acid tolerance)
  before saying anything about function.

**9. Why did you estimate dN/dS within species, not between species?**
Between *S. mutans* and every commensal, the synonymous p-distance is 0.34–0.81 (37 of
42 comparisons above 0.60), close to the 0.75 limit where the Jukes-Cantor correction
breaks down. dS is then not measurable, so ω is meaningless. Within *S. mutans*, pS is
0.008–0.070, so the correction is stable.

**10. Isn't within-species dN/dS biased?**
Yes. Within a species, slightly harmful mutations have not yet been removed by
selection, so ω is inflated compared with long-term divergence (Rocha et al. 2006). That
is why I compare classes measured the same way, and don't read the absolute value as the
strength of selection.

**11. So are virulence genes under positive selection?**
No evidence of it. Every gene has ω < 1, and none of 579 sliding windows has a CI above
1. The difference is driven by dN (more amino-acid polymorphism), with dS similar. The
best description is relaxed purifying selection. It also appears in the commensal
homologs, e.g. *S. gordonii* median 0.138 vs 0.011.

**12. Why Cliff's delta and Mann-Whitney instead of a t-test?**
With 6 vs 8 genes and skewed values (several ω are exactly 0), normality cannot be
assumed. Cliff's delta = P(virulence > control) − P(virulence < control). It is
rank-based, robust and easy to read. I report it with a bootstrap CI because an effect
size with uncertainty says more than a p-value alone.

**13. What does a BH q-value mean?**
When many tests are run, some will be "significant" by chance. Benjamini-Hochberg
controls the false discovery rate: q = 0.012 means that if I call everything with
q ≤ 0.012 a discovery, the expected fraction of false discoveries among them is at most
1.2%. I applied it within each family of tests (4 composition metrics, 20 amino acids,
and so on).

**14. What is the Codon Adaptation Index, and why use ribosomal genes as the reference?**
CAI is the geometric mean of each codon's "relative adaptiveness": how often it is used
compared with the most-used synonymous codon in highly expressed genes. Ribosomal-protein
genes are the classic highly expressed reference set (Sharp & Li 1987). A CAI near 1 means
codon usage like highly expressed genes. The virulence genes' CAI was not different from
the controls' (0.558 vs 0.578, q = 0.41).

**15. How did you identify orthologs without BLAST, and what is the weakness?**
A k-mer prefilter, then Smith-Waterman alignment, then reciprocal best hit (RBH): the best
hit must map back to the original query. RBH is standard and conservative, but with
paralog families it assigns each commensal glucansucrase to only one of *gtfB/C/D*. That
is why gtfB/gtfC look "absent" in commensals even though related GH70 enzymes exist. I
analysed the whole GH70 family separately.

**16. Your multiple alignment isn't MUSCLE or MAFFT. Is it valid?**
It is a center-star progressive alignment: every sequence is aligned to the most central
sequence, and the gaps are merged. It is a known approximation (Gusfield 1993). Gaps
between two non-center sequences are not optimized, which matters most for divergent,
gappy families. I state this as a limitation. The downstream checks (catalytic motifs
aligned and invariant, exact topology recovery on simulated data) suggest it is adequate
here.

**17. Jukes-Cantor vs Kimura two-parameter?**
Both correct the observed p-distance for multiple hits at the same site. JC69 assumes all
substitutions are equally likely: d = −¾ ln(1 − 4p/3). K2P separates transitions (P) from
transversions (Q), which occur at different rates in real DNA:
d = −½ ln(1 − 2P − Q) − ¼ ln(1 − 2Q). I used K2P for trees and JC inside NG86 dN/dS.

**18. One of your predictions failed. Explain.**
P5 predicted the three GtfC catalytic residues would be in the top 10% most conserved
alignment columns. All three are invariant across 11 glucansucrases, the maximum
possible. But 20.4% of columns are also invariant, so their tie-shared rank is 10.2%,
just outside my threshold. The chance of three random columns all being invariant is
0.0085, so the biology is as expected. My percentile test was badly designed for tied
data. I report it as falsified as stated rather than change the rule after seeing the
result.

**19. How do you know the pipeline works?**
Three ways:
- Simulated codon evolution along a known tree with known ω: the pipeline recovers the
  exact topology (RF = 0) and ω (Spearman 0.991, median error 8.8%).
- The NG86 code agrees with Biopython's implementation within 2%.
- 30 unit tests, a clean-clone rebuild that reproduces every result file, and a test
  that checks each README number against `results/`.

**20. So what does this tell us about tooth decay?**
Very little directly, and I would not claim otherwise. It shows that genes associated
with cariogenicity evolve under weaker constraint than housekeeping genes, and that this
is shared with commensal relatives. That is an evolutionary association, not evidence
that any variant causes disease. Linking sequence to cariogenicity would need phenotype
data or experiments.
