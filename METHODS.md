# Methods

This file writes out every formula the pipeline uses, so each number can be traced by
hand. Code references point to `src/cariogenome/`.

## 1. Data and orthology (M1: `m1_retrieval.py`, `homology.py`, `m1_qc.py`)

**Genomes.** 22 complete RefSeq chromosomes: 10 *S. mutans*, 3 each of *S. sanguinis*,
*S. gordonii*, *S. mitis* and *S. salivarius* (`config.yaml`, `ACCESSIONS.md`). Each is
downloaded once with `Bio.Entrez.efetch(rettype="gbwithparts")` with at least 0.5 s
between requests and exponential backoff on HTTP errors.

**Query genes.** Defined by their original locus tag in *S. mutans* UA159 (NC_004350.2).

**Homology search** (a pure-Python replacement for BLAST):
1. *k-mer prefilter*: for every protein in the target genome, count the amino-acid
   3-mers it shares with the query (inverted index), and keep the top 5.
2. *Smith-Waterman local alignment* (BLOSUM62, gap open −11, gap extend −1) of the query
   with each candidate:

   identity = identical aligned pairs / alignment columns (gaps included)

   coverage = aligned span / protein length (computed for both proteins)

3. *Reciprocal best hit (RBH)*: the best hit *h* in genome *X* is an ortholog only if
   searching *h* back against UA159 returns the original query as its best hit, identity
   ≥ 30% and coverage ≥ 60% of both proteins.

**QC exclusion rule.** A record is excluded if it has any non-ACGT base. A coding record
is also excluded if its length is not a multiple of 3, it lacks a start (ATG/GTG/TTG) or
stop codon, it has an internal stop codon (NCBI table 11), it is annotated as a
pseudogene, or it is shorter than 80% of the UA159 query. A 16S record is excluded if it
is shorter than 1400 nt. Records longer than 120% of the query are flagged "extended"
(warning only; see DECISIONS.md).

## 2. Composition (M2: `codon.py`, `m2_composition.py`)

- **GC content** = (G + C) / length. **GC1, GC2, GC3** are the same fraction at codon
  positions 1, 2 and 3.
- **GC skew** = (G − C) / (G + C) on the coding strand.
- **RSCU** (relative synonymous codon usage) of codon *j* for amino acid *i* with *n_i*
  synonymous codons:

  RSCU_ij = X_ij / ( (1/n_i) Σ_k X_ik )

  X_ij is the count of codon *j*. RSCU = 1 means no preference.
- **Codon Adaptation Index** (Sharp & Li 1987). From a reference set of highly expressed
  genes (here, the ribosomal-protein genes of the same genome, with a pseudocount of 0.5
  per codon), the relative adaptiveness of each codon is

  w_ij = RSCU_ij / max_k RSCU_ik

  and for a gene with *L* informative codons (Met, Trp, stops and the start codon
  excluded)

  CAI = exp( (1/L) Σ_l ln w_l )

  CAI is between 0 and 1; values near 1 mean codon usage like the highly expressed
  reference genes.

## 3. Alignment and conservation (M3: `alignment.py`, `m3_conservation.py`)

**Pairwise alignment.** Global Needleman-Wunsch via `Bio.Align.PairwiseAligner`:
proteins with BLOSUM62, gap open −10, extend −0.5; nucleotides with match 2,
mismatch −3, gap open −5, extend −2. Terminal gaps are free.

**Percent identity** = 100 × identical pairs / aligned pairs where neither sequence has a
gap.

**Center-star progressive alignment** (Gusfield 1993):
1. Choose the center *c* that maximizes Σ_j score(c, j) over all pairwise alignments.
2. Align every other sequence to *c*.
3. Merge: for each position of *c*, take the maximum number of residues any sequence
   inserts before it, and pad the others with gaps ("once a gap, always a gap").

This is an approximation. Gaps between two non-center sequences are never optimized.

**Codon alignment.** Each CDS is threaded onto its aligned protein (like PAL2NAL). Each
amino acid becomes its codon and each gap becomes `---`. Every codon is checked against
its amino acid.

**Henikoff sequence weights** (Henikoff & Henikoff 1994). For a column with *k* distinct
residues, where residue *a* appears *n_a* times, a sequence carrying *a* gets
1 / (k · n_a). Weights are summed over columns and normalized to sum to 1, so
near-duplicate sequences share their weight.

**Shannon entropy** of an alignment column, using weighted frequencies *p_a* over non-gap
rows:

  H = − Σ_a p_a log2 p_a (bits)

**Conservation** = 1 − H / log2 K, with K = 20 for proteins and 4 for nucleotides.
1 = invariant column, 0 = all residues equally frequent. Columns with more than 50%
gaps are not scored.

**Most / least conserved regions.** A running mean of conservation over a window of
about 10% of the protein (clamped to 10–30 residues). The top and bottom 3
non-overlapping windows are reported, in UA159 coordinates.

## 4. Phylogenetics (M4: `phylo.py`, `m4_phylogeny.py`)

**p-distance** = mismatches / compared sites (pairwise deletion of gaps).

**Jukes-Cantor (JC69)** correction for multiple substitutions at the same site:

  d = −(3/4) ln(1 − (4/3) p)

**Kimura two-parameter (K2P)** correction. With *P* = proportion of transitions (A↔G,
C↔T) and *Q* = proportion of transversions:

  d = −(1/2) ln(1 − 2P − Q) − (1/4) ln(1 − 2Q)

K2P is used for all trees. It allows transitions and transversions to occur at
different rates.

**Neighbor-joining** (Saitou & Nei 1987). With *n* taxa and distance matrix *d*:
1. For each taxon compute r_i = Σ_k d_ik.
2. Build Q_ij = (n − 2) d_ij − r_i − r_j and join the pair (i, j) with the smallest Q.
3. Branch lengths: d_iu = d_ij/2 + (r_i − r_j) / (2(n − 2)), and d_ju = d_ij − d_iu.
4. Distances to the new node: d_uk = (d_ik + d_jk − d_ij) / 2.
5. Repeat until three nodes remain.

NJ does not assume a molecular clock. The resulting tree is unrooted.

**UPGMA** repeatedly joins the closest pair and averages distances, weighted by cluster
size. It assumes a constant rate (ultrametric tree), so it is shown only for comparison.

**Bootstrap** (Felsenstein 1985). Alignment columns are resampled with replacement
(same length) and the tree is rebuilt, 100 times. The support of a branch is the
percentage of replicate trees that contain the same **split** (bipartition of the taxa).
Support ≥ 70% is called supported.

**Robinson-Foulds distance** (Robinson & Foulds 1981) between trees with split sets A
and B:

  RF = |A Δ B|, normalized RF = |A Δ B| / (|A| + |B|)

**Split compatibility.** Splits X|X' and Y|Y' can coexist in one tree iff at least one of
X∩Y, X∩Y', X'∩Y, X'∩Y' is empty. A gene-tree split is a *conflict* if it is incompatible
with some split of the reference tree (restricted to the gene's taxa).

## 5. Selection (M5: `dnds.py`, `m5_selection.py`)

**Nei-Gojobori (1986) method.**

*Sites.* For a codon, each of its 9 single-nucleotide neighbors is classified as
synonymous or nonsynonymous. A change to a stop codon counts as nonsynonymous (as in
Biopython). For a codon, S = (number of synonymous neighbors) / 3 and N = 3 − S. For a
pair of codons, the sites are the average of the two.

*Differences.* If two codons differ at one position, the change is classified directly.
If they differ at 2 or 3 positions, every mutational pathway (2 or 6 orders) is walked.
The synonymous and nonsynonymous steps are averaged over pathways that do not pass
through a stop codon.

*Proportions and correction.* Summed over codons (and, here, over all sequence pairs):

  p_S = S_d / S, p_N = N_d / N

  d_S = −(3/4) ln(1 − (4/3) p_S), d_N = −(3/4) ln(1 − (4/3) p_N)

  ω = dN/dS = d_N / d_S

**Why a ratio?** Synonymous changes do not alter the protein, so d_S estimates the
neutral mutation rate. It is the internal control for mutation rate and time since
divergence. Dividing by it gives the rate of amino-acid change *relative* to neutral:
ω < 1 purifying selection, ω ≈ 1 neutral, ω > 1 positive selection.

**Pooling.** Within *S. mutans*, all C(10,2) = 45 strain pairs are pooled (sums of S_d,
S, N_d, N) before taking the ratio. This avoids undefined per-pair ratios when a pair has
no synonymous differences.

**Codon bootstrap CI.** Codon columns are resampled with replacement (1000 replicates),
ω is recomputed, and the 2.5th and 97.5th percentiles are reported.

**Saturation.** Between species, p_S is 0.34–0.81, so 37 of 42 comparisons exceed the
0.60 threshold and approach the JC limit of 0.75. There, d_S is not estimable reliably,
so ω is withheld.

**Sliding window.** The same pooled estimate in windows of 60 codons, stepping 20.

## 6. Motifs (M6: `m6_motifs.py`)

**Position frequency matrix** f_{pa} = Henikoff-weighted frequency of amino acid *a* at
position *p*.

**Position weight matrix** (log-odds against a uniform background of 1/20, pseudocount
0.01):

  PWM_{pa} = log2( f'_{pa} / (1/20) )

**Information content** of a position = log2 20 − H_p bits. In a sequence logo, letter
height = f_{pa} × IC_p.

**Conservation percentile of a column** (mid-rank, ties shared) = 100 × (number of
columns with higher conservation + ½ × number with equal conservation) / number of scored
columns.

## 7. Structure (M7: `m7_structure.py`)

- **Mapping.** The PDB chain sequence is globally aligned to UA159 GtfC. Each residue
  gets the GH70 family conservation of its GtfC position, written ×100 into the B-factor
  column for coloring.
- **Catalytic center** = mean coordinate of the side-chain atoms of D477, E515 and D588.
  Distance = Euclidean distance from each residue's Cα atom to that center (Å).
- **Spearman correlation** between conservation and distance, with a 95% CI from 2000
  residue bootstrap resamples.
- **Ramachandran angles.** φ (C−N−Cα−C) and ψ (N−Cα−C−N) dihedrals from
  `Bio.PDB.PPBuilder`. Broad regions: α (φ < 0, −120° < ψ < 50°), β/polyproline (φ < 0,
  ψ ≥ 50° or ψ ≤ −150°), left-handed (φ ≥ 0).

## 8. Statistics (`stats.py`)

- **Unit of replication** = gene (6 virulence vs 8 housekeeping genes). Strain values are
  averaged within a gene first.
- **Mann-Whitney U test**, two-sided (`scipy.stats.mannwhitneyu`).
- **Cliff's delta** (Cliff 1993) = P(X > Y) − P(X < Y) over all virulence/control pairs.
  Range −1 to 1; 0 = no difference. Rank-based, so robust with n = 6 vs 8.
- **Bootstrap CI** of Cliff's delta and of the difference in medians. Each group is
  resampled independently, 5000 times; the 2.5th and 97.5th percentiles are reported.
- **Benjamini-Hochberg** false discovery rate: sort the *m* p-values ascending, set
  q_(i) = min_{j ≥ i} ( m · p_(j) / j ), cap at 1. Applied within each family of tests
  (4 composition metrics, 20 amino acids, the tests of each module, the 6 per-gene dN/dS
  ratios).
- **Seeds.** Every random step uses `numpy.random.default_rng(seed + offset)`, with the
  base seed in `config.yaml`.

## References

- Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.
- Cliff N (1993) Psychol Bull 114:494-509.
- Felsenstein J (1985) Evolution 39:783-791.
- Gusfield D (1993) Bull Math Biol 55:141-154.
- Henikoff S, Henikoff JG (1994) J Mol Biol 243:574-578.
- Ito K et al. (2011) J Mol Biol 408:177-186 (GtfC crystal structures 3AIB/3AIC/3AIE).
- Jukes TH, Cantor CR (1969) In: Mammalian Protein Metabolism, Academic Press.
- Kimura M (1980) J Mol Evol 16:111-120.
- Nei M, Gojobori T (1986) Mol Biol Evol 3:418-426.
- Robinson DF, Foulds LR (1981) Math Biosci 53:131-147.
- Rocha EPC et al. (2006) J Theor Biol 239:226-235.
- Saitou N, Nei M (1987) Mol Biol Evol 4:406-425.
- Sharp PM, Li WH (1987) Nucleic Acids Res 15:1281-1295.
