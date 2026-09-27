# Decisions log

Choices made where the brief left room for interpretation. Each entry says what was
chosen and why.

## Setup

- **Python 3.14 in `.venv`.** The existing PyCharm venv uses Python 3.14.5, and every
  dependency has a wheel for it. Versions are pinned in `requirements.txt`.
- **Removed the PyCharm template `main.py`.** It was the default "Hi, PyCharm" sample and
  not part of the project.
- **Minimal dependencies.** Only biopython, numpy, pandas, scipy, matplotlib, PyYAML,
  py3Dmol, streamlit and pytest. Sequence logos and trees are drawn with matplotlib
  instead of adding extra plotting libraries. Benjamini-Hochberg correction is implemented
  directly (5 lines) instead of pulling in statsmodels.
- **NCBI contact email.** `config.yaml` uses the email already configured for git on this
  machine. Replace it with a reachable address before re-downloading.

## Gene panel

- **Housekeeping controls extended from 4 to 8 coding genes + 16S.** The brief listed
  recA, rpoB, gyrB and 16S as examples. With only 3 coding controls, any virulence vs
  control test would have almost no power, so gyrA, sodA, pheS, atpD and tuf were added.
  These are standard single-copy housekeeping or MLST-type loci in streptococci.
- **Genes are defined by their UA159 locus tag** (`config.yaml`), not by annotation
  names, because gene names are inconsistent between genome annotations (for example,
  `ftf` and `luxS` carry no gene name in the current UA159 RefSeq record).

## M1: retrieval and QC

- **Whole chromosomes, not gene-name searches.** One `efetch` (GenBank with sequence)
  per genome, 22 genomes, is faster and more reliable than hundreds of gene-level queries.
  The download is cached in `data/raw/` (about 176 MB, git-ignored). The small per-gene
  FASTA and GenBank slices derived from it are committed in `data/genes/` and
  `data/genbank/`, so the analysis can be rerun without the network. The first download
  took about 2 minutes; afterwards M1 runs from cache.
- **Orthology by reciprocal best hit (RBH)** with a pure-Python k-mer prefilter plus
  Smith-Waterman (BLOSUM62, gap -11/-1). RBH requires identity >= 30% and coverage
  >= 60% of both proteins. Annotation names are recorded but not used for the call.
- **Paralogs are resolved by RBH.** *S. sanguinis* and *S. gordonii* each have one GH70
  glucansucrase (GtfP / GtfG), which is the reciprocal best hit of UA159 `gtfD`, not
  `gtfB`/`gtfC`. These species are therefore "absent" for `gtfB`/`gtfC`. This is the
  standard orthology definition, but it means "absent" = "no one-to-one ortholog", not
  "no related enzyme". The whole GH70 family is analysed separately in M6.
- **One 16S copy per genome** (the first by coordinate). Streptococcal 16S copies are
  nearly identical within a genome.
- **Species representative** = the first genome of each species in `config.yaml`
  (UA159, SK36, CH1, B6, JIM8777, the classic reference strains). Representatives supply
  the genome-wide CDS background (M2) and the GH70 family set (M6).
- **QC rule changed after first inspection (post hoc, stated openly).** The rule written
  before seeing data excluded records longer than 120% of the UA159 query. The first run
  showed this excluded all three *S. salivarius* `ftf` records (122% of UA159 length),
  which are complete CDSs with valid start and stop codons. The extra length is a known
  biological feature of levansucrases, not a data defect. "Extended" is now a warning
  only. The remaining exclusions (an ambiguous base in LP13 `gtfB`, a missing start codon
  in SK637 `recA`) are unchanged.

## M2: composition

- **Unit of replication = gene.** Per-strain values are averaged per gene within
  *S. mutans* (9-10 strains), giving 6 virulence vs 8 control values. Treating strains as
  independent observations would inflate the sample size (pseudo-replication).
- **CAI reference set = ribosomal-protein genes of the same genome** (the classic Sharp &
  Li choice of highly expressed genes), with a 0.5 pseudocount for codons absent from the
  reference. CAI is computed with our own implementation (see `codon.py`, METHODS.md).
- **16S rRNA is excluded from the coding comparisons** (codon metrics are undefined and
  rRNA GC is not comparable with CDS GC). Its GC is reported in `m2_per_sequence.csv`.
- **Per-gene GC skew is reported but interpreted cautiously.** Within a gene, GC skew mainly
  reflects whether the gene lies on the leading or lagging replication strand, not the
  gene's function.
- **Multiple testing:** BH across the 4 primary metrics, and separately across the
  20 amino acids.
- **Other species** are tested the same way (`m2_tests_other_species.csv`) as a secondary
  analysis. They have only 1-3 virulence orthologs each, so power is very low.

## M3: alignment and conservation

- **Center-star progressive MSA** (Gusfield) built from Biopython global pairwise
  alignments, since MUSCLE/MAFFT are not allowed. It is an approximation: gaps between
  two non-center sequences are not optimised. Terminal gaps are free.
- **Percent identity** = identical pairs / aligned non-gap pairs (gaps ignored), from true
  pairwise global alignments (not read off the MSA).
- **Conservation** = 1 - H/log2(K) with Henikoff position-based sequence weights, so ten
  near-identical *S. mutans* strains do not outweigh one sequence from another species.
  Columns with > 50% gaps are not scored.
- **P2 test uses within-*S. mutans* values** (mean pairwise identity, mean entropy),
  because only this taxon set is shared by every gene. Cross-species identity is
  confounded by which species carry the gene (gtfB/gtfC are *S. mutans*-only).
- **Region window** = about 10% of the protein length, clamped to 10-30 residues (30 for
  the large proteins, 16 for LuxS). Most- and least-conserved windows may not overlap.

## M5: selection

- **Own NG86 implementation**, cross-checked against `Bio.codonalign.cal_dn_ds`
  (identical on luxS/recA pairs, within about 1% on gtfD/rpoB; see the test suite).
  `Bio.codonalign` is labelled experimental by Biopython and does not provide pooled
  estimates or a codon bootstrap, so it is used only as a check.
- **Primary estimate = within *S. mutans*, pooled over all strain pairs** (45 pairs for
  10 strains): pN = sum Nd / sum N, pS = sum Sd / sum S, then JC correction. Pooling avoids
  undefined per-pair ratios when a pair has no synonymous differences.
- **Between-species dN/dS is not used for inference.** The first run showed synonymous
  p-distances of 0.66-0.76 between *S. mutans* and every commensal (37 of 42 comparisons
  above the 0.60 threshold; 0.75 is the JC limit). These are reported and flagged
  "saturated", with omega withheld.
- **Interpretation caveat:** within-species omega measures polymorphism. Slightly
  deleterious mutations inflate it relative to between-species divergence (Rocha et al.
  2006), so omega is compared between gene classes measured the same way, not read as an
  absolute measure of selection.
- **CIs by codon bootstrap** (1000 replicates, codon columns resampled). This captures site
  sampling, not phylogenetic non-independence of strain pairs. Genes with zero
  nonsynonymous differences (recA, tuf) get a degenerate CI of [0, 0].
- **Per-gene test:** ratio of each virulence gene's omega to the pooled control omega
  (all 8 control genes' codons concatenated), with a bootstrap p-value (floor 1/1000) and
  BH across the 6 virulence genes.
- **Sliding window:** 60 codons, step 20, 200 bootstrap replicates per window.
- **Synthetic validation runs every time** (not only in fallback mode): recovery of a
  known 22-taxon topology and of known omega values (0.02-1.0) by the same code.

## M6: motifs and catalytic residues

- **GH70 family set** = every GH70 homolog of GtfB/C/D (identity >= 30%, coverage >= 50%
  of both proteins) in the five species representatives: 11 proteins (UA159 GtfB/C/D,
  SK36 GtfP, CH1 GtfG, six JIM8777 enzymes; *S. mitis* B6 has none). Using all 22 genomes
  would add near-duplicate strains and no new diversity.
- **Catalytic residues are assigned from the canonical GH70 motifs** (D477 nucleophile in
  region II, E515 acid/base in region III, D588 transition-state stabilizer in region IV)
  and the pipeline verifies each motif in UA159 GtfC. UniProt P13470 has no residue-level
  active-site annotation. *Correction made in the polish pass:* earlier documents credited
  these three residue labels directly to Ito et al. 2011. Only the paper's abstract could
  be checked, and it names the adjacent subsite residues (Asn481, Trp517, Ser589), which
  confirms the numbering but does not itself list D477/E515/D588. The wording now says
  exactly that. No number or result changed.
- **Tie handling for P5 was fixed in code before the first M6 run:** the percentile of a
  column is its mid-rank (columns with equal conservation share the average rank).
  Strict and best-case ranks are also written to `m6_catalytic_residues.csv`.
- **Motifs** = the 6 non-overlapping 10-residue windows with the highest mean
  conservation (entropy approach), plus 11-column windows centered on each catalytic
  residue. PWMs are log-odds (bits) against a uniform background with a 0.01 pseudocount.

## M7: structure

- **PDB 3AIE chain A** (GtfC residues 244-1087, 2.1 A, the best resolution of the three
  GtfC entries). 3AIC (acarbose complex, 3.1 A) would mark the active site with a ligand,
  but lower resolution makes the Ramachandran plot noisier. The active site is defined
  instead by the centroid of the catalytic side chains.
- **Conservation is written to the B-factor column** (x100; -1 = not scored), so py3Dmol
  colors the cartoon by conservation. The interactive view is
  `figures/26_m7_structure_3d.html` (self-contained since the polish pass; see below).
- **P6 statistic:** Spearman correlation between residue conservation and C-alpha distance
  to the catalytic center, with a residue bootstrap CI. Neighbouring residues are not
  independent, so the p-value and CI are optimistic. The direction and size of rho are
  the result, not its p-value.
- **Ramachandran regions** are simple phi/psi boxes for description only. They are not a
  MolProbity-style validation.

## Visual pass, app, quality, docs

- **Synthesis forest plot** (`summary.py`) collects every virulence-vs-control test on one
  axis (Cliff's delta with 95% CI). This is the README headline figure.
- **Per-species replication of the dN/dS pattern is exploratory.** It was added after the
  by-species dN/dS figure showed elevated omega in commensal homologs. It was not in
  `HYPOTHESIS.md` and is labelled as exploratory in RESULTS_DISCUSSION.md.
- **Dashboard** is a single Streamlit page with 8 tabs. It only reads `results/` and
  `figures/` and does no heavy computation. The py3Dmol view is embedded with `st.iframe`
  from our own generated HTML file.
- **Reproducibility:** fixed seeds everywhere. Table rows that come from Python sets are
  sorted, gzip files carry no timestamp, and the 3D viewer id is fixed. Two consecutive
  runs, and a clean clone with a new venv running `--offline`, reproduce every results
  file byte for byte (except `results/runtimes.tsv`).
- **Rounding:** README and RESULTS_DISCUSSION numbers are Python-formatted from the CSVs
  (for example, 0.625 prints as 0.62). `tests/test_readme_numbers.py` enforces this.
- **Runtime cap:** every step runs in under 2 minutes; the full run takes about 3 minutes
  (per-step timings are written to `results/runtimes.tsv` on every run). The one-time
  genome download is network-bound (about 106 s for 22 genomes). No bootstrap counts had
  to be reduced below the brief's minimums.
