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
