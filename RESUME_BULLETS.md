# Resume bullets

## Short (one line)

- Built a pure-Python comparative-genomics pipeline (22 streptococcal genomes, 7 analysis
  modules, Streamlit dashboard) showing that caries-associated genes evolve under weaker
  purifying selection than housekeeping genes (dN/dS 0.156 vs 0.019).

## Medium (2-3 bullets)

- Designed and built CARIOGENOME, a reproducible Python pipeline comparing 6
  caries-associated and 8 housekeeping genes across 22 complete genomes of 5 oral
  *Streptococcus* species, with predictions written before analysis.
- Found that virulence-associated genes are under weaker purifying selection
  (dN/dS 0.156 vs 0.019; Cliff's δ = 0.88, FDR q = 0.012) and are less conserved, with no
  compositional or phylogenetic evidence of horizontal transfer. Reported a failed
  prediction openly.
- Delivered 43 automated tests, synthetic-data validation (exact tree recovery,
  ω Spearman 0.99), a clean-clone reproducibility check, and an interactive Streamlit
  dashboard with 3D protein structure views.

## Technical (for bioinformatics roles)

- Implemented ortholog calling (k-mer prefilter + Smith-Waterman + reciprocal best hit),
  center-star progressive multiple alignment, codon alignment, and Henikoff-weighted
  Shannon entropy in Python/Biopython, with no external binaries.
- Built K2P/JC distance phylogenetics (neighbor-joining, UPGMA) with split-based bootstrap
  support (100 replicates), Robinson-Foulds discordance and split-compatibility tests
  against a concatenated housekeeping species tree.
- Implemented pooled Nei-Gojobori dN/dS with codon bootstrap CIs, a synonymous-saturation
  check and sliding-window analysis. Validated it against Biopython (within 2%) and on
  simulated codon evolution with known ω.
- Mapped GH70 family conservation onto the *S. mutans* GtfC crystal structure (PDB 3AIE)
  with py3Dmol. Showed that conservation increases toward the catalytic site (Spearman
  ρ = −0.30) and produced a Ramachandran analysis.
- Applied rigorous statistics: gene-level replication, Cliff's delta with bootstrap CIs,
  Benjamini-Hochberg FDR, and fixed random seeds. Every reported number is traceable to a
  results file and checked by a unit test.
