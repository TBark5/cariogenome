# LinkedIn post (draft)

I just finished a comparative genomics project on the bacteria behind tooth decay.

*Streptococcus mutans* is strongly associated with dental caries, but its close relatives
in the mouth are mostly harmless. I wanted to know whether the genes linked to caries
(the enzymes that build sticky dental plaque, a surface adhesin, and a quorum-sensing
gene) evolve differently from ordinary "housekeeping" genes.

I built a Python pipeline that downloads 22 complete genomes from 5 species, finds the
matching genes, and compares them in four ways: DNA composition, conservation,
evolutionary trees, and selection pressure (dN/dS).

What I found:
- The caries-associated genes are under weaker purifying selection than housekeeping
  genes (median dN/dS 0.156 vs 0.019). They are still constrained, and there is no sign
  of positive selection.
- Their harmless relatives' versions of these genes show the same pattern, so this looks
  like a property of these protein families, not something unique to the caries species.
- I found no evidence that these genes moved between species, and no difference in
  codon usage.
- In the crystal structure of one of the plaque-building enzymes, the most conserved
  amino acids cluster around the active site.

One of my predictions failed, and I kept it in the write-up. That was one of the most
useful parts of the project.

Things I learned: why dN/dS is a ratio, what a bootstrap value does and does not mean,
how easily between-species comparisons saturate, and how much a good control group
matters.

These are evolutionary patterns, not evidence about disease mechanisms, and the gene set
and strain sampling are small. Code, methods, and every number are in the repository.

#bioinformatics #genomics #microbiology #python
