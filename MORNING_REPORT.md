# Morning report

## Status
All 10 phases are complete and committed. Nothing was cut: all seven modules (M1–M7) are
implemented, tested and documented.

## Real or synthetic data?
**Real data.** 22 complete RefSeq genomes were downloaded from NCBI on 2026-09-27
(10 *S. mutans*, 3 each of *S. sanguinis*, *S. gordonii*, *S. mitis*, *S. salivarius*),
plus PDB 3AIE from RCSB. Every accession and access date is in `ACCESSIONS.md`. The
synthetic mode exists as the mandatory fallback. It was run end to end as a check
(`python run_all.py --synthetic --out <folder>`), and its method-validation part runs on
every normal run.

## What works
- `run_all.py` regenerates every result and figure in about 160 s. Every step is under
  2 minutes (the longest is M1 at 74 s).
- 30 tests pass: methods, NG86 vs Biopython, synthetic recovery, the dashboard, and
  README-vs-results number checks.
- A clean clone with a brand-new venv, running `run_all.py --clean --offline`, reproduced
  every results file byte for byte (only the timing file differs).
- The Streamlit dashboard launches (health check OK) and all 8 tabs render for every gene
  (AppTest).
- There is no AI attribution in any file or commit (checked by grep before the final
  commit).

## What does not work / caveats
- One pre-registered prediction (P5) narrowly failed as operationalized. The catalytic
  residues are invariant, but so is 20.4% of the alignment.
- The within-*S. mutans* tree-discordance metric saturates (RF ≈ 1 for 12 of 14 coding
  genes). It cannot separate the gene classes.
- Between-species dN/dS is not estimable: synonymous sites are saturated in 37 of 42
  comparisons.
- The interactive structure (HTML) loads 3Dmol.js from a CDN, so it needs internet access
  when opened.

## Final numbers (from `results/`)
| Result | Value |
|---|---|
| Records passing QC | 284 of 286 |
| dN/dS in *S. mutans*, virulence vs control (median) | 0.156 vs 0.019; Cliff's δ = 0.88 [0.50, 1.00], q = 0.012 |
| Driven by | dN (δ = 0.92, q = 0.012), not dS (δ = 0.42, q = 0.228) |
| Mean aa identity within *S. mutans* | 99.06% vs 99.90%; δ = −0.96 [−1.00, −0.75], q = 0.0055 |
| Composition (GC, GC3, CAI, GC skew) | no difference, all q ≥ 0.33 |
| Supported species-mixing gene-tree conflicts | 0 in all 14 coding genes (1 in 16S) |
| GtfC catalytic residues D477/E515/D588 | invariant across 11 GH70 enzymes; mid-rank top 10.2% |
| Conservation vs distance to active site | Spearman ρ = −0.30 [−0.36, −0.23] |
| Synthetic validation | topology nRF = 0.0; ω Spearman 0.991, median error 8.8% |

## Prediction scorecard
P1 supported (relaxed constraint, not positive selection) · P2 supported · P3 not
supported · P4 partly (distribution) / not supported (discordance) · P5 narrowly
falsified as stated · P6 supported.

## Exact run commands
```
.venv\Scripts\python run_all.py
.venv\Scripts\python -m pytest
.venv\Scripts\streamlit run app.py
```
Optional: `--offline` (no network), `--clean` (delete old results first),
`--synthetic --out synthetic_run` (simulated data in a separate folder).

## Modules cut
None.

## Before you push
- Set a reachable NCBI contact email before re-downloading (see README).

## The 5 things you must understand before discussing this with anyone
1. **What dN/dS is and why it is a ratio.** Synonymous changes estimate the neutral rate,
   so dividing by dS cancels mutation rate and time. Every gene here has ω < 1
   (purifying selection). Virulence genes are *less* constrained, not positively
   selected.
2. **Why the analysis is within *S. mutans*.** Between species, synonymous sites are
   saturated (pS near the 0.75 JC limit), so dS cannot be measured. Within-species ω
   reflects polymorphism and is inflated by slightly deleterious variants. That is
   acceptable only because both gene classes are measured the same way.
3. **The unit of replication is the gene (n = 6 vs 8).** That is why effect sizes with
   CIs (Cliff's δ) and BH-corrected q-values matter, and why "not significant" means "no
   evidence", not "no effect".
4. **Bootstrap support and discordance.** A bootstrap value is the % of resampled trees
   containing a split, not the probability that it is true. Gene-tree discordance can
   suggest HGT but also arises from recombination, paralogy or noise. Here no
   well-supported discordance was found in coding genes.
5. **What the project does not show.** No causal or clinical claims: sequence
   signatures are associations with gene categories. The commensal homologs show the same
   elevated ω, so the pattern is probably about secreted/surface proteins, not
   cariogenicity. And P5 failed as stated; be ready to explain why (ties at maximum
   conservation).
