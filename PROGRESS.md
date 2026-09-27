# PROGRESS

## Done
- Phase 6 (M6 + M7): GH70 family motifs, PWMs and logos, catalytic-residue ranks
  (`results/m6_*.csv`); PDB 3AIE conservation mapping, P6 correlation, Ramachandran,
  interactive `figures/m7_structure_conservation.html` (`results/m7_*.csv`). Runtime ~8 s.
- Phase 5 (M5): pooled NG86 dN/dS within each species with codon-bootstrap CIs,
  between-species saturation check, class test (`results/m5_tests.csv`), per-gene vs
  controls, sliding windows, synthetic recovery validation (`validation.py`), 5 figures.
  Runtime ~8 s + 5 s validation.
- Phase 4 (M4): NJ + UPGMA trees with 100 bootstrap replicates for 15 genes and the
  concatenated housekeeping reference (`results/trees/*.nwk`), discordance table
  (`results/m4_discordance.csv`, `m4_conflicts.csv`), tests, 5 figures. Runtime ~50 s.
- Phase 3 (M3): pairwise identity matrices, center-star MSAs + codon alignments
  (`results/alignments/`), Henikoff-weighted entropy per column (`results/conservation/`),
  extreme regions, GH70 family alignment, P2 test (`results/m3_tests.csv`), 5 figures.
  Runtime ~60 s.
- Phase 2 (M2): GC/GC3/GC skew/CAI per record, gene-level virulence vs control tests
  (`results/m2_tests.csv`), amino-acid tests (BH over 20), RSCU table, genome
  background, 5 figures `figures/m2_*.png`. Runtime ~17 s.
- Phase 1 (M1): 22 complete genomes downloaded (cached in `data/raw`, git-ignored),
  RBH orthologs for the 14 coding genes + first 16S copy, QC (284/286 records pass),
  catalog `results/m1_catalog.csv`, presence matrix, GH70 family set, `ACCESSIONS.md`,
  figures `m1_presence_absence.png`, `m1_length_qc.png`. SYNTHETIC fallback in
  `synthetic.py` (used automatically if NCBI is unreachable and no cache exists).
- Phase 0: folder layout, `.venv` (Python 3.14), pinned `requirements.txt`, git repo,
  `config.yaml`, `HYPOTHESIS.md`, `DECISIONS.md`, shared plotting style
  (`src/cariogenome/plotting.py`).

## In progress
- Phase 7 (visual pass + captions), then Phase 8 (Streamlit app).

## Next
- Phases 2-10 as listed in the project brief.

## Known bugs
- None.

## How to run
```
.venv\Scripts\python -m pip install -r requirements.txt
```

## Key files
- `config.yaml`: genomes, gene panel, all parameters.
- `src/cariogenome/config.py`: paths and seeded RNGs.
- `src/cariogenome/plotting.py`: palette and `save()` at 300 dpi.
