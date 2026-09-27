# PROGRESS

## Done
- Phase 1 (M1): 22 complete genomes downloaded (cached in `data/raw`, git-ignored),
  RBH orthologs for the 14 coding genes + first 16S copy, QC (284/286 records pass),
  catalog `results/m1_catalog.csv`, presence matrix, GH70 family set, `ACCESSIONS.md`,
  figures `m1_presence_absence.png`, `m1_length_qc.png`. SYNTHETIC fallback in
  `synthetic.py` (used automatically if NCBI is unreachable and no cache exists).
- Phase 0: folder layout, `.venv` (Python 3.14), pinned `requirements.txt`, git repo,
  `config.yaml`, `HYPOTHESIS.md`, `DECISIONS.md`, shared plotting style
  (`src/cariogenome/plotting.py`).

## In progress
- Phase 2 (M2 composition).

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
