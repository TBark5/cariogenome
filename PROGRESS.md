# PROGRESS

## Done
- Phase 0: folder layout, `.venv` (Python 3.14), pinned `requirements.txt`, git repo,
  `config.yaml`, `HYPOTHESIS.md`, `DECISIONS.md`, shared plotting style
  (`src/cariogenome/plotting.py`).

## In progress
- Phase 1 (M1 sequence retrieval and QC).

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
