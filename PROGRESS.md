# PROGRESS

## Done
- **Phases 0–10 of the original brief** (M1–M7, figures, dashboard, tests, docs).
  Results are unchanged since then.
- **Polish pass** (bug fixes, consistency, presentation), all committed:
  - Lint and types: ruff and mypy clean; configuration in `pyproject.toml`, tools pinned
    in `requirements-dev.txt`. All tests pass with `-W error`.
  - Privacy: NCBI email is a placeholder (`REPLACE_WITH_YOUR_EMAIL`), overridden by the
    `NCBI_EMAIL` environment variable; no local paths in tracked files (tested).
  - Offline: `--offline` sets `CARIOGENOME_OFFLINE=1` and every network call raises;
    verified in a fresh clone with sockets blocked. The interactive structure page embeds
    3Dmol.js (`data/vendor/`), with a static PNG alongside.
  - Figures: one style module (`plotting.py`), fixed color meanings, numbered file names
    from the registry in `captions.py`, class sample sizes on class-comparison figures,
    opaque backgrounds; every figure was visually reviewed.
  - Docs: README rebuilt (headline, contents, scorecard, "does not show", verified
    literature, per-OS commands, all options, reproducibility, licensing). LICENSE (MIT)
    and LICENSE-DATA added. MORNING_REPORT.md lists every fix.
  - Tests (43): numbers in README.md and MORNING_REPORT.md vs `results/`, figure
    references and captions, links and anchors, table columns, absolute paths, email
    placeholder, offline mode, accessions, PNG opacity, plus the earlier method,
    synthetic-recovery and dashboard tests.

## In progress
- Nothing.

## Next (optional)
- Push to GitHub and view the README there (only checked by tests so far).
- ML trees (IQ-TREE) and codon-model selection tests (PAML/HyPhy); more *S. mutans*
  genomes; recombination detection.

## Known bugs
- None known. `results/runtimes.tsv` changes on every run by design (wall-clock time).

## How to run
```
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python run_all.py --offline      # about 3 minutes, no network
.venv\Scripts\python -m pytest                 # 43 tests
.venv\Scripts\streamlit run app.py
```
Development tools: `.venv\Scripts\python -m pip install -r requirements-dev.txt`, then
`ruff check src tests run_all.py app.py` and `mypy src/cariogenome run_all.py app.py tests`.
To re-download from NCBI, set `NCBI_EMAIL` first.

## Key files
- `config.yaml`: genomes, gene panel, thresholds, seeds (email placeholder).
- `run_all.py`: pipeline order and options (`--offline`, `--clean`, `--synthetic`, `--out`).
- `src/cariogenome/plotting.py`: the only place figure style is defined.
- `src/cariogenome/captions.py`: figure registry (order, file names, captions).
- `src/cariogenome/entrez_client.py`: all network access, offline switch, email handling.
- `tests/test_docs.py`: documentation, figure, privacy and provenance checks.
- `DECISIONS.md` (see "Polish pass"), `MORNING_REPORT.md`, `LICENSE`, `LICENSE-DATA`.

## BLOCKED
- Nothing blocked.
