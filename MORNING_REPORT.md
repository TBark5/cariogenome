# Morning report

## Top line: did any computed number change?

**No.** The polish pass fixed bugs in code quality, documentation, presentation and
reproducibility, but no scientific result changed. After every change the pipeline was
rerun from clean and all of `results/` was compared with the previous commit. The only
file that differs is `results/runtimes.tsv`, which records wall-clock time and changes on
every run.

One documentation claim was corrected (not a number): earlier documents credited the
three catalytic-residue labels (D477, E515, D588) directly to the GtfC structure paper.
Only that paper's abstract could be checked, and it names the neighbouring subsite
residues, not these three. The documents now say the residues are assigned from the
canonical GH70 motifs (which the pipeline verifies) and that the numbering matches the
paper. See DECISIONS.md.

## Status

All phases of the original brief and every item of the polish brief are done. Nothing
was cut; nothing is left undone except the two items under "Not done" below.

## Data

**Real data.** 22 complete RefSeq genomes from NCBI (10 *S. mutans*, 3 each of four
commensal species) and PDB entry 3AIE, all downloaded on 2026-09-27 and listed with
accessions in `ACCESSIONS.md` (22 genomes, 286 gene records, 284 passing QC; a test
checks these counts and the accession format). A simulated SYNTHETIC mode exists as a
fallback and was verified end to end.

## Confirmed figures

| Item | Value |
|---|---|
| Tests (`pytest`) | 43 collected, 43 pass, also with `-W error` |
| Lint / types | ruff: no findings; mypy: no issues |
| Full run (`run_all.py --clean`) | 164 s wall clock (about 3 minutes); every step under 2 minutes |
| Fresh clone, new venv, `--offline` with sockets blocked | 98 s; all results byte-identical |
| Figures | 28 PNG at 300 dpi + 1 self-contained interactive HTML |

## Bugs found and fixed in the polish pass

1. **Inconsistent documents:** the test count (27 vs 30) and runtime (3 minutes vs
   160 s) disagreed between files. All documents now quote the numbers above, and a test
   checks every number in README.md and MORNING_REPORT.md against `results/` (it verifies
   the numeric values found in each file on every run).
2. **Windows commands inside a `bash` block** in the README. Now separate Windows and
   macOS/Linux blocks; the Windows commands were run as written.
3. **Undocumented `--clean` flag.** Every `run_all.py` option is now in one README table,
   and each was run: default, `--offline`, `--clean`, `--synthetic`, `--out`.
4. **Personal data:** an absolute local path and a git identity in MORNING_REPORT.md
   (removed, and the commit that added them was amended; `git log -p` has no hit), and a
   personal address in `config.yaml` (now a placeholder).
5. **`--offline` was not truly offline:** it skipped the NCBI download but still tried the
   network for NCBI reachability and the PDB file. Now every network call raises in
   offline mode.
6. **The structure HTML needed the internet** (3Dmol.js from a CDN). The library is now
   embedded, with a static PNG of the same view.
7. **Hidden warnings:** two blanket `RuntimeWarning` filters were masking empty-window
   and empty-segment cases. The filters are gone and both cases are handled explicitly.
8. **Lint and type findings:** non-strict zips over paired data, lambda assignments,
   missing annotations and docstrings, possible `None` dereferences, and a
   monkey-patched attribute on Biopython tree nodes (replaced with an explicit map).
9. **Color semantics violated:** *S. mutans* and *S. sanguinis* shared the vermillion
   and blue of the gene classes. Fixed across all figures.
10. **Figure layout:** overlapping titles, legends covering data, an unreadable strain
    axis on the identity heatmaps, a colorbar over axis labels. Fixed after viewing every
    figure.
11. **ACCESSIONS.md rendering:** `nan` in empty QC-flag cells, and an unescaped `|` that
    broke the GH70 table on GitHub. Both were found by the new tests.
12. **README screenshot table** mixed header cells into the body. Rebuilt as two proper
    tables with alt text.

## Privacy scan

Scanned every tracked file (results, figures including PNG metadata, test fixtures,
configuration and the vendored library; the project has no notebooks) plus the full git
history for absolute paths, home directories, usernames and email addresses.

- **Removed:** the local path and git identity in MORNING_REPORT.md (also purged from
  history by amending that commit), and the personal address in `config.yaml`.
- **Still in history:** earlier commits of `config.yaml` contain the git noreply address.
  It is the same address as the commit author on every commit, so rewriting history
  would not remove it from the repository; it was left as is.
- **Remaining email-like strings:** only `you@example.org` in the error message that
  explains how to set `NCBI_EMAIL`.
- **Kept by design:** commit author metadata uses the git identity configured on this
  machine, as the brief requires.
- A test now fails if any tracked file contains an absolute or home-directory path.

## Offline verification

A fresh clone into a new virtual environment ran `run_all.py --clean --offline` with
Python's socket functions patched to raise on any connection. It completed in 98 s with
no connection attempted, used the committed sequence and structure caches, and
reproduced every results file byte for byte (except the timing file). The interactive
HTML was opened in headless Chrome with all network requests refused and rendered the
full 3D model.

## Figure-consistency changes

- One style module (`plotting.py`) for fonts, sizes, line widths, widths per figure
  class, 300 dpi and an opaque white background (a test checks opacity for dark mode).
- Fixed color meaning everywhere: vermillion = virulence-associated, blue =
  housekeeping; species in black, sky blue, green, purple and orange.
- Files renamed `01_…` to `28_…` in pipeline order; README, CAPTIONS.md, the dashboard
  and tests all use the registry in `captions.py`.
- n = 6 virulence-associated and n = 8 housekeeping genes printed on every
  class-comparison figure; axis labels carry units.
- Dashboard charts use the same class colors; screenshots retaken at 1500 × 1250.

## Authorship check

Both required checks were run after the final commit and return **zero lines**:

- `git log` over author, committer and full message of every commit, searched for the
  four required attribution terms (AI tool and vendor names, co-author trailers,
  "generated" footers).
- `git grep` over every tracked file for the same four terms.

The exact commands are in the polish brief; they are not reproduced here because the
command text itself contains the search terms and would match itself.

## Not done

- **GitHub rendering was not viewed on github.com** (the repository has not been
  pushed). Instead, tests check table column consistency, internal links and anchors,
  image references and alt text presence; the figures were checked for opacity.
- **Literature comparison is limited to sources whose abstracts could be fetched and
  read** (PubMed). Full texts were not available, so the README cites only claims stated
  in those abstracts.

## Prediction scorecard (unchanged)

P1 supported (relaxed constraint, not positive selection) · P2 supported · P3 not
supported · P4 partly (distribution) / not supported (discordance) · P5 narrowly
falsified as stated (all three residues invariant, but 20.4% of columns are too, so they
rank in the top 10.2%) · P6 supported.

## Run commands

Windows:

```bat
.venv\Scripts\python run_all.py --offline
.venv\Scripts\python -m pytest
.venv\Scripts\streamlit run app.py
```

macOS/Linux: the same commands with `.venv/bin/` in place of `.venv\Scripts\`.

## Before you push

- Set a reachable NCBI contact email through `NCBI_EMAIL` before re-downloading (see
  README); leave the placeholder in `config.yaml`.

## The 5 things to understand before discussing this project

1. **What dN/dS is and why it is a ratio.** Synonymous changes estimate the neutral rate,
   so dividing by dS cancels mutation rate and time. Every gene here has ω < 1 (purifying
   selection); the virulence-associated genes are *less* constrained, not positively
   selected.
2. **Why the analysis is within *S. mutans*.** Between species, synonymous sites are
   saturated, so dS cannot be measured. Within-species ω reflects polymorphism and is
   inflated by slightly deleterious variants. That is acceptable only because both gene
   classes are measured the same way.
3. **The unit of replication is the gene (n = 6 vs 8).** That is why effect sizes with
   CIs and BH-corrected q-values matter, and why "not significant" means "no evidence",
   not "no effect".
4. **Bootstrap support and discordance.** A bootstrap value is the percentage of
   resampled trees containing a split, not the probability that it is true. Discordance
   can suggest transfer but also arises from recombination, paralogy or noise; no
   well-supported species-mixing discordance was found in coding genes.
5. **What the project does not show.** No causal or clinical claims. The commensal
   homologs show the same elevated ω, so the signal is probably about secreted and
   surface proteins, not cariogenicity. And P5 failed as written; be ready to explain
   why (ties at maximum conservation).
