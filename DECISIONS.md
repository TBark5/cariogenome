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
