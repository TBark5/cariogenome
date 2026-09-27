"""Loading cached, QC-filtered sequences for the analysis modules."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .config import GENES_DIR, RESULTS

CATALOG = RESULTS / "m1_catalog.csv"
MODE_FILE = RESULTS / "data_mode.txt"


def read_fasta(path: Path) -> dict[str, str]:
    """FASTA file to an ordered ``{id: sequence}`` dict (sequence upper-cased)."""
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(str(path), "fasta")}


def catalog() -> pd.DataFrame:
    """The M1 catalog of every extracted record with QC flags."""
    return pd.read_csv(CATALOG)


def load_gene(gene: str, kind: str = "nt", included_only: bool = True) -> dict[str, str]:
    """Sequences of one gene keyed by strain label. ``kind`` is ``"nt"`` or ``"aa"``."""
    seqs = read_fasta(GENES_DIR / f"{gene}.{'fna' if kind == 'nt' else 'faa'}")
    if included_only:
        cat = catalog()
        keep = set(cat[(cat["gene"] == gene) & cat["included"]]["label"])
        seqs = {k: v for k, v in seqs.items() if k in keep}
    return seqs


def species_map() -> dict[str, str]:
    """Strain label to species, from the catalog (works for real and synthetic data)."""
    cat = catalog()
    return dict(zip(cat["label"], cat["species"], strict=True))


def data_mode() -> str:
    """``REAL`` or ``SYNTHETIC``, as recorded by M1."""
    return MODE_FILE.read_text().strip() if MODE_FILE.exists() else "UNKNOWN"


def mode_tag() -> str:
    """Suffix for figure titles: empty for real data, a clear label for synthetic data."""
    return "" if data_mode() == "REAL" else "  [SYNTHETIC DATA]"
