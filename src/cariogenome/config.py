"""Project paths, configuration loading and seeded random number generators."""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
GENES_DIR = DATA / "genes"
GENBANK_DIR = DATA / "genbank"
GENOMES_DIR = DATA / "genomes"
STRUCT_DIR = DATA / "structure"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


@dataclass(frozen=True)
class Genome:
    """A complete chromosome used in the study."""

    label: str
    species: str
    accession: str


@lru_cache(maxsize=1)
def load_config(path: str | None = None) -> dict[str, Any]:
    """Read ``config.yaml`` from the project root (cached)."""
    cfg_path = Path(path) if path else ROOT / "config.yaml"
    with open(cfg_path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def genomes() -> list[Genome]:
    """All genomes listed in the configuration, in order."""
    return [Genome(**g) for g in load_config()["genomes"]]


def species_of(label: str) -> str:
    """Species name for a genome/strain label such as ``Smu_UA159``."""
    for g in genomes():
        if g.label == label:
            return g.species
    raise KeyError(label)


def virulence_genes() -> list[str]:
    return list(load_config()["genes"]["virulence"])


def housekeeping_genes() -> list[str]:
    return list(load_config()["genes"]["housekeeping"])


def coding_genes() -> list[str]:
    """Virulence genes followed by housekeeping protein-coding controls."""
    return virulence_genes() + housekeeping_genes()


def all_genes() -> list[str]:
    """All loci, including the 16S rRNA control."""
    return coding_genes() + list(load_config()["genes"]["rrna"])


def gene_class(gene: str) -> str:
    """Return ``"virulence"`` or ``"control"`` for a gene name."""
    return "virulence" if gene in virulence_genes() else "control"


def rng(offset: int = 0) -> np.random.Generator:
    """Seeded generator; ``offset`` gives independent but reproducible streams."""
    return np.random.default_rng(load_config()["seed"] + offset)


def ensure_dirs() -> None:
    """Create every output directory used by the pipeline."""
    for d in (DATA, RAW, GENES_DIR, GENBANK_DIR, GENOMES_DIR, STRUCT_DIR, RESULTS, FIGURES):
        d.mkdir(parents=True, exist_ok=True)
