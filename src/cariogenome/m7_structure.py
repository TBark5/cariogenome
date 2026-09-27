"""M7: structural context of GH70 conservation on the S. mutans GtfC crystal structure.

PDB 3AIE (GtfC catalytic core, residues 244-1087, 2.1 A, Ito et al. 2011) chain A is used.
The GH70 family conservation from M3/M6 is mapped onto residues by aligning the chain
sequence to UA159 GtfC, stored in the B-factor column (x100) for coloring, and rendered
with py3Dmol. P6 is tested with a Spearman correlation between residue conservation and
C-alpha distance to the catalytic center (centroid of D477, E515, D588 side chains).
"""
from __future__ import annotations

import io
import math

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser, PPBuilder, Select
from Bio.SeqUtils import seq1
from scipy.stats import spearmanr

from .alignment import pairwise
from .config import FIGURES, GENES_DIR, RESULTS, STRUCT_DIR, load_config, rng
from .entrez_client import fetch_url
from .m6_motifs import CATALYTIC, GTFC_ID
from .seqio import read_fasta

BACKBONE = {"N", "CA", "C", "O"}
UNSCORED = -1.0


class _ChainA(Select):
    def accept_chain(self, chain):
        return chain.id == "A"


def get_structure() -> tuple[str, object]:
    """Download the first available preferred PDB entry; return (id, chain A)."""
    for pid in load_config()["pdb"]["preferred"]:
        try:
            path = fetch_url(f"https://files.rcsb.org/download/{pid}.pdb", f"{pid}.pdb", "RCSB PDB")
        except Exception as exc:
            committed = STRUCT_DIR / f"{pid}_chainA.pdb"
            if committed.exists():  # offline: use the committed chain-A copy
                print(f"  RCSB unreachable ({exc}); using committed {committed.name}")
                return pid, PDBParser(QUIET=True).get_structure(pid, str(committed))[0]["A"]
            print(f"  could not fetch {pid}: {exc}")
            continue
        struct = PDBParser(QUIET=True).get_structure(pid, str(path))
        STRUCT_DIR.mkdir(parents=True, exist_ok=True)
        io_ = PDBIO()
        io_.set_structure(struct)
        io_.save(str(STRUCT_DIR / f"{pid}_chainA.pdb"), _ChainA())
        return pid, struct[0]["A"]
    raise RuntimeError("no GtfC structure could be retrieved")


def residue_table(chain, cons: pd.DataFrame, gtfc: str) -> pd.DataFrame:
    """Per-residue PDB number, GtfC position, conservation and C-alpha coordinates."""
    res = [r for r in chain if r.id[0] == " " and "CA" in r]
    seq = "".join(seq1(r.get_resname()) for r in res)
    g_chain, g_ref, _ = pairwise(seq, gtfc, "aa")
    mapping, i, j = {}, 0, 0
    for a, b in zip(g_chain, g_ref):
        if a != "-" and b != "-":
            mapping[i] = j + 1
        i += a != "-"
        j += b != "-"
    cmap = cons.dropna(subset=["ref_position"]).set_index("ref_position")["conservation"]
    rows = []
    for k, r in enumerate(res):
        pos = mapping.get(k)
        rows.append({"pdb_resnum": r.id[1], "resname": r.get_resname(), "aa": seq[k],
                     "gtfc_position": pos, "gtfc_aa": gtfc[pos - 1] if pos else "",
                     "conservation": float(cmap.get(pos, np.nan)) if pos else np.nan,
                     "x": r["CA"].coord[0], "y": r["CA"].coord[1], "z": r["CA"].coord[2]})
    return pd.DataFrame(rows)


def active_site_center(chain) -> np.ndarray:
    """Centroid of the side-chain atoms of the three catalytic residues."""
    coords = [a.coord for pos in CATALYTIC for a in chain[pos] if a.get_id() not in BACKBONE]
    return np.mean(coords, axis=0)


def ramachandran(chain) -> pd.DataFrame:
    """phi/psi (degrees) for every residue with both angles defined."""
    rows = []
    for pp in PPBuilder().build_peptides(chain):
        for r, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            if phi is None or psi is None:
                continue
            name = r.get_resname()
            rows.append({"pdb_resnum": r.id[1], "resname": name, "phi": math.degrees(phi),
                         "psi": math.degrees(psi),
                         "type": "Gly" if name == "GLY" else ("Pro" if name == "PRO" else "general")})
    df = pd.DataFrame(rows)
    phi, psi = df["phi"], df["psi"]
    df["region"] = np.select(
        [(phi < 0) & (psi > -120) & (psi < 50), (phi < 0) & ((psi >= 50) | (psi <= -150)), phi >= 0],
        ["alpha (right-handed)", "beta / polyproline", "left-handed (phi > 0)"], "other")
    return df


def write_conservation_pdb(pid: str, table: pd.DataFrame) -> str:
    """Chain A with B-factor = conservation x 100 (-1 where unscored); returns PDB text."""
    struct = PDBParser(QUIET=True).get_structure(pid, str(STRUCT_DIR / f"{pid}_chainA.pdb"))
    score = dict(zip(table["pdb_resnum"], table["conservation"]))
    for r in struct[0]["A"]:
        v = score.get(r.id[1], np.nan)
        for a in r:
            a.set_bfactor(UNSCORED if not np.isfinite(v) else round(100 * v, 2))
    buf = io.StringIO()
    io_ = PDBIO()
    io_.set_structure(struct)
    io_.save(buf)
    text = buf.getvalue()
    (STRUCT_DIR / f"{pid}_chainA_conservation.pdb").write_text(text)
    return text


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 2000) -> tuple[float, float, float, float]:
    """Spearman rho, p-value, and a residue-bootstrap 95% CI."""
    res = spearmanr(x, y)
    gen = rng(8)
    boots = []
    for _ in range(n_boot):
        idx = gen.integers(0, len(x), len(x))
        boots.append(spearmanr(x[idx], y[idx]).statistic)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(res.statistic), float(res.pvalue), float(lo), float(hi)


def run() -> dict:
    """Retrieve the structure, map conservation, test P6, Ramachandran, figures, 3D view."""
    from . import m7_figures
    fam = read_fasta(GENES_DIR / "gtf_family.faa")
    if GTFC_ID not in fam:
        print("  M7 skipped: UA159 GtfC not in the family set (synthetic data)")
        return {}
    cons = pd.read_csv(RESULTS / "conservation" / "gtf_family.csv")
    pid, chain = get_structure()
    table = residue_table(chain, cons, fam[GTFC_ID])
    center = active_site_center(chain)
    table["distance_to_active_site"] = np.linalg.norm(table[["x", "y", "z"]].values - center, axis=1)
    table.to_csv(RESULTS / "m7_residue_conservation.csv", index=False)
    scored = table.dropna(subset=["conservation"])
    rho, p, lo, hi = spearman_ci(scored["conservation"].values, scored["distance_to_active_site"].values)
    near = scored[scored["distance_to_active_site"] <= 12]
    summary = {"pdb_id": pid, "chain": "A", "n_residues_modelled": len(table),
               "n_residues_scored": len(scored),
               "numbering_matches_gtfc": bool((table["pdb_resnum"] == table["gtfc_position"]).mean() > 0.99),
               "spearman_rho_conservation_vs_distance": rho, "rho_ci_low": lo, "rho_ci_high": hi,
               "spearman_p": p, "mean_conservation_within_12A": float(near["conservation"].mean()),
               "n_within_12A": len(near),
               "mean_conservation_beyond_12A": float(scored.loc[scored["distance_to_active_site"] > 12, "conservation"].mean())}
    rama = ramachandran(chain)
    rama.to_csv(RESULTS / "m7_ramachandran.csv", index=False)
    for reg, frac in rama["region"].value_counts(normalize=True).items():
        summary[f"rama_fraction_{reg.split(' ')[0]}"] = float(frac)
    summary["rama_n_residues"] = len(rama)
    pd.DataFrame([summary]).to_csv(RESULTS / "m7_summary.csv", index=False)
    pdb_text = write_conservation_pdb(pid, table)
    m7_figures.plot_all(table, rama, summary, pdb_text, pid)
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in summary.items()})
    return summary
