"""M6: conserved motifs in the glucansucrase (GH70) family and the catalytic residues.

Input: the center-star alignment of every GH70 homolog of GtfB/C/D in the five species
representatives (M1/M3), numbered by S. mutans UA159 GtfC.

Catalytic residues of GtfC (Ito et al. 2011, J. Mol. Biol. 408:177, PDB 3AIC): D477
nucleophile, E515 acid/base, D588 transition-state stabilizer. The run checks that these
positions sit in the canonical GH70 sequence motifs (regions II, III, IV).

Motifs are found by an entropy approach: the non-overlapping windows of MOTIF_LEN columns
with the highest mean conservation. Each motif's position frequency matrix (Henikoff
weighted) gives a position weight matrix and a sequence logo.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .alignment import henikoff_weights
from .config import GENES_DIR, RESULTS
from .seqio import read_fasta

CATALYTIC = {477: ("D", "nucleophile", "region II"), 515: ("E", "acid/base", "region III"),
             588: ("D", "transition-state stabilizer", "region IV")}
MOTIF_PATTERNS = {"region II": r"R.DA.DNV", "region III": r"[ILV]EAW", "region IV": r"F.RAHD"}
GTFC_ID = "Smu_UA159|SMU_RS04625"
AA = "ACDEFGHIKLMNPQRSTVWY"
MOTIF_LEN = 10
N_MOTIFS = 6


def family_alignment() -> dict[str, str]:
    return read_fasta(RESULTS / "alignments" / "gtf_family.aa.fasta")


def family_conservation() -> pd.DataFrame:
    return pd.read_csv(RESULTS / "conservation" / "gtf_family.csv")


def verify_catalytic(gtfc: str) -> pd.DataFrame:
    """Check each catalytic residue's identity and that it lies in its canonical motif."""
    rows = []
    for pos, (aa, role, region) in CATALYTIC.items():
        hits = [m for m in re.finditer(MOTIF_PATTERNS[region], gtfc)
                if m.start() < pos <= m.end()]
        rows.append({"position": pos, "expected": aa, "observed": gtfc[pos - 1], "role": role,
                     "motif": region, "in_motif": bool(hits),
                     "context": gtfc[pos - 6:pos + 5]})
    return pd.DataFrame(rows)


def catalytic_ranks(cons: pd.DataFrame) -> pd.DataFrame:
    """Conservation percentile (mid-rank, ties shared) of each catalytic residue's column."""
    scored = cons.dropna(subset=["conservation"])
    vals = scored["conservation"].values
    n = len(vals)
    frac_invariant = float((vals >= 1 - 1e-9).mean())
    rows = []
    for pos, (aa, role, region) in CATALYTIC.items():
        r = cons[cons["ref_position"] == pos].iloc[0]
        v = r["conservation"]
        greater, equal = (vals > v + 1e-9).sum(), (np.abs(vals - v) <= 1e-9).sum()
        rows.append({"position": pos, "residue": aa, "role": role, "alignment_column": int(r["column"]),
                     "conservation": v, "gap_fraction": r["gap_fraction"],
                     "top_percent_midrank": 100 * (greater + 0.5 * equal) / n,
                     "top_percent_strict": 100 * (greater + equal) / n,
                     "top_percent_best_case": 100 * (greater + 1) / n,
                     "in_top10_midrank": 100 * (greater + 0.5 * equal) / n <= 10})
    out = pd.DataFrame(rows)
    out.attrs["frac_invariant"] = frac_invariant
    out.attrs["n_scored"] = n
    return out


def find_motifs(cons: pd.DataFrame, k: int = MOTIF_LEN, n: int = N_MOTIFS) -> pd.DataFrame:
    """Top-n non-overlapping windows of k reference residues by mean conservation."""
    ref = cons.dropna(subset=["ref_position"]).reset_index(drop=True)
    sm = ref["conservation"].rolling(k).mean()
    order = sm.sort_values(ascending=False).dropna().index
    taken, rows = [], []
    for end in order:
        start = end - k + 1
        if any(abs(start - t) < k for t in taken):
            continue
        taken.append(start)
        block = ref.iloc[start:end + 1]
        cat = [p for p in CATALYTIC if block["ref_position"].min() <= p <= block["ref_position"].max()]
        rows.append({"motif": f"M{len(taken)}", "start": int(block["ref_position"].min()),
                     "end": int(block["ref_position"].max()),
                     "first_column": int(block["column"].min()), "last_column": int(block["column"].max()),
                     "mean_conservation": float(sm[end]), "GtfC_sequence": "".join(block["ref_residue"]),
                     "consensus": "".join(block["consensus"]),
                     "catalytic_residues_inside": ";".join(map(str, cat))})
        if len(taken) == n:
            break
    return pd.DataFrame(rows)


def pfm(msa: dict[str, str], first_col: int, last_col: int) -> pd.DataFrame:
    """Henikoff-weighted position frequency matrix for alignment columns (1-based, inclusive)."""
    w = henikoff_weights(msa)
    rows = []
    for j in range(first_col - 1, last_col):
        col = np.array([s[j] for s in msa.values()])
        f = {a: float(w[col == a].sum()) for a in AA}
        tot = sum(f.values())
        rows.append({a: (v / tot if tot else 0.0) for a, v in f.items()})
    return pd.DataFrame(rows, index=range(first_col, last_col + 1))


def pwm(freq: pd.DataFrame, pseudocount: float = 0.01) -> pd.DataFrame:
    """Log-odds position weight matrix against a uniform background (bits)."""
    f = (freq + pseudocount).div((freq + pseudocount).sum(axis=1), axis=0)
    return np.log2(f * len(AA))


def information_content(freq: pd.DataFrame) -> pd.Series:
    """Bits per position: log2(20) - H."""
    p = freq.values
    with np.errstate(divide="ignore", invalid="ignore"):
        h = -np.nansum(np.where(p > 0, p * np.log2(p), 0), axis=1)
    return pd.Series(np.log2(len(AA)) - h, index=freq.index)


def run() -> dict[str, pd.DataFrame]:
    """Motifs, catalytic-residue ranks and figures for the GH70 family."""
    from . import m6_figures
    fam = read_fasta(GENES_DIR / "gtf_family.faa")
    if GTFC_ID not in fam:  # SYNTHETIC data: no real GtfC, catalytic numbering undefined
        print("  M6 skipped: UA159 GtfC not in the family set (synthetic data)")
        return {}
    msa = family_alignment()
    cons = family_conservation()
    check = verify_catalytic(fam[GTFC_ID])
    check.to_csv(RESULTS / "m6_catalytic_verification.csv", index=False)
    ranks = catalytic_ranks(cons)
    ranks.to_csv(RESULTS / "m6_catalytic_residues.csv", index=False)
    motifs = find_motifs(cons)
    motifs.to_csv(RESULTS / "m6_motifs.csv", index=False)
    windows = [(f"{m.motif} ({m.start}-{m.end})", m.first_column, m.last_column) for m in motifs.itertuples()]
    for pos in CATALYTIC:
        col = int(cons.loc[cons["ref_position"] == pos, "column"].iloc[0])
        windows.append((f"{CATALYTIC[pos][2]} around {CATALYTIC[pos][0]}{pos}", col - 5, col + 5))
    pwms = {}
    for name, a, b in windows:
        f = pfm(msa, a, b)
        pwms[name] = f
        pwm(f).round(3).to_csv(RESULTS / f"m6_pwm_{name.split(' ')[0].replace('region', 'region_')}_{a}.csv")
    summary = {"n_sequences": len(msa), "n_scored_columns": ranks.attrs["n_scored"],
               "fraction_invariant_columns": ranks.attrs["frac_invariant"],
               "p_all_three_invariant_by_chance": ranks.attrs["frac_invariant"] ** 3}
    pd.DataFrame([summary]).to_csv(RESULTS / "m6_summary.csv", index=False)
    n_species = len({k.split("_")[0] for k in msa})
    m6_figures.plot_all(cons, ranks, motifs, pwms, len(msa), n_species)
    print(check.to_string(index=False))
    print(ranks.round(3).to_string(index=False))
    print(motifs[["motif", "start", "end", "mean_conservation", "GtfC_sequence", "catalytic_residues_inside"]].to_string(index=False))
    print(summary)
    return {"ranks": ranks, "motifs": motifs, "check": check}
