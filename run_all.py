"""Regenerate every result, figure and table of CARIOGENOME.

Usage:
    python run_all.py              # real NCBI data (downloads once, then uses the cache)
    python run_all.py --offline    # never touch NCBI; use the committed data/genes cache
    python run_all.py --synthetic  # simulated data only (SYNTHETIC mode)
    python run_all.py --out DIR    # write data/, results/, figures/ under DIR instead
    python run_all.py --clean      # delete results/ and generated figures first
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--synthetic", action="store_true", help="force SYNTHETIC simulated data")
    p.add_argument("--offline", action="store_true", help="use only cached sequences")
    p.add_argument("--out", type=Path, default=None, help="output root (default: project)")
    p.add_argument("--clean", action="store_true", help="remove previous results and figures")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.out is not None:
        os.environ["CARIOGENOME_OUT"] = str(args.out.resolve())
    sys.path.insert(0, str(ROOT / "src"))
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    from cariogenome import (captions, m1_qc, m2_composition, m3_conservation, m4_phylogeny,
                             m5_selection, m6_motifs, m7_structure, summary, validation)
    from cariogenome.config import FIGURES, RESULTS, ensure_dirs

    if args.clean:
        shutil.rmtree(RESULTS, ignore_errors=True)
        for f in FIGURES.glob("*"):
            if f.suffix in (".png", ".html"):
                f.unlink()
    ensure_dirs()
    steps = [
        ("M1 retrieval and QC", lambda: m1_qc.run(force_synthetic=args.synthetic, offline=args.offline)),
        ("M2 composition", m2_composition.run),
        ("M3 alignment and conservation", m3_conservation.run),
        ("M4 phylogenetics", m4_phylogeny.run),
        ("M5 selection (dN/dS)", m5_selection.run),
        ("Synthetic recovery validation", validation.run),
        ("M6 motifs", m6_motifs.run),
        ("M7 structure", m7_structure.run),
        ("Summary", summary.run),
        ("Captions", captions.write_captions),
    ]
    timings = []
    t_all = time.time()
    for name, fn in steps:
        print(f"\n=== {name} ===", flush=True)
        t0 = time.time()
        fn()
        timings.append((name, time.time() - t0))
        print(f"--- {name}: {timings[-1][1]:.1f} s", flush=True)
    lines = [f"{n}\t{t:.1f}" for n, t in timings] + [f"TOTAL\t{time.time() - t_all:.1f}"]
    (RESULTS / "runtimes.tsv").write_text("step\tseconds\n" + "\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nData mode: {(RESULTS / 'data_mode.txt').read_text().strip()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
