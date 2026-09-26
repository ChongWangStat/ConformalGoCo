"""Write every split used in the paper as an explicit table (gene, split, fold), one file per data set.

Module family (Wainberg, Sanger, DRIVE, HAP1): perm = default_rng(20260905 + s).permutation(n), T = first round(0.1 n),
C = next round(0.7 n) - round(0.1 n), E = the rest; global methods calibrate on T ∪ C.
Neighbourhood family (FunMap, STRING view): the same cut with default_rng(20260910 + s).
Development splits: the module-family instantiation was fixed on Wainberg splits 0-9; the neighbourhood family on
FunMap splits 0-9.  Output: splits/<dataset>_splits.csv.gz with columns gene, split, fold in {T, C, E}.
"""
import os, sys, numpy as np, pandas as pd
D = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(D); sys.path.insert(0, D)
out = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/splits"; os.makedirs(out, exist_ok=True)
NREP = 100


def write(name, units, seed0):
    n = len(units); nt, nc = round(0.10 * n), round(0.70 * n)
    fold = np.empty((NREP, n), dtype="U1")
    for s in range(NREP):
        perm = np.random.default_rng(seed0 + s).permutation(n)
        f = np.empty(n, dtype="U1"); f[perm[:nt]] = "T"; f[perm[nt:nc]] = "C"; f[perm[nc:]] = "E"; fold[s] = f
    df = pd.DataFrame({"gene": np.tile(units, NREP), "split": np.repeat(np.arange(NREP), n), "fold": fold.ravel()})
    df.to_csv(f"{out}/{name}_splits.csv.gz", index=False)
    print(name, "n =", n, "(T, C, E) =", (nt, nc - nt, n - nc), flush=True)


import goco_rerun as gr
dset, tset = gr.load_truth()
for name in ["Wainberg", "Sanger", "DRIVE", "HAP1"]:
    ds = gr.Dataset(name, dset, tset); write(name, ds.units, gr.SEED0)

import importlib, goco_second_family as sf
for name, path in [("FunMap", f"{ROOT}/funmap"), ("STRING", f"{ROOT}/string/funmap_view")]:
    os.environ["GOCO_DATA"] = path; importlib.reload(sf)
    ds = sf.FunMap(); write(name, ds.units, sf.SEED0)
