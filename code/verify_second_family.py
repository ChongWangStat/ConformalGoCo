"""Reproducibility gate for the neighbourhood-structure family (FunMap development, STRING validation).

Re-runs the first few splits of every shipped run from clean and requires bit-for-bit agreement with the shipped
split-level files (atol = 0 on every metric of every split), or exits non-zero.

    python verify_second_family.py                 # 3 splits per run, a few minutes
    python verify_second_family.py --repeats 10

Runs checked: run_goco_second_family.py (certified arms), rehearsal_kfold_sf.py (ten-fold selection rule),
generic_sf.py (generic Direct-loss calibrators, coarse and refined grids) and frontier_exhaustive_sf.py.
The STRING view (string/funmap_view) must exist -- python code/build_string_from_shipped_nb.py -- or its gates are skipped.
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ap = argparse.ArgumentParser(); ap.add_argument("--repeats", type=int, default=3); args = ap.parse_args()
DATA = {"funmap": os.path.join(ROOT, "funmap"), "string": os.path.join(ROOT, "string", "funmap_view")}
RES = os.path.join(ROOT, "results", "sf")
TMP = tempfile.mkdtemp(prefix="goco_verify_")
ok_all = True; skipped = []


def same(fresh, stored, keys, cols):
    a = pd.read_csv(fresh); b = pd.read_csv(stored)
    b = b[b.repeat.isin(a.repeat.unique())]
    m = a.merge(b, on=keys, suffixes=("", "_stored"))
    if len(m) != len(a) or len(m) != len(b): return False, f"row mismatch fresh {len(a)} stored {len(b)} merged {len(m)}"
    bad = []
    for c in cols:
        x, y = m[c], m[c + "_stored"]
        if x.dtype.kind in "fc" or y.dtype.kind in "fc":
            if not np.array_equal(x.to_numpy(float), y.to_numpy(float), equal_nan=True): bad.append(c)
        elif not (x.astype(str) == y.astype(str)).all(): bad.append(c)
    return (not bad), ("differences in " + ", ".join(bad) if bad else "identical")


def gate(name, data, script, stored, keys, cols, extra=()):
    global ok_all
    print(f"\n=== {name} ===", flush=True)
    if not os.path.isdir(DATA[data]) or not os.path.exists(stored):
        print("  SKIP: input or stored file missing"); skipped.append(name); return
    out = os.path.join(TMP, os.path.basename(stored))
    env = dict(os.environ); env["GOCO_DATA"] = DATA[data]
    cmd = [sys.executable, os.path.join(HERE, script), "--repeats", str(args.repeats), "--out", out, *extra]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0: print("  FAIL: run error\n" + r.stderr[-2000:]); ok_all = False; return
    good, msg = same(out, stored, keys, cols)
    print(f"  {'PASS' if good else 'FAIL'}: {msg}"); ok_all &= good
    return out


ARM = ["label", "fdp", "pool_risk", "correct", "total", "units", "supp_units", "frontier"]
for data in ["funmap", "string"]:
    fresh = gate(f"{data}: certified arms", data, "run_goco_second_family.py", os.path.join(RES, f"{data}_v4_100.csv"),
                 ["repeat", "alpha", "delta", "arm"], ARM)
    if fresh is not None:
        sel_fresh = fresh.replace(".csv", "_selection.csv")
        good, msg = same(sel_fresh, os.path.join(RES, f"{data}_v4_100_selection.csv"), ["repeat", "alpha"], ["yield_M_T", "yield_N_T", "choose_N", "label_M", "label_N"])
        print(f"  two-half selection rows: {'PASS' if good else 'FAIL'}: {msg}"); ok_all &= good
        # exhaustive frontier on the fresh rows
        fr_out = os.path.join(TMP, f"{data}_frontier.csv")
        env = dict(os.environ); env["GOCO_DATA"] = DATA[data]
        r = subprocess.run([sys.executable, os.path.join(HERE, "frontier_exhaustive_sf.py"), "--run", fresh, "--out", fr_out], env=env, capture_output=True, text=True)
        if r.returncode != 0: print("  frontier: FAIL run error\n" + r.stderr[-1500:]); ok_all = False
        else:
            good, msg = same(fr_out, os.path.join(RES, f"{data}_v4_100.csv"), ["repeat", "alpha", "delta", "arm"], ["frontier_all"])
            print(f"  exhaustive frontier: {'PASS' if good else 'FAIL'}: {msg}"); ok_all &= good
    gate(f"{data}: ten-fold selection rehearsal", data, "rehearsal_kfold_sf.py", os.path.join(RES, f"{data}_v4_100_selection_k10.csv"),
         ["repeat", "alpha"], ["yield_M_T", "yield_N_T", "choose_N", "label_M", "label_N"])
    for grid in ["coarse", "dense"]:
        stored = os.path.join(RES, f"{data}_generic_{grid}.csv")
        if os.path.exists(stored):
            cols = [c for c in pd.read_csv(stored, nrows=1).columns if c not in ("repeat", "alpha", "delta", "rule", "method", "grid")]
            keys = [c for c in ("repeat", "alpha", "delta", "rule", "method") if c in pd.read_csv(stored, nrows=1).columns]
            gate(f"{data}: generic calibrators ({grid} grid)", data, "generic_sf.py", stored, keys, cols, extra=("--grid", grid))

print("\nSKIPPED:", skipped if skipped else "none")
print("ALL GATES PASSED" if ok_all else "GATE FAILURE")
sys.exit(0 if ok_all else 1)
