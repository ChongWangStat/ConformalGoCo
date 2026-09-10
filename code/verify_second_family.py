"""Reproducibility gate for the second predictor family (GoCo-M and GoCo-N on FunMap and STRING).

Re-runs each arm from clean and requires the output to match the shipped split-level results
BIT-FOR-BIT. This is not a comparison of means: every metric of every arm on every re-run split
must be identical at atol = 0, or the gate fails and the script exits non-zero.

    python code/verify_second_family.py            # 3 splits per arm, about 6 minutes
    python code/verify_second_family.py --repeats 10

The STRING arm needs the validation set built first:

    python code/build_string_validation.py
"""
from __future__ import annotations
import argparse, os, subprocess, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results", "second_family")
TMP = os.path.join(ROOT, "results", "second_family", "_verify")
os.makedirs(TMP, exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=3)
args = ap.parse_args()

COLS = ["label", "fdp", "correct", "total", "units"]
ok_all, skipped = True, []


def gate(name, cmd, env_extra, fresh, stored, key, needs=(), hint=""):
    global ok_all
    print("\n=== %s ===" % name, flush=True)
    missing = [q for q in needs if not os.path.exists(q)]
    if missing or not os.path.exists(stored):
        what = os.path.relpath(missing[0], ROOT) if missing else os.path.relpath(stored, ROOT)
        print("  SKIP: %s not present.%s" % (what, (" " + hint) if hint and missing else ""))
        skipped.append(name); return
    env = dict(os.environ); env.update(env_extra)
    r = subprocess.run(cmd, shell=True, env=env, capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        print("  RUN FAILED\n", r.stdout[-1200:], r.stderr[-1200:]); ok_all = False; return
    a, b = pd.read_csv(fresh), pd.read_csv(stored)
    b = b[b.repeat.isin(a.repeat.unique())]
    m = a.merge(b, on=key, suffixes=("_new", "_old"))
    if len(m) != len(a):
        print("  ROW COUNT MISMATCH: %d re-run vs %d matched" % (len(a), len(m))); ok_all = False; return
    bad = []
    for c in COLS:
        x, y = m[c + "_new"], m[c + "_old"]
        same = (x == y).all() if x.dtype == object else np.array_equal(x.to_numpy(), y.to_numpy())
        if not same:
            bad.append(c if x.dtype == object else "%s (max|diff| %.3g)" % (c, np.abs(x - y).max()))
    if bad:
        print("  MISMATCH:", "; ".join(bad)); ok_all = False
    else:
        print("  %d rows x %d metrics identical (atol = 0)" % (len(m), len(COLS)))


N = args.repeats
gate("FunMap: global grid, GoCo-M and GoCo-N",
     'python -u code/run_goco_second_family.py --repeats %d --out "%s/funmap.csv"' % (N, TMP),
     {"GOCO_DATA": os.path.join(ROOT, "funmap")},
     os.path.join(TMP, "funmap.csv"), os.path.join(RES, "goco_funmap_50splits.csv"),
     ["repeat", "alpha", "arm"],
     needs=[os.path.join(ROOT, "funmap", "funmap_gene_go_scores.csv.gz")],
     hint="It is distributed with the Zenodo archive; see README.")

gate("STRING (independent validation)",
     'python -u code/run_goco_second_family.py --repeats %d --out "%s/string.csv"' % (N, TMP),
     {"GOCO_DATA": os.path.join(ROOT, "string", "funmap_view")},
     os.path.join(TMP, "string.csv"), os.path.join(RES, "goco_string_50splits.csv"),
     ["repeat", "alpha", "arm"],
     needs=[os.path.join(ROOT, "string", "funmap_view", "funmap_gene_go_scores.csv.gz")],
     hint="Run: python code/build_string_validation.py")

gate("Wainberg negative control (GoCo-N on module-shared evidence)",
     'python -u code/wainberg_goco_n_control.py --repeats %d --out "%s/wainberg.csv"' % (N, TMP),
     {"GOCO_CODE": HERE},
     os.path.join(TMP, "wainberg.csv"), os.path.join(RES, "wainberg_goco_n_control_50splits.csv"),
     ["repeat", "alpha", "delta", "arm"],
     needs=[os.path.join(ROOT, "w")])

if skipped:
    print("\nSkipped because their inputs are not present:")
    for k in skipped:
        print("  - " + k)
print("\n" + ("REPRODUCIBILITY GATE: PASS" if ok_all else "REPRODUCIBILITY GATE: FAIL"))
sys.exit(0 if ok_all else 1)
