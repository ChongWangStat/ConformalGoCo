#!/usr/bin/env python3
"""Regenerate BMC results, data tables and figures from code and frozen analysis inputs.

All outputs are created locally and are not tracked. The module panels have frozen
per-call prediction scores; FunMap/STRING panel reconstructions use the archived
100-split algorithm outputs because their complete frozen score files are absent.
"""
from pathlib import Path
import argparse, importlib.util, os, subprocess, sys, tarfile, tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
TARGET = (0.10, 0.10)
PANELS = ("Wainberg", "Sanger", "DRIVE", "HAP1")
EXT = ROOT / "bmc_analysis/GoCo_BMC_Reproducibility_Extension_2026-11-15.tar.gz"

def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required frozen input: {path}")
    return pd.read_csv(path)

def subset(df):
    return df[np.isclose(df.alpha, TARGET[0]) & np.isclose(df.delta, TARGET[1])]

def check100(df, split_col, label):
    if len(df) != 100 or df[split_col].nunique() != 100 or df[split_col].duplicated().any():
        raise ValueError(f"{label}: expected exactly 100 unique original split runs; got {len(df)}")

def make_tables(out, staging):
    module = subset(read_csv(ROOT / "results/all_methods_harmonised_100splits.csv"))
    rows = []
    for panel in PANELS:
        for method, display in (("Boger", "Multilabel"), ("GoDag", "GoDag"), ("GoCo", "GoCo-M")):
            x = module[(module.dataset == panel) & (module.method == method)]
            check100(x, "seed", panel + "/" + method)
            rows.append(dict(panel=panel, policy=display, fdp=x.unit_fdp.mean(),
                             calls=x.total_calls.mean(), supported=x.go_yield.mean(),
                             fraction=x.go_yield.sum()/x.total_calls.sum()))
    pd.DataFrame(rows).to_csv(out / "Table1_regenerated.csv", index=False)
    rows = []
    for panel in ("funmap", "string"):
        d = subset(read_csv(ROOT / "results/sf" / (panel + "_v4_100.csv")))
        for method, display in (("global-dense", "refined GoDag"), ("goco-N", "GoCo-N")):
            x = d[d.arm == method]
            check100(x, "repeat", panel + "/" + method)
            rows.append(dict(panel=panel, policy=display, fdp=x.fdp.mean(),
                             calls=x.total.mean(), supported=x.correct.mean(),
                             fraction=x.correct.sum()/x.total.sum()))
    pd.DataFrame(rows).to_csv(out / "Table2_regenerated.csv", index=False)

    source = staging / "bmc_analysis/paired_archive/reproduce_from_repo.py"
    spec = importlib.util.spec_from_file_location("goco_pairs", source)
    pairs = importlib.util.module_from_spec(spec); spec.loader.exec_module(pairs)
    x = pairs.extract(ROOT)
    if len(x) != 800: raise ValueError("Expected 800 paired original split differences")
    tab3 = pairs.summarize(x)
    tab3.to_csv(out / "Table3_regenerated.csv", index=False)

    raw = read_csv(staging / "bmc_analysis/ordering_ablation/raw_100split_4panels.csv")
    if len(raw) != 2800: raise ValueError("Expected 2800 arm/hash results for 400 splits")
    g = raw.groupby(["dataset", "seed", "arm"], as_index=False).agg(
        supported=("go_yield", "mean"), total=("total_calls", "mean"))
    rows = []; ablation = []
    for panel in PANELS:
        learn = g[(g.dataset == panel) & (g.arm == "GoCo-learned")].set_index("seed")
        for control in ("random", "score-only"):
            z = g[(g.dataset == panel) & (g.arm == control)].set_index("seed")
            d = learn.supported - z.supported
            if len(d) != 100: raise ValueError(panel + "/" + control)
            mean = float(d.mean()); se = float(d.std(ddof=1)/np.sqrt(len(d)))
            rows.append(dict(panel=panel, control=control, mean_difference=mean,
                mc_lower=mean-1.96*se, mc_upper=mean+1.96*se, positive_splits=int((d>0).sum()),
                mean_total_difference=float((learn.total-z.total).mean())))
            ablation.append((panel, control, d.to_numpy(float)))
    pd.DataFrame(rows).to_csv(out / "Table4_regenerated.csv", index=False)
    return x, ablation

def plot_new(pairs, ablation, dst):
    dst.mkdir(parents=True, exist_ok=True)
    compare = [
      [("M","Wainberg","GoDag60"),("M","Sanger","GoDag60"),("M","DRIVE","GoDag60"),("M","HAP1","GoDag60")],
      [("N","funmap","GoDagDense"),("N","string","GoDagDense")]]
    fig, axs = plt.subplots(2,1,figsize=(8.1,5.2),layout="constrained")
    for ax, group in zip(axs, compare):
        for i,(family,panel,control) in enumerate(group):
            d = pairs[(pairs.family == family)&(pairs.dataset == panel)&(pairs.control == control)]
            x=d.sort_values("split").delta_supported.to_numpy(float)
            if len(x)!=100: raise ValueError((panel,control,len(x)))
            ax.plot(x,np.full(len(x),i),".",alpha=.25,ms=2)
            ax.errorbar(x.mean(),i,xerr=1.96*x.std(ddof=1)/10,fmt="o",color="black",capsize=3)
        ax.set_yticks(range(len(group)),[p[1] for p in group]); ax.invert_yaxis()
        ax.axvline(0,ls="--",color="0.45"); ax.grid(axis="x",alpha=.15)
        ax.set_xlabel("Paired change in reference-supported released calls")
    fig.savefig(dst/"FigureS3_regenerated.pdf",bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8,4),layout="constrained")
    for i,panel in enumerate(PANELS):
        for shift,control,marker in ((-.15,"random","o"),(.15,"score-only","s")):
            x=next(v for name,c,v in ablation if name==panel and c==control)
            ax.errorbar(x.mean(),i+shift,xerr=1.96*x.std(ddof=1)/10,
                        fmt=marker,capsize=3,label=control if i==0 else None)
    ax.set_yticks(range(4),PANELS);ax.invert_yaxis();ax.axvline(0,ls="--",color="0.45")
    ax.set_xlabel("Learned ordering minus control: reference-supported calls")
    ax.legend();ax.grid(axis="x",alpha=.15)
    fig.savefig(dst/"FigureS4_regenerated.pdf",bbox_inches="tight");plt.close(fig)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("generated/bmc"))
    parser.add_argument("--skip-original-figures",action="store_true")
    args=parser.parse_args()
    dst=args.output.resolve();dst.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="goco_bmc_") as tmp:
        with tarfile.open(EXT,"r:gz") as archive:
            for mem in archive.getmembers():
                pp=Path(mem.name)
                if pp.is_absolute() or ".." in pp.parts: raise ValueError(mem.name)
            archive.extractall(tmp,filter="data")
        pairs, ablation=make_tables(dst,Path(tmp))
        plot_new(pairs,ablation,dst/"figures")
    if not args.skip_original_figures:
        env={**os.environ,"GOCO_FIGURE_DIR":str(dst/"figures")}
        subprocess.run([sys.executable,str(ROOT/"code/make_figures.py")],check=True,cwd=ROOT,env=env)
    for name in ("Table1_regenerated.csv","Table2_regenerated.csv","Table3_regenerated.csv","Table4_regenerated.csv"):
        if not (dst/name).is_file():raise RuntimeError(name)
    print("PASS: generated four numeric BMC tables and supplementary figures S3/S4 from original split-level results")
    print("Input limitation: complete original frozen FunMap/STRING per-call score files are not hosted;")
    print("their results are reconstructed from frozen algorithm split-level outputs, not re-fitted.")
if __name__=="__main__":main()
