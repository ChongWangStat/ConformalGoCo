"""Write CODE_FREEZE.md: sha256 and size of every code file, every split-level result, every paper table, the explicit
splits and the frozen inputs, plus the two inputs that ship with the Zenodo archive rather than with git.
Usage: python code/make_code_freeze.py   (from anywhere)"""
import hashlib, os, datetime
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''): h.update(chunk)
    return h.hexdigest()
def walk(sub, exts=None, exclude=()):
    for root, _, files in sorted(os.walk(os.path.join(ROOT, sub))):
        if any(x in root for x in exclude): continue
        for f in sorted(files):
            if exts is None or f.endswith(exts): yield os.path.join(root, f)
L = [f"# Code and results freeze — release v4.0.0 ({datetime.date.today().isoformat()})", "",
     "Revision after the pre-submission mathematical and reproducibility review of 25 Sep 2026 (see README, 'Revision of 2026-09-25').",
     "Entry points: `code/goco_rerun.py --tag run` or `code/tie_audit.py` (module family, 100 splits), `code/run_goco_second_family.py`",
     "(neighbourhood family), `code/rehearsal_kfold.py` / `code/rehearsal_kfold_sf.py` (selection rule), `code/gocoN_detail.py`,",
     "`code/neard_run.py`, `code/generic_sf.py`, `code/uncalibrated_baselines*.py`, `code/frontier_exhaustive_sf.py`; tables by",
     "`code/make_tables.py`, `code/make_supp_tables.py`, `code/make_paper_tables.py`; gates `code/test_measurability*.py`,",
     "`code/verify_second_family.py`.", ""]
for title, files in [("Code", list(walk("code", (".py", ".txt"), exclude=("__pycache__", "/out")))),
                     ("Split-level results", list(walk("results", (".csv", ".txt", ".npz")))),
                     ("Paper tables and key numbers", list(walk("paper", (".tex", ".csv", ".json", ".md")))),
                     ("Explicit splits", list(walk("splits", (".gz",)))),
                     ("Frozen inputs committed to git", [p for sub in ("w", "ext", "funmap", "string") for p in walk(sub)
                                                         if not p.endswith(("funmap_gene_go_scores.csv.gz",)) and "/funmap_view/" not in p.replace(os.sep, "/") or p.endswith(("SHA256SUMS", "README.txt"))])]:
    L += [f"## {title}", "", "| file | bytes | sha256 |", "|---|---:|---|"]
    for p in files:
        L.append(f"| `{os.path.relpath(p, ROOT)}` | {os.path.getsize(p):,} | `{sha(p)}` |")
    L.append("")
L += ["## Inputs shipped with the Zenodo archive (not committed to git)", "",
      "| file | bytes | sha256 |", "|---|---:|---|"]
for p in [os.path.join(ROOT, "funmap", "funmap_gene_go_scores.csv.gz")] + list(walk("string/funmap_view", (".gz",))):
    if os.path.exists(p): L.append(f"| `{os.path.relpath(p, ROOT)}` | {os.path.getsize(p):,} | `{sha(p)}` |")
    else: L.append(f"| `{os.path.relpath(p, ROOT)}` | (absent in this checkout) | see string/funmap_view/SHA256SUMS |")
open(os.path.join(ROOT, "CODE_FREEZE.md"), "w").write("\n".join(L) + "\n")
print("CODE_FREEZE.md:", len(L), "lines")
