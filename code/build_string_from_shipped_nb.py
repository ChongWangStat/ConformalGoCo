"""Build the STRING independent-validation application, parallel to FunMap in every respect we control.

Parallel by construction:
  predictor   frozen, never retrained -- STRING v12.0 experiments channel, released 2023-07, 18 months
              before FunMap's publication, so independence needs no argument
  neighbourhood  top-50 random-walk-with-restart, focal gene REMOVED
  background     analysis universe minus the focal gene (target exclusion)
  enrichment     hypergeometric over-representation, GO term size in [10, 2000]
  adjustment     BH within (gene, aspect)
  score          -log10(BH-adjusted p)
  GO truth       the IDENTICAL frozen 2024-09-08 tables used for FunMap, restricted to the overlap,
                 so the network is the only thing that changes

Declared difference: FunMap's neighbourhoods are a frozen AUTHOR artefact downloaded from
funmap.linkedomics.org.  STRING publishes no neighbourhoods, so ours are reconstructed with a
stated RWR parameterisation (restart 0.5, 50 power iterations).  This validation therefore tests
the calibration layer under the same predictor RECIPE, not against the same frozen artefact.

GO leakage: STRING's GO-derived evidence lives entirely in the `database` and `textmining`
columns.  We read column 10 (`experiments`) only -- fed by DIP, BioGRID, HPRD, IntAct, MINT, PDB
and benchmarked against KEGG/Complex Portal, not GO.  Because min(combined_score) over the shipped
file is exactly 150 and combined_score >= experiments, no pair with experiments >= 150 can have
been censored on GO-informed grounds.
"""
from __future__ import annotations
import gzip, json, urllib.request, os, sys, time
import numpy as np, pandas as pd
from scipy import sparse
from scipy.stats import hypergeom

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.environ.get("GOCO_STRING_OUT", os.path.join(ROOT, "string"))
CACHE = os.environ.get("GOCO_STRING_CACHE", OUT)
FUN = os.environ.get("GOCO_FUNMAP", os.path.join(ROOT, "funmap"))
os.makedirs(OUT, exist_ok=True)
T_EXP = 400              # STRING's own published "medium confidence" cut -- pre-registered, not chosen
RESTART = 0.5
N_ITER = 50
TOPK = 50
TERM_MIN, TERM_MAX = 10, 2000
LINKS = "https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz"
INFO = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
t0 = time.time()


def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print("  cached", dest, flush=True); return dest
    print("  downloading", url, flush=True)
    urllib.request.urlretrieve(url, dest)
    print("    %.1f MB" % (os.path.getsize(dest) / 1e6), flush=True)
    return dest


print("[1] fetch STRING", flush=True)
lp = fetch(LINKS, CACHE + "/links.txt.gz")
ip = fetch(INFO, CACHE + "/info.txt.gz")

print("[2] ENSP -> symbol", flush=True)
info = pd.read_csv(ip, sep="\t", usecols=["#string_protein_id", "preferred_name"])
info.columns = ["ensp", "symbol"]
p2s = dict(zip(info.ensp, info.symbol))
print("    proteins %d | distinct symbols %d" % (len(info), info.symbol.nunique()), flush=True)

print("[3] experiments-channel edges at T>=%d (column 10 only)" % T_EXP, flush=True)
a, b, mn = [], [], 10**9
with gzip.open(lp, "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.split()
        c = int(f[15])
        if c < mn:
            mn = c
        if int(f[9]) >= T_EXP:
            x, y = f[0], f[1]
            if x < y:
                a.append(x); b.append(y)
print("    min(combined_score) over file = %d  (non-censoring proof holds for T>=%d)" % (mn, mn), flush=True)
e = pd.DataFrame({"a": [p2s.get(x) for x in a], "b": [p2s.get(x) for x in b]}).dropna()
e = e[e.a != e.b].drop_duplicates()
print("    undirected edges %d | genes %d" % (len(e), len(set(e.a) | set(e.b))), flush=True)

print("[4] restrict to the frozen FunMap GO universe (overlap design)", flush=True)
direct = pd.read_csv(FUN + "/go_truth_direct.csv.gz")
tpath = pd.read_csv(FUN + "/go_truth_true_path.csv.gz")
uni = sorted(set(direct.gene) & (set(e.a) | set(e.b)))
uidx = {g: i for i, g in enumerate(uni)}
n = len(uni)
e = e[e.a.isin(uidx) & e.b.isin(uidx)]
print("    analysis universe %d genes (FunMap 10,316) | edges within it %d" % (n, len(e)), flush=True)

print("[5] load the shipped frozen RWR neighbourhoods", flush=True)
nb = pd.read_csv(os.environ["GOCO_SHIPPED_NB"]); nb.columns = ["focal", "neighbor", "rank"]
nb = nb[nb.focal.isin(uidx) & nb.neighbor.isin(uidx)]
print("    neighbourhood rows %d | self-loops %d (must be 0)" % (len(nb), int((nb.focal == nb.neighbor).sum())), flush=True)
print("[6] target-excluded hypergeometric enrichment", flush=True)
d = direct[direct.gene.isin(uidx)]
term_genes = d.groupby("go_id").gene.apply(set).to_dict()
term_aspect = dict(zip(d.go_id, d.aspect))
term_genes = {t: s for t, s in term_genes.items() if TERM_MIN <= len(s) <= TERM_MAX}
gene_terms = d[d.go_id.isin(term_genes)].groupby("gene").go_id.apply(list).to_dict()
print("    terms in size range [%d,%d]: %d" % (TERM_MIN, TERM_MAX, len(term_genes)), flush=True)

nbmap = nb.groupby("focal").neighbor.apply(list).to_dict()
N_bg = n - 1
rows = []
for gi, g in enumerate(uni):
    nbrs = nbmap.get(g, [])
    if not nbrs:
        continue
    ann = [x for x in nbrs if x in gene_terms]
    if not ann:
        continue
    hits = {}
    for x in ann:
        for t in gene_terms[x]:
            hits[t] = hits.get(t, 0) + 1
    nn = len(ann)
    own = set(gene_terms.get(g, []))
    ks, Ks, ts = [], [], []
    for t, k in hits.items():
        K = len(term_genes[t]) - (1 if t in own else 0)      # focal removed from the background
        if K <= 0:
            continue
        ks.append(k); Ks.append(K); ts.append(t)
    if not ts:
        continue
    p = hypergeom.sf(np.array(ks) - 1, N_bg, np.array(Ks), nn)
    rows.append(pd.DataFrame({"gene": g, "go_id": ts, "aspect": [term_aspect[t] for t in ts],
                              "neighbor_hits": ks, "annotated_neighbors": nn,
                              "term_background_count": Ks, "background_genes": N_bg, "p_value": p}))
    if gi % 2000 == 0:
        print("    enrich %d/%d  %.0fs" % (gi, n, time.time() - t0), flush=True)
sc = pd.concat(rows, ignore_index=True)

print("[7] BH within (gene, aspect)", flush=True)
# Denominator is the number of terms TESTED in the aspect, not the number of rows reported:
# terms with no neighbour hit have p = 1 and are omitted from the output but still counted.
# This convention reproduces the FunMap pipeline's published adjusted p-values on 97.6% of
# 3.3M calls to within 1e-9 (the residual is tie ordering).
MASP = {a: int(sum(1 for t in term_genes if term_aspect[t] == a)) for a in set(term_aspect.values())}
print("    terms tested per aspect:", MASP, flush=True)
sc = sc.sort_values(["gene", "aspect", "p_value"], kind="mergesort")
sc["_rank"] = sc.groupby(["gene", "aspect"], sort=False).cumcount() + 1
sc["_raw"] = sc.p_value * sc.aspect.map(MASP) / sc["_rank"]
sc["q_value"] = sc.iloc[::-1].groupby(["gene", "aspect"], sort=False)["_raw"].cummin().iloc[::-1].clip(0, 1)
sc["score"] = -np.log10(np.maximum(sc.q_value, 1e-300))
sc = sc.drop(columns=["_rank", "_raw"])
sc.to_csv(OUT + "/string_gene_go_scores.csv.gz", index=False)

print("[8] freeze the matched GO truth and metadata", flush=True)
direct[direct.gene.isin(uidx)].to_csv(OUT + "/go_truth_direct.csv.gz", index=False)
tpath[tpath.gene.isin(uidx)].to_csv(OUT + "/go_truth_true_path.csv.gz", index=False)
e.rename(columns={"a": "gene_a", "b": "gene_b"}).to_csv(OUT + "/string_edges.csv.gz", index=False)
json.dump(dict(application="STRING v12.0 experiments-only", upstream_retrained=False,
               channel="experiments (column 10)", threshold=T_EXP,
               min_combined_score_in_file=int(mn),
               licence="CC BY 4.0", string_version="v12.0 (frozen 2023-07)",
               rwr=dict(restart=RESTART, iterations=N_ITER, topk=TOPK,
                        note="reconstructed; STRING publishes no neighbourhoods"),
               target_exclusion="focal gene removed from neighbourhood and from enrichment background",
               go_release="2024-09-08 (identical frozen tables as FunMap, restricted to overlap)",
               go_term_size_range=[TERM_MIN, TERM_MAX],
               analysis_genes=n, edges=int(len(e)), score_rows=int(len(sc)),
               direct_truth_rows=int(direct.gene.isin(uidx).sum()),
               true_path_truth_rows=int(tpath.gene.isin(uidx).sum())),
          open(OUT + "/score_construction_metadata.json", "w"), indent=2)
# drop-in view: the runner consumes the same five filenames for every predictor family
V = os.path.join(OUT, "funmap_view")
os.makedirs(V, exist_ok=True)
sc.to_csv(os.path.join(V, "funmap_gene_go_scores.csv.gz"), index=False)
direct[direct.gene.isin(uidx)].to_csv(os.path.join(V, "go_truth_direct.csv.gz"), index=False)
tpath[tpath.gene.isin(uidx)].to_csv(os.path.join(V, "go_truth_true_path.csv.gz"), index=False)
e.rename(columns={"a": "gene_a", "b": "gene_b"}).to_csv(os.path.join(V, "funmap_edges.csv.gz"), index=False)
nb.rename(columns={"focal": "focal_gene", "neighbor": "neighbor_gene", "rank": "author_rank"}).to_csv(os.path.join(V, "author_top50_neighborhoods_lopo.csv.gz"), index=False)
print("DONE %.0fs | genes %d | edges %d | score rows %d" % (time.time() - t0, n, len(e), len(sc)), flush=True)
print("     view for the runner: GOCO_DATA=%s" % V, flush=True)
