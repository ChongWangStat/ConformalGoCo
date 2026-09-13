# ConformalGoCo (GoCo)

Reference implementation, analysis inputs and complete split-level outputs for

> **GoCo: structure-guided conformal calibration for error-controlled Gene Ontology annotation release.**
> Chong Wang, Yongzhao Shao, Peng Liu.

GoCo is a **calibration and release layer** for a frozen Gene Ontology (GO) predictor. It retrains
nothing. Given a gene-by-term score matrix and, for each candidate call, the evidence objects that
produced it, GoCo decides *which predicted annotations to release* under a certified bound on the
gene-level false discovery proportion (FDP).

**The problem.** Module-based predictors give every member of a module the same score for a term, so
thousands of gene-term pairs tie exactly. A single genome-wide cutoff must then take or drop a whole
tied block: it either overshoots the error target or leaves the error budget unspent. In the Wainberg
data the grid step from 800 to 775 moves the full-panel gene-level FDP from 0.092 to 0.140 in one jump,
with no threshold in between.

**The method.** GoCo keeps the global grid and inserts, inside each step, policies that relax only a
fraction of the affected genes, admitted in order of a per-call learned prediction of the loss increment
per supported term. The ordering is fitted on a *ranking fold* whose labels never touch certification.
A fixed-sequence Learn-Then-Test certificate on an independent *certification fold* controls the risk
for **any** ordering computed from the ranking fold, so the biology buys yield, never validity.

## Results

`results/goco_results_100splits.csv` holds every reported number: 3 methods x 4 data sets x 3 targets x 2 delta x 100 splits (7,200 rows).

At alpha = 0.10, over 100 common splits, reference-supported GO-term yield relative to a global
threshold (GoDag): **+8.7% (Wainberg), +7.6% (Sanger), +12.8% (DRIVE), +6.3% (HAP1)**, with mean
held-out gene-level FDP at or below target. Wainberg is the development data set; the other three are
independent dependency resources to which the frozen construction was applied unchanged.

## Layout

    code/      reference implementation and the scripts that build every table and figure
    w/         frozen Wainberg inputs (scores, source attribution, GO truth, modules)
    ext/       frozen Sanger / DRIVE / HAP1 score and loss matrices
    results/   split-level output for every method x data set x alpha x delta (100 splits each)
    funmap/    frozen FunMap objects: per-call sources, GO truth, edges
               (the score table ships with the Zenodo archive, see below)
    string/    derived STRING validation objects (rebuild the rest with build_string_validation.py)
    paper/     tables and figures regenerated from results/

## Quick start

    pip install -r code/requirements.txt
    cd code
    python make_tables.py
    python make_figures.py
    python make_supp_tables.py

To recompute the results themselves from the frozen inputs (hours, not minutes):

    cd code
    python goco_rerun.py --tag run

This writes `results/run_*_100splits.csv` for Boger et al., GoDag and GoCo at alpha in {0.05, 0.10, 0.20} and
delta in {0.10, 0.50}; `results/goco_results_100splits.csv` is the frozen copy of that output.

## The no-leakage check

Validity requires that the ordering be a function of the ranking fold alone. That is testable, and the
test ships with the code:

    cd code && python test_measurability.py Wainberg 0

It scrambles the labels of **every** non-ranking-fold gene and verifies that the ordering of pool genes
is bit-identical:

    block 16 goco       identical under scrambling: True (max |diff| = 0.00e+00)
    PASS

## What maps to what

| Paper item | Built by | From |
|---|---|---|
| Tables 2-4, S1-S2 | `code/make_tables.py` | `results/goco_results_100splits.csv`, `results/first_draft_primary_method_summary_long.csv` |
| Table S4 (risk decomposition) | `code/make_supp_tables.py` | `results/all_methods_harmonised_100splits.csv` |
| Figures 1-4, S1-S2 | `code/make_figures.py` | `results/all_methods_harmonised_100splits.csv`, frozen inputs |
| Table S5 (block diagnostics), Table S6 (case study) | `code/make_block_diagnostics.py`, `code/make_case_study.py` | frozen inputs |

`CODE_FREEZE.md` lists SHA-256 hashes for every code and result file.

## A second predictor family: GoCo-M and GoCo-N

The framework asks a predictor for a score and the evidence behind each call. The only property that
decides the candidate path is whether that evidence is **shared between genes**:

| structure | evidence | score ties | instantiation |
|---|---|---|---|
| **module-shared** | a group many genes belong to (co-essential module, complex, cluster) | exact, by construction | **GoCo-M** |
| **neighbourhood** | assembled per gene (random-walk or interaction neighbourhood) | none material | **GoCo-N** |

Module-shared evidence gives every member the same score, so the budget is trapped inside tied
blocks and the ordering has nothing finer than the gene to act on: GoCo-M relaxes a fraction of the
genes in a block. Neighbourhood evidence is gene-specific, so a refined global grid already recovers
that budget and the ordering is better spent on individual calls: GoCo-N releases the top-*m* calls
under the same statistic. Model, features, cross-fitting and certificate are identical; only the
resolution of the ordering changes. Which applies is read off the predictor by tabulating score
multiplicities, before any calibration.

**Results** (50 splits, delta = 0.10, supported-call yield over the global grid):

| data | alpha | GoCo-M | GoCo-N |
|---|---|---|---|
| FunMap | 0.05 / 0.10 / 0.20 | +29.0% / +9.5% / +15.5% | **+39.5% / +30.1% / +41.1%** |
| STRING | 0.05 / 0.10 / 0.20 | +14.8% / +17.8% / +26.5% | **+48.1% / +41.1% / +64.6%** |

STRING is independent validation: an unrelated network frozen eighteen months before FunMap's
publication, restricted to its **experiments** channel so that no Gene Ontology evidence enters the
graph. The negative control matters as much as the result: run on Wainberg, GoCo-N *loses* 140
supported terms at alpha = 0.10 and is ahead on none of 50 splits, because module-shared evidence
traps the budget where no global path can reach it. The two instantiations are structure-specific,
not ranked.

    python code/run_goco_second_family.py --repeats 50 --out results/second_family/goco_funmap_50splits.csv
    GOCO_DATA=string/funmap_view python code/run_goco_second_family.py --repeats 50 --out results/second_family/goco_string_50splits.csv
    python code/wainberg_goco_n_control.py --repeats 50 --out results/second_family/wainberg_goco_n_control_50splits.csv
    python code/test_measurability_second_family.py          # ranking-fold measurability, both applications

`GOCO_DATA` points at a directory holding the five frozen objects a predictor must supply; `funmap/`
is that directory for FunMap and the default.

### Rebuilding the STRING validation set

STRING's raw download (162 MB) and the derived score table (62 MB) are **not** shipped; they are
rebuilt in about five minutes, which also re-verifies the leakage argument end to end:

    python code/build_string_validation.py

It downloads `9606.protein.links.full.v12.0.txt.gz` (CC BY 4.0), keeps **column 10 only** -- the
experiments channel, fed by physical assays and calibrated against KEGG and Complex Portal, not GO --
thresholds at STRING's published medium-confidence cut, maps ENSP to symbol, builds top-50
random-walk neighbourhoods with the focal gene removed, and runs the same target-excluded
hypergeometric enrichment. It prints `min(combined_score)` over the shipped file, which is exactly
150: since `combined_score >= experiments`, no pair above that floor can have been censored on
GO-derived grounds. `string/` holds the small derived objects (edge list, neighbourhoods, metadata).

GO enters STRING only through the `database` and `textmining` columns, which the script never reads.
Networks whose edge weights are *fitted* to GO co-annotation -- HumanNet, HumanBase/GIANT, GeneMANIA --
were rejected as validation sets for exactly this reason: evaluating GO predictions against a network
built from GO is circular.

### The one file not committed here

`funmap/funmap_gene_go_scores.csv.gz` (48,490,602 bytes, sha256 `b215cf4ca8362efcd94de1e287759882d5d27c4bac932e184eaaf9c239943f6c`) is the FunMap
gene-GO enrichment score table. It is distributed with the Zenodo archive (https://doi.org/10.5281/zenodo.22677009) rather than
committed to git, which keeps this repository to code and small artefacts. Everything else the second family
needs -- per-call source sets, network edges and the frozen GO truth -- is in `funmap/`.

Place it in `funmap/` before running the FunMap arm; `code/verify_second_family.py` skips that arm
cleanly if it is absent. To confirm you have the right file:

    python -c "import hashlib;print(hashlib.sha256(open('funmap/funmap_gene_go_scores.csv.gz','rb').read()).hexdigest())"

STRING needs nothing extra: `code/build_string_validation.py` downloads and rebuilds its inputs,
and that rebuild was checked to reproduce the shipped results exactly (identical integer columns,
p/q to 4e-13, identical certified policies on every re-run split).

### Reproducibility gate

`code/verify_second_family.py` re-runs every arm from clean and requires bit-identical agreement with
the shipped split-level results on every metric of every split (`atol = 0`, no summary comparison).

## Data provenance

The second predictor family: the FunMap network and its top-50 random-walk neighbourhoods are the authors' own release (funmap.linkedomics.org), read with a leave-one-gene-out wrapper that removes the focal gene from its neighbourhood and from the enrichment background; the GO release is the archived 2024-09-08 snapshot, which predates the FunMap paper's acceptance. STRING v12.0 is CC BY 4.0 and is downloaded by `build_string_validation.py` rather than redistributed here. The STRING panel is restricted to the genes covered by the same frozen GO tables used for FunMap, so that the network is the only thing that differs between the two applications.

All underlying resources are public: Wainberg co-essential modules and their CRISPR dependency
profiles; Sanger Project Score; Project DRIVE; and a HAP1 genetic-interaction network. The files in
`w/` and `ext/` are analysis-ready objects derived from those sources. `w/wainberg_gene_go_scores.csv.gz` holds the frozen
gene-GO scores at or above the lowest grid threshold (600), which is the whole range the analysis uses; the
reproduction gate in `goco_rerun.py` checks these rows against the source-attribution file on every run.
`code/goco_rerun.py`
re-verifies the external matrices against the frozen originals on every run and reports the maximum
absolute difference.

## Licence and citation

Code released under the MIT Licence (see `LICENSE`). If you use this software, please cite both the
paper and the archived release (see `CITATION.cff`).
