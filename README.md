# ConformalGoCo (GoCo)

Reference implementation, frozen inputs and complete split-level results for

> **GoCo: co-essentiality-guided conformal calibration for Gene Ontology annotation release.**
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

## Data provenance

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
