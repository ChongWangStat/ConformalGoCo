# ConformalGoCo (GoCo)

Reference implementation and manuscript reproducibility materials for:

> **GoCo: structure-guided conformal calibration for error-controlled Gene Ontology annotation release**  
> Chong Wang, Yongzhao Shao, Peng Liu

GoCo calibrates annotation release from an existing, frozen predictor; it does not retrain the predictor. Its stated asymptotic guarantee concerns the **average per-gene FDP over the target pool**, with abstentions contributing zero. This is not a guarantee for each called gene or for pooled call-level precision.

**GoCo-M** inserts partial-admission policies within score-grid steps, using shared module evidence. **GoCo-N** orders individual calls using gene-specific neighbourhood evidence. The ordering is learned on the ranking fold, then frozen before fixed-sequence certification.

At the **principal BMC setting** (`alpha=delta=0.10`), GoCo-M released **6.3--15.6%** more reference-supported annotations than Multilabel across Wainberg, Sanger, DRIVE and HAP1. GoCo-N gained **34.0%** on FunMap and **35.4%** on STRING. The original `delta=0.50` analyses are retained as supplementary benchmarks.

## BMC Bioinformatics submission (15 November 2026)

The BMC revision adds matched-certification-fold comparisons on six panels and a fixed-resolution ordering ablation on four module panels (100 splits each). Both are documented and verified from archived split results.

- [Updated ablation and paired reproducibility bundle](bmc_analysis/GoCo_BMC_Reproducibility_Extension_2026-11-15.tar.gz)
- [Published BMC numerical Tables 1–4 and Figures S3–S4](bmc_analysis/GoCo_BMC_Results_Tables_Figures_2026-11-15.tar.gz)
- [Verification and reproduction instructions](BMC_RELEASE_2026-11-15.md)

The repository's existing code, frozen predictor inputs, 100-split fold assignments and original benchmark outputs are preserved. Only code, frozen inputs, data, table values, figures, and verification scripts are distributed here. Manuscript and LaTeX sources are **not published in this repository**. The revised data archive contains all 400 ablation splits, 800 paired differences, scripts and verification checks; the tables/figures bundle contains the four main result tables as CSV plus new supplementary figures.

## Current file layout

```text
code/          calibration, analysis, figure/table builders and verification
w/             frozen Wainberg inputs and source reconstruction
ext/           frozen Sanger, DRIVE and HAP1 inputs
funmap/        frozen FunMap network, neighbourhood and GO-reference objects
string/        STRING experiments-channel objects and large-input checksums
splits/        explicit T/C/E assignments for all six panels and 100 streams
results/       current split-level results (sf/ = neighbourhood; wf/ = modules)
paper/tables/  table-value JSON, specifications and locally generated outputs (tracked LaTeX table fragments removed)
paper/figures/ current Figures 1--4 and Supplementary Figure S1
bmc_analysis/   100-split ablation, paired contrasts, result table CSVs, and figures
```

There is one current set of results and publication outputs. Old 50-split results, Near-d sensitivity materials, duplicate table sets, and unused graphics are not retained in the active tree. Git history is preserved. The `v4` component in two result filenames identifies the frozen analysis specification, not a second copy of those results.

## Rebuild the publication outputs

Python 3.13 and the exact analysis dependencies are specified in `code/requirements.txt`.

```bash
pip install -r code/requirements.txt
python code/rebuild_missing_splits.py
python code/make_paper_tables.py
python code/make_figures.py
python code/verify_release.py
```

The table builder recalculates the numeric cells from the stored analysis outputs. The authors' current captions, notes, and expected displayed numbers in `paper/table_specs.json` are used for formatting and verification, not as calculation inputs. Figures are drawn by Python/Matplotlib; no generative-image model or pre-existing illustration is used by their builders. Figure 4 uses **Multilabel** as its percentage-gain denominator. The sensitivity-figure formatting and colours match the submitted version; the workflow remains the original Python-generated figure.

The archived September submission had two main tables and thirteen supplementary tables; the final November BMC numerical tables are shared as CSV in the results/figures bundle. Constants, the method inventory and the delta=0.50 sweep are Tables S11, S12 and S13, respectively. Internal intermediate table filenames in the calculation scripts retain their historical names; paper/table_specs.json maps them to the current publication filenames. Numerical calculations and stored analysis results are unchanged. The manuscript and its LaTeX source are provided directly with the journal submission, not in this software/data repository. Table-building scripts can generate TeX table snippets locally when needed; generated TeX is not tracked. Use the pinned GitHub `v4.1.1-bmc` code/data release for the BMC analysis; no publisher source is included.

## Recompute analyses

```bash
python code/goco_rerun.py --tag run

GOCO_DATA=funmap python code/run_goco_second_family.py --repeats 100 --out results/sf/funmap_v4_100.csv
GOCO_DATA=funmap python code/frontier_exhaustive_sf.py --run results/sf/funmap_v4_100.csv --out results/sf/funmap_v4_100.csv

GOCO_DATA=string/funmap_view python code/run_goco_second_family.py --repeats 100 --out results/sf/string_v4_100.csv
GOCO_DATA=string/funmap_view python code/frontier_exhaustive_sf.py --run results/sf/string_v4_100.csv --out results/sf/string_v4_100.csv

python code/rehearsal_kfold.py --repeats 100 --K 10 --out results/wf/wf_selection_k10.csv
GOCO_DATA=funmap python code/rehearsal_kfold_sf.py --repeats 100 --K 10 --out results/sf/funmap_v4_100_selection_k10.csv
GOCO_DATA=string/funmap_view python code/rehearsal_kfold_sf.py --repeats 100 --K 10 --out results/sf/string_v4_100_selection_k10.csv
```

The neighbourhood commands require the large matching input files described below. The exhaustive-frontier command is required after each raw neighbourhood run: it adds the `frontier_all` values used for the exhaustive matched-risk comparison. The paper table builder deliberately stops with an error if that column is missing rather than silently substituting the refined-grid frontier. The ten-fold commands reproduce the ranking-fold selection analysis reported in the supplement; the two-half companion files are retained because the raw-run verification gate also checks them. The secondary module-calibrator panel is represented by the frozen aggregate summary `results/generic_module_summary.csv`; it is not described as a released split-level primary-analysis output.

```bash
python code/test_measurability.py Wainberg 0
python code/test_measurability_second_family.py
python code/verify_second_family.py --repeats 3
```

A summary/table rebuild is not a full raw-input rerun. The verification tools report failures rather than silently substituting summaries or incomplete runs.

## Large inputs and archive

`funmap/funmap_gene_go_scores.csv.gz` and the large files in `string/funmap_view/` are intentionally not committed to Git. They are needed only for a complete raw-input rerun of the neighbourhood-family analyses; the repository retains the corresponding stored split-level outputs, reconstruction code, checksums, and verification tools used for the manuscript. The small STRING-view checksum and instruction files remain here so the omitted inputs are explicit.

`CODE_FREEZE.md` records checksums for the tracked release files. Original data-source, GO snapshot, and target-exclusion details are documented in the manuscript and source code.

## License and citation

Code is released under the MIT License. Cite the manuscript and the pinned GitHub `v4.1.1-bmc` code-and-results release using `CITATION.cff`.
