# ConformalGoCo (GoCo)

Reference implementation and manuscript reproducibility materials for:

> **GoCo: structure-guided conformal calibration for error-controlled Gene Ontology annotation release**  
> Chong Wang, Yongzhao Shao, Peng Liu

GoCo calibrates annotation release from an existing, frozen predictor; it does not retrain the predictor. Its stated asymptotic guarantee concerns the **average per-gene FDP over the target pool**, with abstentions contributing zero. This is not a guarantee for each called gene or for pooled call-level precision.

**GoCo-M** inserts partial-admission policies within score-grid steps, using shared module evidence. **GoCo-N** orders individual calls using gene-specific neighbourhood evidence. The ordering is learned on the ranking fold, then frozen before fixed-sequence certification.

At the primary operating point (`alpha=0.10`, `delta=0.50`), GoCo-M released 6.4--17.6% more reference-supported calls than Multilabel across Wainberg, Sanger, DRIVE and HAP1. GoCo-N's gains over Multilabel on the same refined grid were 33.8% on FunMap and 36.5% on STRING. These are fixed-panel benchmark comparisons, not universal performance guarantees.

## Current file layout

```text
code/          calibration, analysis, figure/table builders and verification
w/             frozen Wainberg inputs and source reconstruction
ext/           frozen Sanger, DRIVE and HAP1 inputs
funmap/        frozen FunMap network, neighbourhood and GO-reference objects
string/        STRING experiments-channel objects and large-input checksums
splits/        explicit T/C/E assignments for all six panels and 100 streams
results/       current split-level results (sf/ = neighbourhood; wf/ = modules)
paper/tables/  current main Tables 1--5 and Supplementary Tables S1--S10
paper/figures/ current Figures 1--4 and Supplementary Figure S1
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

The table builder recalculates the numeric cells from the stored analysis outputs. The authors' current captions, notes, and expected displayed numbers in `paper/table_specs.json` are used for formatting and verification, not as calculation inputs. Figures are drawn by Python/Matplotlib; no generative-image model or pre-existing illustration is used by their builders. Figure 4 uses **Multilabel** as its percentage-gain denominator.

## Recompute analyses

```bash
python code/goco_rerun.py --tag run
GOCO_DATA=funmap python code/run_goco_second_family.py --repeats 100 --out results/sf/funmap_v4_100.csv
GOCO_DATA=string/funmap_view python code/run_goco_second_family.py --repeats 100 --out results/sf/string_v4_100.csv
```

The neighbourhood commands require the large matching input files described below. Ten-fold ranking-fold selection is implemented in `rehearsal_kfold.py` and `rehearsal_kfold_sf.py`; the two-half companion files are retained because the raw-run verification gate also checks them. The frozen generic module-calibrator summary is `results/generic_module_summary.csv`.

```bash
python code/test_measurability.py Wainberg 0
python code/test_measurability_second_family.py
python code/verify_second_family.py --repeats 3
```

A summary/table rebuild is not a full raw-input rerun. The verification tools report failures rather than silently substituting summaries or incomplete runs.

## Large inputs and archive

`funmap/funmap_gene_go_scores.csv.gz` and the large files in `string/funmap_view/` are intentionally not committed to Git. Use the matching complete reproducibility archive under the stable Zenodo concept DOI **10.5281/zenodo.22677009** and verify the supplied input checksums. Updating this repository does not itself update the Zenodo deposit. The small STRING-view checksum and instruction files remain here so omitted inputs are explicit.

`CODE_FREEZE.md` records checksums for the current tracked files. Original data-source, GO snapshot, and target-exclusion details are documented in the manuscript and source code.

## License and citation

Code is released under the MIT License. Cite the manuscript and archived software using `CITATION.cff`, which retains the stable concept DOI.
