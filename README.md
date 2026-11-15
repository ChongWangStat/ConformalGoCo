# ConformalGoCo: reproducible algorithms, prediction scores and analysis code

**GoCo — structure-guided conformal calibration for error-controlled Gene Ontology annotation release**  
Chong Wang, Yongzhao Shao and Peng Liu

This repository distributes **the code and analysis inputs needed to regenerate** the GoCo results. It is **not** a second copy of the manuscript: no manuscript, publisher LaTeX, pre-rendered figures or publication tables are committed. The BMC manuscript is supplied separately to the journal.

GoCo operates on *frozen predictor scores*; the predictors are not retrained. GoCo-M uses module-shared evidence to construct partial admissions within tied score blocks. GoCo-N constructs call-level policies using gene-specific neighborhood evidence. Both use ranking-fold-only policy construction, fixed-sequence certification, and the gene-averaged GO TruePath loss. See `code/goco_rerun.py`, `code/goco_learned.py`, `code/goco_second_family.py`, and `code/run_goco_second_family.py` for the actual algorithms.

## Inputs: predictions, references and structural evidence — not raw training data

| Analysis | Frozen files required by the algorithm | Scope |
|---|---|---|
| Wainberg | `w/wainberg_gene_go_scores.csv.gz`, `w/wainberg_highscore_source_attribution.csv.gz`, `w/network_communities.csv`, `w/go_truth_direct.csv.gz`, `w/go_truth_true_path.csv.gz` | Predictor scores, the evidence behind each call and frozen GO references; no predictor training data |
| Sanger/DRIVE/HAP1 | `ext/*_external_go_scores_extended100.csv.gz`, frozen verification `.npz` matrices and the same GO reference | External frozen predictor scores and call-source evidence |
| FunMap | `funmap/author_top50_neighborhoods_lopo.csv.gz`, `funmap/funmap_edges.csv.gz`, GO reference files, archived `results/sf/funmap_v4_100.csv` | Neighborhood structure and exact 100-split algorithm output. **The complete original frozen per-call prediction-score file is not publicly tracked**, so a full *from-score* rerun is not available from this checkout alone |
| STRING | `string/string_top50_neighborhoods.csv.gz`, `string/string_edges.csv.gz`, original GO reference, `results/sf/string_v4_100.csv` | Experiments-only network evidence and archived 100-split output. **The complete predictor-score view is not tracked**, so a full *from-score* rerun additionally needs the frozen scores |
| All panels | `splits/*_splits.csv.gz`, `splits/manifest.json` | Original deterministic 100-fold T/C/E assignments |

The GO references, evidence mappings, and split identities **are necessary** to rerun the statistical release algorithms. No original neural-network/AI-model training records or upstream feature matrices are needed when frozen prediction scores are available. The archived split-level results are included as independently checkable outputs and as inputs to summary reconstruction when the two neighborhood score files are unavailable; they should not be described as raw-from-score replications.

## Regenerate numerical results and figures locally

```bash
python -m pip install -r code/requirements.txt
python code/rebuild_missing_splits.py
python code/reproduce_bmc.py --output generated/bmc
```

The last command **computes** the four BMC result-table numeric summaries from original 100-split outputs, recalculates the 800 paired differences and 400-split ordering ablation, and regenerates all original figure graphics plus the two new BMC supplementary plots. Outputs are written only to `generated/bmc/` and are Git-ignored. There are **no stored publication table CSVs or figure PDFs** to substitute for a calculation.

Original numerical methods tables and their independent checks can also be generated locally:

```bash
python code/make_paper_tables.py
python code/verify_release.py
```

These scripts write transient TeX/CSV reports into ignored local output directories; manuscript templates and previously printed tables are not inputs.

## Rerun the method from frozen module scores (slower)

```bash
python code/goco_rerun.py --datasets Wainberg,Sanger,DRIVE,HAP1 \
  --alphas 0.10 --deltas 0.10,0.50 --seeds 100 \
  --methods Boger,GoDag,GoDag-cert60,GoCoGrid:goco --tag fresh_bmc
```

To rerun the *fixed-resolution ordering ablation* from module scores, unpack the code/split-level checkpoint bundle and run:

```bash
mkdir -p generated/extension
python -m tarfile -e bmc_analysis/GoCo_BMC_Reproducibility_Extension_2026-11-15.tar.gz generated/extension
python generated/extension/bmc_analysis/run_matched_ablation.py --repo . \
  --families module --datasets Wainberg,Sanger,DRIVE,HAP1 --repeats 100 \
  --random-reps 3 --out generated/raw_ablation
```

The extension's split-level records are *derived algorithm diagnostics*, not manuscript tables. The three hash controls are averaged within split. Comparisons use the same ranking, certification and evaluation folds, and the common target is `alpha=delta=0.10`. Note that the Wainberg ablation rerun can differ slightly from the archived pre-revision source mean under different numerical environments; the reported paired controls use one internally consistent frozen rerun.

For a full *FunMap/STRING from-score* recomputation, supply the exactly frozen `funmap/funmap_gene_go_scores.csv.gz` and `string/funmap_view/funmap_gene_go_scores.csv.gz` (and the associated frozen GO/evidence view files). `code/run_goco_second_family.py` then performs the algorithm directly on those scores. Absent these inputs, the published neighborhood results remain reproducible at the split-level summary stage, **not** by rerunning the original predictor scores. `code/build_string_from_shipped_nb.py` can rebuild a new score view from network evidence, but it must not silently be substituted for the frozen archived view when testing exact result equality.

## What is in GitHub

- `code/` — loss, certification, module and neighborhood candidate policy algorithms; table and figure generation; original checks and `reproduce_bmc.py`
- `w/`, `ext/`, `funmap/`, `string/` — predictor-score outputs (when available), GO reference and structure evidence; no upstream predictor-training data
- `splits/` — T/C/E assignments and checksums for each 100-split panel
- `results/` — per-split algorithm outputs needed to reconstruct the reported numbers, particularly where large frozen neighborhood scores are not shared
- `bmc_analysis/GoCo_BMC_Reproducibility_Extension_2026-11-15.tar.gz` — raw split-level extension records, ablation/paired scripts, verification and data provenance (no manuscript table or figure PDF)

The pre-existing Git history remains available, but **the current release has no tracked `.tex`, manuscript PDF, table snapshots or rendered figures**. Reproduction outputs are created only when the scripts run.

## Exact software version and citation

Use the pinned [GoCo code-and-data release](https://github.com/ChongWangStat/ConformalGoCo/releases/tag/v4.1.3-bmc), not a moving copy of the main branch. Dependencies are frozen in `code/requirements.txt` (Python 3.13.9). The source code is MIT-licensed. No new Zenodo synchronization is claimed or required for this release.