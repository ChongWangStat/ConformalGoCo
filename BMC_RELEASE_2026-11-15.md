# GoCo BMC Bioinformatics — code, frozen prediction inputs and reproducibility

Release **v4.1.3-bmc** provides *only* analysis algorithms, frozen AI predictor scores where available, GO reference annotations, module/neighborhood source evidence, deterministic split assignments and intermediate split-level algorithm outputs.

**Do not store completed manuscript tables, pre-rendered figures, manuscript PDFs, or LaTeX files on GitHub.**

## Rebuild all published numerical result summaries and figures
From the checked-out tagged release:
```bash
python -m pip install -r code/requirements.txt
python code/rebuild_missing_splits.py
python code/reproduce_bmc.py --output generated/bmc
```

The last command computes BMC numerical Tables 1–4 as local CSV files and regenerates all original plots plus the two new supplementary figures S3/S4 under an ignored directory. It does not read published tables or previously rendered figures. For the original supplementary numerical tables, run:
```bash
python code/make_paper_tables.py
python code/verify_release.py
```

## Regenerate algorithm results directly from frozen input scores

For the four module panels, `code/goco_rerun.py` reruns the GO-hierarchy-aware global threshold and structure-guided calibration algorithms from the frozen predicted scores and evidence at `w/` and `ext/`. The fixed-resolution ablation re-run code is in `bmc_analysis/GoCo_BMC_Reproducibility_Extension_2026-11-15.tar.gz`.

The original archived 100-split algorithm outputs (not publication table copies) remain in `results/` and in the BMC extension. They allow numerical checks and regenerated figures without retaining finished tables/figures.

**Transparency:** The exact frozen per-call FunMap and STRING score matrices are not currently tracked. The neighborhood-family published summaries can therefore be rebuilt from exact archived split-level outputs, but a full rerun of their policy selection *from original scores* requires those missing frozen score files. They must not be replaced silently with newly computed scores.

An upstream predictor's raw training records are **not** needed for GoCo calibration. See README.md for the explicit prediction-score, GO-reference, evidence, and fold schemas.

The scientific manuscript, supplementary PDF and publisher LaTeX are distributed separately with the journal submission. GitHub is the primary code archive. No Zenodo record is assumed to update automatically.