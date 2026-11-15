# GoCo for BMC Bioinformatics — 2026-11-15

The current submitted manuscript uses **alpha = 0.10 and delta = 0.10** as its principal risk setting. The original `delta=0.50` runs remain available as supplementary benchmark results.

The repository's original `code/`, `w/`, `ext/`, `funmap/`, `string/`, `splits/`, and `results/` contain the frozen predictor code, fold assignments, and archived benchmark results. The two new extension artifacts are:

- `bmc_analysis/GoCo_BMC_Reproducibility_Extension_2026-11-15.tar.gz`: exact scripts, 400-split/2,800-row matched-resolution ordering ablation, 800 paired archive differences, detailed protocol and verification scripts.
- `bmc_analysis/GoCo_BMC_Results_Tables_Figures_2026-11-15.tar.gz`: Tables 1–4 as CSV, new Figures S3/S4 as PDF and a plotting script; **no manuscript, TeX, or publisher class**.

## Verify the extension

```bash
mkdir -p bmc_bundled
cd bmc_bundled
tar -xzf ../bmc_analysis/GoCo_BMC_Reproducibility_Extension_2026-11-15.tar.gz
python bmc_analysis/verify_ordering_ablation.py
python bmc_analysis/paired_archive/test_reproduce.py
```

To reproduce raw module ablations from the frozen inputs in the repository root (more time-consuming), run:

```bash
python bmc_bundled/bmc_analysis/run_matched_ablation.py --repo . --families module --datasets Wainberg,Sanger,DRIVE,HAP1 --repeats 100 --random-reps 3 --out bmc_bundled/raw_repeat
```

These analyses use frozen predictor outputs and preserve the same ranking, certification and evaluation folds. The control results do not claim a new independent cohort. Raw FunMap/STRING score inputs needed for a new neighborhood ordering ablation are not bundled; the already archived neighbourhood benchmark results are included in the main repository.

## Citation and exact release

Refer to the repository's **versioned `v4.1.1-bmc` release / commit** when available, and use GitHub as the software source. The BMC submission cites this pinned GitHub version; no Zenodo record is required or claimed for the revision.

This repository contains reproducibility materials only; the author-approved manuscript, LaTeX and supplement are distributed directly with the journal submission. Results, tables and figures are reproducible from the frozen code and inputs in this release.