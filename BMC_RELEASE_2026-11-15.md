# GoCo for BMC Bioinformatics — 2026-11-15

The current submitted manuscript uses **alpha = 0.10 and delta = 0.10** as its principal risk setting. The original `delta=0.50` runs remain available as supplementary benchmark results.

The repository's original `code/`, `w/`, `ext/`, `funmap/`, `string/`, `splits/`, and `results/` contain the frozen predictor code, fold assignments, and archived benchmark results. The two new extension artifacts are:

- `bmc_analysis/GoCo_BMC_Reproducibility_Extension_2026-11-15.tar.gz`: exact scripts, 400-split/2,800-row matched-resolution ordering ablation, 800 paired archive differences, detailed protocol and verification scripts.
- `paper/bmc_2026-11-15/GoCo_BMC_Official_Submission_2026-11-15.tar.gz`: complete current official Springer Nature LaTeX manuscript and supplement, publisher class, editable files, figures, PDFs, and approved changelog.

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

Refer to the repository's **versioned `v4.1.0-bmc` release / commit** when available, and use GitHub as the software source. An earlier Zenodo archive may exist, but Zenodo is **not required** for BMC and automatic release synchronization was not verified for this update. The current manuscript cites GitHub only.

This release preserves the previously approved AI disclosure and original fourfold contributions. The final Discussion removes only the two redundant sentences approved on 2026-11-15.