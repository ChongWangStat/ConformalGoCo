# Code and results freeze

Date: 2026-09-09. This release implements a single method, GoCo, and compares it with two references:
the released Boger et al. calibration on the Direct loss (baseline, ontology ignored) and GoDag
(the same certificate on the GO-hierarchy-aware TruePath loss over one global threshold).

All reported numbers are in `results/goco_results_100splits.csv`: 3 methods x 4 data sets x 3 targets
x 2 delta values x 100 splits = 7,200 rows.

## Code

| file | bytes | sha256 |
|---|---:|---|
| `code/goco_learned.py` | 8,804 | `b55cfd7b021732723e9ec1d3cdc9528f0a0860f3e7b3685117cf2d701cc3be57` |
| `code/goco_rerun.py` | 17,728 | `4fa69e253b29f712b9d752c0802c1dda61b1a042d32387c16fd157e559acde88` |
| `code/make_block_diagnostics.py` | 9,710 | `4a9f377268b6f966605bb4cb371075da603929b088aea339c349f8660134baf2` |
| `code/make_case_study.py` | 5,840 | `4677456ee85130ea2f9340121f8b1c4b0cd5f6a5dcc57026799374a7934a84ae` |
| `code/make_figures.py` | 12,129 | `0e5c86541ee5b02ee044bb1fc3f89cf364b51807c76586719aa354d89d42a40f` |
| `code/make_supp_tables.py` | 5,726 | `b7fcd17f37c8d666a17fea9863cfec4067e82b982d06eda315e6e3d06788bb7d` |
| `code/make_tables.py` | 14,923 | `6ebf42c15959277ebf89c64e15dcede7225bd33e7adfadc24a40940abaedc927` |
| `code/requirements.txt` | 255 | `48e619ccd9c37e53ca901365a416c71d87b7be3a7375dc6a5f6d7ea8cf4e86cd` |
| `code/test_measurability.py` | 2,779 | `dc8eafbda47d2ae3607e697b4164db83af27d4b2842b8f6612bf7f3a589cce4a` |
| `code/wordcount.py` | 2,115 | `b1a7a2bfb99177751e2a94c8b3ea9af7af465e75d0277c9d3b2a079b7d81df92` |

## Results

| file | bytes | sha256 |
|---|---:|---|
| `results/goco_results_100splits.csv` | 1,872,781 | `be8297da81eac1946f1c28928c49e3dd04240183cd309a495926906e2aee7abd` |

## Entry points

- `python code/goco_rerun.py --tag run` reproduces all three methods at alpha in {0.05,0.10,0.20}, delta in {0.10,0.50} over the 100 frozen splits.
- `python code/test_measurability.py <dataset> <seed>` checks that scrambling every non-ranking-fold label leaves the GoCo ordering unchanged.
- `python code/make_tables.py` and `code/make_figures.py` rebuild the paper tables and figures into `paper/`.
