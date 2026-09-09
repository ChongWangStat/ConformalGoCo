# Code and results freeze

Date: 2026-09-09. This release implements a single method, GoCo, and compares it with two references:
the released Boger et al. calibration on the Direct loss (baseline, ontology ignored) and GoDag
(the same certificate on the GO-hierarchy-aware TruePath loss over one global threshold), plus a panel of
generic Direct-loss calibrators whose frozen summary is `results/first_draft_primary_method_summary_long.csv`.

All reported split-level numbers are in `results/goco_results_100splits.csv`: 3 methods x 4 data sets x 3 targets
x 2 delta values x 100 splits = 7,200 rows. The algorithm files `code/goco_rerun.py` and `code/goco_learned.py`
are unchanged since the first freeze (v1.0.0); later commits touched only table, figure and documentation scripts.

## Code

| file | bytes | sha256 |
|---|---:|---|
| `code/conformal_selection.py` | 3,659 | `96f815f1fdc632a1cca970584ac5d51f833b6b68719710720ee2df0a3be5a881` |
| `code/goco_learned.py` | 8,804 | `b55cfd7b021732723e9ec1d3cdc9528f0a0860f3e7b3685117cf2d701cc3be57` |
| `code/goco_rerun.py` | 17,728 | `4fa69e253b29f712b9d752c0802c1dda61b1a042d32387c16fd157e559acde88` |
| `code/make_block_diagnostics.py` | 9,150 | `e9997dd5bbe46700f437e26e1b5493743dc0720bf626eda738150dd552cc50b7` |
| `code/make_case_study.py` | 6,128 | `b08742ccc4edfd726ce7e70447e680526acdfc8cc3322cb681268302a5b5bc63` |
| `code/make_figures.py` | 11,509 | `2f8ff0ed1978e881f08122e31ca519e9b14dfad9c0c5a2e613334f8c4e075e70` |
| `code/make_supp_tables.py` | 2,814 | `d8b7719ef9829a546364740c1afff5a38911c6e4947880c05e1bf26410d75bbb` |
| `code/make_table_s9.py` | 2,957 | `042f70a094c13a8f1e89f4a392dc3ed5ed7a41fb9932a156ecd22b5a38476ea9` |
| `code/make_tables.py` | 14,211 | `826ab92da744ad644ae0cf8ed60330f4c9e4676f8f6301ecca02f93600f3e097` |
| `code/requirements.txt` | 255 | `48e619ccd9c37e53ca901365a416c71d87b7be3a7375dc6a5f6d7ea8cf4e86cd` |
| `code/test_measurability.py` | 2,567 | `753f9a0f8c217bed4173e41e8671ec0659c62e8aa95db024425e0beb1642318b` |
| `code/wordcount.py` | 2,115 | `b1a7a2bfb99177751e2a94c8b3ea9af7af465e75d0277c9d3b2a079b7d81df92` |

## Results

| file | bytes | sha256 |
|---|---:|---|
| `results/all_methods_harmonised_100splits.csv` | 1,962,374 | `9f7033f4abe51b9e875a477c3e3d568ef35748f96a1dd9f26694fe701fe213b5` |
| `results/block_diagnostics_alpha010.csv` | 1,027 | `7f37c2a9268884b4b2dc616c10adc130be3e408640da3b21206087d2c36d57a1` |
| `results/case_study_split0.csv` | 8,059 | `0e967053e841a1b862b174c1e0ef824a2bc94b5605f724fe1e13f465773ae5f4` |
| `results/conformal_selection_alpha0.1.csv` | 24,968 | `c7fa296fedd55df3ecbc126bcbcdb0ebcd91ec220481fee050aacda99fb4fcd4` |
| `results/first_draft_primary_method_summary_long.csv` | 8,768 | `a3fee7c60ef2ea6e35cbfc9f9ec2d35e0cac4cc54c05b872f5814c179abb6394` |
| `results/goco_results_100splits.csv` | 1,872,781 | `be8297da81eac1946f1c28928c49e3dd04240183cd309a495926906e2aee7abd` |
| `results/pool_level_exceedance_alpha010.csv` | 11,967 | `315e485ad4c1fd95ec72f2ce7b9e967a6a34afe3f82141cddef7bc090e4dc527` |
| `results/summary_by_method_harmonised.csv` | 6,439 | `df7db6860ec0dfe47652113fe33a0731bceeb73ae4d1790121470dfa925c9273` |
| `results/tie_counts.csv` | 271 | `424670da41826cb8dada0621172a212ab61e935dde7d4f379b9a8b1959156c89` |

## Entry points

- `python code/goco_rerun.py --tag run` reproduces all three methods at alpha in {0.05,0.10,0.20}, delta in {0.10,0.50} over the 100 frozen splits.
- `python code/test_measurability.py <dataset> <seed>` checks that scrambling every non-ranking-fold label leaves the GoCo ordering unchanged.
- `python code/conformal_selection.py 0.10 100` recomputes the alternative-estimand comparator.
- `python code/make_tables.py`, `code/make_supp_tables.py`, `code/make_figures.py`, `code/make_block_diagnostics.py`, `code/make_case_study.py` and `code/make_table_s9.py` rebuild the paper tables and figures into `paper/`.
