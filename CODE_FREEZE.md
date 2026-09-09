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
| `code/make_case_study.py` | 5,834 | `cb5da0df38b117444ee45fd9d3da13ed38dc715cdc231eca354218de6ab48dc3` |
| `code/make_figures.py` | 11,506 | `a7d8b3663779d2931e6599f270a17fed5218b4ae681e3880823407142bca1a25` |
| `code/make_supp_tables.py` | 2,537 | `c19f7ba7421a1a184de81b4341fedc5f9be75645729bc083eaa248e9ad276300` |
| `code/make_table_s9.py` | 2,784 | `324dd27e08fc30a8bad13564aeccaf63dfb104024b0a47b128e0775186d9d2a5` |
| `code/make_tables.py` | 14,087 | `af09d4fd224b17c251019c9545066257eeb8e8753b1a575ef30b3f95312c510c` |
| `code/requirements.txt` | 255 | `48e619ccd9c37e53ca901365a416c71d87b7be3a7375dc6a5f6d7ea8cf4e86cd` |
| `code/test_measurability.py` | 2,567 | `753f9a0f8c217bed4173e41e8671ec0659c62e8aa95db024425e0beb1642318b` |
| `code/wordcount.py` | 2,115 | `b1a7a2bfb99177751e2a94c8b3ea9af7af465e75d0277c9d3b2a079b7d81df92` |

## Results

| file | bytes | sha256 |
|---|---:|---|
| `results/all_methods_harmonised_100splits.csv` | 1,962,377 | `a26624e24a7a17b1ac99a49880f3c07e00bfbe386e3403af4d255f821866f938` |
| `results/block_diagnostics_alpha010.csv` | 1,027 | `7f37c2a9268884b4b2dc616c10adc130be3e408640da3b21206087d2c36d57a1` |
| `results/case_study_split0.csv` | 2,948 | `228270b69917c9fbb3e731a6419362f0c8641b5ce14ff733a02999a3face4bc6` |
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
