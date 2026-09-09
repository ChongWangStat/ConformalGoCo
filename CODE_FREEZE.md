# Code and results freeze

Date: 2026-09-09 (updated). The algorithms and the reference implementation are frozen; the only run added after the first freeze was the promoted ordering at alpha = 0.05 and 0.20, using the same code.
Manuscript work is limited to writing, structure and journal compliance. Any change to a file below invalidates its hash.

## Code (analysis_rerun/)

| file | bytes | sha256 |
|---|---:|---|
| `goco_rerun.py` | 19,430 | `b661ea89a229528291fd3b168db4fb2fa121162fc8e6017afdb9a16a3ddd3143` |
| `goco_learned.py` | 9,107 | `f213867eecfe61b95706320429d6cbd71cda9abc20ab1d757fed7b3427b72a79` |
| `test_measurability.py` | 2,779 | `dc8eafbda47d2ae3607e697b4164db83af27d4b2842b8f6612bf7f3a589cce4a` |
| `make_tables.py` | 17,739 | `0a8f4d044bbe1a780725132336a1cbea9349ca2801689354518f5ca3a9e1b067` |
| `make_figures.py` | 13,070 | `e614808aad2d9b37441bb57c9b1a69d235b95cfe6c36964af38ff56e1cf4d0a5` |
| `make_block_diagnostics.py` | 9,604 | `063f36f9d35bf8bc48daa25b2c787baa24a16040a3f410c15be397c506c820ac` |
| `make_case_study.py` | 5,592 | `78af40bd7ccff735fc97710ce5117c588061c555707f1f654ce88b8b6a4a5f12` |
| `make_supp_tables.py` | 5,672 | `3b22f4299f63c4e32f9822177cb6426b224276fb06ef94593f3ef81e85a5236b` |
| `wordcount.py` | 2,115 | `b1a7a2bfb99177751e2a94c8b3ea9af7af465e75d0277c9d3b2a079b7d81df92` |
| `requirements.txt` | 255 | `48e619ccd9c37e53ca901365a416c71d87b7be3a7375dc6a5f6d7ea8cf4e86cd` |

## Split-level results and derived tables (analysis_rerun/)

| file | bytes | sha256 |
|---|---:|---|
| `DRIVE_Ldense_100splits.csv` | 44,996 | `52c471e88d01fe76f0759c8c466f41001d2a62ed8c6435f723bd27514d650f6b` |
| `DRIVE_Lsweep_100splits.csv` | 171,076 | `df38f43b738034005ecdbe1d2ea8ef280af6400805ea1db5d33f9fda25436428` |
| `DRIVE_fix1_100splits.csv` | 170,788 | `7742e7e1b230d3788026e7e491a69cb0fefb57e7648aaeeebaba182307d5c08a` |
| `DRIVE_full_100splits.csv` | 1,770,985 | `39770f95797e4bd7d7534b8fb626112fbb5df6d4ab7e6ccb608484431f9c8bf1` |
| `DRIVE_knap2_100splits.csv` | 168,443 | `f7a6839ddb4a4e0bfefc3266e119302c87e787a56f73098967ff25a748d63c94` |
| `DRIVE_knap_100splits.csv` | 252,732 | `0f2c99c1a5b47f1cf43cbe01b2aefbe18e908e861c4e5d16a04ab685d7cfffc9` |
| `DRIVE_learn2_100splits.csv` | 86,264 | `3e93451ba06e066241d909e34a1f6bea3a7e0f484e211ace040ac04a0dc1d500` |
| `DRIVE_learn_100splits.csv` | 212,890 | `a28d0edb64ae1f5fa104c5421fc0df766122b719ff4e32070f0138de95e76b37` |
| `HAP1_Ldense_100splits.csv` | 46,504 | `94ee455ad67fb22fb38b6bf7ee688edfd82dff6398c11aebed7bc78533eac096` |
| `HAP1_Lsweep_100splits.csv` | 178,055 | `48c7f7dfa0686fca11f053bd06accda3e6de887bfeb8c66ce653e6bf797351a4` |
| `HAP1_fix1_100splits.csv` | 170,235 | `65b016f5c0fefb96b9690934b634064d877b2e99432a821fcb2b0c7f8114f99c` |
| `HAP1_full_100splits.csv` | 1,801,450 | `402214af9ef54e1a6b31977a8a5557c2622f413bd3990d0a8f2c659c8f05c3b8` |
| `HAP1_knap2_100splits.csv` | 175,759 | `f05849b759344c8639f24ed4c1cce44479600cd6b5063e2a3859ab5fcec064fd` |
| `HAP1_knap_100splits.csv` | 253,744 | `5e5ea00f869a83acc2def95a6c22d2c4220a4b83fbb0edff8300e035eee49121` |
| `HAP1_learn2_100splits.csv` | 86,186 | `7201b100c04ef2dcc10c6eeb49ebcb67e495c09c7c139b615df72a986944a886` |
| `HAP1_learn_100splits.csv` | 212,356 | `8070c4b6a727959d7aef6dfd62800cf1377a6b6c40123d9e68d099735ad0d603` |
| `Sanger_Ldense_100splits.csv` | 46,024 | `d5628226a008dc552fb82e177846f4ef5fa180f616168698a6f1d82795496ee4` |
| `Sanger_Lsweep_100splits.csv` | 179,290 | `45437a524584b2a02851f51e76363de5ae706006fc464fda87f254a64ce0b83e` |
| `Sanger_fix1_100splits.csv` | 172,173 | `493d08dfe176a4802b9db18378fad7531c00983c2ffe40098d77f700ca4c37cb` |
| `Sanger_full_100splits.csv` | 1,820,559 | `099a129c551a29df04153775e660c2b6b0b09e6b36a8eee9fd6f7bddfeb7b7b5` |
| `Sanger_knap2_100splits.csv` | 176,487 | `c6fbc4cfba5100429c8c8774e9d5959dbf728f31f9c101324066eebd8365fb4a` |
| `Sanger_knap_100splits.csv` | 255,557 | `a97760664fd816147b5e76cf05b5a2b6cf5864f2e50c28f591e4c0256558c4a7` |
| `Sanger_learn2_100splits.csv` | 86,961 | `918ec1545baa0d4fd441ecf0524407a2ed715b192dc179bae3cea43da76caff9` |
| `Sanger_learn_100splits.csv` | 214,547 | `85755f77479aefaee2c4827105f01bf8c85b521bbb7ed244facd891e14c76a5a` |
| `Wainberg_Ldense_100splits.csv` | 47,094 | `3827f3799caf803df82f867e4057e15b12d15da3207dae9b1cf19e3e7df1219f` |
| `Wainberg_Lsweep_100splits.csv` | 189,637 | `f969c0dc408f2c5a97e5550f7f3265acd5184e74ac0025cc0ff90de9bf3e478f` |
| `Wainberg_fix1_100splits.csv` | 187,590 | `686477137228276b34c81a70b7232d993c7ac7ffe0b29664d178667f2c43d04e` |
| `Wainberg_full_100splits.csv` | 1,820,511 | `d4bd4487f5baed85a0cd3a31b084f731bc62051fb752c4eddc8802558a62fbc6` |
| `Wainberg_grid_100splits.csv` | 1,266,779 | `a2e5a7e81322ea3f32e32b1bad458a1a711afab88b8fc1e71f4cfc6ad3c65084` |
| `Wainberg_knap2_100splits.csv` | 175,888 | `850b2820bb75c30b32b889a42462a0bba5eef810e501f02883f798c913345185` |
| `Wainberg_knap_100splits.csv` | 261,391 | `55a075ca652e16acd04b9c57e6183cef942cf10654897c0c1e0ca8856abeb249` |
| `Wainberg_learn2_100splits.csv` | 94,472 | `4053df97ccd1ec6cc547a0396b1116e414546bca61bc2083bea6f12ad7af6122` |
| `Wainberg_learn_100splits.csv` | 233,788 | `0cadc4d40d4c6b422c232b104bcfcd3b6b2c15f8a32e5f76c15cf6b4a1728666` |
| `all_datasets_full_100splits.csv` | 6,605,083 | `6b57cada8abc954b57481271ef990c1601c37fc15d2a11896f03e12a91385759` |
| `all_methods_harmonised_100splits.csv` | 14,339,512 | `c19ff79e517278963d8fa27d9b49c894673c660d84eb73795aa797591070d2c4` |
| `block_diagnostics_alpha010.csv` | 1,347 | `911188e6957ab399e55c02369b67fa6b81f79001a5fcd80dae0cf0cacd2bf11f` |
| `case_study_split0.csv` | 2,948 | `228270b69917c9fbb3e731a6419362f0c8641b5ce14ff733a02999a3face4bc6` |
| `first_draft_primary_method_summary_long.csv` | 8,768 | `a3fee7c60ef2ea6e35cbfc9f9ec2d35e0cac4cc54c05b872f5814c179abb6394` |
| `paired_contrasts_alpha010.csv` | 7,227 | `8ea7f97df1cefe85eced4681ac6e4d321d754e940a7612a3ebe5cfedd9b54312` |
| `pool_level_exceedance_alpha010.csv` | 11,967 | `315e485ad4c1fd95ec72f2ce7b9e967a6a34afe3f82141cddef7bc090e4dc527` |
| `ranking_informativeness_diagnostic.csv` | 1,135 | `7821cc0db4d8b406f4a8779528bcb33a257839c38bb705c9c2962404d278e3b5` |
| `summary_by_method.csv` | 21,637 | `65c9175ec3e2d73ec30da75a7e6f3b27277e938f47aa7a3f751dbcfc3801b41a` |
| `summary_by_method_harmonised.csv` | 49,903 | `674b15e3c1d0b012ac0b071262f7065c566d2ae36435e4334c49bd4934669fdd` |
| `tie_counts.csv` | 271 | `424670da41826cb8dada0621172a212ab61e935dde7d4f379b9a8b1959156c89` |

## Entry points

- `python goco_rerun.py --tag full` reproduces every method at alpha in {0.05,0.10,0.20} and delta in {0.10,0.50} over the 100 frozen splits.
- GoCo (the per-call learned ordering) is `--methods GoCoGrid:learned` on Wainberg and `GoCo:learned` on the external sets; `GoCoDense:learned` is the dense-path variant. Tags `learn2` (alpha=0.10) and `Lsweep`/`Ldense` (alpha=0.05, 0.20) produced the reported runs.
- `python test_measurability.py <dataset> <seed>` runs the ranking-fold-measurability test for every ordering (scramble all non-ranking-fold labels; the pool ordering must not move).
- `python make_tables.py`, `make_figures.py`, `make_supp_tables.py`, `make_block_diagnostics.py`, `make_case_study.py` regenerate all tables and figures from the split-level CSVs.
- `python wordcount.py` reports the body word count per section against the journal limit.
