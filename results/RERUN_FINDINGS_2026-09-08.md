# GoCo re-analysis findings (2026-09-08)

Script: `goco_rerun.py` (self-contained; inputs = raw Wainberg score/truth/attribution files and the frozen external score files). Outputs: `*_full_100splits.csv`, `all_datasets_full_100splits.csv`, `summary_by_method.csv`, `paired_contrasts_alpha010.csv`, `ranking_informativeness_diagnostic.csv`.

## 1. Reproduction gate (alpha = 0.10)
- Recomputed loss/yield/call matrices for Sanger, DRIVE and HAP1 match the frozen `*_extended100.npz` exactly (max abs diff 0).
- Wainberg attribution file (16,480 rows) is identical to the raw score file restricted to score >= 600 and the 15,794 units.
- Boger, GoDag and GoCo rows reproduce the BIB package (`primary_method_summary_long.csv`) to 2 decimals in 22/24 cells; two GoCo cells differ by 0.05-0.3 terms (Sanger delta=.5: 365.88 vs 365.58; DRIVE delta=.1: 255.62 vs 255.67) = one borderline split flipping under float32 storage. Final manuscript numbers must be regenerated from ONE pipeline.
- **C5 resolved:** the Sept-7 CLT package's "GoDag TruePath (CLT)" values for Sanger/DRIVE (358.43 / 282.16) are exactly `GoDag-dense` (global path over every distinct score value); the BIB package's GoDag (353.79 / 258.60) uses the 25-unit grid. The fair global comparator for the externals is GoDag-dense.

## 2. Decomposition of the GoCo gain (GO-term yield, alpha = 0.10; means over 100 splits; paired P = fraction of splits GoCo-source is better)

| Dataset | delta | GoDag (grid) | GoDag-dense | GoCo-random (3 salts) | GoCo-score | GoCo-source (paper) | oracle | source vs random: diff, P(source better) |
|---|---|---|---|---|---|---|---|---|
| Wainberg | .5 | 339.61 | 339.61 | 344.7-347.7 | 341.8 | 353.88 | 381.0 | +6 to +9, 0.90-0.99 |
| Wainberg | .1 | 339.61 | 339.61 | 341.6-342.9 | 340.2 | 349.71 | 379.0 | +7 to +8, 0.94-0.98 |
| Sanger | .5 | 353.79 | 358.43 | 358.0-360.1 | 361.9 | 365.88 | 388.7 | +6 to +8, 0.69-0.80 |
| Sanger | .1 | 319.11 | 326.00 | 330.2-334.0 | 332.2 | 331.36 | 343.7 | -2.6 to +1.2, 0.30-0.56 (none) |
| DRIVE | .5 | 258.60 | 282.16 | 278.6-281.9 | 285.0 | 288.05 | 299.1 | +6 to +9, 0.68-0.79 |
| DRIVE | .1 | 237.23 | 253.06 | 245.6-250.4 | 250.3 | 255.62 | 270.3 | +5 to +10, 0.62-0.81 |
| HAP1 | .5 | 129.45 | 132.11 | 134.4-135.8 | 137.9 | 135.65 | 139.7 | -0.1 to +1.2, 0.47-0.60 (none) |
| HAP1 | .1 | 129.45 | 132.11 | 133.6-134.9 | 137.8 | 134.35 | 139.7 | -0.6 to +0.8, 0.31-0.59 (none) |

Reading: most of the gain over the grid-based GoDag comes from *finer policy granularity* (partial admission of a tie block / grid step), which a random ordering also achieves. Relative to GoDag-dense the source-ordered GoCo adds +14.3 (Wainberg, the tie cannot be split by any threshold), +7.5/+5.4 (Sanger), +5.9/+2.6 (DRIVE), +3.5/+2.2 (HAP1) at delta = .5/.1. The co-regulation ordering beats random ordering clearly in Wainberg (development data), modestly in Sanger and DRIVE at delta = .5, and not at all in HAP1 or in Sanger at delta = .1.

## 3. Why: ranking informativeness vs. yield
Out-of-sample Spearman(d_hat_i, Delta_i) on affected genes (100 splits): Wainberg 0.33, Sanger 0.17-0.24, DRIVE 0.13-0.17, HAP1 0.55 - positive in ~100% of splits. The ranking predicts the CONTROLLED loss everywhere; the yield benefit depends on how the freed budget converts into supported calls. The sqrt(c) utility favours genes with many added calls; in HAP1 (mean Delta 0.76) that buys calls, not supported terms (GoCo-source releases 3,356 calls vs 2,900 for random for the same supported yield). Utility exponent 1 (u = d_hat / c) is equal or slightly better than sqrt in Sanger/DRIVE; exponent 0 is markedly more conservative (does not use the budget). A knapsack-style utility (predicted loss increment per predicted supported added call) is being tested.

### 3b. Knapsack-style utility (run `--tag knap`)
u_i = Delta_hat_i / max{c_i (1 - f_hat_i), 1e-3}, where f_hat_i is the source-smoothed prediction (same estimator as Delta_hat, kappa = 1) of the unsupported fraction of gene i's added calls; i.e. predicted loss increment per predicted supported added call (the greedy knapsack ratio for maximising supported yield under a risk budget). GO-term yield at alpha = 0.10:

| Dataset | delta | GoCo-random1 | GoCo-source (sqrt c) | GoCo-knap | knap vs random: diff, P | knap vs sqrt: diff, P |
|---|---|---|---|---|---|---|
| Wainberg | .5 | 344.73 | 353.88 | 356.19 | +11.5, 0.99 | +2.3, 0.69 |
| Wainberg | .1 | 341.56 | 349.71 | 351.11 | +9.6, 0.93 | +1.4, 0.62 |
| Sanger | .5 | 360.07 | 365.88 | 370.39 | +10.3, 0.92 | +4.5, 0.74 |
| Sanger | .1 | 333.95 | 331.36 | 336.15 | +2.2, 0.63 | +4.8, 0.78 |
| DRIVE | .5 | 281.93 | 288.05 | 292.77 | +10.8, 0.87 | +4.7, 0.74 |
| DRIVE | .1 | 248.78 | 255.62 | 258.13 | +9.4, 0.88 | +2.5, 0.57 |
| HAP1 | .5 | 135.19 | 135.65 | 136.13 | +0.9, 0.59 | +0.5, 0.42 |
| HAP1 | .1 | 134.48 | 134.35 | 135.91 | +1.4, 0.69 | +1.6, 0.67 |

Mean held-out FDP stays <= alpha for all (0.0984/0.0992/0.0994/0.0976 at delta = .5). The knapsack utility is never worse than sqrt(c) and beats random ordering in Wainberg, Sanger and DRIVE; HAP1 remains a null case for ordering. Because it also wins on the development data (Wainberg), it can be adopted as the primary utility provided the utility-family comparison (exp0, sqrt, exp1, knap, random, score, oracle) is reported as a prespecified-family sensitivity analysis and the selection is disclosed.

### 3c. Learned per-call ordering (GoCo-L; final run `--tag learn2`, 2026-09-09)
`goco_learned.py` (audited version: no background-count feature, cross-fitted classifier and term frequency). GO-term yield at alpha = 0.10 (paired difference vs Eq. 6 in brackets; P = fraction of splits ahead):

| Dataset | delta | Eq. 6 | GoCo-L logistic | GoCo-L GB | oracle | held-out FDP L / Eq.6 | exceed L / Eq.6 |
|---|---|---|---|---|---|---|---|
| Wainberg | .5 | 356.2 | 369.0 [+12.9, 1.00] | 365.5 | 381.6 | 0.0964 / 0.0984 | 0.20 / 0.39 |
| Wainberg | .1 | 351.5 | 365.2 [+13.8, 0.98] | 362.1 | 379.0 | 0.0948 / 0.0953 | 0.13 / 0.17 |
| Sanger | .5 | 370.4 | 380.7 [+10.3, 0.70] | 376.7 | 388.5 | 0.0980 / 0.0992 | 0.42 / 0.45 |
| Sanger | .1 | 336.1 | 340.0 [+3.8, 0.73] | 338.4 | 343.5 | 0.0930 / 0.0947 | 0.15 / 0.17 |
| DRIVE | .5 | 292.8 | 291.7 [-1.1, 0.31] | 291.4 | 298.4 | 0.0971 / 0.0994 | 0.33 / 0.47 |
| DRIVE | .1 | 258.1 | 263.5 [+5.3, 0.60] | 263.4 | 269.8 | 0.0931 / 0.0945 | 0.16 / 0.24 |
| HAP1 | .5 | 136.1 | 137.7 [+1.5, 0.69] | 137.8 | 139.7 | 0.0943 / 0.0976 | 0.19 / 0.33 |
| HAP1 | .1 | 135.9 | 137.6 [+1.7, 0.75] | 137.6 | 139.7 | 0.0915 / 0.0940 | 0.06 / 0.13 |

Recovery of the oracle-over-random gap at delta=.5: 64/73/63/55 %. The superseded first version (`--tag learn`, with the background-count feature and in-sample fitting) gave 370.4/384.8/297.6/139.3; the audit showed the removed feature was not exploited (yield unchanged when dropped on 10 dev splits), the differences come from cross-fitting.

T-measurability check (definitive): for Wainberg split 0, block 800->775, permuting the T labels of every non-ranking-fold gene (and scrambling the truth dictionary of pool genes) leaves u^L identical on all affected pool genes (max |diff| = 0.0); the oracle ordering changes under the same permutation. Snippet: see the Bash block in the session log / `make_block_diagnostics.py` uses the same utility() entry point.

Tie anatomy (`tie_counts.csv`): the enrichment score is a function of (k supporting members, n module size ex target, K background). 85% of scored Wainberg calls come from modules with n <= 5. 789.65: 2,470 rows, 1,256 genes, 482 terms, 347 modules; 2,018 rows have (k,n,K) = (1,4,5). External top ties: Sanger 263.217 (1,880 rows, 847 genes, 393 modules), DRIVE 225.614 (723; 218; 88), HAP1 219.347 (14,245; 5,314; 540). HAP1 step 475->450: 2,089 of 2,115 rows share 451.229.

## 4. alpha-sweep
The paper's one-block Wainberg path is stuck at alpha = 0.05 (identical to all globals) and collapses at alpha = 0.20 (381.9 vs GoDag 518.9). A tie-level path (GoCo-dense: partial admission inside every distinct score level, external block rules, no hand-picked 789.65 block) dominates GoDag at every alpha in Wainberg (alpha=.05: 204-210 vs 183; alpha=.10: 353.8 = paper; alpha=.20: 531-533 vs 519) and matches or slightly trails the paper path in the externals. Recommend adopting the tie-level path as the GoCo definition; the 789.65 block becomes a special case.

## 5. Cost metrics (must appear in the paper)
Term-level unsupported fraction 0.79-0.94 across methods/datasets; marginal precision of GoCo's additional calls vs GoDag 5-8% (HAP1 1.4%); ~12% of held-out genes receive any call; P(held-out FDP > alpha) at delta=.5 is 0.36-0.46 for GoCo vs 0.01-0.36 for GoDag, consistent with delta = 0.5 (pool risk certified at the boundary).

## 6. Implications for the manuscript
- Reframe the contribution: (i) module-based predictors produce massive score ties (2,470 rows at 789.65; Sanger/DRIVE top ties 1,880/723), so global thresholds waste risk budget; partial admission within tie blocks recovers it (robust across all four datasets, any ordering); (ii) a source-informed ordering adds a further increment when the source signal is informative, with an explicit diagnostic (nuisance-fold association between d_hat and Delta) and an honest report that the increment is absent in HAP1.
- Report GoDag-dense and GoCo-random as primary baselines; drop the claim that GoCo vs GoDag "isolates the value of co-regulation".
- Include HAP1, delta = 0.10 and the alpha sweep with yields.
