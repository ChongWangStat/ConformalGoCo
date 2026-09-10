# Deriving GoCo for neighbourhood predictors, in parallel with the Wainberg construction

The claim of this note is that the method used on FunMap and STRING is not a new idea. It is the
**same decomposition** that produced GoCo on Wainberg, evaluated in a regime where one of its two
terms is already exhausted. Writing it this way makes the two applications one method with one
proof, and it predicts — rather than rationalises after the fact — which arm should win on which
structure.

Equation numbers in parentheses are the main GoCo paper's.

---

## 1. The predictor interface is identical

Both applications ask the frozen predictor for exactly two objects, and nothing else:

| | Wainberg | FunMap / STRING |
|---|---|---|
| score `s_ig` | module enrichment statistic | −log₁₀(BH-adjusted hypergeometric *p*) |
| source set `R_ig` | annotated members of gene *i*'s winning co-essential module | annotated members of gene *i*'s top-50 RWR neighbourhood carrying term *g* |
| target exclusion | focal gene excluded from module and background | focal gene excluded from neighbourhood and background |
| retrained? | no | no |

Because the interface matches, everything downstream is shared verbatim: the gene-level TruePath
FDP with abstention contributing zero (1), the pool risk `R(π) = N⁻¹ Σ_{i∈P} L_i(π)` (2), the
finite-population normal *p*-value (9), fixed-sequence Learn-Then-Test, Definition 1 (nested and
ranking-fold measurable), Proposition 1 and Corollary 1.

**Nothing in the certificate changes.** That is the point of Corollary 1: validity holds for *any*
ranking-fold-measurable ordering, so a structure-specific path costs nothing in validity and is
spent entirely on power.

One substantive refinement, forced by the structure. Wainberg's module members are unordered, so
(4) spreads each call's unit of evidence uniformly, `a_ir = c_i⁻¹ Σ_{g∈D_i} 1{r∈R_ig}/|R_ig|`. A
random-walk neighbourhood is **ranked by proximity**, so the same attribution is proximity-weighted
with `w(rank) = 1/log₂(1+rank)`. Shrinkage (5) and the statistic (6) are then applied unchanged.

---

## 2. The decomposition that does the work

Fix a block `B` — the set of genes whose released set changes as the threshold falls from `λ_j` to
`λ_{j+1}` — with `M = |B|`, increments `Δ_i = L_i(S^L_i) − L_i(S^H_i)`, and a ranking-fold order
`rk`. Lemma 1 says the risk of admitting the first `k` genes splits into two independent effects:

**(a) Granularity.**
`R(π_{j(k)}) = R(S^H) + N⁻¹ Σ_{rk(i)≤k} Δ_i`, and under a *uniformly random* order
`E[R(π^rand_k)] = R(S^H) + (k/M){R(S^L) − R(S^H)}`.
Partial admission interpolates **linearly** between the two endpoints. This is worth something
only because a global threshold must take the block whole.

**(b) Ordering.**
`E[R(π^rand_k)] − R(π_{j(k)}) = (M/N)·Cov_B(1{rk(i)≤k}, Δ_i)`.
An order that puts small-Δ genes first makes this covariance negative, buying risk back at fixed
`k`. This term has nothing to do with ties.

The whole of GoCo is: **make (a) available by inserting policies inside the block, then make (b)
negative with a ranking-fold statistic.** Everything else is bookkeeping.

---

## 3. Why the same decomposition yields two different paths

The two structures differ in exactly one respect, and it determines which term is available.

**Wainberg — term (a) is large and irrecoverable.** 2,470 gene–GO pairs share the score 789.650;
the largest tie in the operating range holds 3,361. Inside an *exact* tie there is no interior
threshold, so no amount of grid refinement reaches part (a): a global path, however dense, still
takes or leaves the block whole. Partial admission is therefore *necessary*, and the block is the
natural unit for the ordering because the block is where the trapped budget sits. This is GoCo as
published — and it is why its ordering is deliberately block-restricted.

**FunMap / STRING — term (a) is small and free.** BH-adjusted hypergeometric scores carry no
material ties. Every block has an interior, so **densifying the global grid recovers part (a)
entirely, at no cost in validity**, because the grid is a function of frozen scores alone. Part (a)
is thus exhausted by a baseline, not by our method.

The consequence is a derivation, not a design choice:

> If part (a) is exhausted by a free baseline, then *every* remaining gain must come from part (b).
> Part (b) is maximised by letting the ordering act at the finest resolution the loss admits.

For Wainberg that resolution is the gene within a tied block, because within-block genes are score-
indistinguishable and there is nothing finer. For a tie-free score the finest resolution is the
**individual call**. So the ordering should define the path outright: candidate policies are the
top-*m* calls under the ranking-fold statistic, nested in *m*. That is `goco-rank`, and it is
Lemma 1(b) taken to its limit once Lemma 1(a) is spent.

| | Wainberg | FunMap / STRING |
|---|---|---|
| ties in operating range | up to 3,361 pairs at one score | none material |
| part (a) available to a global path | **no** (exact tie has no interior) | **yes** (densification is free) |
| therefore the method must supply | (a) **and** (b) | (b) only |
| finest resolution for the ordering | gene, within a tied block | **call** |
| resulting path | grid + partial admission (7) | top-*m* calls by (6) |
| certificate | (9) + fixed-sequence LTT | **identical** |

---

## 4. The derivation is falsifiable, and it was tested

Read as a prediction rather than a description, §3 says three things. All three were measured on
FunMap over 50 splits, and one of them nearly failed.

1. *Densification should capture most of part (a) and little else.* At α = 0.10 the dense global
   grid is worth **+1.4%** over the coarse grid — part (a) is nearly exhausted there.
2. *A block-restricted ordering should therefore underperform on this structure.* Ported GoCo gains
   only **+78 / +93 / +10** supported calls over the oracle frontier of the global path at α =
   0.05 / 0.10 / 0.20 — at α = 0.20 that is ~3% of its apparent margin, and positive on only 76%
   of splits. GoCo-as-published barely transfers, exactly as predicted.
3. *A call-resolution ordering should capture what is left.* `goco-rank` gains **+166 / +400 /
   +533** on the same frontier (94–100% of splits), at **lower** realised FDP than the block-
   restricted arm, and it returns more supported calls from *fewer* released calls at α = 0.05.

The frontier comparison matters because it is an upper bound on what any global-threshold
certificate can achieve at the realised risk, so it isolates part (b) from "spent more budget".

**Where the derivation stops, and a claim it cost me.**

An earlier version of this note proposed a scope condition: that the call-resolution path dominates
when part (a) is exhausted *and* the base unsupported rate is low enough for spreading to be
affordable. **That condition is refuted.** Measured base unsupported rates in the candidate region
are FunMap 0.855, STRING 0.874, Wainberg 0.887 — effectively identical — yet the call-resolution
path wins decisively on the first two and loses on the third at every target (Wainberg: −40.6,
−140.1, −148.1 supported calls at alpha 0.05 / 0.10 / 0.20, ahead on 12% / 0% / 0% of 50 splits,
and losing even to the plain global baseline at the two larger targets).

A second candidate, candidate-region density (calls per unit: Wainberg 1.0, FunMap 26.4, STRING
61.3), separates those three but is broken by ProteomeHD truepath, which has 120.6 calls per unit
and loses.

So the position is: **part (a)/(b) of Lemma 1 is sound and predicts which of the two paths is even
available, but what determines whether the call-resolution path WINS once available is not yet
established.** The base rate does track the outcome *within* ProteomeHD, where the four loss modes
vary the loss alone on identical data (0.341 win, 0.522 win, 0.708 lose, 0.876 lose) — so it is
plausibly a within-structure predictor whose cross-structure counterpart is confounded by something
else. Settling this needs a designed comparison that varies one factor at a time, not further
storytelling over four datasets that differ in many ways at once.

What IS established, and is enough to act on:

| structure | ties | part (a) reachable by a dense grid | best path | margin |
|---|---|---|---|---|
| Wainberg | exact, up to 3,361 pairs | no | **GoCo-L (block)** | call path loses 140 calls at alpha 0.10 |
| FunMap | none material | yes | **call path** | +400 over the oracle global frontier |
| STRING | none material | yes | **call path** | replicates FunMap |
| ProteomeHD | none material | yes | mode-dependent | wins near1/near2, loses truepath/direct |

The practical rule this supports is narrow and honest: **use the block path where the score has
exact ties, and the call path where it does not — but verify on the ranking fold before committing,
because tie-freeness is necessary and not sufficient.**
