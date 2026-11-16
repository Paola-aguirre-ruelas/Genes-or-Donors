# Methods and limitations

## Question

Does adding biological-individual identity as an auxiliary supervision target change which
64 genes a sparse multi-task selector picks, and how well downstream classifiers built on those
genes label five immune-cell classes in individuals never seen during training?

## Cohort

Processed single-cell object derived from Stephenson et al. (2021) (62,509 cells, 16,299 genes;
see `DATA.md`). Sample labels were mapped to published biological individuals. Cells from
LPS-stimulated samples and one QC-excluded collection were removed. Five annotated classes were
analysed: NK_16hi, CD4.Naive, CD4.CM, CD8.Naive, CD8.TE. An individual was eligible with at
least 20 cells in at least 3 of these classes (counted over all their samples and days):
**88 individuals, 24,456 cells** (sites: 44 / 34 / 10).

## Partition and isolation

One donor-disjoint four-fold partition, stratified by site:
`StratifiedKFold(4, shuffle=True, random_state=20261112)` over the sorted individuals
(22 held-out individuals per fold). Within each fold's training individuals, about 20% of cells
per individual × class were withheld (Stage B; `m = floor(0.2 n + 0.5)`, generator
`[20261112, k]`) and never used. From the remaining *fitting* cells, a selector subsample (≤ 50
cells per individual × class, generator `[20261112, k, 1]`) and a classifier subsample (≤ 200,
`[20261112, k, 2]`) were drawn. Every expression read is restricted to an explicit list of
allowed cells.

## Features

2,000 highly variable genes (Scanpy, `flavor="seurat"`) and per-gene means and standard
deviations (ddof 0) were computed on each fold's fitting cells only; held-out cells were
transformed with those training parameters.

## Panel selection

`MultiTaskLasso` (scikit-learn; intercept, tol 1e-4, max_iter 1000 with one retry at 10,000,
warm starts, cyclic updates) on standardized selector cells. Targets are one-hot blocks,
column-centred and scaled to unit total variance: **w = 0** uses the 5-column class block;
**w = 1** appends the 66-column training-individual block. A per-fold alpha search seeks exactly
64 active genes; if none is found by fit 25 (or 50, when the fallback would discard more than 8
genes), the 64 largest coefficient-row norms of the largest alpha with more than 64 genes are kept.
The alpha sequence differs between fold 0 and folds 1–3 (see `HISTORICAL_FIDELITY.md`, item 1);
the target (64 genes or the same fallback) is identical.

## Classification and evaluation

Multinomial L2 logistic regression (C = 1, lbfgs, max_iter 1000, no class weights) per panel,
trained on the classifier subsample. For each held-out individual, macro-F1 averages per-class F1
over classes with at least 20 true cells (at least 3 required); every test cell of the individual
counts towards true/false positives and negatives, so a cell of a non-evaluable class predicted as
an evaluable class is a false positive. The primary descriptive quantity is the donor-equal mean of
the paired difference Δ = macro-F1(w = 0) − macro-F1(w = 1) over 88 individuals; a negative Δ
means w = 1 scored higher.

## Exploratory controls

* **Random reference:** per fold, 20 random 64-gene panels from the 2,000 candidates, matched to the
  w = 0 panel's counts in 20 bins of training-only mean expression (generator
  `[20261115, k, draw, 444]`); same classifier and evaluation.
* **Fixed-alpha tolerance check (G3):** one cold-start fit per panel at its selected alpha with
  tol 1e-5; stable when the same 64 genes (or the same top-64 for the fallback panel) are recovered.
* **Identity-permutation sanity check (G4):** fold 0, donor labels of the selector rows permuted
  (seed 20261115), cold-start alpha search at tol 1e-5 with the w = 1 target; the resulting panel is
  compared with the w = 0 panel.
* **Fallback diagnosis (H1):** fold 3, w = 1: cold fits at four alphas around the fallback alpha.

## Limitations

1. **One cross-validation repetition** with the pilot seed: four folds are not four independent
   replications, and no confidence intervals or hypothesis tests are reported.
2. **Analysis chronology:** fold 0 was evaluated before the runner for folds 1–3 was written; the
   alpha-search code differs between them. No outcome-driven change was identified, but the
   chronology permits it.
3. **Alpha search** is path-dependent (warm starts, zero-iteration acceptance); one panel uses the
   fallback. The tighter-tolerance checks validate the saved panels at fixed alpha, not a full
   tight-tolerance search.
4. **Labels** are the published annotations from joint clustering that included the held-out cells;
   scores measure agreement with those labels, not with an independent ground truth.
5. **Population:** a selected, processed subset of unknown subsampling provenance; results describe
   these 88 individuals under these technical conditions.
6. **Confounding:** donor identity is entangled with library, site, disease state, sex and cell
   composition; donor-disjoint folds are not library-disjoint.
7. **Controls:** the random reference was introduced after the main held-out results were known;
   stronger selection baselines, a signal-bearing decoy task, displacement controls and fresh-seed
   repetitions, as a confirmatory design would require, were not run, and the exploratory controls are
   not substitutes for them.
8. **Ceiling:** macro-F1 near 0.9 leaves limited room for differences; this does not establish a
   measured performance ceiling.

These limitations rule out confirmatory claims of harm, benefit, equivalence or non-inferiority,
causal attribution to identity, and generalization beyond the evaluated data.
