# Patient-Identity Supervision in Compact Gene Panel Selection: An Exploratory Analysis of Generalization to Unseen Donors in scRNA-seq

**Paola Aguirre** · Research report v2.0 (regenerated from the reproduced computation) · 2026-11-16

*Every number in this report is inserted programmatically from `results/results.json`
(SHA-256 `a0d0dbc5e44626e91d081631fc78a0df00d5a34d53227ece0f8a798a2bff7eef`) by `cgp rebuild-results`; none is typed by hand. The historical computation (report v1.3) with which this reproduction
was compared is kept unchanged; its artifacts are not distributed and are identified by SHA-256 in
`configs/historical_references.toml`.*

## Abstract

**Question.** Does adding biological-individual identity as an auxiliary supervision target alter
compact gene selection and the classification of immune-cell subtypes in previously unseen
individuals? **Methods.** An eligible subset of the Stephenson et al. single-cell peripheral-blood
dataset was analysed after mapping sample labels to biological individuals. Eight 64-gene panels were
selected by multi-task sparse regression across one donor-disjoint four-fold partition
(88 individuals; 24,456 eligible cells): with cell-class labels alone (w=0) or with
class labels plus identity (w=1). L2-regularized classifiers predicted five lymphocyte/NK-cell classes.
The descriptive primary quantity was the donor-equal mean of paired macro-F1 differences,
F1(w=0) − F1(w=1). **Results.** The mean paired Δ was **−0.0037** (39 positive and 49
negative individual-level differences); fold means ranged from −0.0079 to +0.0031.
A negative Δ indicates slightly higher observed macro-F1 for w=1; it does not indicate benefit,
equivalence or non-inferiority. Panel Jaccard similarity ranged from 0.51 to 0.60.
Cell-class-only panels exceeded the mean of 20 expression-matched random 64-gene panels per
fold by **+0.124** macro-F1. Fixed-alpha, tighter-tolerance (tol 10⁻⁵) cold refits recovered the same 64 selected genes for 8 of 8 panels, including the truncated top-64 set of the fallback panel. **Reproduction.** The complete pipeline —
cohort, folds, preprocessing, all eight alpha searches, eight classifiers, 80 random-reference
models, held-out scoring and the exploratory numerical checks — was re-executed from the hash-verified input;
all 54 compared components agreed bit-identically with the historical records. **Interpretation.** Auxiliary identity supervision substantially changed the selected
genes without an observed average loss on this partition and task. This is not a confirmatory test, and
identity may proxy library, site, disease state, sex or cell composition.

## 1. Question

Compact transcriptomic panels can distinguish cell states while remaining interpretable and
inexpensive. Adding donor identity as a second supervised target might favour genes carrying
individual-specific information, displacing cell-type markers and impairing transfer to unseen
individuals; it might also be neutral or beneficial. The exploratory question addressed here is: *how
different are 64-gene panels selected with and without identity supervision, and how do their fixed
downstream classifiers perform on donor-held-out cell-type classification in one four-fold partition?*

## 2. Data and methods

A processed, subsampled object derived from Stephenson et al. (2021) [1] (62,509 cells, 16,299 genes)
was analysed. Sample labels were mapped to published individuals; LPS-stimulated cells and one
QC-excluded collection were removed; individuals with at least 20 cells in at least three of five
classes (NK_16hi, CD4.Naive, CD4.CM, CD8.Naive, CD8.TE) were eligible: **88 individuals and
24,456 cells**. Four site-stratified donor-disjoint folds held out 22 individuals each. Within the
training individuals of each fold, about one fifth of cells were withheld; 2,000 highly variable genes (Scanpy,
seurat flavour [2, 3]) and standardization were estimated on the remaining fitting cells only.
`MultiTaskLasso` [4, 5] selected 64 genes per arm through a per-fold alpha search (fallback: 64 largest
coefficient-row norms); fold 0 followed a different alpha sequence from folds 1–3 towards the same
target. Multinomial L2 logistic regression was trained on each panel. Macro-F1 was computed per
held-out individual over classes with at least 20 true cells, counting every test cell for false
positives. Exploratory controls: a random reference (20 expression-bin-matched 64-gene panels
per fold), fixed-alpha tighter-tolerance refits (G3), an identity-permutation sanity check (G4) and a
fallback diagnosis (H1). Full details: `docs/METHODS.md`.

## 3. Results

### 3.1 Donor-held-out performance

**Table 1. Donor-equal mean macro-F1 by held-out fold.**

| Fold | Individuals | Macro-F1 w=0 | Macro-F1 w=1 | Random-reference mean | Paired Δ (w=0 − w=1) |
|---|---:|---:|---:|---:|---:|
| 0 | 22 | 0.9061 | 0.9140 | 0.8100 | −0.0079 |
| 1 | 22 | 0.9075 | 0.9044 | 0.7936 | +0.0031 |
| 2 | 22 | 0.8997 | 0.9058 | 0.7503 | −0.0061 |
| 3 | 22 | 0.9165 | 0.9201 | 0.7815 | −0.0037 |
| **All folds** | **88** | **0.9074** | **0.9111** | **0.7838** | **−0.0037** |

*Rounded to four decimals from unrounded individual-level values; Δ is computed before rounding.
Unrounded global values: w=0 0.907442, w=1 0.911095, random 0.783828, Δ −0.0036525.*

Of 88 paired differences, 39 were positive, 49 negative and 0 zero; the
median was −0.0019. A negative Δ indicates slightly higher observed macro-F1 for w=1 on this task
and partition; it does not demonstrate benefit, equivalence, non-inferiority or absence of harm.

![Figure 1. Paired macro-F1 differences for the held-out individuals.](../figures/figure_01_paired_donor_deltas.png)

### 3.2 Random reference and panel composition

Cell-class-only panels exceeded the random reference by +0.124 macro-F1 (w=1:
+0.127); w=0 exceeded every random draw on the fold means: yes, in all folds. At the
individual level, w=0 exceeded every corresponding random draw for 73 of 88
individuals. These comparisons show that the selected panels outperform this particular reference, not
that they are optimal.

![Figure 2. w=0, w=1 and expression-matched random panels by fold.](../figures/figure_02_random_baseline.png)

Within each fold, the w=0 and w=1 panels shared 43–48 of 64 genes (Jaccard
0.51–0.60, mean 0.55). The fold 3 w=1 panel used the fallback: 67 active genes reduced to the 64 largest coefficient-row norms.

![Figure 3. Within-fold overlap of the w=0 and w=1 panels.](../figures/figure_03_panel_overlap.png)

### 3.3 Exploratory numerical checks

*Exploratory numerical sensitivity checks (G3 fixed-alpha tolerance refits, G4 identity-permutation sanity check, H1 fallback diagnosis). They are not confirmatory gates and do not replace pre-specified controls such as a decoy task, a displacement control or fresh-seed repetitions.*

Fixed-alpha, tighter-tolerance (tol 10⁻⁵) cold refits recovered the same 64 selected genes for 8 of 8 panels, including the truncated top-64 set of the fallback panel. Per panel — fold 0 w=0: 64 active genes, support Jaccard 1.0000; fold 0 w=1: 64 active genes, support Jaccard 1.0000; fold 1 w=0: 64 active genes, support Jaccard 1.0000; fold 1 w=1: 64 active genes, support Jaccard 1.0000; fold 2 w=0: 64 active genes, support Jaccard 1.0000; fold 2 w=1: 64 active genes, support Jaccard 1.0000; fold 3 w=0: 64 active genes, support Jaccard 1.0000; fold 3 w=1: 67 active genes, support Jaccard 0.9552, top-64 Jaccard 1.0000. The identity-permutation sanity check (fold 0) selected a panel with Jaccard 63/65 = 0.9692 relative to the w=0 panel (PASS at the 0.90 threshold; 6 cold fits). At 4 fixed alphas around the fallback alpha, cold fits had 67, 67, 67, 63 active genes; no exact-64 solution occurred at these alphas.

## 4. Reproducibility

**Level A (verified-input reproduction).** Starting from the processed H5AD (SHA-256 pinned) and two
auxiliary inputs that are re-derivable from public sources, the pipeline re-executed every stage. Each
component was compared with the historical records under tolerances fixed before the comparison
(`configs/acceptance.toml`): fold plans: 4 EXACT; preprocessing: 4 EXACT; alpha searches and panels: 8 EXACT; classifiers: 8 EXACT; random-reference panels and models: 4 EXACT; held-out scores (main arms): 4 EXACT; held-out scores (random reference): 4 EXACT; stability checks: 10 EXACT; results package (tables, summary values, figures): 8 EXACT.

**Level B (source-data reconstruction)** is not claimed: the subsampling and preprocessing that
produced the processed H5AD are not published as code.

## 5. Interpretation and limitations

The principal descriptive observation is a large change in selected gene identity with little change
in unseen-donor classification performance. One explanation is that this five-class task admits several
predictive feature sets; another is that library- and disease-related signatures enter the
identity-supervised objective. Neither is established by these data. Limitations: one cross-validation
repetition; fold 0 was evaluated before the runner for folds 1–3 existed; a path-dependent alpha search
with one fallback panel; dependence on published annotations derived with the held-out cells; unknown
subsampling provenance; confounding of identity with library, site, disease state, sex and composition;
retrospective random controls; no stronger baselines, decoy task, displacement control or fresh-seed
repetitions. These rule out confirmatory claims of harm, benefit or equivalence and causal attribution
to identity.

## 6. Conclusion

In one exploratory donor-disjoint four-fold evaluation of 88 individuals, auxiliary
identity supervision changed 16–21 of the 64 genes in each panel, while the mean paired
individual-level macro-F1 difference was −0.0037. Both arms clearly exceeded an expression-matched
random reference. The full computation is reproducible from the verified input with the code in this
repository.

## AI declaration

AI tools assisted with the project and manuscript preparation. Claude conducted separate computational audits, and ChatGPT assisted with manuscript review. The author retains responsibility for the scientific content and conclusions.

## References

1. Stephenson E, Reynolds G, Botting RA, et al. Single-cell multi-omics analysis of the immune response in COVID-19. *Nature Medicine* 27, 904–916 (2021). doi:10.1038/s41591-021-01329-2
2. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 19, 15 (2018). doi:10.1186/s13059-017-1382-0
3. Satija R, Farrell JA, Gennert D, Schier AF, Regev A. Spatial reconstruction of single-cell gene expression data. *Nature Biotechnology* 33, 495–502 (2015). doi:10.1038/nbt.3192
4. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: Machine Learning in Python. *Journal of Machine Learning Research* 12, 2825–2830 (2011).
5. Obozinski G, Taskar B, Jordan MI. Joint covariate selection and joint subspace selection for multiple classification problems. *Statistics and Computing* 20, 231–252 (2010). doi:10.1007/s11222-008-9111-x
