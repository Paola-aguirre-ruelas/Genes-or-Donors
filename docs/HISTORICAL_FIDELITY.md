# Historical fidelity: from executed scripts to package modules

The package re-implements the analysis **as it was executed**, including its
inconsistencies. This document records the dependency graph of the historical scripts,
the module that now implements each step, and every behaviour that was preserved
deliberately even though a fresh design would differ.

## 1. Dependency graph of the executed analysis

```
identity mapping table (SDRF-derived) ─┐
QC exclusion list (metadata rule) ─────┤
processed H5AD (obs/var) ──────────────┴─► split script ──► cohort, 4 site-stratified donor-disjoint folds
                                                  │
                                                  ├─► Stage B withholding script (fold 0)
                                                  ▼
                                   frozen pilot script (fold 0): HVG, standardization,
                                   selector/classifier subsamples, Selector class, target blocks,
                                   guarded row loader ── frozen parameters (fold 0)
                                                  │ (imported, not executed)
              ┌───────────────────────────────────┼───────────────────────────────┐
              ▼                                   ▼                               ▼
   fold-0 panel recovery script        CV runner (folds 1–3): plans, Stage B,    (both reuse the
   (alpha search, variant A)           HVG, panels (variant B), classifiers,      pilot definitions)
              │                        scoring, 4-fold aggregate
              ▼                                   │
   fold-0 classifier script                       ├─► random-reference script (all folds; imports CV runner)
              ▼                                   ├─► numerical review script: H1 (fold 3), G4 (fold 0)
   fold-0 evaluation script                       └─► fixed-alpha script: G3 (8 panels; imports review)
                                                                  │
                         summary/figure scripts (tables 01–04, descriptive results, figures 1–3)
```

## 2. Step-to-module map

| Historical step | Package module | Verified against history |
|---|---|---|
| Guarded obs/var and CSR row access | `compact_panels.data` | access log ⊆ {obs, var, X parts}; row order |
| patient_id → individual mapping, LPS removal, QC exclusion, eligibility | `compact_panels.identity` | cohort invariants (88 individuals, 24,456 cells, by site and class) |
| StratifiedKFold folds, canonical assignment, redraw rule | `compact_panels.partition` | canonical SHA-256 `88be2045…` |
| Stage B withholding, selector/classifier subsamples | `compact_panels.partition` | all eight lists of all four folds (SHA-256) |
| HVG + standardization (two call conventions) | `compact_panels.preprocessing` | candidate genes, means, SDs bit-identical, all folds |
| Target blocks Y_A, Y_B | `compact_panels.targets` | — (exercised by the panel reproduction) |
| Selector, alpha search variants A/B, fallback, checkpoints | `compact_panels.selection` | all eight recorded alpha traces replay bit-identically |
| Logistic-regression classifiers | `compact_panels.classifier` | M2 comparison |
| Per-individual macro-F1, eligibility, FP accounting | `compact_panels.evaluation` | M3 comparison |
| Expression-matched random reference | `compact_panels.random_reference` | M3 comparison |
| G3 / G4 / H1 numerical checks | `compact_panels.stability` | M3 comparison |
| Summary tables and figures | `compact_panels.reporting` | M4 comparison |

## 3. Behaviours preserved on purpose

| # | Historical behaviour | Where it lives |
|---|---|---|
| 1 | **Two alpha-search variants.** Fold 0 (separate recovery script): w=0 bracket `[0.01·α_max, α_max]`, geometric bisection from fit 2; w=1 halving then bisection; at fit 25 *without* a fit above 64 the run stops. Folds 1–3: bracket `[1e-4·α_max, α_max]`, halving until > 64 genes, then bisection; without a fit above 64 the search continues and stops after fit 50. | `selection.initial_state`, `selection.advance`; `conventions` in `configs/study.toml` |
| 2 | **Two scanpy HVG call conventions** (fold 0 with `uns={'log1p': {'base': None}}`, `inplace=False`, strict tie check; folds 1–3 on a copy, `inplace=True`). They gave identical genes historically and do so here. | `preprocessing.select_hvg` |
| 3 | **Warm-started fits that return after zero iterations** (sklearn's pre-loop duality-gap check) count as converged when the gap is below threshold. This is how fold 3 w=1 reached its 67 → 64 fallback. Not repaired. | `selection.Selector` |
| 4 | **Fallback**: largest tested α with > 64 active genes, keep the 64 largest coefficient-row norms (ties by lower index); used once (fold 3, w=1). | `selection.advance`, `selection.top_k_by_norm` |
| 5 | **w=0 uses the 5-column class block only** (not a 71-column target with a zero block). | `targets.targets_for_arm` |
| 6 | Stage B calls `rng.choice` for every stratum, including m = 0, with one generator consumed sequentially. | `partition.stratified_draw` |
| 7 | The random reference matches the **w=0** panel's expression-bin histogram; w=1 is compared with the same reference. | `random_reference.random_panels` |
| 8 | G4 uses seed 20261115 and tol 1e-5 (not the pilot run's own permutation seed); G3 refits at the *saved* alpha only. | `stability`, `Study.run_g3/run_g4` |

## 4. Historical code paths not re-executed

* The **frozen pilot run itself** (its tolerance and permutation gates on the original 71-column
  w=0 target). Its outputs that the analysis used (frozen parameters, subsample lists) are
  reproduced bit-identically; its gate verdict is not part of the reported study.
* The **full cold-start tight-tolerance alpha searches** for G3 (modes `g3_w0`/`g3_w1` of the
  review script) were implemented historically but never run; only the fixed-alpha G3 check was.
* Identity-mapping **audit reports** (library-sharing documentation, a separate identity review) are documentation,
  not analysis inputs.
