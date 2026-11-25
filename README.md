# Genes or donors? Auxiliary supervision and compact gene selection in single-cell transcriptomics

**Independent Research in Computational Biology**

How does an auxiliary learning objective change the genes selected by a sparse model? This project investigates how donor-identity supervision affects compact gene-panel selection in single-cell transcriptomics.

The initial exploratory study found substantially different panels despite similar observed donor-held-out classification performance. A completed, pre-specified follow-up (ED-01) examines whether gene-selection changes recur across donor-disjoint training sets. Neither study has undergone peer review.

**Research author:** Paola Aguirre · **Code:** MIT · **Report, figures and results:** CC BY 4.0

## ED-01 — How auxiliary supervision reshapes gene selection

**Completed November 2026 · Pre-specified computational evaluation**

Using 24,456 PBMC cells from 88 donors, ED-01 compares 64-gene MultiTaskLasso panels selected for cell-class recognition alone (w = 0) versus cell-class and donor-identity recognition (w = 1). The design uses 16 analytical splits into complementary, site-stratified 44-donor halves. Four pilot splits were excluded from the estimates. The protocol and code were frozen before the experimental seed was drawn from a public randomness beacon.

### Key findings

| Result | Estimate |
|---|---:|
| Mean genes replaced per 64-gene panel | 20.56 (32.1%) |
| Net directional concordance between complementary halves | 19.19 genes |
| Net directional recurrence of displacement | 46.7% |
| Nominal 95% Monte Carlo interval for recurrence | 44.0–49.3% |
| Concordant versus reversed gene changes | 310 vs 3 |
| Relative displacement under within-class identity-label permutation | 6.7% |
| Successful analytical panel searches | 112/112 |

The recurrence estimate measures net agreement in the direction of gene gains and losses between complementary halves, relative to the mean number of gains plus losses per fit. Its interval describes re-splitting of this cohort, not uncertainty for new donors.

### Which genes change?

- **Recurrent gains:** HLA-B, HBB, RPL41 and EEF1G were gained in both halves of all 16 analytical splits. Interferon-response genes and RPS4Y1 also recurred.
- **Changes in support:** 655 of 658 gain events involved genes with zero coefficient norm under the cell-class-only objective.
- **Recurrent losses:** included T-cell and NK-cell markers, concentrated toward the low-norm end of the original panels.
- **Exploratory functional patterns:** all six evaluable pre-specified gene sets were over-represented among gains relative to the candidate background in 29–32 of 32 stream-0 fits. No enrichment significance test was performed.

Auxiliary supervision substantially changed selected gene sets, with partially reproducible directional changes across donor-disjoint training sets. The permutation control supports the importance of the alignment between individual labels and cells, but does not distinguish biological from technical sources of donor-associated variation.

### Robustness and limitations

The recurrence estimate was 46.5% using cold-fit sensitivity panels and 45.2% when restricting to the four splits without stream-0 fallback panels (a non-random subset). An independent computational recomputation agreed with the locked analysis to the reported numerical precision.

ED-01 provides pre-specified evidence conditional on successful selection within one previously analysed cohort. Shared sequencing pools limit donor independence. The experiment is not external validation, does not assess classification performance, and does not establish causal mechanisms or biomarker utility.

**Status:** Completed and publicly released. [Scientific report and reproducibility package (ED-01 v1.2)](https://github.com/Paola-aguirre-ruelas/Genes-or-Donors/releases/tag/ed01-v1.2).

## Historical exploratory study

The following methods, results and reproduction instructions refer to the original study, not ED-01.

### Scientific question

When a sparse multi-task model selects a compact 64-gene panel to recognise immune-cell classes,
does also asking it to recognise *which individual* each cell came from (identity supervision, w = 1)
change the selected genes — and does the panel then transfer worse to individuals never seen in
training than a panel selected for cell class alone (w = 0)?

### Approach

* Processed peripheral-blood scRNA-seq subset derived from Stephenson et al. (2021); sample labels
  mapped to biological individuals; eligible individuals with enough cells in at least three of five
  lymphocyte/NK classes.
* One site-stratified, donor-disjoint four-fold partition; training-only within-fold withholding,
  highly variable genes and standardization.
* `MultiTaskLasso` panels of 64 genes per fold and arm; multinomial logistic-regression classifiers.
* Primary descriptive quantity: donor-equal mean paired macro-F1 difference Δ = F1(w=0) − F1(w=1)
  over held-out individuals.
* Exploratory controls: expression-matched random 64-gene panels (20 per fold) and numerical
  sensitivity checks (fixed-alpha tolerance refits, identity-permutation sanity check, fallback diagnosis).

Details: [`docs/METHODS.md`](docs/METHODS.md) · [regenerated Markdown report](report/compact_gene_panels_report_v2.0.md) · [updated four-page PDF (v2.0)](report/exploratory_report_v2.0.pdf)

### Results

<!-- RESULTS:START -->
*Generated by `cgp rebuild-results` from `results/results.json`; descriptive, one cross-validation repetition.*

| Quantity | Value |
|---|---|
| Held-out individuals (all folds) | 88 |
| Mean macro-F1, w=0 (cell class) | 0.907442 |
| Mean macro-F1, w=1 (cell class + identity) | 0.911095 |
| Mean macro-F1, expression-matched random panels | 0.783828 |
| Donor-equal mean paired Δ = F1(w=0) − F1(w=1) | −0.0036525 |
| Individuals with Δ > 0 / < 0 / = 0 | 39 / 49 / 0 |
| Within-fold panel Jaccard(w=0, w=1) | 0.506–0.600 |

A negative Δ means w=1 scored slightly higher; it does not indicate benefit, equivalence or non-inferiority. Identity supervision changed many selected genes with little change in held-out performance.

![Paired differences](figures/figure_01_paired_donor_deltas.png)

**Reproduction (Level A) against the historical records:** overall `EXACT` — 54 components bit-identical, 0 within tolerance, 0 discrepancies (fold plans: 4 EXACT; preprocessing: 4 EXACT; alpha searches and panels: 8 EXACT; classifiers: 8 EXACT; random-reference panels and models: 4 EXACT; held-out scores (main arms): 4 EXACT; held-out scores (random reference): 4 EXACT; stability checks: 10 EXACT; results package (tables, summary values, figures): 8 EXACT.)
<!-- RESULTS:END -->

## Repository layout

| Path | Content |
|---|---|
| `src/compact_panels/` | Package: data access, identity mapping, partitions, preprocessing, panel selection, classifiers, evaluation, random reference, stability checks, comparison, reporting, CLI |
| `configs/` | `study.toml` (all scientific parameters), `synthetic.toml`, `inputs.toml` (input provenance), `acceptance.toml` (comparison tolerances), `historical_references.toml` (SHA-256 pins of the historical artifacts, directory aliases only), `historical_dirs.example.toml`, `paths.example.toml` |
| `tests/` | Synthetic test-suite (no real data) and an optional historical-trace conformance test |
| `results/`, `figures/` | Tables, figures, `results.json`, comparison with the historical reference, SHA-256 manifest; individual codes in `table_01` are linkable pseudonyms (see [`docs/DATA.md`](docs/DATA.md)) |
| `report/` | Regenerated Markdown report and template; updated four-page research PDF (`exploratory_report_v2.0.pdf`) |
| `docs/` | Methods and limitations, data provenance and access, reproducibility guide, historical fidelity |
| `ledger/` | Execution ledger of the reproduction |

## Quick start (no real data needed)

```bash
python3.13 -m venv .venv && source .venv/bin/activate
pip install -e .
(cd tests && python -m unittest -v)
cgp synthetic-demo --out /tmp/cgp-demo
```

## Full reproduction with the real data

Requires the processed H5AD (SHA-256 `2345f3c5…ab738e`; distributed as a public pertpy example dataset
whose SHA-256 and size were confirmed identical to the verified file by a streamed download on 2026-11-16) and the public
ArrayExpress E-MTAB-10026 SDRF; see [`docs/DATA.md`](docs/DATA.md) and
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for commands, runtime, memory, checkpoint recovery
and comparison with the reference results.

```bash
cgp derive-inputs --h5ad H5AD --sdrf SDRF --out data/
cgp validate --deep
cgp reproduce --through train
cgp score --confirm-held-out-access
cgp stability g3 && cgp stability g4 && cgp stability h1
cgp rebuild-results
```

## Reproducibility status

* **Level A — verified-input reproduction:** see the status block above and
  `results/comparison_with_historical_reference.json`.
* **Level B — source-data reconstruction:** not claimed (the subsampling that produced the processed
  H5AD is not published as code).

## Citation and licences

Cite via [`CITATION.cff`](CITATION.cff). Code: [MIT](LICENSE). Report, figures, tables and documentation
prose: [CC BY 4.0](LICENSE-CONTENT.md). No third-party expression data are redistributed.

## AI assistance

AI tools assisted with the project and manuscript preparation. Claude conducted separate computational audits after the results were obtained, and ChatGPT assisted with manuscript review. The author retains full responsibility for the scientific content and conclusions.
