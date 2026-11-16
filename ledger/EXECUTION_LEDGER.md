# Execution ledger

Persistent record of the reproducibility work. Paths are relative to the repository root;
`LEGACY` denotes the historical project tree (read-only). Times are UTC.

## Environment and resources

| Item | Value |
|---|---|
| Machine | Apple M1, 8 cores, 8 GB RAM; ≈ 18 GB free disk at start |
| Interpreter | Python 3.13.18 (verified conda environment), package installed editable into a `--system-site-packages` venv (`.venv`, no dependency changes) |
| Pinned stack | numpy 2.5.4, scipy 1.18.1, pandas 3.1.1, scikit-learn 1.9.1, anndata 0.13.4, h5py 3.16.0, scanpy 1.12.5, matplotlib 3.10.9 |
| Threads | library defaults (as historically); every training stage run sequentially in one process |
| Integrity baseline | `runs/_integrity/legacy_snapshot_pre_M1.json`: SHA-256 of all 432 LEGACY files; frozen manifest 254/254 OK |

## M1 — Engineering the reproducible pipeline (2026-11-16) — `M1_PASS`

Completed:

* Dependency audit of the executed scripts → `docs/HISTORICAL_FIDELITY.md`.
* Package `compact_panels` (src layout): `data`, `identity`, `partition`, `preprocessing`, `targets`,
  `selection`, `classifier`, `evaluation`, `random_reference`, `stability`, `pipeline`, `validate`,
  `reference`, `comparison`, `reporting`, `synthetic`, `cli`, `provenance`, `config`, `io_utils`.
* Versioned configuration: `configs/study.toml` (SHA-256 `672d0195…5d78`), `configs/synthetic.toml`,
  `configs/inputs.toml`, `configs/historical_references.toml` (156 pinned artifacts),
  `configs/acceptance.toml` (tolerances, SHA-256 `0f5c5cda…fcee`, fixed before any M2 output existed).
* Atomic, write-once artifact writes; resumable per-fit checkpoints; held-out access gated by a
  pre-written attempt record.

Commands and evidence:

| Command | Result |
|---|---|
| `pip install --no-deps --no-build-isolation -e .` | installed `compact-gene-panels 1.0.0` |
| `python -m unittest` (tests/, `CGP_LEGACY_ROOT` set) | 43 tests, OK (35 unit + 8 synthetic end-to-end) |
| historical trace replay (`test_historical_conformance`) | all 8 recorded alpha traces derived bit-identically by the state machine |
| `cgp validate --deep` (real inputs) | **PASS 14/14**; 0 models fitted; no held-out read; peak RSS 1.28 GB; 10 s |
| C7c fold plans vs history | EXACT for all 8 lists × 4 folds |
| C9 preprocessing replay | EXACT (HVG indices; means and SDs max abs diff 0.0) for all 4 folds |
| `cgp derive-inputs` (public SDRF + H5AD) | identity-mapping pairs and QC list match their pins; validation with derived inputs PASS 10/10 |
| `scripts/snapshot_tree.py LEGACY --compare …pre_M1.json` | 432/432 unchanged |
| `scripts/check_public_tree.py .` | PASS (no data, credentials, personal paths) |

Deviations: none from the historical analysis. Engineering additions only (guarded reads, checks,
resumable checkpoints, derivation of auxiliary inputs).

## M2 — Full historical model reproduction (2026-11-16) — `M2_PASS`

Tolerances were fixed in `configs/acceptance.toml` (SHA-256 `0f5c5cda…fcee`) before execution.

| Command | Result |
|---|---|
| `cgp reproduce --through train --max-new-fits 50` | exit 0; 15:11:40 → 16:46:17 UTC; wall clock 5,681 s (94.7 min); peak RSS 1.57 GB |
| `cgp compare` (training components) | **overall EXACT — 28/28 components bit-identical** |

Per component (all `EXACT`, maximum absolute difference 0.0):

| Component | Folds / arms | Evidence |
|---|---|---|
| Fold plans | 4 | all eight lists identical (SHA-256) |
| Preprocessing | 4 | HVG indices identical; means/SDs bit-identical |
| Alpha searches → panels | 8 | identical alpha, nnz and n_iter traces, dual gaps, final coefficient matrices and 64-gene sets; fits per panel 5, 8 / 7, 7 / 12, 9 / 3, 25 (fold 0 → 3, w0, w1); fold 3 w1 reproduces the 67 → 64 fallback |
| Classifiers | 8 | coefficients and intercepts bit-identical |
| Random-reference panels and models | 4 × 20 | panels and bin counts identical; 80 models bit-identical |

Resource note: the alpha searches dominate run time (w=1 fits of 71-column targets, up to ≈ 345 s each);
classifiers take < 1 s each.

## M3 — Evaluation and controls (2026-11-16) — `M3_PASS`

Held-out data were read only after every training component had been compared (M2 `EXACT`), so no
held-out result could influence the reproduction procedure.

| Command | Result |
|---|---|
| first launch of the M3 script | failed at argument parsing before any data access (shell did not split the sub-command); no `evaluation_attempt.json` existed afterwards; log archived as `logs/m3_attempt1_argparse_failure.txt` |
| `cgp score --confirm-held-out-access` | exit 0; each fold scored once (attempt record written before reading); 7 s; peak RSS 1.02 GB |
| `cgp stability g3` | exit 0; 1,792 s; peak RSS 1.58 GB; 7 panels `STABLE_AT_FIXED_ALPHA` (support 64, J = 1.0), fold 3 w=1 `FALLBACK_TOP64_STABLE_AT_FIXED_ALPHA` (support 67, J = 0.9552, top-64 J = 1.0) |
| `cgp stability h1` | exit 0; 1,835 s; peak RSS 1.51 GB; active genes 67, 67, 67, 63 at the four recorded alphas; no exact-64 solution |
| `cgp stability g4` | exit 0; 2,926 s; peak RSS 1.30 GB; 6 cold fits, exact 64 at α = 0.03302943; Jaccard 63/65 = 0.9692 vs the w=0 panel (PASS at 0.90) |
| `cgp compare` / `cgp rebuild-results` | **overall EXACT — 54/54 components** |

Per component (all `EXACT`, maximum absolute difference 0.0):

| Component | Count | Evidence |
|---|---:|---|
| Held-out scores, main arms | 4 folds | 88 per-individual macro-F1 values for w=0 and w=1 and their Δ identical; fold means identical |
| Held-out scores, random reference | 4 folds | all 1,760 random-draw scores identical |
| G3 tolerance refits | 8 | verdict, support, top-64 set, n_iter and dual gap identical |
| G4 permutation check | 1 | alpha/nnz/n_iter trace, selected panel and result fields identical |
| H1 fallback diagnosis | 1 | alpha/nnz/n_iter trace and result fields identical |

G3, G4 and H1 are labelled exploratory everywhere they appear (results, report, README); they are not
confirmatory gates. Historical free-text notes in the stability records are not republished (only the
structured fields that both records share are echoed in the comparison file).

## M4 — Automated reporting (2026-11-16) — `M4_PASS`

| Item | Result |
|---|---|
| `cgp rebuild-results` | writes `results/` (tables 01–04, `results.json`, comparison, build record, SHA-256 manifest), `figures/`, the README results block, `report/compact_gene_panels_report_v2.0.md`, `RELEASE_SHA256SUMS.txt`; 2 s |
| Results package vs historical v1.3 delivery | tables 01–04 byte-identical; three figures byte-identical; all summary values exact (8/8 `EXACT`) |
| Report provenance | `report/report_template.md` contains placeholders only; `check_template_has_no_numbers` rejects any decimal or ≥ 3-digit number outside references and labels; every number is inserted from the newly computed `results/results.json` (its SHA-256 is printed in the report header) |
| Determinism | two consecutive rebuilds give identical release files except the timestamped `build_record.json` and comparison JSON |
| Wording changes during M4 | per-panel labels unified to `w=0`/`w=1`; exploratory label made neutral; README notes that byte identity of the public H5AD copy is not yet confirmed |

## M5 — Release preparation (2026-11-16) — `READY_FOR_PUBLICATION_APPROVAL`

Nothing was pushed, published or uploaded; no remote is configured.

| Check | Evidence |
|---|---|
| Tests | 46 tests OK with the historical conformance test enabled; on a fresh clone of the committed tree 46 run, 45 pass and 1 is skipped (conformance needs the historical files) |
| CI | every step of `.github/workflows/ci.yml` (compileall, public-tree check, TOML parse, unit tests, synthetic demo) passes locally on a fresh clone; the workflow has not run on GitHub (no remote) |
| Licences | `LICENSE` MIT (2026, Paola Aguirre) for code; `LICENSE-CONTENT.md` CC BY 4.0 for report, figures, tables and documentation prose; third-party data excluded |
| Data | no expression matrices, cell-level records or raw donor identifiers tracked; `.gitignore` excludes `data/`, `runs/`, `*.h5ad`, local paths; per-individual results use anonymized codes |
| Sensitive content, working tree | `scripts/check_public_tree.py`: PASS (65 files) with a private deny-list kept outside the repository |
| Sensitive content, Git history | all 5 commits authored and committed by `Paola Aguirre <271601722+koalaresearch@users.noreply.github.com>`; four unreachable commits from the earlier identity rewrite (placeholder `author@example.invalid`, file trees identical to the rewritten commits) purged; an earlier version of the public-tree script stored its deny-list as escaped bytes, so every commit was rewritten to carry the current script (final tree unchanged, dates preserved), backup refs deleted, reflog expired, unreachable objects 0; `check_public_tree.py --history`: PASS (127 objects); an additional scan found no home paths, personal e-mail addresses, credentials or institution identifiers; the only sample labels are synthetic test fixtures and the public ArrayExpress label used by the QC rule |
| Integrity | historical project 432/432 files unchanged against the pre-M1 snapshot; frozen delivery manifest 254/254 OK; `RELEASE_SHA256SUMS.txt` 15/15 and `results/SHA256SUMS.txt` 9/9 OK on the fresh clone |
| Regeneration | `cgp rebuild-results` deterministic apart from two timestamped records |
| Claims | README, report v2.0, `docs/` and `CITATION.cff` agree: Level A demonstrated (54/54 `EXACT`), Level B not claimed, G3/G4/H1 exploratory, public H5AD byte identity not yet confirmed by download |

Release-candidate review: `ledger/RELEASE_CANDIDATE_REVIEW.md`.

## Public source verification (2026-11-16) — `PUBLIC_SOURCE_SHA256_VERIFIED`

Authorized by the research owner as one download transfer, solely to compare the publicly hosted H5AD with the
verified local input. This entry supersedes the statements above (M4 wording, M5 claims) that byte identity of the
public copy was not yet confirmed.

| Item | Result |
|---|---|
| Source | `https://exampledata.scverse.org/pertpy/stephenson_2021_subsampled.h5ad` |
| Transfer | exactly one HTTP `GET` (19:21:30–19:21:40 UTC, 10.2 s, about 68 MB/s); HTTP 200; no redirect, no retry |
| Size | expected 700,549,709 bytes; `content-length` 700,549,709; observed 700,549,709 |
| SHA-256 | expected `2345f3c512552892bb8327f4779dc69972faccac3b00dcf90c0e9d9ae5ab738e`; observed identical |
| Method | `scripts/verify_public_h5ad.py` (standard-library `http.client` and `hashlib`); body read into one fixed 1 MiB buffer and discarded after each digest update; nothing written to disk; dataset never loaded |
| Record | `ledger/public_source_verification.json` (copy of the script's JSON output) |
| No second copy | the only H5AD on disk is the original input; a search for files over 100 MB created during the transfer found none |
| Original input | unchanged: same size, modification time and inode before and after; re-hashed to the pinned SHA-256 afterwards |
| Scientific outputs | no model, evaluation or rebuild was run; `results/`, `figures/` and report v2.0 unchanged |
| Historical files | 432/432 unchanged against the pre-M1 snapshot; frozen delivery manifest 254/254 OK; report manifests v1.2 8/8, v1.3 9/9, release 14/14 OK |
| Documentation | `docs/DATA.md` (verification table and streaming command), `README.md` (wording), `ledger/RELEASE_CANDIDATE_REVIEW.md` (verification section; owner decision 3 resolved) |

Status remains `READY_FOR_PUBLICATION_APPROVAL`. Nothing was pushed, published or uploaded.

## Publication candidate (2026-11-16) — `READY_FOR_INDEPENDENT_RELEASE_REVIEW`

Packaging, documentation and verification only: no model was trained, no evaluation rerun and no scientific
result changed.

| Item | Result |
|---|---|
| Owner decisions | destination `koalaresearch/compact-gene-panels-donor-generalization`; public contact in `CITATION.cff`; no-reply commit identity kept; `CV0201` retained; release date kept; AI declaration restored to the approved v1.3 wording; legacy folder names removed |
| Historical references | `configs/historical_references.toml` uses directory aliases; the git-ignored `configs/historical_dirs.local.toml` maps them locally (template `configs/historical_dirs.example.toml`); 156 pins and semantics unchanged; `validate.py` takes its protected directories from the same mapping |
| Rebuild | `cgp rebuild-results` from the existing run outputs: comparison 54/54 `EXACT`; `results.json`, tables and figures byte-identical; report changed only in the AI declaration |
| History transformation | one `git filter-branch --tree-filter` pass over all 8 commits replaced 7 legacy folder names with the aliases and 2 legacy file-name stems with neutral names; final tree identical; authors, dates and messages preserved; `refs/original` deleted, reflog expired, `git gc --prune=now`, 0 unreachable objects; private pre-rewrite bundle and commit map kept under the git-ignored `runs/_integrity/` |
| Checks | 46 tests OK (45 + 1 skipped without historical artifacts); synthetic demo OK; `cgp validate --deep` 14/14; public-tree and full-history scan PASS; historical project 432/432 unchanged; frozen manifest 254/254 |
| Export | review ZIP made with `git archive` from the final commit; private history bundle kept separately; checksums reported with the export |

Nothing was pushed, published or uploaded; no remote is configured.

## Documentation revision r2 (2026-11-16) — `READY_FOR_INDEPENDENT_RELEASE_REVIEW`

Documentation only: no model was trained, rescored or rebuilt, and no file in `results/`, `figures/` or
`report/` changed.

| Item | Result |
|---|---|
| Individual codes | reviewed `results/table_01_individuals_anonymized.csv` (88 rows, columns `donor_code`, fold, cell and class counts, scores). The `F<k>_D<nn>` codes follow the sorted source identifier within each fold, and the fold assignment and cell counts can be recomputed from the public data, so the codes are linkable pseudonyms, not anonymization. Documented in `docs/DATA.md` (new section), the README layout table and the `reporting.py` docstring. No source identifier or mapping was added; the file name is kept because it is part of the byte-identical historical comparison |
| `cgp rebuild-results` | `docs/REPRODUCIBILITY.md` and the README now state that it overwrites `results/`, `figures/`, report v2.0, the README results block and `RELEASE_SHA256SUMS.txt` without confirmation, and recommend a clean working copy (fresh clone or extracted archive) for independent reruns; `cgp report --out DIR` is given as the non-overwriting alternative |
| Checks | `RELEASE_SHA256SUMS.txt` 15/15 and `results/SHA256SUMS.txt` 9/9 OK; `results.json` SHA-256 `a0d0dbc5…eef7eef`; committed release files unchanged; comparison 54/54 `EXACT`; 46 tests OK (45 + 1 skipped without historical artifacts); synthetic demo OK; public-tree and full-history scan PASS; historical project 432/432 unchanged |

Nothing was pushed, published or uploaded; no remote is configured.
