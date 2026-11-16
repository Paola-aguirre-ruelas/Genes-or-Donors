# Release-candidate review — compact-gene-panels 1.0.0

**Status: `READY_FOR_INDEPENDENT_RELEASE_REVIEW`** (2026-11-16). Local repository only; no remote, nothing pushed or published.

## What the release contains

Pipeline package `compact_panels` with the `cgp` command, versioned configuration, a synthetic test-suite,
regenerated results (`results/`, `figures/`), research report v2.0 rendered from those results, data
provenance and reproducibility documentation, and the execution ledger. No third-party data.

## Scientific status

| Milestone | Status | Evidence |
|---|---|---|
| M1 pipeline engineering | PASS | 46 tests; `cgp validate --deep` 14/14; historical files unchanged |
| M2 panels and classifiers | PASS | 28/28 training components bit-identical (4 plans, 4 preprocessing, 8 alpha searches and panels, 8 classifiers, 4 × 20 random panels and models) |
| M3 evaluation and controls | PASS | held-out scores for 88 individuals and 1,760 random draws bit-identical; G3 ×8, G4, H1 bit-identical (exploratory) |
| M4 reporting | PASS | tables 01–04 and three figures byte-identical to the historical delivery; report numbers inserted from newly computed results |
| M5 release preparation | READY | tests, CI steps, licences, sensitive-content and history audit, integrity and claims checks in the ledger |

Classification: **exact reproduction** of every compared component (54/54 `EXACT`, maximum difference 0.0);
no component needed the numerical-agreement tolerances; no discrepancies. Level A (verified-input) only;
Level B (source-data reconstruction) is not claimed.

## Public source verification

**`PUBLIC_SOURCE_SHA256_VERIFIED`** (2026-11-16, authorized by the research owner as a one-time transfer).
One HTTP `GET` of `https://exampledata.scverse.org/pertpy/stephenson_2021_subsampled.h5ad` returned
700,549,709 bytes (expected 700,549,709) with SHA-256
`2345f3c512552892bb8327f4779dc69972faccac3b00dcf90c0e9d9ae5ab738e` (identical to the pinned digest).
The body was streamed through a fixed 1 MiB buffer into `hashlib.sha256` and discarded; no copy was written
to disk and the dataset was not loaded. The existing verified H5AD was not touched (same size, modification
time and inode before and after), no model or evaluation was rerun, and the historical project is still
432/432 unchanged. Method and record: `scripts/verify_public_h5ad.py`, `ledger/public_source_verification.json`,
`docs/DATA.md`.

## Limitations to keep visible

* One cross-validation repetition; descriptive, not confirmatory; identity is confounded with library,
  site, disease state, sex and composition.
* The historical alpha-search difference between fold 0 and folds 1–3, zero-iteration warm-start fits and
  the single fallback panel are reproduced as they were, not repaired.
* Bit-identity was shown on one machine (Apple M1, pinned stack); other BLAS builds or versions may change
  floating-point paths.
* The public pertpy copy of the H5AD was confirmed byte-identical on 2026-11-16 only; a file served later at
  the same address could differ (the pipeline refuses any file whose SHA-256 differs from the pin).

* `results/table_01_individuals_anonymized.csv` uses fold-local codes that are **linkable pseudonyms**:
  anyone with the public data can recompute which source-dataset record each code denotes. No source
  identifiers or mapping are published (see `docs/DATA.md`).
* `cgp rebuild-results` overwrites the committed release files in place; independent reruns should use a
  clean working copy (`docs/REPRODUCIBILITY.md`).

## Owner decisions applied (2026-11-16)

1. **Destination:** GitHub account `koalaresearch`, repository `compact-gene-panels-donor-generalization`
   (`repository-code` in `CITATION.cff`). The remote is not created and nothing is pushed until the
   independent review of the exported candidate is complete.
2. **Authorship and contact:** Paola Aguirre. The approved public research contact
   `paolajanetaguirre@gmail.com` is in `CITATION.cff` (`authors`, `contact`). Commits keep the GitHub no-reply
   identity `Paola Aguirre <271601722+koalaresearch@users.noreply.github.com>`.
3. **Public H5AD:** `PUBLIC_SOURCE_SHA256_VERIFIED` (see above).
4. **Public sample label:** `CV0201` is retained in `configs/inputs.toml` and `docs/DATA.md`; it is needed to
   re-derive the QC exclusion list from public records.
5. **Historical directory names:** removed from the public repository and its history.
   `configs/historical_references.toml` now refers to the historical artifacts through neutral directory
   aliases (`@splits`, `@fold0_partition`, `@fold0_run`, `@cv_runs`, `@random_reference`,
   `@numerical_checks`, `@delivery`). The optional local file `configs/historical_dirs.local.toml`
   (git-ignored; template `configs/historical_dirs.example.toml`) maps them to local directories. All 156
   SHA-256 pins, every parameter and the comparison semantics are unchanged; the comparison recomputed
   through the aliases is again 54/54 `EXACT`. The private historical project was not modified.
6. **AI declaration:** report v2.0 now carries the approved v1.3 wording exactly; the sentence added during
   the reproduction was removed.
7. **Release date:** `CITATION.cff` keeps `date-released: 2026-11-16`.
8. **Licences:** MIT for original code (`LICENSE`); CC BY 4.0 for the report, figures and scientific content
   (`LICENSE-CONTENT.md`); no third-party expression data are redistributed.

## Remaining decision

Publication (create the remote, push `main`, publish a release) after the independent review of the exported
publication candidate.

## Prepublication checks for the publication candidate (2026-11-16)

| Check | Result |
|---|---|
| Working tree | clean on `main`; every release file committed |
| Tests | 46 tests OK with the historical artifacts available; without them 45 pass and 1 is skipped |
| Synthetic demonstration | `cgp synthetic-demo` completes |
| Real-data validation | `cgp validate --deep` 14/14 (models fitted: 0; held-out expression read: no) |
| Historical comparison | recomputed through the directory aliases: 54/54 `EXACT`, 0 numerical agreement, 0 discrepancies |
| Scientific numbers | `results/results.json` (SHA-256 `a0d0dbc5…eef7eef`), tables 01–04 and the three figures byte-identical to the previous commit; only the report's AI declaration, the timestamped records and the path labels in the comparison file changed |
| Licences and metadata | MIT `LICENSE` and CC BY 4.0 `LICENSE-CONTENT.md` (Paola Aguirre, 2026); `CITATION.cff` with the approved contact, repository URL and `date-released: 2026-11-16` |
| README | approved research positioning present verbatim |
| Data instructions | public URL, size and SHA-256 in `docs/DATA.md`, verified by the one-time streamed download |
| Sensitive content | `scripts/check_public_tree.py --history` PASS on the working tree and every reachable object, with a private deny-list of institution, supervisor, project-code, local-account and legacy folder names; no H5AD, local paths, credentials or raw donor identifiers |
| Git history | the 8 local commits were rewritten once to replace legacy internal folder and file names with the neutral aliases (final tree unchanged; authors, committers, dates and messages preserved); backup refs deleted, reflog expired, 0 unreachable objects. The pre-rewrite history and the old-to-new commit map are kept privately outside the repository. |
| Historical project | 432/432 files unchanged against the pre-M1 snapshot; frozen delivery manifest 254/254 OK |
| Documentation revision r2 | pseudonym linkability of the `table_01` codes and the overwrite behaviour of `cgp rebuild-results` documented; release files unchanged (`RELEASE_SHA256SUMS.txt` 15/15, `results/SHA256SUMS.txt` 9/9); tests, scan and historical checks rerun and passing |
