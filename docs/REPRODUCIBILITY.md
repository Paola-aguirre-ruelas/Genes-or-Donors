# Reproducing the study

## Reproducibility levels

| Level | Meaning | Status |
|---|---|---|
| **A — verified-input reproduction** | Re-execute the whole analysis from the exact, hash-verified processed H5AD | **Demonstrated** (see `results/comparison_with_historical_reference.json` and the ledger) |
| **B — source-data reconstruction** | Rebuild that H5AD from publicly distributed original data | **Not claimed**: the subsampling/preprocessing provenance of the processed object is undocumented |

Data-access prerequisite for Level A: the processed H5AD with SHA-256
`2345f3c512552892bb8327f4779dc69972faccac3b00dcf90c0e9d9ae5ab738e`, plus the public ArrayExpress
E-MTAB-10026 SDRF from which the two auxiliary inputs are derived (`docs/DATA.md`). Without the
H5AD, everything except the real-data run can still be executed: the test-suite and the synthetic
end-to-end demonstration.

## 1. Environment

The analysis was verified with Python 3.13.18 and the exact versions in `pyproject.toml`
(complete package set in `requirements-lock.txt`). Other versions can change floating-point paths
in Scanpy and scikit-learn and therefore panel membership; `cgp validate` reports any mismatch.

```bash
python3.13 -m venv .venv && source .venv/bin/activate
pip install -e .            # exact pins from pyproject.toml
```

## 2. Without real data (minutes)

```bash
cd tests && python -m unittest -v && cd ..
cgp synthetic-demo --out /tmp/cgp-demo
```

The demonstration generates a synthetic dataset with the same structure (sites, multi-sample
individuals, an excluded-status individual, a QC list, ineligible individuals), then runs validation,
all four folds of panel search, classifiers, random reference, scoring, the G3/G4 checks and the
results package.

## 3. Real-data reproduction

```bash
cgp derive-inputs --h5ad /path/to/stephenson_2021_subsampled.h5ad \
                  --sdrf /path/to/E-MTAB-10026.sdrf.txt --out data/
cp configs/paths.example.toml configs/paths.local.toml   # edit the four paths
cgp validate --deep                                        # read-only; fits nothing
cgp reproduce --through train                              # 8 panels, 8 classifiers, 80 random models
cgp score --confirm-held-out-access                        # reads each held-out fold once
cgp stability g3 && cgp stability g4 && cgp stability h1   # numerical checks (training data only)
cgp rebuild-results                                        # tables, figures, report, manifests
```

With `legacy_root` set to a copy of the historical project and the optional local file
`configs/historical_dirs.local.toml` (git-ignored; template `configs/historical_dirs.example.toml`)
mapping each `@alias` of `configs/historical_references.toml` to a directory of that project,
`cgp rebuild-results` also compares every component with the historical records (`cgp compare` does
this alone). Every historical file is verified against its SHA-256 pin before use; without the local
mapping the comparison is skipped. The historical artifacts themselves are not distributed.

### Checkpoints and recovery

* Every alpha-search fit writes `fold<k>/checkpoint_<arm>.npz` atomically (state, coefficients, best
  norms). Rerunning the same command resumes after the last completed fit; the warm-start path is
  bit-identical to an uninterrupted run (tested).
* Completed stages are verified and skipped; run-directory artifacts are never overwritten.
* **`cgp rebuild-results` overwrites the repository's committed release files in place**, in the
  directory given by `--repo-root` (default: the current directory): everything in `results/` and
  `figures/`, `report/compact_gene_panels_report_v2.0.md`, the results block of `README.md` and
  `RELEASE_SHA256SUMS.txt`. It does not ask for confirmation. All of these are deterministic except
  `results/build_record.json` and `results/comparison_with_historical_reference.json`, which carry
  timestamps.
* **For an independent rerun, use a clean working copy** — a fresh clone or the extracted release
  archive — not a checkout you need to keep. After the rerun, `git status` and `git diff --stat` show
  exactly which release files changed (normally only the two timestamped records). Without Git, copy
  `RELEASE_SHA256SUMS.txt` aside before the rebuild and run `shasum -a 256 -c` on that copy afterwards. To build the tables, figures and `results.json` without touching the repository, use
  `cgp report --out DIR` instead.
* A checkpoint whose status is `NONCONVERGED_REVIEW_REQUIRED`, `MEMORY_GATE_REVIEW_REQUIRED` or
  `FIT_LIMIT_REVIEW_REQUIRED` stops the run for review; nothing is retried automatically.
* Scoring writes `evaluation_attempt.json` before reading held-out expression and refuses to run
  twice; training stages refuse to run once a fold has been opened for evaluation.
* A run directory is bound to one configuration hash and one set of input hashes.

### Resources (measured on an Apple M1 with 8 GB RAM; see ledger)

Each stage ran sequentially in one process with library-default threading; time is wall clock, memory
is the maximum resident set size reported by `/usr/bin/time -l`.

| Stage | Wall clock | Peak memory |
|---|---:|---:|
| `cgp validate --deep` (read-only) | 10 s | 1.28 GB |
| `cgp reproduce --through train` (8 alpha searches, 8 classifiers, 80 random models) | 5,681 s (94.7 min) | 1.57 GB |
| `cgp score --confirm-held-out-access` (4 folds) | 7 s | 1.02 GB |
| `cgp stability g3` (8 cold refits at tol 10⁻⁵) | 1,792 s (29.9 min) | 1.58 GB |
| `cgp stability h1` (4 cold fits) | 1,835 s (30.6 min) | 1.51 GB |
| `cgp stability g4` (6 cold fits) | 2,926 s (48.8 min) | 1.30 GB |
| `cgp rebuild-results` (comparison, tables, figures, report) | 2 s | 0.18 GB |
| **Total real-data computation** | **≈ 3.4 h** | **< 1.6 GB** |
| Test-suite (46 tests) / `cgp synthetic-demo` (no real data) | 7 s / 5 s | 0.56 GB |

The w=1 alpha-search fits (71-column targets) dominate: up to about 345 s per fit. Disk: the input
H5AD is 700.5 MB; the complete run directory (plans, preprocessing, checkpoints, panels, 88 models,
scores, stability records, comparisons, logs) occupies about 1.6 MB.

## 4. Comparing with the historical reference

`configs/acceptance.toml` fixes the tolerances (set before any reproduced output existed).
Components are classed `EXACT` (identical discrete outcomes and bit-identical numbers),
`NUMERICAL_AGREEMENT` (identical discrete outcomes, continuous values within tolerance) or
`DISCREPANCY`. The comparison covers fold plans, preprocessing, all alpha traces and panels, classifier
parameters, random panels and their 80 models, per-individual held-out scores (main arms and every
random draw), and the G3/G4/H1 checks.
