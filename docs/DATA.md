# Data provenance and access

No expression data, cell-level records or donor-level source identifiers are distributed in
this repository.

## Source study

Stephenson E, Reynolds G, Botting RA, et al. *Single-cell multi-omics analysis of the immune
response in COVID-19.* Nature Medicine 27, 904–916 (2021). doi:10.1038/s41591-021-01329-2.
Raw and processed data: ArrayExpress **E-MTAB-10026**. Cell-type labels (`author_cell_type`) are the
published annotations; they were not re-derived.

## Inputs required for the real-data reproduction (Level A)

| Input | What it is | How to obtain | Verification |
|---|---|---|---|
| Processed H5AD | Subsampled object, 62,509 cells × 16,299 genes, CSR float32 `X`; obs columns `patient_id`, `Site`, `Status`, `author_cell_type`, `Collection_Day` | Public example dataset of the scverse **pertpy** library (see "Independent access" below) | SHA-256 `2345f3c512552892bb8327f4779dc69972faccac3b00dcf90c0e9d9ae5ab738e`, 700,549,709 bytes |
| Identity mapping | H5AD `patient_id` (sample label) → published biological individual (119 rows) | `cgp derive-inputs` from the public E-MTAB-10026 **SDRF** file and the H5AD | pair SHA-256 `98fef9e2…22aaf` (`configs/inputs.toml`) |
| QC exclusion list | 498 cells of one collection whose library token disagrees with the published sample records | `cgp derive-inputs` (rule: not LPS, `patient_id = CV0201`, `Collection_Day = D28`) | SHA-256 `48c6bf67…b69a3c` |

```bash
cgp derive-inputs --h5ad /path/to/stephenson_2021_subsampled.h5ad \
                  --sdrf /path/to/E-MTAB-10026.sdrf.txt --out data/
```

The command refuses an H5AD with a different hash and reports whether both derived files match
their pins.

## Independent access to the exact H5AD

The file is distributed as the example dataset `pertpy.data.stephenson_2021_subsampled()` of the
scverse **pertpy** library, whose loader downloads

```
https://exampledata.scverse.org/pertpy/stephenson_2021_subsampled.h5ad
```

(loader source: `src/pertpy/data/_datasets.py` in github.com/scverse/pertpy; described there as
"downsampled to approx. 500 cells per donor" and mapped to a healthy-PBMC reference with scArches).

```bash
# download, then hash
curl -L -o stephenson_2021_subsampled.h5ad \
     https://exampledata.scverse.org/pertpy/stephenson_2021_subsampled.h5ad
shasum -a 256 stephenson_2021_subsampled.h5ad   # must print 2345f3c5…ab738e

# or stream and hash without saving the file (Python standard library only)
python -I scripts/verify_public_h5ad.py --record verification.json
```

**Verification status: `PUBLIC_SOURCE_SHA256_VERIFIED` (2026-11-16).**

| Item | Value |
|---|---|
| URL | `https://exampledata.scverse.org/pertpy/stephenson_2021_subsampled.h5ad` |
| Transfer | one HTTP `GET`, 2026-11-16 19:21:30–19:21:40 UTC, HTTP 200, `last-modified: Sun, 14 Sep 2025 15:05:44 GMT` |
| Size | expected 700,549,709 bytes; `content-length` 700,549,709; received 700,549,709 |
| SHA-256 | expected and observed `2345f3c512552892bb8327f4779dc69972faccac3b00dcf90c0e9d9ae5ab738e` |
| Method | `scripts/verify_public_h5ad.py`: standard-library `http.client`, a single request with no redirect followed, no retry and `Accept-Encoding: identity`; the body was read into one fixed 1 MiB buffer, `hashlib.sha256` was updated for every chunk and each chunk was discarded. The file was neither written to disk nor loaded as a dataset. |
| Record | `ledger/public_source_verification.json` |

The publicly hosted file is therefore byte-identical to the input analysed here, as of the
verification date. A file served later at the same address could differ; `cgp validate` and
`cgp derive-inputs` refuse any file whose SHA-256 differs from the pin, so a changed file cannot
silently enter the analysis. If the public file ever differs, the analysis is not reproducible from
it at Level A.

## Why identity mapping matters

Several biological individuals contributed more than one sample under different `patient_id`
labels (three such individuals are eligible). Treating sample labels as donors would place one
person on both sides of a split. All isolation in this pipeline is enforced on the mapped
individual identifier and additionally checked at the sample-label and cell levels.

## Individual codes in the public tables are pseudonyms, not anonymization

`results/table_01_individuals_anonymized.csv` (88 rows) reports per-individual held-out results under
fold-local codes `F<k>_D<nn>`. Despite the historical file name, these codes are **pseudonyms**:

* within each fold, codes are numbered in the sorted order of the source-dataset individual identifier
  (`compact_panels.reporting`);
* the fold assignment is deterministic and reproducible from the public H5AD, the public SDRF, the
  configuration and this code (`cgp validate` re-derives it);
* the table also gives each individual's held-out cell count and number of evaluable classes, which
  can be computed from the public H5AD.

Anyone holding the public data can therefore link every code to the corresponding individual record of
the source dataset (E-MTAB-10026). What the table adds to those public records is the model
performance of each individual in this study. The table contains no names, no direct identifiers and no
source identifiers, and this repository does not publish the code-to-source mapping or any other
source-dataset identifier beyond the single sample label `CV0201` needed to re-derive the QC exclusion
list. Re-identifying the people behind the source records is not possible from this repository alone;
it depends on what the source study itself has published, which this project does not extend. The file
name is kept unchanged because it is part of the byte-identical comparison with the historical
delivery.

## Known provenance limits

* **Subsampling of the processed H5AD is not documented as code.** The distributor describes it only
  as "approx. 500 cells per donor" after scArches mapping; no seed or script is published.
  Reconstructing the file from the original data (Level B) is therefore not possible, and not claimed.
* The QC-excluded collection shows a library-token disagreement with the published records,
  consistent with (but not proof of) a donor-label exchange with another collection; that other
  collection belongs to a non-eligible individual and enters no analysis set.
* Donor-disjoint is not library-disjoint: some held-out individuals share sequencing libraries with
  training individuals (pooled designs at two sites).
