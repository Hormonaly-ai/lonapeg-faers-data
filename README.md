# LONAPEG-FAERS: aggregate outputs and analysis code

Supporting data and code for:

**Post-marketing safety of lonapegsomatropin: a disproportionality analysis of US Food and Drug Administration adverse event reports with active comparators and a systematic review** (manuscript under review).

Authors: Fady Hannah-Shmouni (first and corresponding author; Division of Endocrinology, University of British Columbia, and Hormonaly Inc.; ORCID 0000-0001-6829-6205) and Constantine A. Stratakis (Hospital for Endocrine Surgery Research Center, HCA Florida Healthcare, HCA Healthcare Research Institute, Tampa, FL, USA; ORCID 0000-0002-4058-5520).

Prospectively registered protocol: Open Science Framework, doi:[10.17605/OSF.IO/UB6W2](https://doi.org/10.17605/OSF.IO/UB6W2) (registered 23 November 2026, before any lonapegsomatropin-specific FAERS query). A copy of the protocol is in `protocol/`.

## Data source

US FDA Adverse Event Reporting System (FAERS) quarterly ASCII extracts, publicly available at https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html.

- Quarters downloaded for deduplication: 2010Q4–2026Q3 (legacy AERS and FAERS); URL list in `code/upstream_faers_download_dedup/urls_2010Q4_2026Q3.txt`.
- Deduplication: FDA-recommended rule (latest CASEVERSION per CASEID, deleted cases removed).
- Analysis window: 2021Q3–2026Q3 (21 quarters; lonapegsomatropin first US approval August 2021).

## Contents

| Folder | What it holds |
|---|---|
| `outputs/tables/` | Aggregate result tables: report characteristics (Table 1), PT-level disproportionality (ROR, PRR, chi-square) for lonapegsomatropin primary-suspect, suspect and any-role reports and for comparators, prespecified event groups with BCPNN IC and MGPS EBGM, active-comparator RORs (daily somatropin, somatrogon), paediatric vs adult, sensitivity analyses S1–S12, signal robustness, 2026Q3 reporting surge, quarterly counts |
| `outputs/figures/` | Figures 1–4 (PNG) |
| `outputs/aggregate_background/` | Whole-window database totals, quarterly report counts, PT background counts across all retained reports, event-group background counts, fitted MGPS prior |
| `outputs/*.json` | `results.json` (all aggregate numbers used in the manuscript), `derived_numbers.json`, `review_results.json` (systematic-review counts), `S11_openfda_crosscheck.json`; `run_log.txt` (timestamped run log) |
| `review_extraction/` | Systematic review: search results (PubMed, Europe PMC, ClinicalTrials.gov), full-text decisions, case-report extraction (Naranjo, JBI) and trial/registry safety extraction. All from published literature |
| `code/upstream_faers_download_dedup/` | Download, reduction and deduplication of the FAERS quarterly files (shared with the companion tesamorelin FAERS study, OSF H5PBN) |
| `code/lonapeg/` | Study code: `lonaconfig.py` (prespecified drug patterns, event groups, non-clinical PT rules), `lonastats.py` (ROR, PRR, IC, EBGM), `03`–`04` window aggregation, `06` MGPS prior, `07` main analysis, `09`/`12` literature search and synthesis, `13` figures, `14` derived numbers, `16` number audit |

## Reproducing

Python 3.12 with pandas, NumPy, SciPy, pyarrow and Matplotlib.

1. Download and reduce the FAERS quarters: `bash code/upstream_faers_download_dedup/00_pipeline.sh START END` (line numbers of the URL list), which calls `01_process_quarter.py`; then `python 02_deduplicate.py`. The study scripts expect the reduced files at `../../tesamorelin-faers-2026-11/data/processed/` relative to `code/lonapeg/`; edit `TESA` in `03_aggregate_window.py` and `07_analysis.py` to point elsewhere.
2. `python code/lonapeg/03_aggregate_window.py 2021Q3 ... 2026Q3`, then `04_combine.py`, `06_mgps_prior.py`, `07_analysis.py`.
3. `12_review_synthesis.py`, `13_figures.py`, `14_derived_numbers.py`; `16_audit_numbers.py` traces manuscript numbers to the outputs.

## What is not included

Raw FAERS files, reduced case-level parquet files, report-level extracts and any list of FAERS primaryid/caseid values (including the per-report time-to-onset values and the case-level listing used for the supplement) are not redistributed; they can be regenerated from the public FDA files with the code above. Manuscript drafts and correspondence are not included.

## Licence

Code (`code/`): MIT. Data, outputs and documentation: CC BY 4.0. See `LICENSE`.
