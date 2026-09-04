# Peptide T-Cell Epitope Immunogenicity Prediction Pipeline

A six-stage pipeline that curates T-cell assay data from IEDB/CEDAR, engineers
sequence + physicochemical + MHC-binding + protein-language-model (ESM-2)
features, trains an explainable ensemble classifier that predicts whether a
peptide is immunogenic (elicits a T-cell response), and validates that model
on an independent, non-overlapping external peptide set.

The six notebooks in this repository are meant to be run **in order**; each
one consumes the CSV produced by the previous stage.

## Pipeline overview

| # | Notebook | Purpose | Key input(s) | Key output(s) |
|---|----------|---------|---------------|----------------|
| 1 | `Data_collection_ML_Model.ipynb` | Collect labeled peptides from the IEDB/CEDAR APIs for 12 named antigens (SARS-CoV-2, Influenza, HIV-1, *M. tuberculosis*, EBV) plus 5 organism-level external pathogens (Dengue, Zika, HCV, RSV, HKU15-CoV); verify PRIMARY-target peptides against the UniProt reference sequence; annotate every peptide with predicted MHC class I/II binding (IEDB `netmhcpan`/`netmhciipan` tools) | IEDB/CEDAR/UniProt REST APIs | `dataset_outputs/raw_evidence.csv`, `numbered_pool.csv`, and a balanced `Data_collection_*.csv` (equal Positive/Negative) |
| 2 | `External_Raw_Data_collection.ipynb` | Re-run the same collection logic without class balancing (keep every labeled peptide) to build a large reference pool; cross-check a curated training CSV against a full IEDB export and strip out peptides that already appear in it, so a downstream external/negative corpus doesn't leak into training | Same APIs, plus a full IEDB export CSV | Full unbalanced pool CSV, `iedb_full_dataset_without_duplicates.csv` |
| 3 | `Dataset_Clean_main_Model.ipynb` | Re-query IEDB `tcell_search` per peptide to gather every independent literature record; compute an evidence-agreement "confidence" tier (`high` / `low_single_record` / `conflicting` / `not_found_in_iedb`); relabel high-confidence disagreements to the literature majority, drop unreliable rows, then flag/​resolve near-duplicate peptides (Hamming distance ≤ 2) and add an organism-based CV grouping column | `main_dataset_with_mhc_scores.csv` (output of stage 1) | `main_ml_dataset_clean.csv` — the training-ready dataset |
| 4 | `ML_Full_Model.ipynb` | **Core modeling notebook.** Feature engineering (AAC, DPC, CTD, DHKR, autocorrelation, BLOSUM62, N/C-terminal anchor, MHC-derived, and optional ESM-2 embeddings) → similarity-aware grouped/nested cross-validation → hyperparameter search over Logistic Regression, SVC/NuSVC, Random Forest, Extra Trees, HistGradientBoosting, and XGBoost → weighted ensemble with an MCC-optimized decision threshold → SHAP + LIME explainability → publication figures → statistical robustness checks (bootstrap CIs, Wilcoxon + Benjamini-Hochberg model comparison, class-imbalance robustness, k-mer/motif enrichment, per-feature-family ablation) | `main_dataset_cleaned_2.csv` (output of stage 3) | `final_model.pkl` (deployable model artifact), `figures/`, and ~15 CSV/JSON diagnostic reports under `ml_report_v5/` |
| 5 | `Dataset_Clean_external_validation.ipynb` | Build the external validation set: repeat the stage-3 confidence-scoring steps on a held-out peptide pool, then remove near-duplicates against the training set at the sequence level (RapidFuzz, ≥ 85% identity) | `validation_dataset.csv` (drawn from stage 1/2 with `--exclude-from`) | `external_validation_clean.csv` |
| 6 | `External_Validations.ipynb` | Load `final_model.pkl` and the external peptide CSV, reconstruct the exact training-time feature matrix (including the custom variance/correlation-filter transformer), audit feature-column compatibility, generate predictions, compute external metrics (ROC-AUC, PR-AUC, accuracy, balanced accuracy, sensitivity, specificity, MCC) with bootstrap 95% CIs, plot ROC/PR/confusion-matrix figures, and break performance down by max-sequence-identity-to-training-set and by MHC class/organism subgroup | `final_model.pkl` + `external_validation_peptides.csv` | `external_validation_results/` (predictions, figures, JSON report, zipped for download) |

**Training branch:**
`Data_collection_ML_Model` (collect + MHC-annotate) → `Dataset_Clean_main_Model` (confidence-filter) → `ML_Full_Model` (train + explain)

**External validation branch:**
`External_Raw_Data_collection` (full pool + dedup) → `Dataset_Clean_external_validation` (confidence-filter + fuzzy dedup) → `External_Validations` (score external set, report metrics)

## Repository contents

```
Data_collection_ML_Model.ipynb              Stage 1 — IEDB/CEDAR data collection + MHC binding annotation
External_Raw_Data_collection.ipynb          Stage 2 — Full unbalanced pool + IEDB duplicate removal
Dataset_Clean_main_Model.ipynb              Stage 3 — Literature-confidence filtering & near-duplicate cleanup
ML_Full_Model.ipynb                         Stage 4 — Feature engineering, model training, explainability, robustness checks
Dataset_Clean_external_validation.ipynb     Stage 5 — External validation set construction
External_Validations.ipynb                  Stage 6 — External validation & scoring of the trained model
README.md
requirements.txt
CITATION.cff
```

## Setup

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Notes on dependencies:
- **PyTorch** is CPU-only by default (`requirements.txt` points at the PyTorch
  CPU wheel index). If you have a CUDA GPU and want faster ESM-2 embedding
  in notebooks 04/06, install a CUDA build of `torch` instead before running
  `pip install -r requirements.txt`.
- **ESM-2** (`fair-esm` + `torch`) is only needed if `use_esm_embeddings=True`
  in the `Config` inside notebook 04/06. Set it to `False` to skip the
  protein-language-model features and drop this dependency.
- **XGBoost, SHAP, LIME** are optional at import time in notebook 04 (the
  notebook detects missing packages and skips the corresponding step with a
  warning) but are required to reproduce the full reported pipeline.

## Running the pipeline

1. Run `Data_collection_ML_Model.ipynb` to collect and MHC-annotate the primary training pool.
2. Run `External_Raw_Data_collection.ipynb` against a full IEDB export to build/clean the larger
   reference pool used later for the external split.
3. Run `Dataset_Clean_main_Model.ipynb` on stage 1's output to filter by evidence confidence and
   remove near-duplicate/label-noisy peptides. This produces the CSV that
   `ML_Full_Model.ipynb` trains on.
4. Run `ML_Full_Model.ipynb` end-to-end (each cell is numbered and meant to run in
   order; the first two cells install packages and set thread-limiting
   environment variables and should be run first). This produces
   `final_model.pkl` plus the full diagnostic report.
5. Run `Data_collection_ML_Model.ipynb` again with `--exclude-from <stage-1-output>.csv` (or reuse
   `External_Raw_Data_collection.ipynb`'s pool) to draw a non-overlapping validation pool, then run
   `Dataset_Clean_external_validation.ipynb` on it to build a clean external validation CSV.
6. Run `External_Validations.ipynb`, pointing it at `final_model.pkl` and the
   CSV produced by `Dataset_Clean_external_validation.ipynb`, to score the
   external set and generate the validation report.

### Google Colab cells

A few cells (file upload/download in notebooks 05 and 06) call
`google.colab.files`. These are optional convenience hooks for running in
Google Colab — if you're running locally, replace them with a direct file
path (e.g. `CSV_PATH = "external_validation_peptides.csv"`) and skip the
`from google.colab import files` cell.

## Methodology summary

**Features** (notebook 04): amino-acid composition (AAC), dipeptide
composition (DPC), CTD physicochemical group descriptors, charge/DHKR
fractions, BLOSUM62-based autocorrelation, N/C-terminal anchor residues,
MHC restriction features (percentile rank, IC50, class, allele), and
optional ESM-2 protein language model embeddings.

**Modeling**: nested, repeated, group-aware cross-validation
(`StratifiedGroupKFold`/`GroupKFold`, organism- and similarity-cluster-based
grouping to prevent leakage) with `RandomizedSearchCV` over six model
families, combined into a weighted ensemble whose blend weight and decision
threshold are frozen by optimizing the Matthews correlation coefficient
(MCC) on out-of-fold predictions.

**Explainability**: SHAP (global + per-sample) and LIME (per-sample local
explanations) on the final ensemble.

**Robustness checks**: bootstrap 95% confidence intervals on all metrics,
pairwise model-comparison significance testing (Wilcoxon signed-rank with
Benjamini-Hochberg FDR correction), a class-imbalance robustness sweep, a
k-mer/dipeptide motif-enrichment scan independent of SHAP, and a
per-feature-family ablation (how much AAC/DPC/CTD/ESM-2/etc. contribute on
their own).

**External validation** (notebook 06): the packaged model is re-applied to
an independent peptide set with zero training overlap, evaluated with
bootstrap CIs, and broken down by sequence similarity to the training set
and by MHC class/organism subgroup — the standard way to check whether
performance holds up outside the training distribution.

## Reproducibility

All notebooks default to `random_seed = 42`. Notebook 04 also pins
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
`VECLIB_MAXIMUM_THREADS`, and `NUMEXPR_NUM_THREADS` to `1` at the top of the
notebook to keep BLAS/OpenMP thread counts from interfering with
`n_jobs`-based parallelism in scikit-learn.

## Caveats

- Notebooks 01–03 and 05 depend on live third-party APIs (IEDB, CEDAR,
  UniProt, IEDB's MHC prediction tools). Query behavior, field names, and
  rate limits are outside this repository's control and may change; the
  scripts include schema-discovery and fallback logic for this reason.
- MHC binding predictions in notebook 01 restrict to a fixed 6-allele
  class I / 7-allele class II reference panel; results are specific to
  that panel.

## Citation

If you use this pipeline, please cite it — see `CITATION.cff` (fill in your
name/affiliation and repository URL before publishing).

## License

This repository includes a `LICENSE` file. The code is released under the MIT License.

## Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22301492.svg)](https://doi.org/10.5281/zenodo.22301492)
