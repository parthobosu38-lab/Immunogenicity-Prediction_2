# Peptide T-Cell Epitope Immunogenicity Prediction Pipeline

A six-stage pipeline that curates T-cell assay data from IEDB/CEDAR, engineers
sequence + physicochemical + MHC-binding + protein-language-model (ESM-2)
features, trains an explainable ensemble classifier that predicts whether a
peptide is immunogenic (elicits a T-cell response), and validates that model
on an independent, non-overlapping external peptide set.

The repository is organised as nine numbered folders that follow the order of the
work. Folders `01`–`05` hold the six core pipeline notebooks, which are meant to be
run **in order**; each one consumes the CSV produced by the previous step. Folders
`06`–`09` hold four supplementary analyses that build on the outputs of the core
pipeline (see [Supplementary analyses](#supplementary-analyses)).

Folder names contain spaces (for example `03. ML_Model`), so quote them in a shell:
`cd "03. ML_Model"`.

## Pipeline overview

| Folder | Notebook | Purpose | Key input(s) | Key output(s) |
|---|----------|---------|---------------|----------------|
| `01. Raw_Data_collection` | `Data_collection_ML_Model.ipynb` | Collect labeled peptides from the IEDB/CEDAR APIs for 12 named antigens (SARS-CoV-2, Influenza, HIV-1, *M. tuberculosis*, EBV) plus 5 organism-level external pathogens (Dengue, Zika, HCV, RSV, HKU15-CoV); verify PRIMARY-target peptides against the UniProt reference sequence; annotate every peptide with predicted MHC class I/II binding (IEDB `netmhcpan`/`netmhciipan` tools) | IEDB/CEDAR/UniProt REST APIs | `dataset_outputs/raw_evidence.csv`, `numbered_pool.csv`, and a balanced `Data_collection_*.csv` (equal Positive/Negative) |
| `01. Raw_Data_collection` | `External_Raw_Data_collection.ipynb` | Re-run the same collection logic without class balancing (keep every labeled peptide) to build a large reference pool; cross-check a curated training CSV against a full IEDB export and strip out peptides that already appear in it, so a downstream external/negative corpus doesn't leak into training | Same APIs, plus a full IEDB export CSV | Full unbalanced pool CSV, `iedb_full_dataset_without_duplicates.csv` |
| `02. ML_Dataset_Clean` | `Dataset_Clean_main_Model.ipynb` | Re-query IEDB `tcell_search` per peptide to gather every independent literature record; compute an evidence-agreement "confidence" tier (`high` / `low_single_record` / `conflicting` / `not_found_in_iedb`); relabel high-confidence disagreements to the literature majority, drop unreliable rows, then flag/​resolve near-duplicate peptides (Hamming distance ≤ 2) and add an organism-based CV grouping column | `main_dataset_with_mhc_scores.csv` (output of folder `01`) | `main_ml_dataset_clean.csv` — the training-ready dataset |
| `03. ML_Model` | `ML_Full_Model.ipynb` | **Core modeling notebook.** Feature engineering (AAC, DPC, CTD, DHKR, autocorrelation, BLOSUM62, N/C-terminal anchor, MHC-derived, and optional ESM-2 embeddings) → similarity-aware grouped/nested cross-validation → hyperparameter search over Logistic Regression, SVC/NuSVC, Random Forest, Extra Trees, HistGradientBoosting, and XGBoost → weighted ensemble with an MCC-optimized decision threshold → SHAP + LIME explainability → publication figures → statistical robustness checks (bootstrap CIs, Wilcoxon + Benjamini-Hochberg model comparison, class-imbalance robustness, k-mer/motif enrichment, per-feature-family ablation) | `main_dataset_cleaned_2.csv` (output of folder `02`) | `final_model.pkl` (deployable model artifact), `figures/`, and ~15 CSV/JSON diagnostic reports under `ml_report_v5/` |
| `04. External_Dataset_Clean` | `Dataset_Clean_external_validation.ipynb` | Build the external validation set: repeat the folder-`02` confidence-scoring steps on a held-out peptide pool, then remove near-duplicates against the training set at the sequence level (RapidFuzz, ≥ 85% identity) | `validation_dataset.csv` (drawn from folder `01` with `--exclude-from`) | `external_validation_clean.csv` |
| `05. External_Validation_Result` | `External_Validations.ipynb` | Load `final_model.pkl` and the external peptide CSV, reconstruct the exact training-time feature matrix (including the custom variance/correlation-filter transformer), audit feature-column compatibility, generate predictions, compute external metrics (ROC-AUC, PR-AUC, accuracy, balanced accuracy, sensitivity, specificity, MCC) with bootstrap 95% CIs, plot ROC/PR/confusion-matrix figures, and break performance down by max-sequence-identity-to-training-set and by MHC class/organism subgroup | `final_model.pkl` + `external_validation_peptides.csv` | `external_validation_results/` (predictions, figures, JSON report, zipped for download) |

**Training branch:**
`Data_collection_ML_Model` (collect + MHC-annotate) → `Dataset_Clean_main_Model` (confidence-filter) → `ML_Full_Model` (train + explain)

**External validation branch:**
`External_Raw_Data_collection` (full pool + dedup) → `Dataset_Clean_external_validation` (confidence-filter + fuzzy dedup) → `External_Validations` (score external set, report metrics)

In folder terms: the training branch is `01` → `02` → `03`; the external validation
branch is `01` → `04` → `05`. Folder `01` holds two notebooks (the balanced primary
collection and the larger unbalanced pool); every other core folder holds one.

## Supplementary analyses

Four additional notebooks (folders `06`–`09`) test the trained model and the curation
choices more closely. They are **not** part of the core chain: each one reads outputs
of the core pipeline and can be run on its own, after folders `02`, `03` and `05` are done.

| Folder | Notebook | Purpose | Key input(s) | Key output(s) |
|---|----------|---------------|----------------|----------------|
| `06. Benchmark` | `benchmark_colab_v2.ipynb` | Benchmark the ensemble against external tools and simple baselines on the **same** external peptides: NetMHCpan/NetMHCIIpan percentile rank, PRIME 2.1 (`Score_bestAllele` and `%Rank_bestAllele`), IEDB class I immunogenicity (Calis, 9-mers only), and baselines (shorter peptide, predicted class, rank + class + length fitted on the training set only). Confidence intervals use a cluster bootstrap (similar peptides are resampled together) and the ensemble-minus-score differences are paired. The helper script `benchmark_tools.py` is written by the notebook's first code cell | `external_validation_clean.csv`, `external_validation_predictions.csv` (output of folder `05`, with a `predicted_probability` column), `train718_with_cluster_groups.csv` | `bench/` folder with `benchmark_table.csv`, `benchmark_summary.txt`, `benchmark_scores_per_peptide.csv` |
| `07. ESM2_ablation_redundancy` | `esm2_ablation_redundancy_local.ipynb` | Feature-block ablation (ESM2, MHC_ANCHOR, CLASSICAL; 7 configurations, fast fixed Extra Trees mode or the full tuned ensemble) with Nadeau–Bengio corrected CIs and paired tests, plus a redundancy analysis of ESM-2 against the hand-engineered features (linear CKA with permutation null, effective dimensionality, out-of-fold ridge R², CCA, Spearman, residualised-ESM test) | `main_dataset_for_ml_model.csv` (the cleaned training dataset) | `ablation_*_table.csv`, `ablation_*_increments.csv`, `redundancy_report.txt`, `redundancy_*.csv/.json`, `figures/` |
| `08. Similarity_Threshold_Sensitivity` | `similarity_threshold_sensitivity.ipynb` | Sensitivity of the 80% sequence-identity cutoff used for leakage control. Diagnostics (cluster/cohort counts, residual close neighbours across clusters, external cutoff sweep); **Design A** keeps the cohort fixed and changes only the CV grouping (`naive`, 0.90 … 0.50, `core9` = identity ≥ 0.80 **or** a shared 9-aa stretch); **Design B** re-runs the full pipeline at other thresholds, including external scoring | `main_dataset_for_ml_model.csv`, `external_validation_clean.csv`, optional `ext_pred_merged.csv` (baseline external predictions) | `sens_r24/` folder with `designA_table.csv`, `designB_table.csv`, `diag_*.csv`, `figR24_designA.png` |
| `09. Consensus_Threshold_Sensitivity` | `consensus_threshold_sensitivity.ipynb` | Sensitivity of all results to the 70% label-confidence threshold and the 75% cluster-consensus threshold. Scenarios: `base`, `cons_0.60`, `conf_0.80`, `conf_1.00`, optionally original (un-relabeled) labels and other similarity thresholds. Each scenario runs the full pipeline and scores the same external cohort; differences are tested with a paired bootstrap on shared external peptides | `main_dataset_for_ml_model.csv`, `external_validation_clean.csv` | `sens/` folder with `sensitivity_table.csv`, `cohort_table.csv`, `sensitivity_plot.png` |

## Repository contents

```
01. Raw_Data_collection/
    Data_collection_ML_Model.ipynb           IEDB/CEDAR data collection + MHC binding annotation
    External_Raw_Data_collection.ipynb       Full unbalanced pool + IEDB duplicate removal
02. ML_Dataset_Clean/
    Dataset_Clean_main_Model.ipynb           Literature-confidence filtering & near-duplicate cleanup
03. ML_Model/
    ML_Full_Model.ipynb                      Feature engineering, model training, explainability, robustness checks
04. External_Dataset_Clean/
    Dataset_Clean_external_validation.ipynb  External validation set construction
05. External_Validation_Result/
    External_Validations.ipynb               External validation & scoring of the trained model
06. Benchmark/
    benchmark_colab_v2.ipynb                 Benchmark vs NetMHCpan rank / PRIME / IEDB and simple baselines
07. ESM2_ablation_redundancy/
    esm2_ablation_redundancy_local.ipynb     Feature-block ablation + ESM-2 redundancy analysis
08. Similarity_Threshold_Sensitivity/
    similarity_threshold_sensitivity.ipynb   Sensitivity of the sequence-identity cutoff
09. Consensus_Threshold_Sensitivity/
    consensus_threshold_sensitivity.ipynb    Sensitivity to label-confidence and cluster-consensus thresholds
README.md
requirements.txt
CITATION.cff
LICENSE
```

The notebook in `06. Benchmark` also uses a helper script, `benchmark_tools.py`, which the
notebook writes itself in its first code cell (keeping a copy next to it is optional).

## Setup

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Notes on dependencies:
- **PyTorch** is CPU-only by default (`requirements.txt` points at the PyTorch
  CPU wheel index). If you have a CUDA GPU and want faster ESM-2 embedding
  in folders `03` and `05`, install a CUDA build of `torch` instead before running
  `pip install -r requirements.txt`.
- **ESM-2** (`fair-esm` + `torch`) is only needed if `use_esm_embeddings=True`
  in the `Config` inside the notebooks of folders `03` and `05`. Set it to `False` to skip the
  protein-language-model features and drop this dependency.
- **XGBoost, SHAP, LIME** are optional at import time in folder `03` (the
  notebook detects missing packages and skips the corresponding step with a
  warning) but are required to reproduce the full reported pipeline.
- **scikit-learn is pinned below 1.8** (`scikit-learn>=1.3,<1.8`). Newer versions
  produce different `StratifiedGroupKFold` folds, so cross-validation results would
  not match the reported run.
- **RapidFuzz** is also used by folders `06`, `08` and `09` (edit-identity
  similarity) and **LIME/SHAP** by folders `03` and `08`.
- **PRIME and MixMHCpred** (used only by folder `06`) are not pip packages; see
  "Running the supplementary analyses" below.

## Running the pipeline

1. Run `Data_collection_ML_Model.ipynb` (folder `01. Raw_Data_collection`) to collect and MHC-annotate the primary training pool.
2. Run `External_Raw_Data_collection.ipynb` (folder `01. Raw_Data_collection`) against a full IEDB export to build/clean the larger
   reference pool used later for the external split.
3. Run `Dataset_Clean_main_Model.ipynb` (folder `02. ML_Dataset_Clean`) on the output of step 1 to filter by evidence confidence and
   remove near-duplicate/label-noisy peptides. This produces the CSV that
   `ML_Full_Model.ipynb` trains on.
4. Run `ML_Full_Model.ipynb` (folder `03. ML_Model`) end-to-end (each cell is numbered and meant to run in
   order; the first two cells install packages and set thread-limiting
   environment variables and should be run first). This produces
   `final_model.pkl` plus the full diagnostic report.
5. Run `Data_collection_ML_Model.ipynb` again (folder `01. Raw_Data_collection`) with `--exclude-from <stage-1-output>.csv` (or reuse
   `External_Raw_Data_collection.ipynb`'s pool) to draw a non-overlapping validation pool, then run
   `Dataset_Clean_external_validation.ipynb` (folder `04. External_Dataset_Clean`) on it to build a clean external validation CSV.
6. Run `External_Validations.ipynb` (folder `05. External_Validation_Result`), pointing it at `final_model.pkl` and the
   CSV produced by `Dataset_Clean_external_validation.ipynb`, to score the
   external set and generate the validation report.

### Running the supplementary analyses

The inputs are CSVs produced by the core pipeline (see the table above); put them in
the folder the notebook uses as its working directory.

- **`06. Benchmark`** is written for Google Colab. Run it cell by cell (not "Run all"):
  create the script → upload the three CSVs → `prepare` → PRIME → IEDB → `analyze` →
  download. PRIME 2.1 and MixMHCpred 3.x are cloned from GitHub by the `setup-prime`
  step only if you pass `--accept-license`.
- **`07. ESM2_ablation_redundancy`, `08. Similarity_Threshold_Sensitivity` and `09. Consensus_Threshold_Sensitivity`** are written for local Jupyter. Set the path near
  the top (`BASE_DIR` or the data path in the settings cell) to your own working
  folder before running. The `%pip` cell installs the packages; restart the kernel
  once after installing, then continue.
- **Run time:** folders `08` and `09` run the full nested cross-validation for every
  scenario and can take many hours; folder `07` is much faster in its default fixed
  Extra Trees mode. Finished scenarios and ablation configurations are saved as they
  complete and are skipped when you run the notebook again.
- **`QUICK_TEST = True`** (folders `08`, `09`) and `FAST = True` (folder `07`) are only
  for checking that the flow runs. Never report results from these modes.
- The stability / robustness criteria (`TOL`) in folders `08` and `09` are fixed before
  running and should not be changed after looking at the results.

### Google Colab cells

A few cells (file upload/download in folders `04` and `05`) call
`google.colab.files`. These are optional convenience hooks for running in
Google Colab — if you're running locally, replace them with a direct file
path (e.g. `CSV_PATH = "external_validation_peptides.csv"`) and skip the
`from google.colab import files` cell. Folder `06` is written for Colab (it uses the upload
dialog on purpose); folders `07`–`09` run locally in Jupyter.

## Methodology summary

**Features** (folder `03`): amino-acid composition (AAC), dipeptide
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

**External validation** (folder `05`): the packaged model is re-applied to
an independent peptide set with zero training overlap, evaluated with
bootstrap CIs, and broken down by sequence similarity to the training set
and by MHC class/organism subgroup — the standard way to check whether
performance holds up outside the training distribution.

**Supplementary analyses** (folders `06`–`09`): benchmarking against NetMHCpan/NetMHCIIpan
rank, PRIME 2.1 and the IEDB immunogenicity tool with cluster-bootstrap confidence
intervals; a feature-block ablation and ESM-2 redundancy analysis; and sensitivity
analyses of the similarity (80%), label-confidence (70%) and cluster-consensus (75%)
thresholds, including a leakage check that also catches short peptides lying inside
longer ones (`core9`).

## Reproducibility

All notebooks default to `random_seed = 42`. The notebook in folder `03` also pins
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
`VECLIB_MAXIMUM_THREADS`, and `NUMEXPR_NUM_THREADS` to `1` at the top of the
notebook to keep BLAS/OpenMP thread counts from interfering with
`n_jobs`-based parallelism in scikit-learn.

The notebooks in folders `07`–`09` use the same thread-limiting environment variables and
`random_seed = 42`, and require `scikit-learn<1.8` (see Setup). Their baseline
scenarios (`cv_0.80` in folder `08`, `base` in folder `09`) are meant to reproduce the
main pipeline's cross-validation AUC; a mismatch means a package version or a setting
has changed.

## Caveats

- The notebooks in folders `01`, `02` and `04` depend on live third-party APIs (IEDB, CEDAR,
  UniProt, IEDB's MHC prediction tools). Query behavior, field names, and
  rate limits are outside this repository's control and may change; the
  scripts include schema-discovery and fallback logic for this reason.
- MHC binding predictions in folder `01` restrict to a fixed 6-allele
  class I / 7-allele class II reference panel; results are specific to
  that panel.
- The notebook in folder `06` calls the IEDB Next-Generation Tools API and, at run time,
  clones PRIME 2.1 and MixMHCpred 3.x (Gfeller Lab). These tools are **free for
  academic, non-commercial use only**, are not part of this repository and are not
  covered by its MIT License. PRIME and the IEDB tool only score class I peptides, so
  for class II the comparison is against the NetMHCIIpan rank and the baselines.
- The rank used in the benchmark is also an input of the ensemble, and the training
  data of PRIME and the IEDB tool also come from IEDB, so the external peptides may
  appear in their training sets. The tools also predict different things (binding /
  CD8 immunogenicity versus T-cell assay labels).

## Citation

If you use this pipeline, please cite it — see `CITATION.cff` (fill in your
name/affiliation and repository URL before publishing).

## License

This repository includes a `LICENSE` file. The code is released under the MIT License. Third-party tools that the notebook in
folder `06` can download at run time (PRIME, MixMHCpred) are subject to their own
licenses.

## Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22301492.svg)](https://doi.org/10.5281/zenodo.22301492)
