# Methylation Deconvolution Benchmark

This repository contains the data, scripts, and results accompanying the manuscript:

> **Cross-platform benchmarking of DNA methylation deconvolution methods across methylation array and whole-genome bisulfite sequencing data**

---

## Overview

DNA methylation deconvolution methods estimate cell-type composition from bulk methylation profiles, yet systematic comparisons across sequencing platforms remain limited. This study provides a comprehensive benchmark of 18 deconvolution methods evaluated across three platforms — Illumina 450k array, EPIC (850k) array, and whole-genome bisulfite sequencing (WGBS) — using both simulated and real blood methylation data spanning six immune cell types.

This repository documents the complete benchmarking pipeline as executed in the study, including reference data construction, simulated data generation, per-method execution scripts, and aggregated results. It is intended as a transparency and reproducibility resource accompanying the manuscript, not as a general-purpose deconvolution toolkit.

---

## Repository Structure

```text
methylation_deconvolution_benchmark/
├── data/
│   ├── reference_data/          # Purified cell-type reference profiles (450k / EPIC / WGBS)
│   ├── real_data/               # Real bulk blood methylation datasets
│   └── test_data/               # Simulated bulk methylation mixtures used for benchmarking
├── scripts/
│   ├── array_process/           # Preprocessing pipelines for 450k and EPIC array data
│   ├── wgbs_process/            # Preprocessing pipeline for WGBS data
│   ├── simulated_data_generate/ # Scripts for generating simulated mixture datasets
│   ├── deconvolution/           # Per-method execution scripts (18 methods)
│   └── requirements/            # R and Python package dependencies
└── results/
    ├── real_data/               # Deconvolution results on real datasets
    └── simulated_data/          # Deconvolution results on simulated datasets (by platform)
```

> **Note:** Some large files (reference profiles, pre-trained model weights) are managed via [Git LFS](https://git-lfs.github.com/). Run `git lfs pull` after cloning to retrieve these files.

---

## Datasets

### Reference Data

Purified cell-type methylation profiles were compiled from publicly available datasets for six immune cell types across all three platforms. Raw per-sample profiles are provided in `data/reference_data/<platform>/raw_data/`.

| Cell Type     | Directory name |
|---------------|----------------|
| B cell        | `bcell`        |
| CD4+ T cell   | `cd4+Tcell`    |
| CD8+ T cell   | `cd8+Tcell`    |
| Monocyte      | `monocyte`     |
| Neutrophil    | `neutrophil`   |
| NK cell       | `nkcell`       |

### Real Benchmark Datasets

Five real bulk blood methylation datasets with known cell-type composition were used to evaluate deconvolution performance under real-world conditions. Processed beta-value matrices are provided in `data/real_data/`.

| Dataset   | Platform | Location                      |
|-----------|----------|-------------------------------|
| GSE127824 | 450k     | `data/real_data/450k/`        |
| GSE77797  | 450k     | `data/real_data/450k/`        |
| GSE58888  | 450k     | —                             |
| GSE110554 | EPIC     | `data/real_data/epic/`        |
| GSE112618 | EPIC     | `data/real_data/epic/`        |

### Simulated Datasets

Synthetic bulk mixtures were generated from purified reference profiles under multiple biological and technical conditions to systematically assess method robustness. Scripts for simulation are provided in `scripts/simulated_data_generate/`. Simulated test matrices are stored in `data/test_data/`.

| Scenario         | Platform          | Description                                              |
|------------------|-------------------|----------------------------------------------------------|
| `random`         | 450k, EPIC, WGBS  | Fully random Dirichlet-sampled cell-type proportions     |
| `simulated_real` | 450k, EPIC, WGBS  | Proportions mimicking real blood composition             |
| `low`            | 450k, EPIC, WGBS  | Uniformly low cell-type fractions                        |
| `less_onetype`   | 450k, EPIC, WGBS  | One cell type present at very low abundance              |
| `more_onetype`   | 450k, EPIC, WGBS  | One cell type at dominant abundance                      |
| `noisy`          | 450k, EPIC        | Technical noise added to beta values                     |
| `sparsity`       | 450k, EPIC        | Sparse CpG feature coverage                             |
| `depth`          | WGBS              | Variable sequencing depth                                |
| `CpGcoverage`    | WGBS              | Variable CpG site coverage across samples                |

---

## Benchmarked Methods

All 18 methods were implemented as closely as possible to their original publications and software documentation to ensure a fair comparison. Each method directory under `scripts/deconvolution/` contains the execution scripts and method-specific reference inputs used in this study. Detailed workflow descriptions are provided in the `readme.md` within each method directory.

| Method            | Platform          | Language     | Reference                                                                 | Directory                                                             |
|-------------------|-------------------|--------------|---------------------------------------------------------------------------|-----------------------------------------------------------------------|
| ARIC              | 450k, EPIC, WGBS  | Python       | [GitHub](https://github.com/XWangLabTHU/ARIC)                            | [`scripts/deconvolution/ARIC`](scripts/deconvolution/ARIC)           |
| CelFEER           | WGBS              | Python       | [GitHub](https://github.com/pi-zz-a/CelFEER)                             | [`scripts/deconvolution/CelFEER`](scripts/deconvolution/CelFEER)     |
| CelFiE            | WGBS              | Python       | [GitHub](https://github.com/christacaggiano/celfie)                       | [`scripts/deconvolution/CelFiE`](scripts/deconvolution/CelFiE)       |
| EDec              | 450k, EPIC, WGBS  | R / Python   | [GitHub](https://github.com/BRL-BCM/EDec)                                | [`scripts/deconvolution/EDec`](scripts/deconvolution/EDec)           |
| EMeth             | 450k, EPIC, WGBS  | R            | [GitHub](https://github.com/Sun-lab/dMeth)                               | [`scripts/deconvolution/EMeth`](scripts/deconvolution/EMeth)         |
| EpiDISH           | 450k, EPIC, WGBS  | R            | [GitHub](https://github.com/sjczheng/EpiDISH)                            | [`scripts/deconvolution/EpiDISH`](scripts/deconvolution/EpiDISH)     |
| EpiSCORE          | 450k, EPIC, WGBS  | R / Python   | [GitHub](https://github.com/aet21/EpiSCORE)                              | [`scripts/deconvolution/EpiSCORE`](scripts/deconvolution/EpiSCORE)   |
| Houseman's QP     | 450k, EPIC, WGBS  | R            | [Houseman et al. 2012](https://doi.org/10.1186/1471-2105-13-86)           | [`scripts/deconvolution/Houseman's QC_QP`](scripts/deconvolution/Houseman's%20QC_QP) |
| MEnet             | 450k, EPIC, WGBS  | Python       | [GitHub](https://github.com/yyoshiaki/MEnet)                             | [`scripts/deconvolution/MEnet`](scripts/deconvolution/MEnet)         |
| MeDeCom           | 450k, EPIC, WGBS  | R / Python   | [GitHub](https://github.com/CompEpigen/MeDeCom)                          | [`scripts/deconvolution/MeDeCom`](scripts/deconvolution/MeDeCom)     |
| MetDecode         | WGBS              | Python       | [GitHub](https://github.com/JorisVermeeschLab/MetDecode)                 | [`scripts/deconvolution/MetDecode`](scripts/deconvolution/MetDecode) |
| MethAtlas         | 450k, EPIC, WGBS  | Python       | [GitHub](https://github.com/nloyfer/meth_atlas)                          | [`scripts/deconvolution/MethAtlas`](scripts/deconvolution/MethAtlas) |
| MethylBERT        | WGBS              | Python       | [GitHub](https://github.com/CompEpigen/methylbert)                       | [`scripts/deconvolution/MethylBERT`](scripts/deconvolution/MethylBERT) |
| MethylCIBERSORT   | 450k, EPIC, WGBS  | R            | [GitHub](https://github.com/WonyoungCho/MethylCIBERSORT)                 | [`scripts/deconvolution/MethylCIBERSORT`](scripts/deconvolution/MethylCIBERSORT) |
| PRMeth            | 450k, EPIC, WGBS  | R            | [GitHub](https://github.com/hedingqin/PRMeth)                            | [`scripts/deconvolution/PRMeth`](scripts/deconvolution/PRMeth)       |
| RefFreeEWAS       | 450k, EPIC, WGBS  | R            | [GitHub](https://github.com/eahouseman/RefFreeEWAS)                      | [`scripts/deconvolution/RefFreeEWAS`](scripts/deconvolution/RefFreeEWAS) |
| Tsisal            | 450k, EPIC, WGBS  | R            | [GitHub](https://github.com/ziyili20/TOAST)                              | [`scripts/deconvolution/Tsisal`](scripts/deconvolution/Tsisal)       |
| UXM               | WGBS              | Shell        | [GitHub](https://github.com/nloyfer/UXM_deconv)                          | [`scripts/deconvolution/UXM`](scripts/deconvolution/UXM)             |

---

## Results

Aggregated deconvolution results are provided in `results/` as CSV files. Each file contains the predicted cell-type proportions for all evaluated methods and samples within a given dataset or simulation scenario, combined into a single table for downstream analysis.

### Real Data (`results/real_data/`)

| File            | Dataset   | Platform |
|-----------------|-----------|----------|
| `GSE127824.csv` | GSE127824 | 450k     |
| `GSE77797.csv`  | GSE77797  | 450k     |
| `GSE58888.csv`  | GSE58888  | 450k     |
| `GSE110554.csv` | GSE110554 | EPIC     |
| `GSE112618.csv` | GSE112618 | EPIC     |

### Simulated Data (`results/simulated_data/`)

Results are organized by platform into three subdirectories: `450k/`, `epic/`, and `wgbs/`. Each CSV file corresponds to one simulation scenario and contains predictions from all applicable methods.

---

## Software Environment

All analyses were performed under the following environments:

- **R**: v4.5.1
- **Python**: v3.10.19

Package dependencies are listed in `scripts/requirements/`:

- `python.txt` — Python environment specification (pip)
- `R.csv` — R package list with version information

---

## Citation

If you use the data, scripts, or results from this repository, please cite:

> *Cross-platform benchmarking of DNA methylation deconvolution methods across methylation array and whole-genome bisulfite sequencing data*  
> [Authors] · [Journal] · [Year] · DOI: [DOI]

---

## Contact

For questions regarding this repository or the benchmarking framework, please contact the corresponding authors of the manuscript.
