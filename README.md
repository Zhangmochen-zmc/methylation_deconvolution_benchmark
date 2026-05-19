# Methylation Deconvolution Benchmark

This repository accompanies the manuscript:

> **Cross-platform benchmarking of DNA methylation deconvolution methods across methylation array and whole-genome bisulfite sequencing data**

It provides the complete data, scripts, and aggregated results for a systematic evaluation of 21 DNA methylation deconvolution methods across the Illumina 450K array, EPIC (850K) array, and whole-genome bisulfite sequencing (WGBS) platforms, using both simulated mixtures and real bulk blood methylation datasets spanning six immune cell types.

This repository documents the benchmarking pipeline as executed in the study and is intended as a transparency and reproducibility resource accompanying the manuscript.

---

## Repository Structure

```text
methylation_deconvolution_benchmark/
├── data/
│   ├── reference_data/          # Purified cell-type reference profiles (450K / EPIC / WGBS)
│   ├── real_data/               # Real bulk blood methylation datasets (450K / EPIC)
│   └── test_data/               # Simulated bulk methylation mixtures used for benchmarking
├── scripts/
│   ├── array_process/           # Preprocessing pipelines for 450K and EPIC array data
│   ├── wgbs_process/            # Preprocessing pipeline for WGBS data
│   ├── simulated_data_generate/ # Scripts for generating simulated mixture datasets
│   ├── deconvolution/           # Per-method execution scripts (21 methods)
│   └── requirements/            # Software dependency lists (R and Python)
└── results/
    ├── real_data/               # Deconvolution results on real datasets
    └── simulated_data/          # Deconvolution results on simulated datasets (by platform)
```

> **Note:** Large files (reference profiles, pre-trained model weights) are managed via [Git LFS](https://git-lfs.github.com/). Run `git lfs pull` after cloning to retrieve these files.

---

## Datasets

### Reference Data

Purified cell-type methylation profiles for six immune cell types were compiled from publicly available datasets across all three platforms. Raw per-sample profiles are provided in `data/reference_data/<platform>/raw_data/`.

| Cell Type           | Directory       |
|---------------------|-----------------|
| B cell              | `bcell`         |
| CD4+ T cell         | `cd4+Tcell`     |
| CD8+ T cell         | `cd8+Tcell`     |
| Monocyte            | `monocyte`      |
| Neutrophil          | `neutrophil`    |
| NK cell             | `nkcell`        |

### Real Benchmark Datasets

Five publicly available bulk blood methylation datasets with known cell-type composition were used to evaluate method performance under real-world conditions. Processed beta-value matrices are provided in `data/real_data/`.

| Dataset    | Platform | Accession                                                          |
|------------|----------|--------------------------------------------------------------------|
| GSE127824  | 450K     | [GSE127824](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127824) |
| GSE77797   | 450K     | [GSE77797](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77797)   |
| GSE58888   | 450K     | [GSE58888](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58888)   |
| GSE110554  | EPIC     | [GSE110554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554) |
| GSE112618  | EPIC     | [GSE112618](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112618) |

### Simulated Datasets

Synthetic bulk mixtures were generated from purified reference profiles under multiple biological and technical conditions. Mixing proportions were drawn from Dirichlet distributions and are provided alongside the resulting beta-value matrices. Scripts are in `scripts/simulated_data_generate/`; test matrices are in `data/test_data/`.

| Scenario         | Platform          | Description                                                        |
|------------------|-------------------|--------------------------------------------------------------------|
| `random`         | 450K, EPIC, WGBS  | Fully random Dirichlet-sampled cell-type proportions               |
| `simulated_real` | 450K, EPIC, WGBS  | Proportions mimicking typical peripheral blood composition         |
| `low`            | 450K, EPIC, WGBS  | All cell types present at uniformly low fractions                  |
| `less_onetype`   | 450K, EPIC, WGBS  | One cell type present at very low abundance                        |
| `more_onetype`   | 450K, EPIC, WGBS  | One cell type present at dominant abundance                        |
| `noisy`          | 450K, EPIC        | Gaussian noise added in M-value space to simulate technical variation |
| `sparsity`       | 450K, EPIC        | Sparse CpG feature coverage                                        |
| `depth`          | WGBS              | Variable sequencing depth across samples                           |
| `CpGcoverage`    | WGBS              | Variable CpG site coverage across samples                          |

---

## Benchmarked Methods

All 21 methods were implemented as closely as possible to their original publications and software documentation. Per-method execution scripts, reference construction pipelines, and any method-specific preprocessing steps are provided in `scripts/deconvolution/<method>/`.

| Method              | Platform          | Language    |
|---------------------|-------------------|-------------|
| ARIC                | 450K, EPIC, WGBS  | Python      |
| CelFEER             | WGBS              | Python      |
| CelFiE              | WGBS              | Python      |
| EDec                | 450K, EPIC, WGBS  | R / Python  |
| EMeth-noraml/laplace| 450K, EPIC, WGBS  | R           |
| EpiDISH-RPC/CP/CBS  | 450K, EPIC, WGBS  | R           |
| EpiSCORE            | 450K, EPIC, WGBS  | R / Python  |
| Houseman's QP       | 450K, EPIC, WGBS  | R           |
| MEnet               | 450K, EPIC, WGBS  | Python      |
| MeDeCom             | 450K, EPIC, WGBS  | R / Python  |
| MetDecode           | WGBS              | Python      |
| MethAtlas           | 450K, EPIC, WGBS  | Python      |
| MethylBERT          | WGBS              | Python      |
| MethylCIBERSORT     | 450K, EPIC, WGBS  | R           |
| PRMeth              | 450K, EPIC, WGBS  | R           |
| RefFreeEWAS         | 450K, EPIC, WGBS  | R           |
| Tsisal              | 450K, EPIC, WGBS  | R           |
| UXM                 | WGBS              | Shell       |

EpiDISH was evaluated with three internal algorithms (RPC, CBS, CP) and is reported as three entries in the results.
EMeth was evaluated with three internal algorithms (noraml, laplace) and is reported as three entries in the results.

---

## Results

Aggregated deconvolution results are provided in `results/` as CSV files. Each file contains predicted cell-type proportions and performance metrics (RMSE, PCC, R², JSD) for all applicable methods within a given dataset or simulation scenario.

### Real Data (`results/real_data/`)

| File            | Dataset    | Platform |
|-----------------|------------|----------|
| `GSE127824.csv` | GSE127824  | 450K     |
| `GSE77797.csv`  | GSE77797   | 450K     |
| `GSE58888.csv`  | GSE58888   | 450K     |
| `GSE110554.csv` | GSE110554  | EPIC     |
| `GSE112618.csv` | GSE112618  | EPIC     |

### Simulated Data (`results/simulated_data/`)

Results are organized into three subdirectories by platform (`450k/`, `epic/`, `wgbs/`). Each CSV corresponds to one simulation scenario and reports aggregated metrics across all 100 simulated samples.

---

## Software Environment

All analyses were performed under the following environments:

- **R**: v4.5.1  
- **Python**: v3.10.19

Package dependency lists are provided in `scripts/requirements/`:

- `python.txt` — Python environment (pip)
- `R.csv` — R package list with version information

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
