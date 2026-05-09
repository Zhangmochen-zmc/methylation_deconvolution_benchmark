# Methylation_deconvolution_benchmark

This repository contains the processed data, benchmarking pipelines, and analysis results used in the study:

> **Cross-platform benchmarking of DNA methylation deconvolution methods across methylation array and whole-genome bisulfite sequencing data**

The project provides a comprehensive benchmarking framework for evaluating DNA methylation deconvolution methods across both methylation array and whole-genome bisulfite sequencing (WGBS) platforms.

---

## Repository Structure

```text
methylation_deconvolution_benchmark/
│
├── data/                  
│   ├── reference_data/          # Reference data used for benchmarking
│   ├── real_data/               # Real methylation data used for benchmarking
│   ├── test_data/               # Simulated bulk methylation data for benchmarking
├── scripts/         
│   ├── array_process/           # Array data preprocessing scripts
│   ├── wgbs_process/            # WGBS data preprocessing scripts
│   ├── simulated_data_generate/ # Simulated bulk methylation data generation scripts
│   ├── deconvolution/           # Deconvolution methods execution scripts
│   ├── requirements/            # R and Python package dependencies
├── results/                
│   ├── real_data/               # Benchmark results on real data
│   ├── simulated_data/          # Benchmark results on simulated data
└── README.md
```

## Installation

```text
git clone https://github.com/Zhangmochen-zmc/methylation_deconvolution_benchmark.git
cd methylation_deconvolution_benchmark
```

---

## Software Environment

The analyses were performed under the following environments:

- **R**: v4.5.1
- **Python**: v3.10.19

Additional package dependencies are described in the corresponding folders.

---

## Included Methods

This repository contains execution pipelines for all deconvolution methods included in the manuscript, with method-specific configurations and running instructions provided within the corresponding directories, and all methods implemented as closely as possible to their original publications and software documentation.

---

## Usage

Please refer to the README or documentation within each subdirectory for detailed execution instructions.

---

## Citation

If you use this repository, please cite:

> *Cross-platform benchmarking of DNA methylation deconvolution methods across methylation array and whole-genome bisulfite sequencing data*

---

## Contact

For questions regarding this repository or the benchmarking framework, please contact the corresponding authors of the study.
