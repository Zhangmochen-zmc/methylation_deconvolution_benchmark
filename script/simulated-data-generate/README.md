# Simulation Data Generation

This directory contains the pipeline and scripts for generating synthetic datasets for both **DNA Methylation Array** and **WGBS (Whole Genome Bisulfite Sequencing)** platforms. The primary goal is to create realistic DNA methylation mixtures to benchmark deconvolution algorithms.

## Directory Structure

### 1. `dirichlet_proportions/` (Proportion Generation)
This module handles the generation of cell-type mixing proportions using Dirichlet distributions.
* `random.py`: Generates fully random proportions for general testing.
* `simulated_real.py`: Generates proportions that mimic real-world biological scenarios.
* `ultra_low.py`: Specifically generates proportions with **ultra-low abundance** cell types to test detection sensitivity.

### 2. `array_data/` (Infinium Array Simulation)
Scripts for simulating Illumina Methylation Array data.
* `mixdata.py`: Combines purified cell-type profiles based on generated proportions to generate the mixed array matrix.
* `noisy.py`: Adds various levels of technical or biological noise to the synthetic array samples.

### 3. `wgbs_data/` (WGBS Simulation)
Tools for generating simulated sequencing data in PAT and BED formats.
* `mix_pat.py`: Mixes purified PAT files to create simulation samples at various sequencing depths.
* `bed_generate.sh`: Automates the workflow from mixing to BED generation.

## Workflow Summary

1.  **Define Proportions**: Use scripts in `dirichlet_proportions/` to create a CSV file of target cell fractions.
2.  **Generate Mixtures**:
    * **For Array**: Run `array_data/mixdata.py` to produce beta-value matrices.
    * **For WGBS**: Use `wgbs_data/mixdata.py` to generate depth-specific PAT files.
3.  **Process/Refine**: Add noise to array data using `noisy.py` or convert WGBS files to BED format using `bed_generate.sh`.

## Requirements
* Python 3.10
* Pandas
* wgbstools
* numpy/scipy (for Dirichlet distribution and noise injection)
