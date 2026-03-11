# Workflow

The execution is divided into three main steps. Please follow them in order:

## Step 0: Data Preparation

Prepare your input data in `.csv` format. It is recommended to place them in a `test_data/` folder.
Prepare the ref_raw file
Prepare the EMeth

## Step 1: Feature Alignment (`ref.py`)

This step aligns the features (e.g., Gene Symbols or Probe IDs) between the reference and the mixture data, keeping only the intersection of the two.

```bash
Rscript ref.R
```
## Step 2: Deconvolution (`decon.py`)

This step performs the core ARIC algorithm, including marker gene selection and weighted support vector regression.

```bash
Rscript decon.R
```
