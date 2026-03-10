# Workflow

The execution is divided into three main steps. Please follow them in order:

## Step 0: Data Preparation

Prepare your input data in `.csv` format. It is recommended to place them in a `test_data/` folder.
Prepare the ref_raw file

## Step 1: Feature Alignment (`ref.py`)

This step aligns the features (e.g., Gene Symbols or Probe IDs) between the reference and the mixture data, keeping only the intersection of the two.

1.  **Reference Matrix (`ref.csv`)**: A signature matrix where rows are features (Genes/CpG sites) and columns are known cell types.
2.  **Mixture Matrix (`mix.csv`)**: The bulk data matrix where rows are features and columns are samples to be deconvolved.

```bash
python ref.py
```
## Step 2: Deconvolution (`decon.py`)

This step performs the core ARIC algorithm, including marker gene selection and weighted support vector regression.

```bash
python decon.py
```
