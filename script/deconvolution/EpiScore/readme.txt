# Workflow

The execution is divided into three main steps. Please follow them in order:

## Step 0: Data Preparation

Prepare your input data in `.csv` format. It is recommended to place them in a `test_data/` folder.
Prepare the ref_raw file
Prepare the EMeth

## Step 1: Feature Alignment (`ref.R`)

This step aligns the features (e.g., Gene Symbols or Probe IDs) between the reference and the mixture data, keeping only the intersection of the two.

```bash
Rscript ref.R
```
## Step 2: Deconvolution (`decon.R`)

This step performs the core ARIC algorithm, including marker gene selection and weighted support vector regression.

```bash
Rscript decon.R
```

We provide pre-processed data files in the test_data that can be used directly for testing the pipeline. Users can point the input path to these files to immediately evaluate the performance of EpiScore. Specifically, the filtered reference files used for deconvolution are located in the episcore_ref, ensuring you have a streamlined starting point for your analysis.
