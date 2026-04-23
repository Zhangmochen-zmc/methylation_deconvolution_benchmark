## Workflow

The execution is divided into three main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/` folder.

*   **Reference Matrix (`ref_raw.RData`)**: 
    *   **Rows**: Features (Probe IDs).
    *   **Columns**: Known cell types.
*   **Mixture Matrix (`test.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.  

### Step 2: Marker Selection

Run `ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix.

```bash
Rscript ref.R
```

**Input:** `.txt`(`ref_data`) and `test_data/`
**Output:** `test_data_ref.csv` (`methylcibersort_ref`)


### Step 3: Deconvolution

Run `decon.R` to perform the core deconvolution. This process includes:

```bash
Rscript decon.R
```

**Input:** `test_data_ref.csv` and `test_data/` from Step 2.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/WonyoungCho/MethylCIBERSORT.git](https://github.com/WonyoungCho/MethylCIBERSORT.git)

