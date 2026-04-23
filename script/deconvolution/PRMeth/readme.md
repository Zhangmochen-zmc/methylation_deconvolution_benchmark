## Workflow

The execution is divided into three main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/` folder.

*   **Reference Matrix (`.txt`)**: 
    *   **Rows**: Features (Probe IDs).
    *   **Columns**: Known cell types.
    *   Organization: Data for different cell types are organized as individual subfolders within the reference directory.
*   **Mixture Matrix (`test.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.  

### Step 2: Marker Selection

Run `ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix.

```bash
Rscript ref.R
```

**Input:** `.txt`(`ref_data`)    
**Output:** `reference_output_PRMeth.csv` (`prmeth_ref`)


### Step 3: Deconvolution

Run `decon.R` to perform the core deconvolution. This process includes:

```bash
Rscript decon.R
```

**Input:** `reference_output_PRMeth.csv` and test_data/` from Step 2.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/hedingqin/PRMeth.git](https://github.com/hedingqin/PRMeth.git)
*   **Environment**: Users need to download ‘PEMeth’ before proceeding with the workflow.

