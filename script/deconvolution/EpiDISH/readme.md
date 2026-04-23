## Workflow

The execution is divided into three main steps. Please follow them in order.Since Episdh handles 450k and 850k methylation arrays differently, this workflow and the provided examples are based on the 450k array.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/`folder.

*   **Reference Matrix (`.txt`)**: 
    *   **Rows**: Features (Probe IDs).
    *   **Columns**: Known cell types.
    *   Organization: Data for different cell types are organized as individual subfolders within the reference directory.
*   **Mixture Matrix (`test.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.  

### Step 2: Marker Selection

Run `450k_ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix.

```bash
Rscript 450k_ref.R
```

**Input:** `.txt`(`ref_data`)    
**Output:** `EpiDISH_450k_reference_result.csv` (`epidish_ref`)


### Step 3: Deconvolution

Run `450k_decon.R` to perform the core deconvolution. This process includes:

```bash
Rscript 450k_decon.R
```

**Input:** `EpiDISH_450k_reference_result.csv` from Step 2.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/sjczheng/EpiDISH.git](https://github.com/sjczheng/EpiDISH.git)

