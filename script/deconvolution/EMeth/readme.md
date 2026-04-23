## Workflow

The execution is divided into three main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/`folder.

*   **Reference Matrix (`.txt`)**: 
    *   **Rows**: Features (Gene Symbols, Probe IDs, or CpG sites).
    *   **Columns**: Known cell types.
    *   Organization: Data for different cell types are organized as individual subfolders within the reference directory.
*   **Mixture Matrix (`test.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.

### Step 2: Rreference

Run `ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix.

```bash
Rscript ref.R
```

**Input:** `.txt`(`ref_data`)    
**Output:** `edec_stage0_markers.rds` (`edec_ref`)


### Step 3: Deconvolution

Run `decon.R` to perform the core ARIC algorithm. This process includes:

```bash
Rscript decon.R
```

**Input:** `edec_stage0_markers.rds` from Step 2.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/BRL-BCM/EDec.git](https://github.com/BRL-BCM/EDec.git)

