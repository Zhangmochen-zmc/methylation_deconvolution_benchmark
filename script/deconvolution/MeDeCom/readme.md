## Workflow

The execution is divided into three main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/`folder.

*   **Reference Matrix (`refdata.txt`)**: A signature matrix where:
    *   **Rows**: Features (Probe IDs).
    *   **Columns**: Known cell types.
*   **Reference Metadata Matrix (`refmeta.csv`)**:
    *   **Rows**:Features.
    *   **Columns**:Cell types information.
*   **Mixture Matrix (`test.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.

### Step 2: Deconvolution

Run `decon.R` to perform the deconvolution. This process includes:

```bash
Rscript decon.R
```

**Input:** `test.csv` from Step 1.  
**Output:** ``, ``.

### Step 3: Deconvolution Result Processing

```bash
Rscript decon_process.R
```

**Input:** ``, ``. from Step 2.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/CompEpigen/MeDeCom.git](https://github.com/CompEpigen/MeDeCom.git)
