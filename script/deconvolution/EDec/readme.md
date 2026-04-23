## Workflow

The execution is divided into three main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/`folder.

*   **Reference Matrix (`refdata.txt`)**: A signature matrix where:
    *   **Rows**: Features (Gene Symbols, Probe IDs, or CpG sites).
    *   **Columns**: Known cell types.
*   **Reference Metadata Matrix (`refmeta.csv`)**:
    *   **Rows**:Features.
    *   **Columns**:Cell types information.
*   **Mixture Matrix (`test.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.

### Step 2: Feature Alignment

Run `ref.R` to align the features between the reference and the mixture data. This script identifies the intersection of features present in both datasets to ensure compatibility.

```bash
Rscript ref.R
```

**Input:** `ref_data.txt`, `refmeta.csv`(`ref_data`)
**Output:** `` (`edec_ref`).


### Step 3: Deconvolution

Run `decon.py` to perform the core ARIC algorithm. This process includes:
*   **Marker Selection**: Identifying highly informative features for each cell type.
*   **Deconvolution**: Performing deconvolution using a weighted Support Vector Regression approach.

```bash
python decon.py
```

**Input:** Aligned files from Step 2.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/XWangLabTHU/ARIC.git](https://github.com/XWangLabTHU/ARIC.git)

