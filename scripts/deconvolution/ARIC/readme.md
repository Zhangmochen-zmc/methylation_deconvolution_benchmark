## Workflow

The execution is divided into three main steps. Please follow them in order. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in  `methylation_deconvolution_benchmark/data/reference_data/450k/` and `methylation_deconvolution_benchmark/data/test_data/450k/`folder.

*  **Reference Matrix (`Ravg.csv`)**: A signature matrix where:
    *   **Rows**: Features (Probe IDs).
    *   **Columns**: Known cell types.
*  **Mixture Matrix (`simulated_real.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.

### Step 2: Feature Alignment

Run `ref.py` to align the features between the reference and the mixture data. This script identifies the intersection of features present in both datasets to ensure compatibility.

```bash
python ref.py
```

**Input:** `Ravg.csv` 
**Output:** `ref.csv`, `mix.csv`

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

