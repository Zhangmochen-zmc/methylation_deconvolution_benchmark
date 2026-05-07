## Workflow

The execution is divided into two main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/` folder.

*   **Reference Matrix (`ref_raw.csv`)**: 
    *   **Rows**: Features (Probe IDs).
    *   **Columns**: Known cell types.
*   **Mixture Matrix (`test.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.  

### Step 2: Deconvolution

Run `decon.py` to perform the core deconvolution. This process includes:

```bash
python decon.py
```

**Input:** `ref_raw.csv` and `test_data/` from Step 1.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/nloyfer/meth_atlas.git](https://github.com/nloyfer/meth_atlas.git)

