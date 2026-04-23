## Workflow

The execution is divided into three main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `test_data/` folder.

*   **Mixture Matrix (`test.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference matrix).
    *   **Columns**: Samples.  

### Step 2: Deconvolution

Run `array_decon.py`/`wgbs_decon.py` to perform the core deconvolution. The pre-trained `.pkl` models are provided. This process includes:

```bash
python array_decon.py <model>
```

**Input:** `test_data` from Step 1.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/yyoshiaki/MEnet.git](https://github.com/yyoshiaki/MEnet.git)
