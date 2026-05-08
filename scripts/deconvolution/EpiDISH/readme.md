## Workflow

The execution is divided into three main steps. Please follow them in order.Since Episdh handles 450k and 850k methylation arrays differently, this workflow and the provided examples are based on the 450k array.

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.

### Step 2: Marker Selection

Run `450k_ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix.

```bash
Rscript 450k_ref.R
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref/`    
**Output:** `EpiDISH_450k_reference_result.csv` (marker_ref/)


### Step 3: Deconvolution

Run `450k_decon.R` to perform the core deconvolution. This process includes:

```bash
Rscript 450k_decon.R
```

**Input:** `EpiDISH_450k_reference_result.csv` from Step 2, `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/).       
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/sjczheng/EpiDISH.git](https://github.com/sjczheng/EpiDISH.git)

