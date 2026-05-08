
## Workflow

The execution is divided into three main steps. Please follow them in order. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.


### Step 2: Marker Selection

Run `ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix.

```bash
Rscript ref.R
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref/`     
**Output:** `450k_reference_output_Tsisal.csv` (marker_ref/)


### Step 3: Deconvolution

Run `decon.R` to perform the core deconvolution. This process includes:

```bash
Rscript decon.R
```

**Input:** `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/), `450k_reference_output_Tsisal.csv` from Step 2.     
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/ziyili20/TOAST.git](https://github.com/ziyili20/TOAST.git)

