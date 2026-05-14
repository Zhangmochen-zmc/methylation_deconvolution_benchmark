## Workflow

The execution is divided into three main steps. Please follow them in order. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.

Run `data_processing.R` to generate a methylation matrix by integrating the `raw_ref/` for each sample.

```bash
Rscript data_processing.R
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref/`    
**Output:** `ref_data.RData`   


### Step 2: Marker Selection

Run `ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix.

```bash
Rscript ref.R
```

**Input:** `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/), `ref_data.RData` from Step 1.      
**Output:** `test_ref_Signature.txt` (marker_ref/)


### Step 3: Deconvolution

Run `decon.R` to perform the core deconvolution. This process includes:

```bash
Rscript decon.R
```

**Input:** `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/), `test_ref_Signature.txt` from Step 2.     
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/WonyoungCho/MethylCIBERSORT.git](https://github.com/WonyoungCho/MethylCIBERSORT.git)

