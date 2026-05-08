## Workflow

The execution is divided into three main steps. Please follow them in order. Houseman's QC_QP handles 450k and 850k methylation arrays differently.(Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.

### Step 2: Marker Selection

Run `ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix. 

```bash
Rscript ref.R
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref`      
**Output:** `450k_reference_output_houseman.csv` (`marker_ref`)


### Step 3: Deconvolution

Run `450k_decon.R` to perform the core deconvolution. (Taking 450k as an example)

```bash
Rscript 450k_decon.R
```

**Input:** `450k_reference_output_houseman.csv` from Step 2, `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/).     
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://link.springer.com/article/10.1186/1471-2105-13-86](https://link.springer.com/article/10.1186/1471-2105-13-86)


