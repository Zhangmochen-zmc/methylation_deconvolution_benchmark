## Workflow

The execution is divided into three main steps. Please follow them in order. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.

Run `data_processing.py` to generate a methylation matrix by integrating the `raw_ref` for each sample.

```bash
python data_processing.py
```
**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref/`    
**Output:** `ref_data.txt` 

*   **Reference Metadata Matrix (`refmeta.csv`)**:
    *   **Rows**:Features.
    *   **Columns**:Cell types information.


### Step 2: Deconvolution

Run `decon.R` to perform the deconvolution. This process includes:

```bash
Rscript decon.R
```

**Input:** `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/)
**Output:** `medecom_rds/simulated_real.rds`, `medecom_pdf/simulated_real.pdf`.

### Step 3: Deconvolution Result Processing

```bash
Rscript decon_process.R
```

**Input:** `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/), `ref_data.txt` from Step 1, `refmeta.csv`, `simulated_real.rds`(medecom_rds/), `simulated_real.pdf`(medecom_pdf/) from Step 2.
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/CompEpigen/MeDeCom.git](https://github.com/CompEpigen/MeDeCom.git)
