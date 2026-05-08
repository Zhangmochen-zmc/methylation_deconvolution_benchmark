## Workflow

The execution is divided into three main steps. Please follow them in order. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data and metadata in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.

Run `data_processing.py` to generate a methylation matrix by integrating the `raw_ref/` for each sample.

```bash
python data_processing.py
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref/`    
**Output:** `ref_data.txt`   

*   **Reference Metadata Matrix (`refmeta.csv`)**:
    *   **Rows**:Features.
    *   **Columns**:Cell types information.

### Step 2: Marker Selection

Run `ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix.

```bash
Rscript ref.R
```

**Input:** `ref_data.txt` from Step 1, `refmeta.csv`(`methylation_deconvolution_benchmark/data/reference_data/450k/)       
**Output:** `edec_stage0_markers.rds` (marker_ref/)


### Step 3: Deconvolution

Run `decon.R` to perform the deconvolution. This process includes:

```bash
Rscript decon.R
```

**Input:** `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/), `edec_stage0_markers.rds` from Step 2.
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/BRL-BCM/EDec.git](https://github.com/BRL-BCM/EDec.git)

