## Workflow

The execution is divided into three main steps. Please follow them in order. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.

Run `data_processing.py` to generate a  methylation matrix by averaging the `raw_ref` for each cell type. 

```bash
python data_processing.py
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref/`    
**Output:** `ref_data.csv`   

### Step 2: Marker Selection and data integration

Run `ref.R` using reference data to extract cell type specific marker genes and generate the signature matrix. The processed reference matrix and the list of all project datasets are finally saved into a single `.RData` file.

```bash
Rscript ref.R
```

**Input:** `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/),`ref_data.csv` from Step 1.     
**Output:** `episcore.RData` (marker_ref/)


### Step 3: Deconvolution

Run `decon.R` to perform the core deconvolution. This process includes:

```bash
Rscript decon.R
```

**Input:** `episcore.RData` from Step 2.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/aet21/EpiSCORE.git](https://github.com/aet21/EpiSCORE.git)


