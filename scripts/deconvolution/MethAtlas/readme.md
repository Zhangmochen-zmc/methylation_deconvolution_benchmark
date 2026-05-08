## Workflow

Episdh handles 450k and 850k methylation arrays differently. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.

Run `data_processing.py` to generate a  methylation matrix by averaging the `raw_ref` for each cell type. 

```bash
python data_processing.py
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref/`    
**Output:** `ref_data.csv`   


### Step 2: Deconvolution

Run `decon.py` to perform the core deconvolution. This process includes:

```bash
python decon.py
```

**Input:** `simulated_real.csv`(methylation_deconvolution_benchmark/data/test_data/450k/), `ref_data.csv` from Step 1.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/nloyfer/meth_atlas.git](https://github.com/nloyfer/meth_atlas.git)

