## Workflow

The execution is divided into three main steps. Please follow them in order. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.

Run `data_processing.py` to generate a  methylation matrix by averaging the `original_data` for each cell type. 

```bash
python data_processing.py
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/original_data`    
**Output:** `ref_data.csv`   


### Step 2: Feature Alignment

Run `ref.py` to align the features between the reference and the mixture data. This script identifies the intersection of features present in both datasets to ensure compatibility.

```bash
python ref.py
```

**Input:** `ref_data.csv` from Step 1, `simulated real.csv'(methylation_deconvolution_benchmark/data/test_data/450k/)         
**Output:** `ref.csv`, `mix.csv` (marker_ref)

### Step 3: Deconvolution

Run `decon.py` to perform the core ARIC algorithm. This process includes:
*   **Marker Selection**: Identifying highly informative features for each cell type.
*   **Deconvolution**: Performing deconvolution using a weighted Support Vector Regression approach.

```bash
python decon.py
```

**Input:** Aligned files from Step 2.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/XWangLabTHU/ARIC.git](https://github.com/XWangLabTHU/ARIC.git)

