## Workflow 1 (450K and EPIC)

The execution is divided into two main steps. Please follow them in order. (Taking 450k as an example)

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `methylation_deconvolution_benchmark/data/test_data/450k/` folder.


### Step 2: Deconvolution

Run `array_decon.py` to perform the core deconvolution. The pre-trained `.pkl` models are provided. This process includes:

```bash
python array_decon.py <model>
```

**Input:** `methylation_deconvolution_benchmark/data/reference_data/450k/raw_ref/`   
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

## Workflow 2 (WGBS)

The execution is divided into two main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `test_data/wgbs/` folder.


### Step 2: Deconvolution

Run `wgbs_decon.py` to perform the core deconvolution. The pre-trained `.pkl` models are provided. This process includes:

```bash
python wgbs_decon.py <model>
```

**Input:** `sample.bismark.cov.gz`(methylation_deconvolution_benchmark/data/test_data/wgbs/menet/).  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **More Information**: [https://github.com/yyoshiaki/MEnet.git](https://github.com/yyoshiaki/MEnet.git)
