## Workflow

The execution is divided into three main steps. Please follow them in order.

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/`folder.

*   **Reference Matrix (`atlas.tsv`)**:
    *   **Rows**: Regions.
    *   **Columns**: Each cell type has two dedicated columns, namely the number of methylated CpG sites spanned in the marker region, and the total number of CpG sites (both methylation and unmethylated).  
    *   Notes: For the atlas file format, please refer to `ref_data/atlas.tsv` for example. The first 3 columns contain respectively the chromosome, start position and end position of each marker region. The file must contain a header of the form: CHROM START END CELL1_METH CELL1_DEPTH CELL2_METH ...
*   **Mixture Matrix (`test.tsv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Regions (must use the same naming convention as the reference matrix).
    *   **Columns**: Each sample has two dedicated columns, namely the number of methylated CpG sites spanned in the marker region, and the total number of CpG sites (both methylation and unmethylated).  


### Step 2: Deconvolution

Run `decon.py` to perform the core deconvolution. This process includes:

```bash
python decon.py
```

**Input:** `atlas.tsv` and `test_data/` from Step 1.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **Atlas**: [https://github.com/JorisVermeeschLab/MetDecode.git](https://github.com/JorisVermeeschLab/MetDecode.git)
*   **More Information**: [https://github.com/JorisVermeeschLab/MetDecode.git](https://github.com/JorisVermeeschLab/MetDecode.git)
