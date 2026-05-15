## Workflow

The execution is divided into two main steps. Please follow them in order. Users need to download `MetDecode` to the project.

### Step 1: Data Preparation

Before running the scripts, organize your tested data in `test_data` folder.

*   **Reference Matrix (`atlas.tsv`)**:
    *   **Rows**: Regions.
    *   **Columns**: Each cell type has two dedicated columns, namely the number of methylated CpG sites spanned in the marker region, and the total number of CpG sites (both methylation and unmethylated).  
    *   Notes: For the atlas file format, please refer to `marker_ref/atlas.tsv` for example. The first 3 columns contain respectively the chromosome, start position and end position of each marker region. The file must contain a header of the form: CHROM START END CELL1_METH CELL1_DEPTH CELL2_METH ...


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
