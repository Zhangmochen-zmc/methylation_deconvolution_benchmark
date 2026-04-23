## Workflow

The execution is divided into four main steps. Please follow them in order:

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in  `ref_data/` and `test_data/`folder.

*  **Reference Matrix (`ref.txt`)**: A signature matrix where:
    *   **Rows**: Defined by chrom, chrom_start, and chrom_end.
    *   **Columns**: The model expects the methylation data in the form: # of methylated reads (METH) and # of total reads (DEPTH) for each cell type. For example it could look like:
```text
CHR   START END METH DEPTH
chr1	10	11	44.0	63.0
chr1	50	51	71.0	133.0
chr1	60	61	89.0	115.0
```
      
*  **Mixture Matrix (`test_raw.txt`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Defined by chrom, chrom_start, and chrom_end (match the Reference Matrix).
    *   **Columns**: The model expects the methylation data in the form: # of methylated reads (METH) and # of total reads (DEPTH) for each sample.


### Step 2: Marker Selection (`ref.py`)

Run `ref.py` (tim.py in CelFiE) using your processed scRNA-seq reference data (`ref.txt` in `ref_data/`) to extract cell-type-specific marker genes and generate the signature matrix.

```bash
python tim.py <input file> <output file> <num of tim/tissue> <num of tissues> <depth filter> <nan filter>
```

**Input:** `ref.txt`(ref_data)  
**Output:** `marker.txt`(celfie_ref)

### Step 3: Integration

Integrate `marker.txt` and `test_raw.txt`, and sort them according to chromosome order.

A single input file may look as follows:

```text
CHROM START END SAMPLE1_METH SAMPLE1_DEPTH CHROM START END CELL1_METH CELL1_DEPTH
chr1	10	11	44.0	63.0  chr1	10	11	25.0	29.0
chr1	50	51	71.0	133.0 chr1	50	51	85.0	99.0
chr1	60	61	89.0	115.0 chr1	60	61	92.0	117.0
```

This file contains one different individuals, and the reference data of one different cell types. 

### Step 4: Deconvolution (`decon.sh`)

Run `decon.sh` (celfie.py in CelFiE) using your prepared bulk RNA-seq data and the signature matrix generated in the previous step.

```bash
sh decon.sh
```

**Input:** Integrated files from Step 3.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **Data processing**：[https://github.com/christacaggiano/celfie.git](https://github.com/christacaggiano/celfie.git)
*   **More Information**:  [https://github.com/christacaggiano/celfie.git](https://github.com/christacaggiano/celfie.git)
