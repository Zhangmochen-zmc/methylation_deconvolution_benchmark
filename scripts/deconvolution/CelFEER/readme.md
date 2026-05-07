## Workflow

The execution is divided into four main steps. Please follow them in order:

### Step 1: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in `ref_data/` and `test_data/`folder.

*  **Reference Matrix (`ref.txt`)**: A signature matrix where:
    *   **Rows**: Defined by chrom, chrom_start, and chrom_end.
    *   **Columns**: Each cell type consists of 5 space-separated integers representing the read counts for the same five methylation levels (0%, 25%, 50%, 75%, and 100%).
*  **Mixture Matrix (`test_raw.txt`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Defined by chrom, chrom_start, and chrom_end (match the Reference Matrix).
    *   **Columns**: Each sample consists of 5 tab-separated groups of space-separated integers, representing the number of reads at 0%, 25%, 50%, 75%, and 100% methylation levels.


### Step 2: Marker Selection (`ref.py`)

Run `ref.py` (markers.py in CelFEER) using data to extract cell-type-specific marker genes and generate the signature matrix.

```bash
python markers.py <input_file> <output_file> <num_values> <tissues> <depth_filter> <nan_filter> <extra_filter> <variant>
```

**Input:** `ref.txt`(ref_data)  
**Output:** `marker.txt`(marker_ref)

### Step 3: Integration

Integrate `marker.txt` and `test_raw.txt`, and sort them according to chromosome order.

A single input line may look as follows:

`chr1 <tab> 1 <tab> 500 <tab> 0 1 0 2 34 <tab> 12 8 0 0 0 <tab> chr1 <tab> 1 <tab> 500 <tab>  12 5 1 1 1 <tab> 0 1 1 5 41  <tab> 2 2 5 2 3`

This line contains two different individuals, and the reference data of three different cell types. 

### Step 4: Deconvolution (`decon.py`)

Run `decon.py` (celfeer.py in CelFEER) using methylation data and the signature matrix generated in the previous step.

```bash
python decon.py
```

**Input:** Integrated files from Step 3.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

### Notes
*   **Data processing**：[https://github.com/pi-zz-a/CelFEER.git](https://github.com/pi-zz-a/CelFEER.git)
*   **More Information**:  [https://github.com/pi-zz-a/CelFEER.git](https://github.com/pi-zz-a/CelFEER.git)
