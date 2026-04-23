## Workflow

The execution is divided into three main steps. Please follow them in order.


- Sample data for Step 4: `*.pat.gz` files
- Metadata:
  - Step 2: `uxm_wgbs_meta.csv`
  - Step 3: `wgbstools_pat_meta.csv`
 
### Step 1: Genome Segmentation

Segment the genome into blocks based on CpG density:

```bash
wgbstools segment \
    --betas *beta \
    --min_cpg 5 \
    --max_bp 2000 \
    -o ../blocks_blood.bed
```

*  **betas (`*.beta`)**: Reference data used for genome segmentation.

### Step 2: Marker Selection

Identify informative methylation markers:

```bash
wgbstools find_markers \
    --blocks_path ../blocks_blood.bed \
    --groups_file ../uxm_wgbs_meta.csv \
    --betas *beta \
    --sort_by delta_maxmin \
    -o ../markers/ \
    --delta_quants .3 \
    --pval 1
```

Merge all marker BED files:

```bash
awk 'FNR==1 && NR!=1 {next} 1' *.bed > all_markers_merged.bed
```

*  **blocks_path (blocks_blood.bed)**: Segmented genomic regions.
*  **betas (`*.beta`)**: Reference data used for marker selection.
*  **groups_file (uxm_wgbs_meta.csv)**: Metadata of the reference data.

### Step 3: Atlas Construction

Build a reference atlas using selected markers:

```bash
uxm build \
    -m ../25markers/all_25markers_merged.bed \
    --pats *.pat.gz \
    -o ../25markers/atlas/25marker_atlas.csv \
    --groups ../wgbstools_pat_meta.csv \
    --threads 1 \
    --rlen 4
```

*  **pats (`*.pat.gz`)**: Reference data used for atlas construction.
*  **groups (wgbstools_pat_meta.csv)**: Metadata of the reference data.
*  **m (all_25markers_merged.bed)**: Selected specific markers.

### Step 4: Deconvolution

Run batch deconvolution using the provided script:

```bash
bash UXM_deconv.sh <input_dir> <output_dir>
```

### Notes

*   **More Information**: [https://github.com/nloyfer/wgbs_tools.git](https://github.com/nloyfer/wgbs_tools.git)
*   **More Information**: [https://github.com/nloyfer/UXM_deconv.git](https://github.com/nloyfer/UXM_deconv.git)

