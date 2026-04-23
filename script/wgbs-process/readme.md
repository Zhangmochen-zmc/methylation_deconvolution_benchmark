# WGBS Data Processing for BED Generation

This project provides a complete **pipeline** for processing Whole-Genome Bisulfite Sequencing (WGBS) data. It automatically detects input file types (`.bam` or `.pat.gz`) and converts them into high-quality **BED files** containing methylation information.

---

## Project Overview

This pipeline is designed to streamline WGBS data processing by integrating multiple steps into a single workflow. It supports both raw alignment files and processed methylation pattern files, making it flexible for different stages of analysis.

The output files include:
- PAT files: Read-level methylation data
- BED files: Methylation level at each CpG genomic coordinate

---

## Workflow Summary

| Step | Purpose | Key Tool |
|------|--------|----------|
| **Input Detection** | Automatically detects whether input files are BAM or PAT | bash |
| **Sorting (BAM only)** | Sorts BAM files for downstream processing | `samtools sort` |
| **Conversion** | Converts BAM to PAT format | `wgbstools bam2pat` |
| **Indexing** | Generates `.beta` files from PAT | `wgbstools index` |
| **BED Generation** | Converts beta values to BED format and merges statistics | `wgbstools beta2bed`, `awk` |

---

## Requirements

- `wgbstools`: [https://github.com/nloyfer/wgbs_tools.git](https://github.com/nloyfer/wgbs_tools.git)
- `samtools` : [https://github.com/samtools/samtools.git](https://github.com/samtools/samtools.git)

---

## Usage

```bash
bash auto_detect_wgbs.sh <input_folder> [threads]
