## reference_data

This directory contains the 450K, EPIC, WGBS reference inputs used for benchmarking the DNA methylation deconvolution methods.

For all array-based deconvolution methods, a unified reference input was used and is provided in this directory to ensure consistency and comparability across methods.

Because sequencing-based deconvolution methods require different WGBS input formats and reference contents depending on the specific implementation, the corresponding reference inputs for each method are provided within the respective method directories under:

```text
/scripts/deconvolution/<method_name>/ref_data
```

Each method-specific folder contains the required reference files and detailed instructions for reproducing the corresponding analyses.

