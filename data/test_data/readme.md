## test_data

This directory contains the simulated 450K, EPIC, and WGBS datasets used as deconvolution test inputs for benchmarking DNA methylation deconvolution methods.

For all array-based deconvolution methods, a unified deconvolution input was used and is provided in this directory to ensure consistency and comparability across methods.

Because sequencing-based deconvolution methods require different WGBS input formats depending on the specific implementation, the corresponding deconvolution test inputs for each method are provided within the respective method directories under:

```text
/scripts/deconvolution/<method_name>/test_data
```

Each method-specific folder contains the required deconvolution test files and detailed instructions for reproducing the corresponding analyses.

