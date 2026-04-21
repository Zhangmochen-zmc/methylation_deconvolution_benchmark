Step 0: Data Preparation
Before running the scripts, organize your input data. It is recommended to place all files in a test_data/ folder.

Reference Matrix (ref_raw.csv): A signature matrix where:

Rows: Features (Gene Symbols, Probe IDs, or CpG sites).

Columns: Known cell types.

Mixture Matrix (mix.csv): The bulk data matrix to be deconvolved where:

Rows: Features (must use the same naming convention as the reference).

Columns: Samples.

Step 1: Feature Alignment
Run ref.py to align the features between the reference and the mixture data. This script identifies the intersection of features present in both datasets to ensure compatibility.

code
Bash
python ref.py
Input: ref_raw.csv, mix.csv
Output: Aligned matrices (e.g., ref.csv) containing only shared features.

Step 2: Deconvolution
Run decon.py to perform the core ARIC algorithm. This process includes:

Marker Selection: Identifying highly informative features for each cell type.

Weighted SVR: Performing deconvolution using a weighted Support Vector Regression approach.

code
Bash
python decon.py
Input: Aligned files from Step 1.
Output: Predicted cell type proportions for each sample in the mixture matrix.

File Requirements
File Name	Format	Description
ref_raw.csv	CSV	Raw signature matrix (Features x Cell Types)
mix.csv	CSV	Bulk mixture data (Features x Samples)
