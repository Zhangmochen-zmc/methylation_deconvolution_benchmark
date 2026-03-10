# Workflow

The execution is divided into three main steps. Please follow them in order:

## Step 0: Data Preparation
Prepare your input data in `.csv` format. It is recommended to place them in a `data/` folder.

1.  **Reference Matrix (`ref.csv`)**: A signature matrix where rows are features (Genes/CpG sites) and columns are known cell types.
2.  **Mixture Matrix (`mix.csv`)**: The bulk data matrix where rows are features and columns are samples to be deconvolved.

**Project Structure:**
```text
.
├── data/
│   ├── ref.csv          # Reference signature matrix
│   └── mix.csv          # Bulk mixture matrix
├── ref.py               # Step 1: Preprocessing script
├── decon.py             # Step 2: Deconvolution script
└── README.md
```

---

## Step 1: Feature Alignment (`ref.py`)
This step aligns the features (e.g., Gene Symbols or Probe IDs) between the reference and the mixture data, keeping only the intersection of the two.

**Run command:**
```bash
python ref.py
```

*   **Function**: Reads `ref.csv` and `mix.csv`, identifies common features, and generates normalized intermediate files.
*   **Output**: Prepared matrices ready for the deconvolution algorithm.

---

## Step 2: Deconvolution (`decon.py`)
This step performs the core ARIC algorithm, including marker gene selection and weighted support vector regression.

**Run command:**
```bash
python decon.py
```

*   **Function**: 
    1.  Selects optimal marker genes using ARIC's adaptive strategy.
    2.  Calculates the proportions of each cell type for every sample.
    3.  Outputs the final proportion matrix.

---

# Output Description

Upon completion, the script will generate a result file (e.g., `cell_proportions.csv`):

*   **Rows**: Sample IDs from your mixture file.
*   **Columns**: Predicted proportions for each cell type.
*   **Constraint**: The sum of proportions for each sample will be scaled to 1 (100%).

| SampleID | CellType_A | CellType_B | CellType_C | ... |
| :--- | :--- | :--- | :--- | :--- |
| Sample_01 | 0.25 | 0.40 | 0.35 | ... |
| Sample_02 | 0.10 | 0.80 | 0.10 | ... |

