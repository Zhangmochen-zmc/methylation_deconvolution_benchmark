## Workflow

The execution is divided into three main steps. Please follow them in order.

### Step 0: Data Preparation

Before running the scripts, organize your input data. It is recommended to place all files in a `test_data/` folder.

1.  **Reference Matrix (`ref_raw.csv`)**: A signature matrix where:
    *   **Rows**: Features (Gene Symbols, Probe IDs, or CpG sites).
    *   **Columns**: Known cell types.
2.  **Mixture Matrix (`mix.csv`)**: The bulk data matrix to be deconvolved where:
    *   **Rows**: Features (must use the same naming convention as the reference).
    *   **Columns**: Samples.

### Step 1: Feature Alignment

Run `ref.py` to align the features between the reference and the mixture data. This script identifies the intersection of features present in both datasets to ensure compatibility.

```bash
python ref.py
```

**Input:** `ref_raw.csv`, `mix.csv`  
**Output:** Aligned matrices (e.g., `ref.csv`) containing only shared features.

---

### Step 2: Deconvolution

Run `decon.py` to perform the core ARIC algorithm. This process includes:
*   **Marker Selection**: Identifying highly informative features for each cell type.
*   **Weighted SVR**: Performing deconvolution using a weighted Support Vector Regression approach.

```bash
python decon.py
```

**Input:** Aligned files from Step 1.  
**Output:** Predicted cell type proportions for each sample in the mixture matrix.

---

## File Requirements

| File Name | Format | Description |
| :--- | :--- | :--- |
| `ref_raw.csv` | CSV | Raw signature matrix (Features x Cell Types) |
| `mix.csv` | CSV | Bulk mixture data (Features x Samples) |

## Environment

Ensure you have the following Python libraries installed:
*   `numpy`
*   `pandas`
*   `scikit-learn`

```bash
pip install numpy pandas scikit-learn
```

---

### 修改建议说明：
1.  **增加表格**：使用 Markdown 表格清晰展示输入文件要求。
2.  **环境说明**：增加了 `Environment` 章节，提示用户需要安装的依赖库（如 `pandas`, `sklearn`）。
3.  **结构化**：将“数据准备”、“对齐”和“反褶积”分段，并明确了每一步的输入输出。
4.  **专业用语**：使用了更通用的术语，如 "Features x Cell Types"，方便其他开发者理解。

您可以直接将上述内容复制并保存为 `README.md`。
