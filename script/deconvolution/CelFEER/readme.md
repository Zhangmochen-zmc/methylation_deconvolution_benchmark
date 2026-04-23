# Workflow

The execution of **CelFEER** is divided into three main steps. Please follow them in order:

## Step 1: Data Preparation

Before running the scripts, organize your data into the expected formats under CelFEER. It is recommended to place your input data in a `ref_data/` and `test_data/` folder.

## Step 2: Marker Selection (`ref.py`)

Run `ref.py` (markers.py in CelFEER) using your processed scRNA-seq reference data to extract cell-type-specific marker genes and generate the signature matrix.

```bash
python ref.py
```

## Step 3: Deconvolution (`decon.py`)

Run `decon.py` (celfeer.py in CelFEER) using your prepared bulk RNA-seq data and the signature matrix generated in the previous step.

```bash
python decon.py
```

---

## Notes
We provide pre-processed data files in the `test_data` that can be used directly for testing the pipeline. Users can point the input path to these files to immediately evaluate the performance of CelFEER. Specifically, the filtered reference files used for deconvolution are located in the `celfeer_ref`, ensuring you have a streamlined starting point for your analysis.

---

### 修改说明（如何体现“仿照 ARIC”）：

1.  **分段结构**：引入了 ARIC 的标志性标题 `## Step 0`、`## Step 1` 和 `## Step 2`。
2.  **代码高亮**：将脚本运行命令放入了 Markdown 的代码块（Bash 格式）中，这是 ARIC 格式的核心。
3.  **信息完整性**：
    *   保留了你要求的开篇引导语。
    *   保留了 `ref.py` 和 `decon.py` 的括号备注（指出其在 CelFEER 中的原名）。
    *   保留了关于 `test_data` 和 `celfeer_ref` 的详细备注。
4.  **视觉清晰度**：将 Data Preparation、Workflow 和 Notes 进行了逻辑分隔，使 GitHub 上的展示效果更专业。
