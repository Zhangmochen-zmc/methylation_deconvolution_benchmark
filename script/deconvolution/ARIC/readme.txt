这是一个专门为你的 **ARIC** 项目定制的 `README.md` 完整模板。你可以直接将其中的内容复制到你的 GitHub 仓库中。

---

# ARIC: Accurate and Robust Inference of Cell Type Proportions

ARIC 是一种用于大体转录组（Bulk RNA-seq）或 DNA 甲基化数据（DNA Methylation）解卷积的计算框架。它通过两步标记基因选择策略（消除特征共线性与自适应异常值剔除），实现基于加权 υ-SVR 的高精度细胞比例预测。

## 目录
- [环境配置](#环境配置)
- [快速开始](#快速开始)
  - [Step 0: 准备输入文件](#step-0-准备输入文件)
  - [Step 1: 特征对齐与提取交集](#step-1-特征对齐与提取交集)
  - [Step 2: 执行解卷积计算](#step-2-执行解卷积计算)
- [输出说明](#输出说明)

---

## 环境配置

在运行脚本之前，请确保您的 Python 环境中已安装以下依赖库：

```bash
# 克隆仓库
git clone https://github.com/YourUsername/ARIC.git
cd ARIC

# 安装必要的 Python 包
pip install numpy pandas scikit-learn
```

---

## 快速开始

项目的核心流程分为三步，请按顺序执行：

### Step 0: 准备输入文件

准备两个 `.csv` 格式的矩阵文件，建议存放在 `data/` 目录下：

1.  **Reference Matrix (`ref.csv`)**: 参考数据集。
    *   行：特征（如 Gene Symbol 或 CpG 位点）
    *   列：已知的细胞类型
2.  **Mixture Matrix (`mix.csv`)**: 待解卷积的 Bulk 数据。
    *   行：特征（必须包含参考文件中的特征）
    *   列：待测样本名

**文件结构示例：**
```text
ARIC/
├── data/
│   ├── ref.csv          # 参考矩阵
│   └── mix.csv          # 待解卷积矩阵
├── ref.py               # 第一步脚本：取交集
├── decon.py             # 第二步脚本：解卷积
└── README.md
```

---

### Step 1: 特征对齐与提取交集

运行 `ref.py` 将参考文件和混合文件中的共有特征（基因/位点）提取出来，并生成预处理后的中间文件。

```bash
# 运行预处理脚本
python ref.py
```

*   **输入**: `data/ref.csv`, `data/mix.csv`
*   **输出**: 自动生成对齐后的中间矩阵文件（如 `ref_intersect.csv` 和 `mix_intersect.csv`）。

---

### Step 2: 执行解卷积计算

运行 `decon.py` 开始执行 ARIC 核心算法（包括标记基因选择和加权 SVR 计算）。

```bash
# 运行解卷积主程序
python decon.py
```

*   **主要逻辑**: 
    1.  自动加载第一步生成的中间文件。
    2.  利用 ARIC 的自适应策略筛选最优 Marker Genes。
    3.  输出每个样本的细胞类型比例。

---

## 输出说明

脚本运行完成后，结果将保存在 `results/` 文件夹（或当前目录）下：

*   **`final_proportions.csv`**: 这是最终的结果文件。
    *   **行**: 样本名称（Sample ID）
    *   **列**: 各种细胞类型的占比（Sum to 1）

| SampleID | CellType_A | CellType_B | CellType_C |
| :--- | :--- | :--- | :--- |
| Sample_01 | 0.45 | 0.30 | 0.25 |
| Sample_02 | 0.12 | 0.68 | 0.20 |

---

## 常见问题 (FAQ)

1.  **特征不匹配怎么办？**
    `ref.py` 会自动处理不匹配的情况，只保留两个矩阵共有的特征。请确保两者的 ID 类型一致（例如同为 Gene Symbol 或同为 Ensembl ID）。
2.  **运行速度慢？**
    如果特征（基因）数量超过 20,000，计算量会增加。建议在高性能服务器上运行或适当预过滤非表达基因。

---

## 引用 (Citation)

如果您在学术工作中使用了 ARIC，请引用：
> *Author, et al. "ARIC: Accurate and robust inference of cell type proportions from bulk gene expression or DNA methylation data." (202X).*

---

### 如何自定义此文件？
1.  **替换链接**：将 `YourUsername` 替换为你真实的 GitHub 用户名。
2.  **细化参数**：如果你的脚本支持命令行参数（例如 `python decon.py --input mydata.csv`），请在 Bash 命令中补充完整。
3.  **添加联系方式**：可以在末尾留下你的邮箱。
