# Workflow

The execution of **CelFEER** is divided into three main steps. Please follow them in order:

## Step 1: Data Preparation

Before running the scripts, organize your data into the expected formats under CelFEER. It is recommended to place your input data in a `ref_data/` and `test_data/` folder.

## Step 2: Marker Selection (`ref.py`)

Run `ref.py` (markers.py in CelFEER) using your processed scRNA-seq reference data (`ref.txt` in `ref_data/`) to extract cell-type-specific marker genes and generate the signature matrix.

```bash
python markers.py <input_file> <output_file> <num_values> <tissues> <depth_filter> <nan_filter> <extra_filter> <variant>
```

## Step 3: Deconvolution (`decon.py`)

Run `decon.py` (celfeer.py in CelFEER) using your prepared bulk RNA-seq data and the signature matrix generated in the previous step.

```bash
python decon.py
```

## Notes
*   **More Information**:  [https://github.com/pi-zz-a/CelFEER.git](https://github.com/pi-zz-a/CelFEER.git)
