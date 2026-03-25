To run CelFEER effectively, please ensure your data is prepared and processed according to the following steps:
1. Data Preparation
   Before running the scripts, organize your data into the expected formats under CelFEER.
2. Workflow
   Run ref.py (markers.py in CelFEER) using your processed scRNA-seq reference data to extract cell-type-specific marker genes and generate the signature matrix.
   Run decon.py (celfeer.py in CelFEER) using your prepared bulk RNA-seq data and the signature matrix generated in the previous step.
