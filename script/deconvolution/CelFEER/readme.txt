To run CelFEER effectively, please ensure your data is prepared and processed according to the following steps:
1. Data Preparation
   Before running the scripts, organize your data into the expected formats under CelFEER.
2. Workflow
   Run ref.py (markers.py in CelFEER) using your processed scRNA-seq reference data to extract cell-type-specific marker genes and generate the signature matrix.
   Run decon.py (celfeer.py in CelFEER) using your prepared bulk RNA-seq data and the signature matrix generated in the previous step.

Notes: We provide pre-processed data files in the test_data that can be used directly for testing the pipeline. Users can point the input path to these files to immediately evaluate the performance of CelFEER.Specifically, the filtered reference files used for deconvolution are located in the celfeer_ref, ensuring you have a streamlined starting point for your analysis.
