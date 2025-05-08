# FinalProject
This project aims to detect DNA copy number alterations (CNAs) in scRNA-seq data, with per‑cluster reference selection and CNAs detection.
## Features
-   Identifies reference cells within clusters.
-   Detects CNA events (gains/losses) per cell.
## Core Functions
1.  `find_cluster_references`: Identifies reference (e.g., diploid) cells within each cluster based on low gene expression variance.
2.  `detect_cnas_per_cluster`: Calculates per-window log2 fold-changes against reference cells and calls CNA gains/losses per chromosome for each cell.
3.  `extract_cna_segments`: Converts the per-cell, per-chromosome CNA calls into a table of genomic segments with their associated CNA status.
4.  `compute_mean_f1`: Computed precision, recall, and F1 scores of different input parameters comparing to provided predefined CNAs. 
## Dependencies
-   pandas
-   numpy
-   scipy
-   scikit-learn
-   scanpy
-   matplotlib
## Getting Started

1. Cloning the Repository

```bash
git clone [https://github.com/3536qiu/FinalProject.git](https://github.com/3536qiu/FinalProject.git)
import os
os.chdir('FinalProject')
```
2. Once you have the files, you'll need to install the necessary Python packages.
```bash
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import scipy.sparse
from sklearn.mixture import GaussianMixture
from sklearn.metrics import precision_recall_fscore_support
from hmmlearn import hmm
from scipy import stats
from scipy.ndimage import uniform_filter1d
from scipy.ndimage import gaussian_filter1d
from scipy import stats
from scipy.signal import medfilt, savgol_filter
import re
import itertools
```
3. After obtaining the code and installing dependencies, you can import and use the functions from core_functions.py in your Python scripts or Jupyter notebooks:
```bash
import core_functions as cf
```
4. Now you can utilize the defined functions to analyze your scRNA-seq data! Enjoy! Below is an example usage
```bash
# Example: Load your AnnData object
# adata = sc.read_h5ad("path/to/your/data.h5ad")
# Important: remember to process your data using pca, neighboring, umap and leiden

# 1. Find reference cells
# cluster_references = cf.find_cluster_references(adata, cluster_key=cluster_key_name)

# 2. Detect CNAs
# calls_df, windows_data, var_sorted_info = cf.detect_cnas_per_cluster(
#     adata,
#     cluster_key=cluster_key_name,
#     cluster_refs=cluster_references
# )

# 3. Extract CNA segments
# segments_df = cf.extract_cna_segments(calls_df, windows_data, var_sorted_info)

# print(segments_df.head())
```

