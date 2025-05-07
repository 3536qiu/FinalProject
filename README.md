# FinalProject
This project aims to detect DNA copy number alterations (CNAs) in scRNA-seq data, with per‑cluster reference selection and CNAs detection.
## Features
-   Identifies reference cells within clusters.
-   Detects CNV events (gains/losses) per cell.
## Core Functions
1.  `find_cluster_references`: Identifies reference (e.g., diploid) cells within each cluster based on low gene expression variance.
2.  `detect_cnas_per_cluster`: Calculates per-window log2 fold-changes against reference cells and calls CNV gains/losses per chromosome for each cell.
3.  `extract_cna_segments`: Converts the per-cell, per-chromosome CNV calls into a table of genomic segments with their associated CNV status.
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
