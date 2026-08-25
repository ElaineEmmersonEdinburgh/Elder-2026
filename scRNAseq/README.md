# Single-cell RNA-seq analysis

## Overview

This repository contains the R code used for the single-cell RNA-seq analysis
and visualization presented in the manuscript. This code was used to produce the data in figures S4B.

The analysis focuses on the identification and characterization of epithelial
ductal cell populations across the M (control), Day 3 (D3), and Day 28 (D28) timepoints.

The workflow includes:

- SoupX correction for ambient RNA contamination
- Quality control and filtering
- SCTransform normalization
- PCA and UMAP dimensionality reduction
- Harmony integration across timepoints
- Unsupervised clustering
- Cell-type annotation
- Ductal epithelial cell reclustering
- Differential gene expression analysis
- Comparison of gene expression between M and D3
- Analysis of D3 vs M within Krt14+ progenitor cells
- Generation of figures and marker-gene tables

## Dataset

Three timepoints are analyzed:

| Timepoint | Description |
|-----------|-------------|
| M | Control |
| D3 | Day 3 |
| D28 | Day 28 |

The input data consist of 10X Genomics filtered and raw feature-barcode
matrices for each timepoint.

Expected input directory structure:

data/
├── M/
│   ├── filtered_feature_bc_matrix/
│   └── raw_feature_bc_matrix/
├── D3/
│   ├── filtered_feature_bc_matrix/
│   └── raw_feature_bc_matrix/
└── D28/
    ├── filtered_feature_bc_matrix/
    └── raw_feature_bc_matrix/

## Analysis workflow

### 1. SoupX correction

Ambient RNA contamination was estimated and corrected independently for
each timepoint using SoupX.

The corrected count matrices are generated as:

- `M_out`
- `D3_out`
- `D28_out`

### 2. Quality control and preprocessing

The SoupX-corrected matrices are converted into Seurat objects.

Cells are filtered using:

- `nFeature_RNA > 200`
- `percent.mt < 20`

SCTransform normalization is subsequently performed.

### 3. Dataset integration and clustering

The three timepoints are merged and integrated using Harmony, with
`timepoint` used as the batch/integration variable.

PCA, Harmony, UMAP, nearest-neighbor analysis, and graph-based clustering
are used to identify cell populations.

### 4. Cell-type annotation

Clusters are annotated based on marker-gene expression.

The major populations identified include:

- Ductal epithelial cells
- Acinar cells
- Endothelial cells
- Macrophages
- Myofibroblasts
- Krt14+ progenitor cells

Clusters identified as contaminants are removed before downstream
ductal-cell analysis.

### 5. Ductal epithelial cell analysis

Ductal epithelial populations are isolated and reclustered independently.

The resulting populations include:

- Ascl3+ duct
- Krt18+ duct
- Wfdc18+ duct
- Krt14+ progenitors
- cKit+ duct
- Granular duct
- Myoepithelial cells
- Myofibroblasts

Marker genes are identified using Seurat's `FindAllMarkers()` function.

### 6. Differential gene expression

Differential expression is performed using Seurat's `FindMarkers()`.

Two comparisons are performed:

#### All ductal epithelial cells

D3 vs M:

```text
D3_vs_M_all_ductal_cells.xlsx
