# Single-cell RNA-seq and CellChat Analysis of Ductal Epithelial Cells and Macrophages

## Overview

This repository contains the R scripts used for the single-cell RNA-seq analysis of ductal epithelial cell populations and macrophages. This code was used to generate the data found in figure 3.

The analysis compares:

- **M** — control condition
- **D3** — day 3 post-irradiation (IR)

The analysis focuses on the interaction between ductal epithelial cells, including Krt14+ progenitor cells, and macrophage populations.

The workflow includes:

1. SoupX ambient RNA correction
2. Quality control and preprocessing
3. Single-cell clustering and UMAP visualization
4. Identification of epithelial and macrophage populations
5. Reclustering of ductal epithelial/macrophage populations
6. Cell-type annotation
7. Differential gene expression analysis
8. CellChat cell-cell communication analysis
9. Comparison of cell-cell communication between M and D3
10. Ligand-receptor interaction visualization
11. Manuscript figure generation

---

# Experimental design

The analysis contains two single-cell RNA-seq datasets:

| Condition | Description |
|---|---|
| M | Control |
| D3 | Day 3 post-irradiation |

The primary cell populations of interest are:

- Ductal epithelial cells
- Krt14+ progenitors
- Macrophages

Cell-cell communication is subsequently compared between:

```text
M control
vs
D3 post-irradiation
