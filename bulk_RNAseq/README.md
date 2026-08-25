# Krt14+ RNA-seq Differential Expression Analysis

## Overview

This repository contains the analysis workflow used to compare the transcriptomic profiles of **Krt14+ cells** under two macrophage conditions:

- **Mac-Def** — macrophage-deficient condition
- **Mac-Prof** — macrophage-proficient condition

The analysis is performed on aligned RNA-seq BAM files and includes:

1. Gene-level read counting
2. Low-expression filtering
3. TMM normalization
4. Differential expression analysis using edgeR
5. Gene annotation using Ensembl/BioMart
6. Identification of differentially expressed protein-coding genes
7. Volcano plot generation
8. Top differentially expressed gene visualization
9. Reactome pathway GSEA
10. Redundant pathway collapsing
11. Leading-edge gene identification
12. Leading-edge heatmap visualization

---

## Experimental comparison

The primary comparison is:

**Mac-Def vs Mac-Prof**

The edgeR design defines:

- `Mac-Prof` as the reference condition
- `Mac-Def` as the comparison condition

Therefore:

- **Positive logFC** = higher expression in Mac-Def
- **Negative logFC** = higher expression in Mac-Prof

---

## Important statistical limitation

### No biological replicates

This dataset does **not contain biological replicates** for the two conditions.

Because biological replicates are unavailable, biological dispersion cannot be estimated empirically from replicate-to-replicate variation.

The differential expression analysis therefore uses an **assumed biological coefficient of variation (BCV)**:

```r
bcv <- 0.4
