# Biological question

The overall aim of this project is to investigate how irradiation affects the
transcriptional state and cellular interactions of ductal epithelial cells,
with a particular focus on Krt14+ progenitor populations.

In addition, the study aims to identify cell-cell communication between
Krt14+ progenitors and macrophage populations to investigate how interactions
between these cell populations may contribute to salivary gland regeneration.

The single-cell analysis examines cellular composition and transcriptional
states of ductal epithelial cells at different timepoints:

- **M** — control / baseline
- **D3** — day 3 post-irradiation
- **D28** — day 28 post-irradiation

A separate single-cell analysis focuses on ductal epithelial cells and
macrophages to investigate potential ligand-receptor interactions between
these populations, with particular attention to Krt14+ progenitor cells.

A separate bulk RNA-seq analysis compares the transcriptional state of Krt14+
cells under the following conditions:

- **Mac-Def** — macrophage-deficient
- **Mac-Prof** — macrophage-proficient

These analyses are complementary. The single-cell analysis is used to identify
cellular populations, characterize their transcriptional states, and
investigate cell-cell communication. The bulk RNA-seq analysis is used to
evaluate transcriptional differences in Krt14+ cells associated with
macrophage presence or absence.

# Repository structure

.
├── README.md
│
├── scRNAseq/
│   ├── README.md
│   ├── scRNAseq_ductal_analysis.R
│   └── scRNAseq_cellchat_analysis.R
│
├── bulk_RNAseq/
│   ├── README.md
│   └── Krt14_bulk_RNAseq_analysis.R
│
├── input/
│   └── KRT14-mac_pathways.xlsx
│
└── .gitignore

# Data availability

Raw sequencing data and large intermediate files are not included in this
repository due to file size and data-sharing restrictions.

The single-cell RNA-seq analysis requires the following input datasets:

- M (control)
- D3 (day 3 post-irradiation)
- D28 (day 28 post-irradiation)

Each dataset consists of raw and filtered 10x Genomics expression matrices.

The bulk RNA-seq analysis requires aligned BAM files for:

- Mac-Def
- Mac-Prof

A corresponding gene annotation GTF file is also required.

These files should be placed in the appropriate `data/` directory before
running the analysis scripts.

For example:

```text
data/
├── M/
│   ├── filtered_feature_bc_matrix/
│   └── raw_feature_bc_matrix/
│
├── D3/
│   ├── filtered_feature_bc_matrix/
│   └── raw_feature_bc_matrix/
│
├── D28/
│   ├── filtered_feature_bc_matrix/
│   └── raw_feature_bc_matrix/
│
├── genes.gtf
├── Mac-Def_Aligned.sortedByCoord.out.bam
└── Mac-Prof_Aligned.sortedByCoord.out.bam


## Data availability

The raw sequencing data are not included in this repository.

The raw scRNA-seq and bulk RNA-seq data are available via the European Molecular Biology Laboratory (EMBL) - European Bioinformatics Institute (EBI) public database ArrayExpress (https://www.ebi.ac.uk/biostudies/arrayexpress): scRNA-seq - accession code: E-MTAB-13374; bulk-seq - accession code: E-MTAB-17509. 

