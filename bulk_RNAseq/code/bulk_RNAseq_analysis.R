
# Krt14+ RNA-seq analysis
# ============================================================
#
# Comparison:
#   Macrophage-Deficient vs Macrophage-Proficient Krt14+ cells
#
#
# Analysis:
#   1. Read alignment/counting
#   2. Low-expression filtering
#   3. TMM normalization
#   4. Differential expression with edgeR
#   5. Gene annotation
#   6. Volcano plot
#   7. Top DEG visualization
#   8. Reactome GSEA
#   9. Leading-edge heatmaps
#
# Important:
#   There are no biological replicates in this dataset.
#   Differential expression therefore uses an assumed BCV.
#
# ============================================================

# ============================================================
# 1. Load required packages
# ============================================================


library(Rsubread)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(fgsea)
library(msigdbr)
library(tibble)
library(openxlsx)
library(edgeR)
library(biomaRt)
library(tibble)


# Input BAM files

bam_files <- c("data/Mac-Def_Aligned.sortedByCoord.out.bam",
"data/Mac-Prof_Aligned.sortedByCoord.out.bam")

# Sample names

sample_names <- c("Mac-Def", "Mac-Prof")


# Reference annotation

gtf_file <- "data/genes.gtf"


# ============================================================
# 3. Generate gene-level count matrix
# ============================================================


fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,
  nthreads = 8
)

# Count matrix
count_matrix <- fc$counts


# ============================================================
# 4. Create sample metadata
# ============================================================

sample_names <- sub("_Aligned.*", "", basename(bam_files))

condition <- ifelse(grepl("Mac-Def", sample_names), "Def", "Prof")

colData <- data.frame(
  row.names = sample_names,
  condition = factor(condition)
)


# ============================================================
# 5. Differential expression analysis with edgeR
# ============================================================

# ------------------------------------------------------------
# 5.1 Create DGEList
# ------------------------------------------------------------

group <- factor(colData$condition, levels = c("Prof", "Def"))
dge <- DGEList(counts = count_matrix, group = group)


# ------------------------------------------------------------
# 5.2 TMM normalization
# ------------------------------------------------------------
dge <- calcNormFactors(dge)




# ------------------------------------------------------------
# 5.3 Filter low-expression genes
# ------------------------------------------------------------


keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Recalculate normalization after filtering
dge <- calcNormFactors(dge)

# QC: see how many genes were removed
cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(dge), "\n")


# ------------------------------------------------------------
# 5.4 Generate normalized logCPM values
# ------------------------------------------------------------
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

# Clean Ensembl IDs once
rownames(logCPM) <- sub("\\..*", "", rownames(logCPM))






# ------------------------------------------------------------
# 5.5 Differential expression testing
# ------------------------------------------------------------

# No biological replicates are available.
# Therefore, an assumed biological coefficient of variation
# is used to estimate dispersion.


bcv <- 0.4
# Perform exact test with assumed BCV

et <- exactTest(dge, dispersion = bcv^2)

# Extract top DE genes
res_edgeR <- topTags(et, n = Inf)
res_edgeR_df <- as.data.frame(res_edgeR)

write.csv(res_edgeR_df, "edgeR_results_no_replicates.csv")



# ============================================================
# 5.6 Gene Annotation and Differential Expression Results
# ============================================================


# Add Ensembl gene IDs as a column
res_edgeR_df <- rownames_to_column(res_edgeR_df, var = "ensembl_gene_id")


# Connect to Ensembl mouse annotation
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Extract Ensembl gene IDs for annotation
gene_list <- res_edgeR_df$ensembl_gene_id

# Retrieve mouse gene annotations from Ensembl
annot <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol", "description", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = gene_list,
  mart = mart
)


# Add gene annotation to the edgeR results
res_annot_final <- merge(
  res_edgeR_df,
  annot,
  by = "ensembl_gene_id", # Use the shared column for merging
  all.x = TRUE            # Keep all rows from your results, even if annotation is missing
)

# Export annotated differential expression result
write.xlsx(res_annot_final, file = "res_annot_final.xlsx")



# ============================================================
# 5.7 Identify Genes with Large Fold Changes
# ============================================================


# Genes with log2 fold change > 2
# Positive logFC indicates higher expression in Mac-Def
up_def <- res_annot_final %>% 
  filter(logFC > 2)

# Genes with log2 fold change < -2
# Negative logFC indicates higher expression in Mac-Prof
up_prof<- res_annot_final %>% 
  filter(logFC < -2)

# Select the 20 genes with the largest fold changes

top_prof <- up_prof %>%
  arrange(logFC) %>%
  head(20)
top_def <- up_def %>%
  arrange(desc(logFC)) %>%
  head(20)

# Export top genes

write.xlsx(top_prof, file = "top_prof.xlsx")
write.xlsx(top_def, file = "top_def.xlsx")


# ============================================================
# 5.8 Filter for Protein-Coding Genes
# ============================================================

# Restrict the differential expression results to
# protein-coding genes
res_protein_coding <- res_annot_final %>%
  filter(gene_biotype == "protein_coding")

# Report the number of genes removed
cat("Total genes in full results:", nrow(res_annot_final), "\n")
cat("Protein-coding genes remaining:", nrow(res_protein_coding), "\n")


# ============================================================
# 5.9 Identify Differentially Expressed Protein-Coding Genes
# ============================================================

# Protein-coding genes upregulated in Mac-Def
up_def_protein <- res_protein_coding %>% 
  filter(logFC > 2)

# Protein-coding genes upregulated in Mac-Prof
up_prof_protein<- res_protein_coding %>% 
  filter(logFC < -2)


# Select the top 20 protein-coding genes by fold chang
# Genes Upregulated in Mac-def (logFC > +2)
top_def_protein <- up_def_protein %>% 
  arrange(desc(logFC)) %>% 
  head(20)

# Genes Upregulated in Mac-Prof (logFC < -2)
top_prof_protein <- up_prof_protein %>% 
  arrange(logFC) %>% 
  head(20)


top_prof_protein  <- up_prof_protein %>% arrange(desc(logFC)) %>% head(20)
top_def_protein <- up_def_protein %>% arrange(logFC) %>% head(20)

# Export top protein-coding genes
write.xlsx(top_prof_protein, file = "top_prof_protein.xlsx")
write.xlsx(top_def_protein, file = "top_def_protein.xlsx")


# ============================================================
# 5.10 Volcano Plot
# ============================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

# Calculate -log10(P-value) and assign regulation groups

df_volcano <- res_protein_coding %>%
  mutate(
    negLogP = -log10(PValue),
    Regulation = case_when(
      logFC > 2 & FDR < 0.05  ~ "Up in Mac-Def",
      logFC < -2 & FDR < 0.05 ~ "Up in Mac-Prof",
      TRUE                    ~ "Not significant"
    )
  )


# Count significant genes in each direction

n_up_def <- nrow(filter(df_volcano, Regulation == "Up in Mac-Def"))
n_up_prof <- nrow(filter(df_volcano, Regulation == "Up in Mac-Prof"))

# Identify the top 10 genes on each side of the volcano plot
# for gene-labeling purposes

top_10_def <- df_volcano %>% 
  filter(Regulation == "Up in Mac-Def") %>% 
  arrange(desc(logFC)) %>% 
  head(10)

top_10_prof <- df_volcano %>% 
  filter(Regulation == "Up in Mac-Prof") %>% 
  arrange(logFC) %>% 
  head(10)

genes_to_label <- bind_rows(top_10_def, top_10_prof)


# Create volcano plot
ggplot(df_volcano, aes(x = logFC, y = negLogP)) +
  
  
# Plot individual genes 
geom_jitter(aes(color = Regulation), 
              alpha = 0.5, 
              size = 1.2, 
              width = 0.1, 
              height = 0.1) +
  
  
  # Label Mac-Prof side
  annotate("text", x = -7.5, y = 7.5, 
           label = "Upregulated in Mac-Prof", 
           color = "#4682B4", fontface = "bold", size = 3.5) +
  annotate("text", x = -7.5, y = 7.1, 
           label = paste0("(Total: ", n_up_prof, ")"), 
           color = "#4682B4", size = 3) +
  
  # Label Mac-Def side
  annotate("text", x = 7.5, y = 7.5, 
           label = "Upregulated in Mac-Def", 
           color = "#CD5C5C", fontface = "bold", size = 3.5) +
  annotate("text", x = 7.5, y = 7.1, 
           label = paste0("(Total: ", n_up_def, ")"), 
           color = "#CD5C5C", size = 3) +
  
  # Label selected genes
  
  geom_text_repel(
    data = genes_to_label,
    aes(label = mgi_symbol),
    size = 2.8,
    fontface = "italic",
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = "grey50",
    force = 2 # Increases repulsion to keep names away from the center
  ) +
  
  # Define colours for each regulation category
  scale_color_manual(values = c(
    "Up in Mac-Def" = "#CD5C5C", 
    "Up in Mac-Prof" = "#4682B4", 
    "Not Significant" = "grey92"
  )) +
  
  # Add fold-change thresholds
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "grey40") +
  
  # Add significance threshold
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40") +
  
  labs(
    title = "Macrophage Transcriptomic Profile",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~italic(P))
  ) +
  theme_classic() +
  theme(legend.position = "none")



# ============================================================
# 5.11 Top 20 Differentially Expressed Protein-Coding Genes
# ============================================================

library(tidyr)

# Making the bar plot of the top 20 up and down-regulated protein coding genes

# Combine the top Mac-Prof and Mac-Def genes
top_genes_combined <- bind_rows(top_prof_protein, top_def_protein)

# Reorder genes according to log2 fold change
top_genes_combined <- top_genes_combined %>%
  mutate(mgi_symbol = reorder(mgi_symbol, logFC))


# Generate bar plot
ggplot(top_genes_combined, aes(x = mgi_symbol, y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +  # Flip to make gene names readable
  scale_fill_manual(values = c("skyblue", "tomato"), 
                    labels = c("Deficient (Down)", "Proficient (Up)"),
                    name = "Direction") +
  theme_minimal() +
  labs(
    title = "Top 20 Up- and Down-regulated Proteins",
    subtitle = "Based on log2 Fold Change (Pooled Samples)",
    x = "Protein/Gene Symbol",
    y = "log2 Fold Change"
  ) +
  theme(axis.text.y = element_text(size = 8)) # Shrink text if names overlap





# ============================================================
# 5.12 Gene Ranking for GSEA
# ============================================================


library(GSVA)
library(clusterProfiler)

# Prepare protein-coding genes for pathway analysis
gsea_df <- res_protein_coding %>%
  filter(
    !is.na(mgi_symbol),
    mgi_symbol != ""
  ) %>%
  distinct(mgi_symbol, .keep_all = TRUE)

# Create ranked gene list using log2 fold change
ranked_genes <- gsea_df$logFC
names(ranked_genes) <- gsea_df$mgi_symbol

# Remove genes without valid gene symbols
ranked_genes <- ranked_genes[!is.na(names(ranked_genes))]
ranked_genes <- ranked_genes[names(ranked_genes) != ""]

# Sort genes from highest to lowest log2 fold change
ranked_genes <- sort(ranked_genes, decreasing = TRUE)


# Check for missing or empty gene names
any(names(ranked_genes) == "")
any(is.na(names(ranked_genes)))






# ============================================================
# 5.13 Reactome Pathway Enrichment Analysis
# ============================================================

# Retrieve Reactome gene sets for Mus musculus
msig_reactome <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "REACTOME"
)

# Convert Reactome gene sets into a list
reactome_list <- split(msig_reactome$gene_symbol, msig_reactome$gs_name)


# Run FGSEA
fgsea_reactome <- fgsea(
  pathways = reactome_list,
  stats    = ranked_genes,
  minSize  = 15,
  maxSize  = 500,
  nproc    = 1,
  nPermSimple = 10000
)


# Remove pathways with missing enrichment statistics
fgsea_reactome_clean <- fgsea_reactome %>%
  filter(!is.na(NES))


# Export Reactome enrichment results
write.xlsx(fgsea_reactome_clean, file = "fgsea_reactome_clean.xlsx")


# ============================================================
# 5.14 Reactome Pathway Analysis and Pathway Collapsing
# ============================================================

# Load packages required for pathway analysis

library(ReactomePA)
library(enrichplot)
library(msigdbr)
library(fgsea)
library(dplyr)
library(openxlsx)
library(ReactomePA)
library(enrichplot)
library(ggplot2)


# Getting reactome gene sets


msig_reactome <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "REACTOME"
)

reactome_list <- split(
  msig_reactome$gene_symbol,
  msig_reactome$gs_name
)



# Running FSEA


fgsea_reactome <- fgsea(
  pathways = reactome_list,
  stats = ranked_genes,
  minSize = 15,
  maxSize = 500,
  nproc = 1,
  nPermSimple = 10000
)



# Clean results


fgsea_reactome_clean <- fgsea_reactome %>%
  filter(
    !is.na(NES),
    !is.na(padj)
  )

write.xlsx(
  fgsea_reactome_clean,
  file = "fgsea_reactome_clean.xlsx"
)



# Ordering pathways by statistical significance

fgsea_ordered <- fgsea_reactome_clean %>%
  arrange(padj)


# ------------------------------------------------------------
# 5. Collapse redundant pathways
#
# Use the top 200 significant pathways as input to
# collapsePathways().
# ------------------------------------------------------------

top200_for_collapse <- fgsea_ordered %>%
  slice_head(n = 200)

collapsed_pathways <- collapsePathways(
  fgseaRes = top200_for_collapse,
  pathways = reactome_list,
  stats = ranked_genes
)


# ------------------------------------------------------------
# 6. Keep only the independent/main pathways
# ------------------------------------------------------------

fgsea_collapsed <- fgsea_reactome_clean %>%
  filter(pathway %in% collapsed_pathways$mainPathways) %>%
  arrange(padj)


# Save collapsed results
write.xlsx(
  fgsea_collapsed,
  file = "fgsea_collapsed.xlsx"
)


# ------------------------------------------------------------
# 7. Select the 20 MOST SIGNIFICANT non-redundant pathways
# ------------------------------------------------------------

top20_reactome <- fgsea_collapsed %>%
  arrange(padj) %>%
  slice_head(n = 20)


# Saving top 20
write.xlsx(
  top20_reactome,
  file = "fgsea_top20_reactome.xlsx"
)



# ============================================================
# 5.15 Reactome Pathway Visualisation
# ============================================================

# Plotting the top 20


top20_reactome %>%
  ggplot(aes(
    x = NES,
    y = reorder(pathway, NES),
    size = size,
    color = padj
  )) +
  geom_point() +
  scale_color_viridis_c(
    option = "plasma",
    direction = -1
  ) +
  labs(
    title = "Reactome pathway enrichment (Mac-Def vs Mac-Prof)",
    x = "Normalized Enrichment Score (NES)",
    y = "Reactome Pathway",
    size = "Gene set size",
    color = "Adjusted p-value"
  ) +
  theme_bw()




# Identify the 20 pathways with the largest absolute NES

top_20_pathways <- fgsea_collapsed %>%
  arrange(desc(abs(NES))) %>%
  head(20)


ggplot(
  top_20_pathways,
  aes(
    x = NES,
    y = reorder(pathway, NES),
    size = size,
    color = padj
  )
) +
  geom_point() +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  scale_color_viridis_c(
    option = "plasma",
    direction = -1
  ) +
  
  labs(
    title = "Top 20 Reactome Pathways",
    subtitle = "Mac-Def vs Mac-Prof",
    x = "Normalized Enrichment Score (NES)",
    y = "Reactome Pathway",
    size = "Gene count",
    color = "Adjusted p-value"
  ) +
  
  theme_bw()



# ============================================================
# 5.16 Leading-Edge Gene Heatmaps
# ============================================================


library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)
library(msigdbr)
library(fgsea)


# ------------------------------------------------------------
# Define Reactome pathways of biological interest
# ------------------------------------------------------------


# Define Pathway Strings
pathway_list <- list(
  "b-Catenin"   = "REACTOME_DEGRADATION_OF_BETA_CATENIN_BY_THE_DESTRUCTION_COMPLEX",
  "PTEN"        = "REACTOME_REGULATION_OF_PTEN_STABILITY_AND_ACTIVITY",
  "Starvation"  = "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
  "Hypoxia"     = "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
  "Unfolded"    = "REACTOME_MITOCHONDRIAL_UNFOLDED_PROTEIN_RESPONSE",
  "Mito import" = "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
  "Mito deg"    = "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION"
)


# ------------------------------------------------------------
# Function to extract top leading-edge genes
# ------------------------------------------------------------


get_top_genes <- function(pw_name, pw_vec, n = 10) {
  fgsea_reactome_clean %>%
    filter(pathway %in% pw_vec) %>%
    unnest(leadingEdge) %>%
    left_join(res_annot_final, by = c("leadingEdge" = "mgi_symbol")) %>%
    filter(!is.na(logFC)) %>%
    group_by(pathway) %>%
    slice_max(abs(logFC), n = n) %>%
    ungroup() %>%
    transmute(
      Gene = leadingEdge,
      Pathway = pw_name,
      Reactome = pathway
    )
}


# ------------------------------------------------------------
# Generate individual heatmaps for each selected pathway
# ------------------------------------------------------------ 

# Clear open graphical devices
while(!is.null(dev.list())) dev.off()

for (name in names(pathway_list)) {
  
  # Identify the top leading-edge genes
  current_gene_info <- get_top_genes(pw_name = name, pw_vec = pathway_list[[name]], n = 10)
  
  
  # Skip pathways if no genes are found
  
  if(nrow(current_gene_info) == 0) {
    message(paste("Skipping - No genes found for:", name))
    next
  }
  
# ----------------------------------------------------------
# Extract expression matrix for selected genes
# ----------------------------------------------------------
  gene_subset <- res_annot_final %>%
    filter(mgi_symbol %in% current_gene_info$Gene) %>%
    distinct(mgi_symbol, .keep_all = TRUE)
  
  sub_mat <- logCPM[rownames(logCPM) %in% gene_subset$ensembl_gene_id, , drop = FALSE]
  rownames(sub_mat) <- gene_subset$mgi_symbol[match(rownames(sub_mat), gene_subset$ensembl_gene_id)]
  colnames(sub_mat) <- c("Mac-Def", "Mac-Prof")
  
  # Manual Z-score Calculation (in order to prevent division-by-zero errors)
  row_means <- rowMeans(sub_mat)
  row_sds   <- apply(sub_mat, 1, sd)
  row_sds[row_sds == 0] <- 1
  sub_scaled <- (sub_mat - row_means) / row_sds
  
  # Choosing the colours
  limit <- max(abs(sub_scaled), na.rm = TRUE)
  if(limit == 0) limit <- 1
  col_fun <- colorRamp2(c(-limit, 0, limit), c("#2166AC", "white", "#B2182B"))
  
  # Annotating the columns
  ann_col <- data.frame(Condition = colData$condition)
  rownames(ann_col) <- colnames(sub_scaled)
  ann_colors <- list(Condition = c("Def" = "#D3D3D3", "Prof" = "#4F4F4F"))
  
  
  
  # ----------------------------------------------------------
  # Generate heatmap
  # ----------------------------------------------------------
  
  grid.newpage()
  p <- pheatmap(
    as.matrix(sub_scaled),
    color = col_fun,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = paste("Top 10 Genes:", name),
    border_color = "white",
    cellwidth = 40,
    cellheight = 18,
    heatmap_legend_param = list(title = "Z-score")
  )
  draw(p)
  message(paste("Successfully drew:", name))
}













