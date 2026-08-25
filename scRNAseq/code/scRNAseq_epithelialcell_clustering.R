
############################################################
# Single-cell RNA-seq analysis
#
# This script contains the analysis and visualization
# code used for the manuscript. This code was used for the clustering and analysis of epithelial ductal cells. This script includes: 
#   -  SoupX correction
#   - Quality control and preprocessing
#   - Clustering and UMAP analysis
#   - Cell-type annotation
#   - Differential gene expression analysis
#   - Comparison of cell-cell communication between timepoints
#   - Manuscript figure generation

############################################################

library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(BiocManager)
library(multtest)
library (writexl)
library (biomaRt)
library (clustree)
library(ReactomeGSA)
library("GEOquery")
library(data.table)
library(SoupX)
library(DoubletFinder)
library(CellChat)
library(patchwork)
library(RColorBrewer)
library(EnhancedVolcano)
library(scCustomize)

############################################################
# 1. M (control) dataset - SoupX correction
############################################################

toc = Seurat::Read10X(data.dir = "data/M/filtered_feature_bc_matrix")
tod = Seurat::Read10X(data.dir = "data/M/raw_feature_bc_matrix")
sc = SoupChannel(tod, toc)

srat    <- CreateSeuratObject(toc)
srat    <- SCTransform(srat, verbose = F, return.only.var.genes=F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, verbose = T)

meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings

sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
sc  <- setDR(sc, umap)
sc = setContaminationFraction(sc, 0.45)
M_out = adjustCounts(sc, roundToInt = T)


############################################################
# 2. D3 dataset - SoupX correction
############################################################


toc = Seurat::Read10X(data.dir = "data/D3/filtered_feature_bc_matrix")
tod = Seurat::Read10X(data.dir = "data/D3/raw_feature_bc_matrix")
sc = SoupChannel(tod, toc)

srat    <- CreateSeuratObject(toc)
srat    <- SCTransform(srat, verbose = F, return.only.var.genes=F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, verbose = T)

meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
sc  <- setDR(sc, umap)
sc = setContaminationFraction(sc, 0.45)
D3_out = adjustCounts(sc, roundToInt = T)

############################################################
# 3. D28 dataset - SoupX correction
############################################################

toc = Seurat::Read10X(data.dir = "data/D28/filtered_feature_bc_matrix")
tod = Seurat::Read10X(data.dir = "data/D28/raw_feature_bc_matrix")
sc = SoupChannel(tod, toc)

srat    <- CreateSeuratObject(toc)
srat    <- SCTransform(srat, verbose = F, return.only.var.genes=F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, verbose = T)

meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
sc  <- setDR(sc, umap)
sc = setContaminationFraction(sc, 0.45)
D28_out = adjustCounts(sc, roundToInt = T)


############################################################
# 4. Dataset integration, QC and clustering
############################################################


################### Importing the adjusted matrix files

M <- CreateSeuratObject(M_out)
M$timepoint <- "M"
M[["percent.mt"]] <- PercentageFeatureSet(M, pattern = "^mt-")
M[["percent.ribo"]] <- PercentageFeatureSet(M, pattern = "^Rp[Sl]")
M <- subset(M, subset = nFeature_RNA > 200 & percent.mt < 20)

D3 <- CreateSeuratObject(D3_out)
D3$timepoint <- "D3"
D3[["percent.mt"]] <- PercentageFeatureSet(D3, pattern = "^mt-")
D3[["percent.ribo"]] <- PercentageFeatureSet(D3, pattern = "^Rp[Sl]")
D3 <- subset(D3, subset = nFeature_RNA > 200 & percent.mt < 20)

D28 <- CreateSeuratObject(D28_out)
D28$timepoint <- "D28"
D28[["percent.mt"]] <- PercentageFeatureSet(D28, pattern = "^mt-")
D28[["percent.ribo"]] <- PercentageFeatureSet(D28, pattern = "^Rp[Sl]")
D28 <- subset(D28, subset = nFeature_RNA > 200 & percent.mt < 20)


################### integrating using harmony
M  <- SCTransform(M, verbose = F, return.only.var.genes=F)
D3  <- SCTransform(D3, verbose = F, return.only.var.genes=F)
D28  <- SCTransform(D28, verbose = F, return.only.var.genes=F)

Combined1  <- merge(D3, y = c(D28, M), add.cell.ids = c("D3", "D28", "M"), project = "COMBINED")
Combined1  <- SCTransform(Combined1, verbose = F, return.only.var.genes=F, vars.to.regress = "percent.mt")
Combined1  <- RunPCA(Combined1, verbose = F)

Combined1 <- RunHarmony(Combined1, assay.use = "SCT", group.by.vars = "timepoint")
Combined1 <- RunUMAP(Combined1, reduction = "harmony", dims = 1:30)
Combined1 <- FindNeighbors(Combined1, reduction = "harmony", dims = 1:30)
Combined1 <- FindClusters(Combined1, resolution = 0.3)

DimPlot(Combined1, reduction = "umap", label = TRUE, pt.size = 1)
FeaturePlot(Combined1, features =c("Csf1r"))
FeaturePlot(Combined1, features =c("Epcam"))
FeaturePlot(Combined1, features =c("Vim"))


DefaultAssay(Combined1) <- "SCT"
Epi.markers1 <- FindAllMarkers(Combined1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
Combined1.markers.top10 <- Epi.markers1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)



####################  Making a dotplot to check the validity of the combined clusters

DotPlot(object = Combined1, features = c(unique(Combined1.markers.top10$gene)),  cols = "Spectral", col.min = 0, dot.scale = 2) + theme(axis.text.x = element_text(angle = 90, size = 6, face = "italic", hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.2), axis.title = element_blank()) + NoLegend() 

dev.off()


DotPlot(object = Combined1, 
        features = c(unique(Combined1.markers.top10$gene)), 
        cols = "Spectral", 
        col.min = 0, 
        dot.scale = 2) + 
  theme(
    axis.text.x = element_text(angle = 90, size = 6, face = "italic", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5),
    axis.title = element_blank()
  ) + 
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
  NoLegend()


####################  Refining and renaming the clusters for the combined object
####################  Removing contaminant clusters 17 and 18

Combined1_2 <- subset(Combined1, idents=c(0,1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15, 16))
DimPlot(Combined1_2, reduction = "umap", label = TRUE, pt.size = 1)

FeaturePlot(Combined1_2, features =c("Vim"))

Combined1_2 <- RenameIdents(object = Combined1_2, `0` = "Macrophage zero", `1` = "Acinar 1", `2` = "Endothelial 1", `3` = "Macrophage 1", `4`= "Duct 3" , `5` = "Macrophage 4", `6` = "Duct 1", `7` = "Duct 2", `8` = "Endothelial 2", `9` = "Acinar 4", `10` = "Acinar 3", `11` = "Krt14+ progenitors", `12` = "Endothelial 3", `13` = "Acinar 2", `14` = "Macrophage 2", `15` = "Macrophage 3", `16` = "Myofibroblasts")
DimPlot(Combined1_2, reduction = "umap", label = TRUE, pt.size = 1)



#################### Re-order identities for easier visualization of sub clusters from similar populations
target <- c("Krt14+ progenitors","Duct 1","Duct 2","Duct 3","Acinar 1","Acinar 2","Acinar 3","Acinar 4","Myofibroblasts", "Endothelial 1", "Endothelial 2", "Endothelial 3", "Macrophage zero", "Macrophage 1", "Macrophage 2", "Macrophage 3", "Macrophage 4")
Idents(Combined1_2) <- factor(Idents(Combined1_2), levels = target)
DimPlot(Combined1_2, reduction = "umap", label = TRUE, pt.size = 1)

DefaultAssay(Combined1_2) <- "SCT"
Combined1_2.markers <- FindAllMarkers(Combined1_2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
Combined1_2.markers.top10 <- Combined1_2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


Combined1_2.markers$cluster <- factor(Combined1_2.markers$cluster, levels = target)
Combined1_2.markers.top10 <- Combined1_2.markers.top10[order(Combined1_2.markers.top10$cluster), ]



#################### Making a dotplot with the top 10 genes in each cluster  


DotPlot(object = Combined1_2, features = c(unique(Combined1_2.markers.top10$gene)),  cols = "Spectral", col.min = 0, dot.scale = 2) + theme(axis.text.x = element_text(angle = 90, size = 6, face = "italic", hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.2), axis.title = element_blank()) + NoLegend() 

dev.off()

DotPlot(object = Combined1_2, 
        features = c(unique(Combined1_2.markers.top10$gene)), 
        cols = "Spectral", 
        col.min = 0, 
        dot.scale = 2) + 
  theme(
    axis.text.x = element_text(angle = 90, size = 6, face = "italic", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5),
    axis.title = element_blank()
  ) + 
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) + # Adjust this for spacing
  NoLegend()




############################################################
# 5. Ductal epithelial cell reclustering 
############################################################


################# Subsetting and reclustering the ductal epithelial cells only

EpCAM_subset <- subset(Combined1_2, idents=c(4,6,7,11))

EpCAM_subset    <- SCTransform(EpCAM_subset, verbose = F, return.only.var.genes=F)
EpCAM_subset    <- RunPCA(EpCAM_subset, verbose = F)
EpCAM_subset    <- RunUMAP(EpCAM_subset, dims = 1:30, verbose = F)
EpCAM_subset    <- FindNeighbors(EpCAM_subset, dims = 1:30, verbose = F)
EpCAM_subset    <- FindClusters(EpCAM_subset, verbose = T, resolution = 0.2)

##############UMAP visualisation
DimPlot(EpCAM_subset, reduction = "umap", label = TRUE, pt.size = 1)
FeaturePlot(EpCAM_subset, features =c("Krt14"))
FeaturePlot(EpCAM_subset, features =c("Krt5"))
FeaturePlot(EpCAM_subset, features =c("Cd83"))

DimPlot(EpCAM_subset, reduction = "umap", label = F, pt.size = 2, split.by = "timepoint")

EpCAM_subset$timepoint <- factor(EpCAM_subset$timepoint, levels = c("M", "D3", "D28"))
DimPlot(EpCAM_subset, split.by = "timepoint")

DefaultAssay(EpCAM_subset) <- "SCT"
Epi.markers <- FindAllMarkers(EpCAM_subset, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)

VlnPlot(EpCAM_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

EpCAM_subset <- FindVariableFeatures(EpCAM_subset, selection.method = "vst", nfeatures = 2000)


##############UMAP visualisation
top10 <- head(VariableFeatures(EpCAM_subset ), 10)

############## Plot variable features with and without labels
plot1 <- VariableFeaturePlot(EpCAM_subset )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


############## Refining the Epcam subset to remove any non-ductal cell contaminants 

EpCAM_subset_2 <- subset(EpCAM_subset, idents=c(0,1,2,3,4,5,8,9))

EpCAM_subset_2    <- SCTransform(EpCAM_subset_2, verbose = F, return.only.var.genes=F)
EpCAM_subset_2    <- RunPCA(EpCAM_subset_2, verbose = F)
EpCAM_subset_2    <- RunUMAP(EpCAM_subset_2, dims = 1:30, verbose = F)
EpCAM_subset_2    <- FindNeighbors(EpCAM_subset_2, dims = 1:30, verbose = F)
EpCAM_subset_2    <- FindClusters(EpCAM_subset_2, verbose = T, resolution = 0.2)

DimPlot(EpCAM_subset_2, reduction = "umap", label = F, pt.size = 2)


dim_plot <- DimPlot(EpCAM_subset_2, reduction = "umap", split.by = "timepoint", pt.size = 2)

dim_plot$data$timepoint <- factor(x = dim_plot$data$timepoint, levels = c("M", "D3", "D28"))
dim_plot

DimPlot(EpCAM_subset_2, split.by = "timepoint")
FeaturePlot(EpCAM_subset_2, features =c("Cdkn1a"), split.by = "timepoint")


##############  Labelling the clusters and changing their colours

EpCAM_subset_2 <- RenameIdents(object = EpCAM_subset_2, `0` = "Ascl3+ duct", `1` = "Krt18+ duct", `2` = "Wfdc18+ duct", `3` = "Krt14+ progenitors", `4` = "cKit+ duct", `5` = "Granular duct", `6` = "Myofibroblasts", `7` = "Myoepithelial cells")
DimPlot(EpCAM_subset_2, reduction = "umap", label = TRUE, pt.size = 1)

DimPlot(EpCAM_subset_2, reduction = "umap", label = FALSE, split.by = "timepoint", pt.size = 1,  cols = c(`Ascl3+ duct`='#B6E6BD',`Krt18+ duct`='#BAE3F2',`Wfdc18+ duct`= '#F7C8EE',`Krt14+ progenitors`='#CAC3F7',`cKit+ duct`='#FFCA96',
                                                                                                          `Granular duct`='#F59A8E',`Myoepithelial cells`='#FEC8D8',`Myofibroblasts`='#F7D05E'))
EpCAM_subset_2 <- RenameIdents(object = EpCAM_subset_2, `0` = "Ascl3+ duct", `1` = "Krt18+ duct", `2` = "uknown", `3` = "Krt14+ progenitors", `4` = "cKit+ duct", `5` = "Granular duct", `6` = "Myoepithelial cells", `7` = "Myofibroblasts")
DimPlot(EpCAM_subset_2, reduction = "umap", label = TRUE, pt.size = 1)

EpCAM_subset_2$timepoint <- factor(EpCAM_subset$timepoint, levels = c("M", "D3", "D28"))
DimPlot(EpCAM_subset_2, split.by = "timepoint")



############################################################
# 6. Differential gene expression analysis for the final ductal epithelial cell cluster
############################################################


############## Finding conserved markers for the Epcam_2 subset


DefaultAssay(EpCAM_subset_2) <- "SCT"
Epi.markers <- FindAllMarkers(EpCAM_subset_2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
top10 <- Epi.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(EpCAM_subset_2, features = top10$gene) + NoLegend()

xlstop10<- cbind(" "=rownames(top10), top10)
writexl::write_xlsx(xlstop10, "xlstop10.xlsx")


##############  Finding differentially expressed genes between D3 vs M

############################################################
# Differential gene expression: D3 vs M
############################################################

DefaultAssay(EpCAM_subset_2) <- "SCT"

# Set identities to timepoint
Idents(EpCAM_subset_2) <- "timepoint"

# Compare D3 vs M
D3_vs_M <- FindMarkers(
  EpCAM_subset_2,
  ident.1 = "D3",
  ident.2 = "M",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = FALSE
)

D3_up <- D3_vs_M %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)

writexl::write_xlsx(
  D3_up,
  "D3_upregulated_vs_M.xlsx"
)


DefaultAssay(EpCAM_subset_2) <- "SCT"

# Keep only Krt14+ progenitors
Krt14_subset <- subset(
  EpCAM_subset_2,
  idents = "Krt14+ progenitors"
)

# Set identity to timepoint
Idents(Krt14_subset) <- "timepoint"

# D3 vs M
Krt14_D3_vs_M <- FindMarkers(
  Krt14_subset,
  ident.1 = "D3",
  ident.2 = "M",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = FALSE
)

Krt14_D3_vs_M <- Krt14_D3_vs_M %>%
  tibble::rownames_to_column("gene")

writexl::write_xlsx(
  Krt14_D3_vs_M,
  "Krt14_progenitors_D3_vs_M.xlsx"
)


##############  Violin plots for the expression of Krt14, Krt5, Trp53 and Cdk1na split by timepoint 

VlnPlot(
  EpCAM_subset_2,
  features = "Krt14",
  group.by = "ident",
  split.by = "timepoint",
  pt.size = 0
)

VlnPlot(
  EpCAM_subset_2,
  features = "Krt5",
  group.by = "ident",
  split.by = "timepoint",
  pt.size = 0
)

VlnPlot(
  EpCAM_subset_2,
  features = "Trp53",
  group.by = "ident",
  split.by = "timepoint",
  pt.size = 0
)

VlnPlot(
  EpCAM_subset_2,
  features = "Cdkn1a",
  group.by = "ident",
  split.by = "timepoint",
  pt.size = 0
)


