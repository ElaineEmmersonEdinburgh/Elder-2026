############################################################
# Single-cell RNA-seq and CellChat analysis
#
# This script contains the analysis and visualization code
# used for the manuscript. This script includes the analysis of ductal epithelial cell populations and macrophages including:
#   - SoupX correction
#   - Quality control and preprocessing
#   - Clustering and UMAP analysis
#   - Cell-type annotation
#   - Differential gene expression analysis
#   - CellChat analysis
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
library(NMF)
library(ComplexHeatmap)
library(readxl)

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
# 2. M dataset - QC and clustering
############################################################


################### Importing the adjusted matrix files

M <- CreateSeuratObject(M_out)
M$timepoint <- "M"
M[["percent.mt"]] <- PercentageFeatureSet(M, pattern = "^mt-")
M[["percent.ribo"]] <- PercentageFeatureSet(M, pattern = "^Rp[Sl]")
M <- subset(M, subset = nFeature_RNA > 200 & percent.mt < 20)


###################integrating using harmony

M <- SCTransform(M, verbose = F, return.only.var.genes=F, vars.to.regress = "percent.mt")
M  <- RunPCA(M, verbose = F)

M <- RunUMAP(M, dims = 1:30)
M <- FindNeighbors(M, dims = 1:30)
M <- FindClusters(M, resolution = 0.3)



DimPlot(M, reduction = "umap", label = T)
DimPlot(M, reduction = "umap", label = F)

ElbowPlot(M)

FeaturePlot(M, features = c("Epcam"))
FeaturePlot(M, features = c("Aqp5"))
FeaturePlot(M, features = c("Adgre1"))
FeaturePlot(M, features = c("Cx3cr1"))
FeaturePlot(M, features = c("Csf1r"))
FeaturePlot(M, features = c("Cd163"))



############################################################
# 3. M dataset - epithelial/macrophage reclustering
###########################################################


## Ductal epithelial cells and macrophages reclustering 
EpCAM_subset_M <- subset(M, idents=c(6,7,10,2,5,13,12))

EpCAM_subset_M    <- SCTransform(EpCAM_subset_M, verbose = F, return.only.var.genes=F)
EpCAM_subset_M    <- RunPCA(EpCAM_subset_M, verbose = F)
EpCAM_subset_M    <- RunUMAP(EpCAM_subset_M, dims = 1:12, verbose = F)
EpCAM_subset_M    <- FindNeighbors(EpCAM_subset_M, dims = 1:12, verbose = F)
EpCAM_subset_M    <- FindClusters(EpCAM_subset_M, verbose = T, resolution = 0.3)

DimPlot(EpCAM_subset_M, reduction = "umap", label = TRUE, pt.size = 1)
ElbowPlot(EpCAM_subset_M)

##############UMAP visualisation
DimPlot(EpCAM_subset_M, reduction = "umap", label = TRUE, pt.size = 1)
FeaturePlot(EpCAM_subset_M, features =c("Krt14"))
FeaturePlot(EpCAM_subset_M, features =c("Cd163"))
FeaturePlot(EpCAM_subset_M, features =c("Ascl3"))
FeaturePlot(EpCAM_subset_M, features =c("Top2a"))


EpCAM_subset_F <- subset(EpCAM_subset_M, idents=c(0,1,2,3,4,5,7,8))
DimPlot(EpCAM_subset_F, reduction = "umap", label = TRUE, pt.size = 1)

FeaturePlot(EpCAM_subset_F, features =c("Csf1r"))
FeaturePlot(EpCAM_subset_F, features =c("Hexb"))
FeaturePlot(EpCAM_subset_F, features =c("Ngf"))

############################################################
# 4. M dataset - differential expression analysis
############################################################

#############Differential Gene expression analysis (ALL)
DefaultAssay(EpCAM_subset_F) <- "SCT"
xlsEpi.markers_Epcam_F <- FindAllMarkers(EpCAM_subset_F, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)

xlsEpi.markers_Epcam_F<- cbind(" "=rownames(xlsEpi.markers_Epcam_F), xlsEpi.markers_Epcam_F)





### top 10 differentially expressd genes per cluster, IR and no IR 
top10 <- xlsEpi.markers_Epcam_F %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)




##identify the clusters of M
DefaultAssay(EpCAM_subset_M) <- "SCT"
xlsEpi.markers_Epcam_M <- FindAllMarkers(EpCAM_subset_M, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)

top10_M <- xlsEpi.markers_Epcam_M %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

##comparing mac C0 to mac C3

xlsC0vsC3 <- FindMarkers(EpCAM_subset_F, assay = "SCT", ident.1 = "0", ident.2 = "3", min.pct = 0.25, verbose = FALSE, recorrect_umi = F)


############################################################
# 5. M dataset - cell type annotation
############################################################

## Naming the clusters

EpCAM_subset_F <- RenameIdents(object = EpCAM_subset_F, `0` = "Wfdc18+ duct", `1` = "John Mac 0", `2` = "Ascl3+ duct", `3` = "John mac 1", `4` = "Krt18+ duct", `5` = "K14 progenitors", `7` = "Granular duct", `8` = "John mac 3 cd163")
DimPlot(EpCAM_subset_F, reduction = "umap", label = TRUE, pt.size = 1)

############################################################
# 6. M dataset - CellChat analysis
############################################################

##cellchat on M

data.input = EpCAM_subset_F

data.input <- GetAssayData(EpCAM_subset_F, assay = "RNA", slot ="data")
labels <- Idents(EpCAM_subset_F)
meta <- data.frame(group = labels, row.names = names(labels))

EpCAM_subset_cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
EpCAM_subset_cellchat<- updateCellChat(EpCAM_subset_cellchat)

CellChatDB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
EpCAM_subset_cellchat@DB <- CellChatDB.use
EpCAM_subset_cellchat <- subsetData(EpCAM_subset_cellchat)

EpCAM_subset_cellchat <- identifyOverExpressedGenes(EpCAM_subset_cellchat)
EpCAM_subset_cellchat <- identifyOverExpressedInteractions(EpCAM_subset_cellchat)
EpCAM_subset_cellchat <- projectData(EpCAM_subset_cellchat, PPI.mouse)
EpCAM_subset_cellchat <- computeCommunProb(EpCAM_subset_cellchat)
EpCAM_subset_cellchat <- filterCommunication(EpCAM_subset_cellchat, min.cells = 10)
EpCAM_subset_cellchat <- netAnalysis_computeCentrality(EpCAM_subset_cellchat)

df.net.M <- subsetCommunication(EpCAM_subset_cellchat)

writexl::write_xlsx(df.net.M, "output/df.net.M.xlsx")


EpCAM_subset_cellchat <- computeCommunProbPathway(EpCAM_subset_cellchat)
EpCAM_subset_cellchat <- aggregateNet(EpCAM_subset_cellchat)
groupSize <- as.numeric(table(EpCAM_subset_cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(EpCAM_subset_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(EpCAM_subset_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



############################################################
# 7. D3 dataset - SoupX correction
############################################################

################# D3 soup

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



D3 <- CreateSeuratObject(D3_out)
D3$timepoint <- "D3"
D3[["percent.mt"]] <- PercentageFeatureSet(D3, pattern = "^mt-")
D3[["percent.ribo"]] <- PercentageFeatureSet(D3, pattern = "^Rp[Sl]")
D3 <- subset(D3, subset = nFeature_RNA > 200 & percent.mt < 20)

############################################################
# 8. D3 dataset - QC and clustering
############################################################

D3 <- SCTransform(D3, verbose = F, return.only.var.genes=F, vars.to.regress = "percent.mt")
D3  <- RunPCA(D3, verbose = F)

D3 <- RunUMAP(D3, dims = 1:30)
D3 <- FindNeighbors(D3, dims = 1:30)
D3 <- FindClusters(D3, resolution = 0.3)

DimPlot(D3, reduction = "umap", label = T)
DimPlot(D3, reduction = "umap", label = F)

FeaturePlot(D3, features = c(""))
FeaturePlot(D3, features = c("Adgre1"))

FeaturePlot(D3, features = c("Epcam"))

DefaultAssay(D3) <- "SCT"
xlsEpi.markers_D3 <- FindAllMarkers(D3, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)

############################################################
# 9. D3 dataset - epithelial/macrophage reclustering
############################################################


##Epcam macrophage reclustering
EpCAM_subset_D3 <- subset(D3, idents=c(0,6,4,10,9,8,11))

EpCAM_subset_D3    <- SCTransform(EpCAM_subset_D3, verbose = F, return.only.var.genes=F)
EpCAM_subset_D3    <- RunPCA(EpCAM_subset_D3, verbose = F)
EpCAM_subset_D3   <- RunUMAP(EpCAM_subset_D3, dims = 1:20, verbose = F)
EpCAM_subset_D3   <- FindNeighbors(EpCAM_subset_D3, dims = 1:20, verbose = F)
EpCAM_subset_D3 <- FindClusters(EpCAM_subset_D3, verbose = T, resolution = 0.3)

ElbowPlot(EpCAM_subset_D3)

##############UMAP visualisation
DimPlot(EpCAM_subset_D3, reduction = "umap", label = TRUE, pt.size = 1)
FeaturePlot(EpCAM_subset_D3, features =c("Hexb"))
FeaturePlot(EpCAM_subset_D3, features =c("Mki67"))
FeaturePlot(EpCAM_subset_D3, features =c("Prkdc"))


EpCAM_subset_F_D3 <- subset(EpCAM_subset_D3, idents=c(0,1,2,3,4,5,6,7))
DimPlot(EpCAM_subset_F_D3, reduction = "umap", label = TRUE, pt.size = 1)

FeaturePlot(EpCAM_subset_F_D3, features =c("Ascl3"))


############################################################
# 10. D3 dataset - differential expression analysis
############################################################


#############Differential Gene expression analysis (ALL)
DefaultAssay(EpCAM_subset_D3) <- "SCT"
Epi.markers_D3_Epcam <- FindAllMarkers(EpCAM_subset_D3, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)


### top 10 differentially expressd genes per cluster, IR and no IR 
top10_d3 <- Epi.markers_D3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


############################################################
# 11. D3 dataset - cell type annotation
############################################################

##naming clusters
EpCAM_subset_F_D3 <- RenameIdents(object = EpCAM_subset_F_D3, `0` = "John mac 0", `1` = "Ascl3+ duct", `2` = "Wfdc18+ duct", `3` = "John mac 1", `4` = "Granular duct", `5` = "Krt18+ duct", `6` = "K14 progenitors", `7` = "John mac 3 cd163")
DimPlot(EpCAM_subset_F_D3, reduction = "umap", label = TRUE, pt.size = 1)

############################################################
# 12. D3 dataset - CellChat analysis
############################################################


data.input = EpCAM_subset_F_D3

data.input <- GetAssayData(EpCAM_subset_F_D3, assay = "RNA", slot ="data")
labels <- Idents(EpCAM_subset_F_D3)
meta <- data.frame(group = labels, row.names = names(labels))

EpCAM_subset_cellchat_D3 <- createCellChat(object = data.input, meta = meta, group.by = "group")
EpCAM_subset_cellchat_D3 <- updateCellChat(EpCAM_subset_cellchat_D3)

CellChatDB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
EpCAM_subset_cellchat_D3@DB <- CellChatDB.use
EpCAM_subset_cellchat_D3 <- subsetData(EpCAM_subset_cellchat_D3)

EpCAM_subset_cellchat_D3 <- identifyOverExpressedGenes(EpCAM_subset_cellchat_D3)
EpCAM_subset_cellchat_D3 <- identifyOverExpressedInteractions(EpCAM_subset_cellchat_D3)
EpCAM_subset_cellchat_D3 <- projectData(EpCAM_subset_cellchat_D3, PPI.mouse)
EpCAM_subset_cellchat_D3 <- computeCommunProb(EpCAM_subset_cellchat_D3)
EpCAM_subset_cellchat_D3 <- filterCommunication(EpCAM_subset_cellchat_D3, min.cells = 10)
EpCAM_subset_cellchat_D3 <- netAnalysis_computeCentrality(EpCAM_subset_cellchat_D3)

df.net_D3 <- subsetCommunication(EpCAM_subset_cellchat_D3)

writexl::write_xlsx(df.net_D3, "output/df.net.D3.xlsx")


EpCAM_subset_cellchat_D3 <- computeCommunProbPathway(EpCAM_subset_cellchat_D3)
EpCAM_subset_cellchat_D3 <- aggregateNet(EpCAM_subset_cellchat_D3)
groupSize <- as.numeric(table(EpCAM_subset_cellchat_D3@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(EpCAM_subset_cellchat_D3@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(EpCAM_subset_cellchat_D3@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


############################################################
# 13. CellChat comparison between timepoints
############################################################

#### Integrating the cell chat objects

cellchat.control <- EpCAM_subset_cellchat
cellchat.D3irradiated <- EpCAM_subset_cellchat_D3
object.list <- list(control = cellchat.control, D3irradiated = cellchat.D3irradiated)

cellchat_merged <- mergeCellChat(object.list, cell.prefix = TRUE, add.names = names(object.list))



############################################################
# 14. Manuscript figure generation
############################################################


### Creating the cell chat plots

## Plot showing differential number of interactions or interaction strengths among different cell populations


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = T)
netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight")

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}



## Plot showing differential number of interactions or interaction strengths using a heatmap 

gg1 <- netVisual_heatmap(cellchat_merged)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


## Plot comparing the overall information flow of each signalling pathway 

gg1 <- rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2



# Chord diagram showing CX3CL1-CX3CR1 signalling pathways in control (M) and day 3 post IR
pathways.show <- c("CX3C") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


# Diagram showing the ligand-receptor pairs between Krt14 cells nad macrophage populations CD14+ and CD81+


# Read the Excel file
df <- read_excel("input/KRT14-mac_pathways.xlsx")

# Create a ligand-receptor pair name
df <- df %>%
  mutate(Pair = paste(Ligand, Receptor, sep = "_"))


# Continuous heatmap (p-values)

ggplot(df, aes(x = Condition, y = Pair, fill = -log10(Probability))) +
  geom_tile(color = "grey70") +
  scale_fill_gradient(low = "white", high = "red",
                      name = "-log10(p)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Ligand–Receptor Interactions",
       x = "Condition",
       y = "Ligand–Receptor Pair")


# Binary heatmap (p < 0.05)

df <- df %>%
  mutate(Significant = ifelse(Probability < 0.05, 1, 0))

ggplot(df, aes(x = Condition, y = Pair, fill = factor(Significant))) +
  geom_tile(color = "grey70") +
  scale_fill_manual(values = c("white", "blue"),
                    name = "Significance",
                    labels = c("Not Sig", "p < 0.05")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Binary Heatmap of Significant Interactions",
       x = "Condition",



