rm(list = ls())
##RunningEnvironment
library(future)
library(data.table)
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)
library(stringr)
library(tidyverse)
library(tidyr)
library(harmony)
library(Rcpp)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggrastr)
library(ggpubr)
library(gridExtra)
options(future.globals.maxSize = 8000 * 8000^2)
plan(multisession, workers=6)
setwd("E:/Github Data")

##Import Exercise Data
base_dir <- "E:/Github Data/GEO_HumanMacrophage_Aging"
groups <- c("YG", "OG1","OG2")
seurat_list <- list()
for (grp in groups) {
  message("Importing group: ", grp)
  dir_path <- file.path(base_dir, grp)
  matrix_file <- file.path(dir_path, "matrix.mtx")
  barcodes_file <- file.path(dir_path, "barcodes.tsv")
  features_file <- file.path(dir_path, "features.tsv")
  if (!file.exists(matrix_file) || !file.exists(barcodes_file) || !file.exists(features_file)) {
    warning("Missing files in group: ", grp)
    next
  }
  counts <- readMM(matrix_file)
  barcodes <- readLines(barcodes_file)
  features <- read.delim(features_file, header = FALSE, stringsAsFactors = FALSE)
  rownames(counts) <- features[,1]
  colnames(counts) <- barcodes
  seu <- CreateSeuratObject(counts = counts, project = grp)
  seu$group <- grp
  seurat_list[[grp]] <- seu
}
if (length(seurat_list) >= 2) {
  seurat_merge <- merge(seurat_list[[1]], y = seurat_list[-1], 
                        add.cell.ids = groups[1:length(seurat_list)], 
                        project = "MacrophageExercise")
} else if (length(seurat_list) == 1) {
  seurat_merge <- seurat_list[[1]]
} else {
  stop("No Seurat objects loaded.")
}
seurat_merge

##MergingDataQC
seurat_merge$log10GenesPerUMI <- log10(seurat_merge$nFeature_RNA) / log10(seurat_merge$nCount_RNA)
seurat_merge[["percent.mt"]] <- PercentageFeatureSet(seurat_merge, pattern = "^mt-")
VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(seurat_merge@meta.data$orig.ident)
seurat_merge <- subset(seurat_merge, subset = nFeature_RNA > 200  & percent.mt < 50)
table(seurat_merge@meta.data$orig.ident)

##NormalizingData
seurat_merge<- NormalizeData(seurat_merge,verbose = F)
seurat_merge<- FindVariableFeatures(seurat_merge, nfeatures = 2000)
seurat_merge <- ScaleData(seurat_merge)

##PCA
seurat_merge <- RunPCA(seurat_merge,npcs = 30)
seurat_merge <- RunHarmony(seurat_merge, group.by.vars = "orig.ident")
seurat_merge <- seurat_merge %>% 
  RunTSNE(reduction = "harmony", dims = 1:20, verbose = F)
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:20)
seurat_merge <- FindClusters(seurat_merge, resolution = 0.04, verbose = FALSE)
seurat_merge <- RunTSNE(seurat_merge, dims = 1:20)
DimPlot(seurat_merge, reduction = 'tsne', group.by = 'seurat_clusters',
                     label = F, pt.size = 1.5, alpha= 1)
Main_Macrophage <- subset(seurat_merge, ident=c("0"))

##MainMacrophageSubclustering
Main_Macrophage <- NormalizeData(Main_Macrophage,verbose = F)
Main_Macrophage <- FindVariableFeatures(Main_Macrophage, nfeatures = 2000)
Main_Macrophage <- ScaleData(Main_Macrophage)
Main_Macrophage <- RunPCA(Main_Macrophage)
ElbowPlot(Main_Macrophage, ndims = ncol(Embeddings(Main_Macrophage, "pca")))
Main_Macrophage <- RunPCA(Main_Macrophage, npcs = 30, verbose = F)
Main_Macrophage <- FindNeighbors(Main_Macrophage, dims = 1:20, verbose = FALSE)
Main_Macrophage <- FindClusters(Main_Macrophage, resolution = 0.05, verbose = FALSE)
Main_Macrophage <- RunTSNE(Main_Macrophage, dims = 1:20)
DimPlot(Main_Macrophage, reduction = "tsne", pt.size = 1) 
plot1<- DimPlot(Main_Macrophage, reduction = 'tsne', group.by = 'orig.ident',
                label = F, pt.size = 1.5, alpha= 0.7, cols = c("#a6caeb","#cc3300","#0070c0")) 
plot1
table(Main_Macrophage$orig.ident)

ggsave("Human_Aging_Testis_Macrophage_Subclustering_Plots.tif", plot = plot1, width = 6, height = 5, dpi = 600)
Macrophage_Markers <- c("CD68","MRC1","HLA-DOA")
plot2 <- FeaturePlot(Main_Macrophage,
                     features = Macrophage_Markers,
                     reduction = "tsne",       
                     cols = c("lightgrey", "red"),  
                     pt.size = 1.5,
                     order = TRUE)       
plot2
ggsave("Human_Testis_CanonicalMacrophage_Markers.tif", plot = plot2, width = 8, height = 8, dpi = 600)

plot3 <- VlnPlot(Main_Macrophage,features = "MRC1",
                 cols = c("0" = "lightcoral", "1" = "mediumaquamarine"),pt.size=0)  + 
  geom_boxplot(fill="white", alpha = 0.9, width=0.1)
plot3
plot4 <- VlnPlot(Main_Macrophage,features = "HLA-DOA",
                 cols = c("0" = "lightcoral", "1" = "mediumaquamarine"),pt.size=0)  + 
  geom_boxplot(fill="white", alpha = 0.9, width=0.1)
plot4
plot5 <- grid.arrange(plot3,plot4,nrow = 1)
ggsave("Human_Aging_Testis_Macrophage_Markers_VlnPlot.tif", plot = plot5, width = 6, height = 5, dpi = 600)

##Find Average Expression Genes
Main_Macrophage <- JoinLayers(Main_Macrophage)
MHCIIHi_Macrophage <- subset(Main_Macrophage, ident=c("0"))
CD206Hi_Macrophage <- subset(Main_Macrophage, ident=c("1"))

CD206HiOG1vsYG <- FindMarkers(CD206Hi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG1",ident.2 = "YG")
CD206HiOG2vsYG <- FindMarkers(CD206Hi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG2",ident.2 = "YG")
CD206HiOG2vsOG1 <- FindMarkers(CD206Hi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG2",ident.2 = "OG1")

MHCIIHiOG1vsYG <- FindMarkers(MHCIIHi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG1",ident.2 = "YG")
MHCIIHiOG2vsYG <- FindMarkers(MHCIIHi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG2",ident.2 = "YG")
MHCIIHiOG2vsOG1 <- FindMarkers(MHCIIHi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG2",ident.2 = "OG1")

##Write CSV
write.csv(MHCIIHiOG1vsYG,"MHCIIHi_Macrophage_OG1vsYG_DEG.csv")
write.csv(MHCIIHiOG2vsYG,"MHCIIHi_Macrophage_OG2vsYG_DEG.csv")
write.csv(MHCIIHiOG2vsOG1,"MHCIIHi_Macrophage_OG2vsOG1_DEG.csv")

##Write CSV
write.csv(CD206HiOG1vsYG,"CD206Hi_Macrophage_OG1vsYG_DEG.csv")
write.csv(CD206HiOG2vsYG,"CD206Hi_Macrophage_OG2vsYG_DEG.csv")
write.csv(CD206HiOG2vsOG1,"CD206Hi_Macrophage_OG2vsOG1_DEG.csv")

##DEG Analysis
head(MHCIIHiOG1vsYG)
head(MHCIIHiOG2vsYG)
cut_off_P =0.05
cut_off_log2FC =0.585
MHCIIHiOG1vsYG$Sig = ifelse(MHCIIHiOG1vsYG$p_val_adj < cut_off_P &
                              abs(MHCIIHiOG1vsYG$avg_log2FC) >= cut_off_log2FC,  
                            ifelse(MHCIIHiOG1vsYG$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
MHCIIHiOG1vsYG = data.frame(MHCIIHiOG1vsYG)
table(MHCIIHiOG1vsYG$Sig)
MHCIIHiOG2vsYG$Sig = ifelse(MHCIIHiOG2vsYG$p_val_adj < cut_off_P &
                              abs(MHCIIHiOG2vsYG$avg_log2FC) >= cut_off_log2FC,  
                            ifelse(MHCIIHiOG2vsYG$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
MHCIIHiOG2vsYG = data.frame(MHCIIHiOG2vsYG)
table(MHCIIHiOG2vsYG$Sig)
CD206HiOG1vsYG$Sig = ifelse(CD206HiOG1vsYG$p_val_adj < cut_off_P &
                              abs(CD206HiOG1vsYG$avg_log2FC) >= cut_off_log2FC,  
                            ifelse(CD206HiOG1vsYG$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
CD206HiOG1vsYG = data.frame(CD206HiOG1vsYG)
table(CD206HiOG1vsYG$Sig)
CD206HiOG2vsYG$Sig = ifelse(CD206HiOG2vsYG$p_val_adj < cut_off_P &
                              abs(CD206HiOG2vsYG$avg_log2FC) >= cut_off_log2FC,  
                            ifelse(CD206HiOG2vsYG$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
CD206HiOG2vsYG = data.frame(CD206HiOG2vsYG)
table(CD206HiOG2vsYG$Sig)

a1 <- EnhancedVolcano(MHCIIHiOG1vsYG,
                      lab = rownames(MHCIIHiOG1vsYG),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 0.05,
                      FCcutoff = 0.585,
                      title = 'MHCIIHiOG1vsYG',
                      subtitle = NULL,
                      col = c(
                        'grey70',
                        'lightblue3',
                        'blue3',
                        'firebrick'))+ 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA))
a1
a2 <- EnhancedVolcano(MHCIIHiOG2vsYG,
                      lab = rownames(MHCIIHiOG2vsYG),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 0.05,
                      FCcutoff = 0.585,
                      title = 'MHCIIHiOG2vsYG',
                      subtitle = NULL,
                      max.overlaps = 1000,
                      col = c(
                        'grey70',
                        'lightblue3',
                        'blue3',
                        'firebrick'))+ 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA))
a2
a3 <- EnhancedVolcano(CD206HiOG1vsYG,
                      lab = rownames(CD206HiOG1vsYG),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 0.05,
                      FCcutoff = 0.585,
                      title = 'CD206HiOG1vsYG',
                      subtitle = NULL,
                      col = c(
                        'grey70',
                        'lightblue3',
                        'blue3',
                        'firebrick'))+ 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA))
a3
a4 <- EnhancedVolcano(CD206HiOG2vsYG,
                      lab = rownames(CD206HiOG2vsYG),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 0.05,
                      FCcutoff = 0.585,
                      title = 'CD206HiOG2vsYG',
                      subtitle = NULL,
                      col = c(
                        'grey70',
                        'lightblue3',
                        'blue3',
                        'firebrick'))+ 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA))
a4
ggsave("Human_Aging_Testis_CanonicalMacrophage_MHCIIHiOG1vsYG_VolcanoPlots.tif", plot = a1, width = 10, height = 10, dpi = 300)
ggsave("Human_Aging_Testis_CanonicalMacrophage_MHCIIHiOG2vsYG_VolcanoPlots.tif", plot = a2, width = 10, height = 10, dpi = 300)
ggsave("Human_Aging_Testis_CanonicalMacrophage_CD206HiOG1vsYG_VolcanoPlots.tif", plot = a3, width = 10, height = 10, dpi = 300)
ggsave("Human_Aging_Testis_CanonicalMacrophage_CD206HiOG2vsYG_VolcanoPlots.tif", plot = a4, width = 10, height = 10, dpi = 300)

my_comparisons <- list( c("YG", "OG1"), c("YG", "OG2"),c("OG1","OG2"))
#Polarization
plot10 <- VlnPlot(MHCIIHi_Macrophage,"STAT3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot10
plot11 <- VlnPlot(MHCIIHi_Macrophage,"FCGR1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot11
plot12 <- VlnPlot(MHCIIHi_Macrophage,"CD163", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot12
plot13 <- VlnPlot(CD206Hi_Macrophage,"STAT3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot13
plot14 <- VlnPlot(CD206Hi_Macrophage,"FCGR1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot14
plot15 <- VlnPlot(CD206Hi_Macrophage,"CD163", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot15
merged_plot2 <- grid.arrange(plot10,plot11,plot12,plot13,plot14,plot15,nrow=1)
ggsave("Human_Aging_Testis_CanonicalMacrophage_Polarization_Vlnplot.tif", plot = merged_plot2, width = 12, height = 4, dpi = 300)

##MHCIIHi Macrophage OXPHOS
plot16 <- VlnPlot(MHCIIHi_Macrophage,"HIF1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot16
plot17 <- VlnPlot(CD206Hi_Macrophage,"HIF1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot17
merged_plot3 <- grid.arrange(plot16,plot17,nrow=1)
ggsave("Human_Aging_Testis_CanonicalMacrophage_HIF1A.tif", plot = merged_plot3, width = 4, height = 4, dpi = 300)

##MHCIIHiGlYGolysis
plot18 <- VlnPlot(MHCIIHi_Macrophage,"SLC2A3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot18
plot19 <- VlnPlot(MHCIIHi_Macrophage,"HK1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot19
plot20 <- VlnPlot(MHCIIHi_Macrophage,"PFKL", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot20
plot21 <- VlnPlot(MHCIIHi_Macrophage,"PKM", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot21
plot22 <- VlnPlot(MHCIIHi_Macrophage,"LDHA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot22
merged_plot4 <- grid.arrange(plot18,plot19,plot20,plot21,plot22,nrow=1)
ggsave("Human_Aging_Testis_CanonicalMacrophage_MHCIIHiGlYGolysis.tif", plot = merged_plot4, width = 10, height = 4, dpi = 300)

##CD206HiGlYGolysis
plot23 <- VlnPlot(CD206Hi_Macrophage,"SLC2A3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot23
plot24 <- VlnPlot(CD206Hi_Macrophage,"HK1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot24
plot25 <- VlnPlot(CD206Hi_Macrophage,"PFKL", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot25
plot26 <- VlnPlot(CD206Hi_Macrophage,"PKM", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot26
plot27 <- VlnPlot(CD206Hi_Macrophage,"LDHA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot27
merged_plot5 <- grid.arrange(plot23,plot24,plot25,plot26,plot27,nrow=1)
ggsave("Human_Aging_Testis_CanonicalMacrophage_CD206HiGlYGolysis.tif", plot = merged_plot5, width = 10, height = 4, dpi = 300)

###MHCIIHiLipoMetabolic
plot28 <- VlnPlot(MHCIIHi_Macrophage,"CD36", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 5) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot28
plot29 <- VlnPlot(MHCIIHi_Macrophage,"LIPA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 5) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot29
plot30 <- VlnPlot(MHCIIHi_Macrophage,"LIPE", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 5) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot30
plot31 <- VlnPlot(MHCIIHi_Macrophage,"CPT1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot31
plot32 <- VlnPlot(MHCIIHi_Macrophage,"CPT1B", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot32
plot33 <- VlnPlot(MHCIIHi_Macrophage,"CPT2", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot33
plot34 <- VlnPlot(MHCIIHi_Macrophage,"ACLY", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot34
plot35 <- VlnPlot(MHCIIHi_Macrophage,"ACACA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot35
plot36 <- VlnPlot(MHCIIHi_Macrophage,"FASN", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot36
plot37 <- VlnPlot(MHCIIHi_Macrophage,"DGAT1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot37
merged_plot6 <- grid.arrange(plot28,plot29,plot30,plot31,plot32,plot33,plot34,plot35,plot36,plot37,nrow=2)
ggsave("Human_Aging_Testis_CanonicalMacrophage_MHCIIHi_LipidMetabolism.tif", plot = merged_plot6, width = 10, height = 8, dpi = 300)

###CD206HiLipoMetabolic
plot38 <- VlnPlot(CD206Hi_Macrophage,"CD36", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 5) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot38
plot39 <- VlnPlot(CD206Hi_Macrophage,"LIPA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 5) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot39
plot40 <- VlnPlot(CD206Hi_Macrophage,"LIPE", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 5) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot40
plot41 <- VlnPlot(CD206Hi_Macrophage,"CPT1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot41
plot42 <- VlnPlot(CD206Hi_Macrophage,"CPT1B", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot42
plot43 <- VlnPlot(CD206Hi_Macrophage,"CPT2", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot43
plot44 <- VlnPlot(CD206Hi_Macrophage,"ACLY", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot44
plot45 <- VlnPlot(CD206Hi_Macrophage,"ACACA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot45
plot46 <- VlnPlot(CD206Hi_Macrophage,"FASN", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot46
plot47 <- VlnPlot(CD206Hi_Macrophage,"DGAT1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot47
merged_plot7 <- grid.arrange(plot38,plot39,plot40,plot41,plot42,plot43,plot44,plot45,plot46,plot47,nrow=2)
ggsave("Human_Aging_Testis_CanonicalMacrophage_CD206Hi_LipidMetabolism.tif", plot = merged_plot7, width = 10, height = 8, dpi = 300)
AverageExpression(CD206Hi_Macrophage,features = "FASN", group.by = "orig.ident")