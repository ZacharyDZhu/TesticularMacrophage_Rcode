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
plot1 <- DimPlot(object = seurat_merge, reduction = "harmony", 
                 pt.size = .1, group.by = "orig.ident") + NoLegend()
plot2 <- VlnPlot(object = seurat_merge, features = "harmony_1", 
                 group.by = "orig.ident", pt.size = .1) + NoLegend()
plot1 + plot2
plot3 <- DimPlot(seurat_merge,reduction = "tsne") 
plot4 <- DimPlot(seurat_merge, reduction = "tsne", 
                 group.by = "orig.ident", pt.size = .1, 
                 split.by = 'orig.ident') + NoLegend()
plot3 + plot4

Human_Macrophage <- subset(seurat_merge,ident=c("0","3"))
DimPlot(Human_Macrophage, reduction = "tsne", pt.size = 1) 
plot5 <- FeaturePlot(Human_Macrophage,features = c("CD68","CD67","CMKLR1","CCRL2","GPR1"),cols = c("lightgrey", "red"),pt.size = 2)
plot5

##DistinguishMacrophage
plot6 <- VlnPlot(Human_Macrophage,features = "FCGR1A",
                 cols = c("0" = "#E78AC3", "3" = "#4B9CD3"),pt.size=0)  + 
  geom_boxplot(fill="white", alpha = 0.9, width=0.1)
plot6
plot7 <- VlnPlot(Human_Macrophage,features = "CD14",
                 cols = c("0" = "#E78AC3", "3" = "#4B9CD3"),pt.size=0)  + 
  geom_boxplot(fill="white", alpha = 0.9, width=0.1)
plot7
plot8 <- VlnPlot(Human_Macrophage,features = "CEACAM8",
                 cols = c("0" = "#E78AC3", "3" = "#4B9CD3"),pt.size=0)  + 
  geom_boxplot(fill="white", alpha = 0.9, width=0.1)
plot8
plot9 <- grid.arrange(plot6,plot7,plot8,nrow = 2)
ggsave("Human_Aging_Testis_CMKLR1Macrophage_Markers_Plots.tif.tif", plot = plot5, width = 8, height = 4, dpi = 300)

define_color1 <- c("#E78AC3","#4B9CD3")  
plot10<- DimPlot(Human_Macrophage, reduction = 'tsne', group.by = 'seurat_clusters',
                label = F, pt.size = 1.5, alpha= 1) + scale_color_manual(values = define_color1) 
plot10
ggsave("Human_Aging_Testis_CMKLR1Macrophage_Barcode_Plots.tif", plot = plot6, width = 10, height = 5, dpi = 300)

##SubclusteringData
Human_Macrophage <- JoinLayers(Human_Macrophage)
CPos_Macrophage <- subset(Human_Macrophage, ident=c("0"))
CNeg_Macrophage <- subset(Human_Macrophage, ident=c("3"))

##Polarization Markers Comparison
my_comparisons <- list( c("YG", "OG1"), c("YG", "OG2"),c("OG1","OG2"))
plot11 <- VlnPlot(CPos_Macrophage,"STAT3", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot11
plot12 <- VlnPlot(CPos_Macrophage,"FCGR1A", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot12
plot13 <- VlnPlot(CPos_Macrophage,"CD163", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot13
plot14 <- VlnPlot(CNeg_Macrophage,"STAT3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot14
plot15 <- VlnPlot(CNeg_Macrophage,"FCGR1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot15
plot16 <- VlnPlot(CNeg_Macrophage,"CD163", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot16
merged_plot2 <- grid.arrange(plot11,plot12,plot13,plot14,plot15,plot16,nrow=1)
ggsave("Human_Aging_Testis_CPosMacrophage_Polarization.tif", plot = merged_plot2, width = 12, height = 4, dpi = 300)

##Recruiting Neutrophil Markers Comparison
plot17 <- VlnPlot(CPos_Macrophage,"CCL2", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot17
plot18 <- VlnPlot(CNeg_Macrophage,"CCL2", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot18
merged_plot3 <- grid.arrange(plot17,plot18,nrow=1)
ggsave("Human_Aging_Testis_CPosMacrophage_RecruitingNeutrophil.tif", plot = merged_plot3, width = 4, height = 4, dpi = 300)

##Glycolysis Markers Comparison
plot19 <- VlnPlot(CPos_Macrophage,"LDHA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot19
plot20 <- VlnPlot(CNeg_Macrophage,"LDHA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot20
merged_plot4 <- grid.arrange(plot19,plot20,nrow=1)
ggsave("Human_Aging_Testis_CPosMacrophage_Glycolysis.tif", plot = merged_plot4, width = 4, height = 4, dpi = 300)

##OXPHOS Markers Comparison
plot21 <- VlnPlot(CPos_Macrophage,"NDUFA3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot21
plot22 <- VlnPlot(CNeg_Macrophage,"NDUFA3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot22
merged_plot5 <- grid.arrange(plot21,plot22,nrow=1)
ggsave("Human_Aging_Testis_CPosMacrophage_OXPHOS.tif", plot = merged_plot5, width = 4, height = 4, dpi = 300)

##Cmklr1MacrophageSubclustering
CPos_Macrophage <- NormalizeData(CPos_Macrophage,verbose = F)
CPos_Macrophage <- FindVariableFeatures(CPos_Macrophage, nfeatures = 2000)
CPos_Macrophage <- ScaleData(CPos_Macrophage)
CPos_Macrophage <- RunPCA(CPos_Macrophage)
ElbowPlot(CPos_Macrophage, ndims = ncol(Embeddings(CPos_Macrophage, "pca")))
CPos_Macrophage <- RunPCA(CPos_Macrophage, npcs = 30, verbose = F)
CPos_Macrophage <- FindNeighbors(CPos_Macrophage, dims = 1:20, verbose = FALSE)
CPos_Macrophage <- FindClusters(CPos_Macrophage, resolution = 0.05, verbose = FALSE)
CPos_Macrophage <- RunTSNE(CPos_Macrophage, dims = 1:20)
DimPlot(CPos_Macrophage, reduction = "tsne", pt.size = 1) 
CPos_Macrophage <- subset(CPos_Macrophage,ident=c("0","1"))
DimPlot(CPos_Macrophage, reduction = "tsne", pt.size = 1) 
plot23<- DimPlot(CPos_Macrophage, reduction = 'tsne',
                label = F, pt.size = 1.5, alpha= 0.7, cols = c("#a6caeb","#cc3300","#0070c0")) 
plot23
table(CPos_Macrophage$orig.ident)
ggsave("Human_Aging_Testis_CPosMacrophage_BarcodePlots.tif", plot = plot23, width = 6, height = 5, dpi = 300)
getwd()
Macrophage_Markers <- c("CMKLR1","MRC1","HLA-DOA")
plot24 <- FeaturePlot(CPos_Macrophage,
                     features = Macrophage_Markers,
                     reduction = "tsne",       
                     cols = c("lightgrey", "red"),  
                     pt.size = 2,
                     order = TRUE)       
plot24
ggsave("Human_Aging_Testis_CPosMacrophage_Markers_Featureplots.tif", plot = plot24, width = 8, height = 8, dpi = 300)

plot25 <- VlnPlot(CPos_Macrophage,features = "MRC1",
                 cols = c("0" = "#E31A1C", "1" = "#33A02C"),pt.size=0)  + 
  geom_boxplot(fill="white", alpha = 0.9, width=0.1)
plot25
plot26 <- VlnPlot(CPos_Macrophage,features = "HLA-DOA",
                 cols = c("0" = "#E31A1C", "1" = "#33A02C"),pt.size=0)  + 
  geom_boxplot(fill="white", alpha = 0.9, width=0.1)
plot26
plot6 <- grid.arrange(plot25,plot26,nrow = 1)
ggsave("Human_Aging_Testis_CPosMacrophage_Markers_VlnPlot.tif", plot = plot6, width = 4, height = 4, dpi = 300)

CPos_Macrophage <- JoinLayers(CPos_Macrophage)
CPosMHCIIHi_Macrophage <- subset(CPos_Macrophage, ident=c("0"))
CPosCD206Hi_Macrophage <- subset(CPos_Macrophage, ident=c("1"))

CPosMHCIIHiOG1vsYG <- FindMarkers(CPosMHCIIHi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG1",ident.2 = "YG")
CPosMHCIIHiOG2vsYG <- FindMarkers(CPosMHCIIHi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG2",ident.2 = "YG")
CPosMHCIIHiOG2vsOG1 <- FindMarkers(CPosMHCIIHi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG2",ident.2 = "OG1")

CPosCD206HiOG1vsYG <- FindMarkers(CPosCD206Hi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG1",ident.2 = "YG")
CPosCD206HiOG2vsYG <- FindMarkers(CPosCD206Hi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG2",ident.2 = "YG")
CPosCD206HiOG2vsOG1 <- FindMarkers(CPosCD206Hi_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "OG2",ident.2 = "OG1")

##Write CSV
write.csv(CPosMHCIIHiOG1vsYG,"CPosMHCIIHi_Macrophage_OG1vsYG_DEG.csv")
write.csv(CPosMHCIIHiOG2vsYG,"CPosMHCIIHi_Macrophage_OG2vsYG_DEG.csv")
write.csv(CPosMHCIIHiOG2vsOG1,"CPosMHCIIHi_Macrophage_OG2vsOG1_DEG.csv")

##Write CSV
write.csv(CPosCD206HiOG1vsYG,"CPosCD206Hi_Macrophage_OG1vsYG_DEG.csv")
write.csv(CPosCD206HiOG2vsYG,"CPosCD206Hi_Macrophage_OG2vsYG_DEG.csv")
write.csv(CPosCD206HiOG2vsOG1,"CPosCD206Hi_Macrophage_OG2vsOG1_DEG.csv")

cut_off_P=0.05
cut_off_log2FC=0.585

CPosMHCIIHiOG1vsYG$Sig = ifelse(CPosMHCIIHiOG1vsYG$p_val_adj < cut_off_P &
                                  abs(CPosMHCIIHiOG1vsYG$avg_log2FC) >= cut_off_log2FC,  
                                ifelse(CPosMHCIIHiOG1vsYG$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
CPosMHCIIHiOG1vsYG = data.frame(CPosMHCIIHiOG1vsYG)
table(CPosMHCIIHiOG1vsYG$Sig)

CPosMHCIIHiOG2vsYG$Sig = ifelse(CPosMHCIIHiOG2vsYG$p_val_adj < cut_off_P &
                                  abs(CPosMHCIIHiOG2vsYG$avg_log2FC) >= cut_off_log2FC,  
                                ifelse(CPosMHCIIHiOG2vsYG$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
CPosMHCIIHiOG2vsYG = data.frame(CPosMHCIIHiOG2vsYG)
table(CPosMHCIIHiOG2vsYG$Sig)

CPosCD206HiOG1vsYG$Sig = ifelse(CPosCD206HiOG1vsYG$p_val_adj < cut_off_P &
                                  abs(CPosCD206HiOG1vsYG$avg_log2FC) >= cut_off_log2FC,  
                                ifelse(CPosCD206HiOG1vsYG$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
CPosCD206HiOG1vsYG = data.frame(CPosCD206HiOG1vsYG)
table(CPosCD206HiOG1vsYG$Sig)

CPosCD206HiOG2vsYG$Sig = ifelse(CPosCD206HiOG2vsYG$p_val_adj < cut_off_P &
                                  abs(CPosCD206HiOG2vsYG$avg_log2FC) >= cut_off_log2FC,  
                                ifelse(CPosCD206HiOG2vsYG$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
CPosCD206HiOG2vsYG = data.frame(CPosCD206HiOG2vsYG)
table(CPosCD206HiOG2vsYG$Sig)


library(EnhancedVolcano)
a1 <- EnhancedVolcano(CPosMHCIIHiOG1vsYG,
                      lab = rownames(CPosMHCIIHiOG1vsYG),
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      boxedLabels = FALSE,
                      drawConnectors = TRUE,
                      pCutoff = 0.05,
                      FCcutoff = 0.585,
                      title = "CPosMHCIIHiOG1vsYG",
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

a2 <- EnhancedVolcano(CPosMHCIIHiOG2vsYG,
                      lab = rownames(CPosMHCIIHiOG2vsYG),
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      boxedLabels = FALSE,
                      drawConnectors = TRUE,
                      pCutoff = 0.05,
                      FCcutoff = 0.585,
                      title = "CPosMHCIIHiOG2vsYG",
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
a2

a3 <- EnhancedVolcano(CPosCD206HiOG1vsYG,
                      lab = rownames(CPosCD206HiOG1vsYG),
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      boxedLabels = FALSE,
                      drawConnectors = TRUE,
                      pCutoff = 0.05,
                      FCcutoff = 0.585,
                      title = "CPosCD206HiOG1vsYG",
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

a4 <- EnhancedVolcano(CPosCD206HiOG2vsYG,
                      lab = rownames(CPosCD206HiOG2vsYG),
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      boxedLabels = FALSE,
                      drawConnectors = TRUE,
                      pCutoff = 0.05,
                      FCcutoff = 0.585,
                      title = "CPosCD206HiOG2vsYG",
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
ggsave("Human_Aging_Testis_CPosMHCIIHiOG1vsYG_VolcanoPlots.tif", plot = a1, width = 10, height = 10, dpi = 300)
ggsave("Human_Aging_Testis_CPosMHCIIHiOG2vsYG_VolcanoPlots.tif", plot = a2, width = 10, height = 10, dpi = 300)
ggsave("Human_Aging_Testis_CPosCD206HiOG1vsYG_VolcanoPlots.tif", plot = a3, width = 10, height = 10, dpi = 300)
ggsave("Human_Aging_Testis_CPosCD206HiOG2vsYG_VolcanoPlots.tif", plot = a4, width = 10, height = 10, dpi = 300)

#Inflammation
my_comparisons <- list( c("YG", "OG1"), c("YG", "OG2"),c("OG1","OG2"))
plot27 <- VlnPlot(CPosMHCIIHi_Macrophage,"IL1B", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot27
plot28 <- VlnPlot(CPosMHCIIHi_Macrophage,"TGFB1", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot28
merged_plot <- grid.arrange(plot27,plot28,nrow=1)
ggsave("Human_Aging_Testis_CPosMHCIIHi_Inflammation.tif", plot = merged_plot, width = 4, height = 4, dpi = 300)

#Inflammation
plot29 <- VlnPlot(CPosCD206Hi_Macrophage,"IL1B", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot29
plot30 <- VlnPlot(CPosCD206Hi_Macrophage,"TGFB1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot30
merged_plot2 <- grid.arrange(plot29,plot30,nrow=1)
ggsave("Human_Aging_Testis_CPosCD206Hi_Inflammation.tif", plot = merged_plot2, width = 4, height = 4, dpi = 300)

#MHCIIHiGlycolysis
plot31 <- VlnPlot(CPosMHCIIHi_Macrophage,"SLC2A3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot31
plot32 <- VlnPlot(CPosMHCIIHi_Macrophage,"HK1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot32
plot33 <- VlnPlot(CPosMHCIIHi_Macrophage,"PFKL", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot33
plot34 <- VlnPlot(CPosMHCIIHi_Macrophage,"PKM", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot34
plot35 <- VlnPlot(CPosMHCIIHi_Macrophage,"LDHA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot35
merged_plot3 <- grid.arrange(plot31,plot32,plot33,plot34,plot35,nrow=1)
ggsave("Human_Aging_Testis_CPosMHCIIHi_Glycolysis.tif", plot = merged_plot3, width = 10, height = 4, dpi = 300)

#CD206HiGlycolysis
plot36 <- VlnPlot(CPosCD206Hi_Macrophage,"SLC2A3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("YG", "OG1","OG2")) 
plot36
plot37 <- VlnPlot(CPosCD206Hi_Macrophage,"HK1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("YG", "OG1","OG2")) 
plot37
plot38 <- VlnPlot(CPosCD206Hi_Macrophage,"PFKL", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("YG", "OG1","OG2")) 
plot38
plot39 <- VlnPlot(CPosCD206Hi_Macrophage,"PKM", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("YG", "OG1","OG2")) 
plot39
plot40 <- VlnPlot(CPosCD206Hi_Macrophage,"LDHA", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("YG", "OG1","OG2")) 
plot40
merged_plot4 <- grid.arrange(plot36,plot37,plot38,plot39,plot40,nrow=1)
ggsave("Human_Aging_Testis_CPosCD206Hi_Glycolysis.tif", plot = merged_plot4, width = 10, height = 4, dpi = 300)

##Polarization
plot41 <- VlnPlot(CPosMHCIIHi_Macrophage,"FCGR1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot41
plot42 <- VlnPlot(CPosMHCIIHi_Macrophage,"CD163", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot42
merged_plot5 <- grid.arrange(plot41,plot42,nrow=1)
ggsave("Human_Aging_Testis_CPosMHCIIHi_Polarization.tif", plot = merged_plot5, width = 4, height = 4, dpi = 300)

##Polarization
plot43 <- VlnPlot(CPosCD206Hi_Macrophage,"FCGR1A", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot43
plot44 <- VlnPlot(CPosCD206Hi_Macrophage,"CD163", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot44
merged_plot6 <- grid.arrange(plot43,plot44,nrow=1)
ggsave("Human_Aging_Testis_CPosCD206Hi_Polarization.tif", plot = merged_plot6, width = 4, height = 4, dpi = 300)

##FAO
plot45 <- VlnPlot(CPosMHCIIHi_Macrophage,"CD36", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot45
plot46 <- VlnPlot(CPosCD206Hi_Macrophage,"CD36", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#a6caeb","#cc3300","#0070c0")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c( "YG", "OG1","OG2")) 
plot46
merged_plot7 <- grid.arrange(plot45,plot46,nrow=1)
ggsave("Human_Aging_Testis_CPos_FAO.tif", plot = merged_plot7, width = 4, height = 4, dpi = 300)

