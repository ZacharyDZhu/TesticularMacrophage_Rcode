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
base_dir <- "E:/Github Data/GEO_MiceMacrophage_P12C5"
groups <- c("P12C5", "Scramble")
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

##Merge_QC
seurat_merge$log10GenesPerUMI <- log10(seurat_merge$nFeature_RNA) / log10(seurat_merge$nCount_RNA)
seurat_merge[["percent.mt"]] <- PercentageFeatureSet(seurat_merge, pattern = "^mt-")
VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(seurat_merge@meta.data$orig.ident)
seurat_merge <- subset(seurat_merge,  subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & nCount_RNA > 500 & log10GenesPerUMI > 0.7 & percent.mt < 20)
table(seurat_merge@meta.data$orig.ident)

##NormalizingDataSUM
seurat_merge <- NormalizeData(seurat_merge)
seurat_merge <- FindVariableFeatures(seurat_merge, nfeatures = 2000)
seurat_merge <- ScaleData(seurat_merge)
##RunningtSNE
seurat_merge <- RunPCA(seurat_merge,npcs = 40)
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:20)
seurat_merge <- FindClusters(seurat_merge, resolution = 0.5, verbose = FALSE)
seurat_merge <- RunTSNE(seurat_merge, dims = 1:20)
mycol<- brewer.pal(20, 'Dark2')
mycol0<- colorRampPalette(mycol)(5)
DimPlot(seurat_merge, reduction = "tsne", pt.size = 1, label = T) 
DimPlot(seurat_merge, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 1)
DimPlot(seurat_merge, reduction = "tsne", label = FALSE, repel = TRUE, pt.size = 2, alpha=1) + 
  scale_color_manual(values = mycol0) + 
  scale_fill_manual(values = mycol0)

plot3 <- DimPlot(seurat_merge, reduction = "tsne", label = FALSE, repel = TRUE, pt.size = 2, alpha=1) + 
  scale_color_manual(values = mycol0) + 
  scale_fill_manual(values = mycol0)
plot3
ggsave("C57BL6_Mice_12C5_Testis_Macrophage_Barcode_Plots.tif",plot = plot3, width = 5, height = 4, dpi = 300)
plot4 <- FeaturePlot(seurat_merge, features = c("Adgre1","H2-Aa","Mrc1","Cmklr1"), cols =c("lightgray","red"),pt.size = 2,ncol=2)
plot4
ggsave("C57BL6_Mice_12C5_Testis_Macrophage_Feature_Plots.tif",plot = plot4, width = 10, height = 10, dpi = 300)

##
MHCIIHiMacrophage <- subset(seurat_merge, idents = ("1"))
CD206HiMacrophage <- subset(seurat_merge, idents = ("0"))
my_comparisons <- list(c("P12C5","Scramble"))
a1 <- AverageExpression(MHCIIHiMacrophage, group.by = "orig.ident")
a2 <- AverageExpression(CD206HiMacrophage, group.by = "orig.ident")
write.csv(a1,"C57BL6Mice_Testis_P12C5_MHCIIHi_AverageExpression.csv")
write.csv(a2,"C57BL6Mice_Testis_P12C5_CD206Hi_AverageExpression.csv")

##MacrophageInflammation
plot6 <- VlnPlot(MHCIIHiMacrophage,"Cd163", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+ 
  ylim(0, 6) + scale_fill_manual(values=c("#19a0ea","#bbc6ce")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("Scramble","P12C5")) 
plot6
plot7 <- VlnPlot(CD206HiMacrophage,"Stat3", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+ 
  ylim(0, 6) + scale_fill_manual(values=c("#19a0ea","#bbc6ce")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("Scramble","P12C5")) 
plot7
merged_plot1 <- grid.arrange(plot6,plot7,nrow=1)
merged_plot1
ggsave("C57BL6_Mice_12C5_Testis_Macrophage_Il1b.tif",plot = merged_plot1, width = 4, height = 4, dpi = 300)

AverageExpression(MHCIIHiMacrophage,features = "Pfkl", group.by = "orig.ident")