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
base_dir <- "E:/Github Data/GEO_MiceMacrophage_CMRKO"
groups <- c("WT", "KO")
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

##QualityControl
seurat_merge[["percent.mt"]] <- PercentageFeatureSet(seurat_merge, pattern = "^mt-")
seurat_merge$log10GenesPerUMI <- log10(seurat_merge$nFeature_RNA) / log10(seurat_merge$nCount_RNA)
VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(seurat_merge@meta.data$orig.ident)
seurat_merge <- subset(seurat_merge, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & log10GenesPerUMI > 0.7 & percent.mt < 20)
table(seurat_merge@meta.data$orig.ident)

##Normalization_Merge
seurat_merge <- NormalizeData(seurat_merge)
seurat_merge <- FindVariableFeatures(seurat_merge, nfeatures = 2000)
seurat_merge <- ScaleData(seurat_merge)
seurat_merge <- RunPCA(seurat_merge,npcs = 40)
seurat_merge <- RunHarmony(seurat_merge, group.by.vars = "orig.ident")
seurat_merge <- seurat_merge %>% 
  RunTSNE(reduction = "harmony", dims = 1:20, verbose = F)
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:20)
seurat_merge <- FindClusters(seurat_merge, resolution = 0.08, verbose = FALSE)
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

##Clustering
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:20)
seurat_merge <- FindClusters(seurat_merge, resolution = 0.5, verbose = FALSE)
seurat_merge <- RunTSNE(seurat_merge, dims = 1:20)
mycol<- brewer.pal(12, 'Set3')
mycol0<- colorRampPalette(mycol)(5)
plot5 <- DimPlot(seurat_merge, reduction = "tsne", label = FALSE, repel = TRUE, pt.size = 2, alpha=1) + 
  scale_color_manual(values = mycol0) + 
  scale_fill_manual(values = mycol0)
ggsave("C57BL6_Adult_Testis_Barcodeplot.tif", plot = plot5, width = 5, height = 4, dpi = 300)

Cell_Markers <- c("Col1a1","Dcn","Adgre1","Mrc1","H2-Aa")
DotPlot(seurat_merge, group.by = 'seurat_clusters',
        features = unique(Cell_Markers),cols = c("#ffffff", "#448444"), dot.scale=6) + RotatedAxis() + 
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20)) + labs(x='',y='')
Seurat_Macrophage <- subset(seurat_merge, ident=c("0"))
Macrophage_Markers <- c("Adgre1","Mrc1","H2-Aa")
plot6 <- DotPlot(Seurat_Macrophage, group.by = 'seurat_clusters',
        features = unique(Macrophage_Markers),cols = c("#ffffff", "#448444"), dot.scale=6) + RotatedAxis() + 
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20)) + labs(x='',y='')
plot6
ggsave("C57BL6_Adult_Testis_Featureplot.tif", plot = plot6, width = 5, height = 4, dpi = 300)

FeaturePlot(seurat_merge, features = "Cmklr1", cols =c("lightgrey","red"),pt.size = 2)

## Violin Plot and Group Comparison (e.g. Lipa)
A <- AverageExpression(Seurat_Macrophage, group.by = "orig.ident")
write.csv(A,"Macrophage Average Expression.csv")
my_comparisons <- list(c("WT", "KO"))
plot6 <- VlnPlot(Seurat_Macrophage,"Lipe", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#1B9E77","#D95F02")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("WT", "KO")) 
plot6
ggsave("C57BL6_Adult_Testis_Cmklr1.tif", plot = plot6, width = 3, height = 4, dpi = 300)
