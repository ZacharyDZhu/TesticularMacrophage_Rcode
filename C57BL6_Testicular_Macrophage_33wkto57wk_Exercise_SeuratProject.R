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
base_dir <- "E:/GEO_MiceMacrophage_Exercise"
groups <- c("SDO", "MICTO", "HIITO")
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
seurat_merge$log10GenesPerUMI <- log10(seurat_merge$nFeature_RNA) / log10(seurat_merge$nCount_RNA)
seurat_merge[["percent.mt"]] <- PercentageFeatureSet(seurat_merge, pattern = "^mt-")
VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(seurat_merge@meta.data$orig.ident)
seurat_merge <- subset(seurat_merge, subset = nFeature_RNA > 500  & percent.mt < 20)
table(seurat_merge@meta.data$orig.ident)

##NormalizingData
seurat_merge <- NormalizeData(seurat_merge)
seurat_merge <- ScaleData(seurat_merge)
seurat_merge <- FindVariableFeatures(seurat_merge, nfeatures = 1000)

##RunningtSNE
seurat_merge <- RunPCA(seurat_merge,npcs = 20)
seurat_merge <- RunHarmony(seurat_merge, group.by.vars = "orig.ident")
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:20)
seurat_merge <- FindClusters(seurat_merge, resolution = 0.05, verbose = FALSE)
seurat_merge <- RunTSNE(seurat_merge, dims = 1:20,reduction = "harmony")
mycol<- brewer.pal(15, 'PiYG')
mycol0<- colorRampPalette(mycol)(5)
DimPlot(seurat_merge, reduction = "tsne", pt.size = 1, label = T) 
DimPlot(seurat_merge, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 1)
plot0 <- DimPlot(seurat_merge, reduction = "tsne", label = FALSE, repel = TRUE, pt.size = 2, alpha=1) + 
  scale_color_manual(values = mycol0) + 
  scale_fill_manual(values = mycol0)
plot0
plot1 <- FeaturePlot(seurat_merge, features = c("Adgre1","H2-Aa","Mrc1","Cmklr1"), cols =c("lightgray","red"),pt.size = 2,ncol=2)
plot1
plot2 <- DotPlot(seurat_merge,features = "Cmklr1")
plot2
ggsave("C57BL6Mice_Testis_Exercise_Macrophage_BarcodePlots_VlnComparison.tif", plot = plot0, width = 6, height = 5, dpi = 300)
ggsave("C57BL6Mice_Testis_Exercise_Macrophage_Barcode_FeaturePlots_VlnComparison.tif", plot = plot1, width = 10, height = 10, dpi = 300)
ggsave("C57BL6Mice_Testis_Exercise_Macrophage_Barcode_DotPlots_VlnComparison.tif", plot = plot2, width = 5, height = 5, dpi = 300)


my_comparisons <- list(c("SDO", "MICTO"),c("SDO","HIITO"))
##Inflammation
plot3 <- VlnPlot(seurat_merge,"Il1b", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+ 
  ylim(0, 8) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot3
plot4 <- VlnPlot(seurat_merge,"Tgfb1", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 8) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot4
merged_plot1 <- grid.arrange(plot3,plot4,nrow=1)
ggsave("C57BL6Mice_Testis_Exercise_Macrophage_Inflammation_VlnComparison.tif", plot = merged_plot1, width = 6, height = 4, dpi = 300)

##Plarization
plot5 <- VlnPlot(seurat_merge,"H2-Aa", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+ 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot5
plot6 <- VlnPlot(seurat_merge,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot6
merged_plot2 <- grid.arrange(plot5,plot6,nrow=1)
ggsave("C57BL6Mice_Testis_Exercise_Macrophage_CellPolarization_VlnComparison.tif", plot = merged_plot2, width = 6, height = 4, dpi = 300)

##
Cmklr1Pos_Macrophage <- subset(seurat_merge, ident=c("0"))
Cmklr1Neg_Macrophage <- subset(seurat_merge, ident=c("1"))
##Cmklr1CoupledMacrophageGlycolysis
plot7 <- VlnPlot(Cmklr1Pos_Macrophage,"Il1b", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot7
plot8 <- VlnPlot(Cmklr1Neg_Macrophage,"Il1b", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot8
merged_plot3 <- grid.arrange(plot7,plot8,nrow=1)
ggsave("C57BL6Mice_Testis_Exercise_Cmklr1CoupledMacrophage_Glycolysis_VlnComparison.tif", plot = merged_plot3, width = 4, height = 4, dpi = 300)

##Cmklr1CoupledMacrophageM2Polarization
plot9 <- VlnPlot(Cmklr1Pos_Macrophage,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot9
plot10 <- VlnPlot(Cmklr1Neg_Macrophage,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot10
merged_plot4 <- grid.arrange(plot9,plot10,nrow=1)
ggsave("C57BL6Mice_Testis_Exercise_Cmklr1CoupledMacrophage_M2Polarization_VlnComparison.tif", plot = merged_plot4, width = 4, height = 4, dpi = 300)

##Cmklr1PosMacrophageOXPHOS
plot11 <- VlnPlot(Cmklr1Pos_Macrophage,"Sdha", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot11 
plot12 <- VlnPlot(Cmklr1Pos_Macrophage,"Sdhb", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO"))  
plot12
plot13 <- VlnPlot(Cmklr1Pos_Macrophage,"Sdhc", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot13
plot14 <- VlnPlot(Cmklr1Pos_Macrophage,"Sdhd", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot14
merged_plot5 <- grid.arrange(plot11,plot12,plot13,plot14,nrow=1)
ggsave("C57BL6Mice_Testis_Exercise_Cmklr1PosMacrophage_OXPHOSCII_VlnComparison.tif", plot = merged_plot5, width = 6, height = 4, dpi = 300)

##Cmklr1NegMacrophageOXPHOS
plot15 <- VlnPlot(Cmklr1Neg_Macrophage,"Sdha", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot15
plot16 <- VlnPlot(Cmklr1Neg_Macrophage,"Sdhb", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO"))  
plot16
plot17 <- VlnPlot(Cmklr1Neg_Macrophage,"Sdhc", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot17
plot18 <- VlnPlot(Cmklr1Neg_Macrophage,"Sdhd", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot18
merged_plot6 <- grid.arrange(plot15,plot16,plot17,plot18,nrow=1)
ggsave("C57BL6Mice_Testis_Exercise_Cmklr1NegMacrophage_OXPHOSCII_VlnComparison.tif", plot = merged_plot6, width = 6, height = 4, dpi = 600)

##Subclustering
CPosMacrophage <- subset(seurat_merge, ident=c("0"))
CPosMacrophage <- SCTransform(CPosMacrophage)
CPosMacrophage <- FindVariableFeatures(CPosMacrophage, nfeatures = 1000)
CPosMacrophage <- RunPCA(CPosMacrophage,npcs = 20)
CPosMacrophage <- FindNeighbors(CPosMacrophage, dims = 1:20)
CPosMacrophage <- FindClusters(CPosMacrophage, resolution = 0.09, verbose = FALSE)
CPosMacrophage <- RunTSNE(CPosMacrophage, dims = 1:20)
mycol2<- brewer.pal(20, 'Paired')
mycol3<- colorRampPalette(mycol2)(10)
plot19 <- DimPlot(CPosMacrophage, reduction = "tsne", pt.size = 2, label = F) + 
  scale_color_manual(values = mycol3) + 
  scale_fill_manual(values = mycol3)
plot19
ggsave("C57BL6Mice_Testis_Exercise_Cmklr1PosMacrophage_BarcodePlots.tif", plot = plot19, width = 5, height = 4, dpi = 300)

plot20 <- FeaturePlot(CPosMacrophage, features = c("Adgre1","H2-Aa","Mrc1"), cols =c("lightgrey","red"),pt.size = 2)
plot20
ggsave("C57BL6Mice_Testis_Exercise_Cmklr1PosMacrophage_Barcode_FeaturePlots.tif", plot = plot20, width = 8, height = 8, dpi = 300)


##Subclustering
M1_Macrophage <- subset(CPosMacrophage, ident=c("0"))
M2_Macrophage <- subset(CPosMacrophage, ident=c("1"))
Norm_Macrophage <- subset(CPosMacrophage, ident=c("2"))
my_comparisons <- list(c("SDO", "MICTO"),c("SDO","HIITO"))

##Cmklr1PosMacrophageCII
plot21 <- VlnPlot(M1_Macrophage,"Sdhb", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot21
plot22 <- VlnPlot(M2_Macrophage,"Sdhb", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO"))  
plot22
plot23 <- VlnPlot(Norm_Macrophage,"Sdhb", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot23
merged_plot7 <- grid.arrange(plot21,plot22,plot23,nrow=1)
ggsave("C57BL6Mice_Testis_Exercise_Cmklr1PosMacrophage_Sdhb_VlnComparison.tif", plot = merged_plot7, width = 6, height = 4, dpi = 300)


plot24 <- VlnPlot(M1_Macrophage,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                 combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot24
plot25 <- VlnPlot(M2_Macrophage,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO"))  
plot25
plot26 <- VlnPlot(Norm_Macrophage,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#d50d8c","#00b3b0","#c2d0ea")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDO", "MICTO","HIITO")) 
plot26
merged_plot8 <- grid.arrange(plot24,plot25,plot26,nrow=1)
ggsave("C57BL6Mice_Testis_Exercise_Cmklr1PosMacrophage_CD206_VlnComparison.tif", plot = merged_plot8, width = 6, height = 4, dpi = 300)

