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
base_dir <- "E:/Github Data/GEO_MiceMacrophage_Aging"
groups <- c("SDY", "SDO")
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
seurat_merge <- SCTransform(seurat_merge, vars.to.regress = "percent.mt", verbose = FALSE)
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

##BarcodePlotting
mycol<- brewer.pal(12, 'PiYG')
mycol0<- colorRampPalette(mycol)(5)
plot5 <- DimPlot(seurat_merge, reduction = "tsne", label = FALSE, repel = TRUE, pt.size = 2, alpha=1) + 
  scale_color_manual(values = mycol0) + 
  scale_fill_manual(values = mycol0)
plot5
ggsave("C57BL6Mice_Aging_Testis_Cmklr1macrophage_BarcodePlot.tif", plot = plot5, width = 6, height = 5, dpi = 300)
plot6 <- FeaturePlot(seurat_merge, features = c("Adgre1","Mrc1","H2-Aa","Cmklr1"), cols =c("lightgrey","red"),pt.size = 2)
plot6
ggsave("C57BL6Mice_Aging_Testis_Cmklr1macrophage_Featureplots.tif", plot = plot6, width = 10, height = 8, dpi = 300)


##CMKLR1MonocyteDEGAnalysis
CPos_Macrophage <- subset(seurat_merge, ident=c("0"))
CPos_Macrophage <- PrepSCTFindMarkers(CPos_Macrophage)
Idents(CPos_Macrophage) <- CPos_Macrophage$orig.ident
table(Idents(CPos_Macrophage))
cut_off_log2FC=0.585
cut_off_P=0.05
CPosSDOvsSDY <- FindMarkers(CPos_Macrophage, min.pct=0.25, logfc.threshold=0.25,
                            ident.1 = "SDO", ident.2 = "SDY", assay = "SCT")
CPosSDOvsSDY$Sig = ifelse(CPosSDOvsSDY$p_val_adj < cut_off_P &
                            abs(CPosSDOvsSDY$avg_log2FC) >= cut_off_log2FC,  
                          ifelse(CPosSDOvsSDY$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
CPosSDOvsSDY = data.frame(CPosSDOvsSDY)
table(CPosSDOvsSDY$Sig)
plot7 <- EnhancedVolcano(CPosSDOvsSDY,
                         lab = rownames(CPosSDOvsSDY),
                         x = "avg_log2FC",
                         y = "p_val_adj",
                         boxedLabels = FALSE,
                         drawConnectors = FALSE,
                         pCutoff = 0.05,
                         FCcutoff = 0.585,
                         title = "CPosSDOvsSDY",
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
plot7
ggsave("C57BL6Mice_Aging_Testis_Cmklr1Monocyte_SDOvsSDY_VolcanoPlot.tif", plot = plot7, width = 10, height = 10, dpi = 300)

CPosSDOvsSDY <- FindMarkers(CPos_Macrophage, min.pct=0.25, logfc.threshold=0.25, group.by = "orig.ident", ident.1 = "SDO",ident.2 = "SDY")
write.csv(CPosSDOvsSDY,"C57BL6Mice_Aging_Testis_Cmklr1Monocyte_SDOvsSDY_DEG.csv")

plot8 <- DimPlot(seurat_merge, reduction = 'tsne', group.by = 'orig.ident',
                 label = F, pt.size = 1.5, alpha= 0.7, cols = c("#c2d0ea","#0070b8")) 
plot8
ggsave("C57BL6Mice_Aging_Testis_Cmklr1Monocyte_Splitbarcodeplots.tif", plot = plot8, width = 5, height = 4, dpi = 300)
table(seurat_merge$orig.ident)

##Freq Analysis
prop.table(table(Idents(seurat_merge)))
table(Idents(seurat_merge), seurat_merge$orig.ident)
Cellratio <- prop.table(table(Idents(seurat_merge), seurat_merge$orig.ident), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
plot9 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = NA)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA, color= NA, size=0.5, linetype="solid"))+
  scale_color_manual(values = mycol0) + 
  scale_fill_manual(values = mycol0) +
  scale_x_discrete(limits = c("SDY", "SDO")) 
plot9
ggsave("C57BL6Mice_Aging_Testis_Cmklr1Monocyte_SDOvsSDY_CMKLRL1_Ratio.tif", plot = plot9, width = 4, height = 4, dpi = 300)

##Cmklr1+MacrophageGlycolysis
my_comparisons <- list( c("SDY","SDO"))
plot10 <- VlnPlot(CPos_Macrophage,"Slc2a3", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot10
plot11 <- VlnPlot(CPos_Macrophage,"Hk1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot11
plot12 <- VlnPlot(CPos_Macrophage,"Pfkp", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot12
plot13 <- VlnPlot(CPos_Macrophage,"Pkm", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot13
plot14 <- VlnPlot(CPos_Macrophage,"Ldha", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 4) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot14
merged_plot1 <- grid.arrange(plot10,plot11,plot12,plot13,plot14,nrow=1)
ggsave("C57BL6Mice_Aging_Testis_Cmklr1PosMacrophage_Glycolysis.tif", plot = merged_plot1, width = 8, height = 4, dpi = 300)

##Cmklr1+MacrophagePolarization
plot15 <- VlnPlot(CPos_Macrophage,"H2-Aa", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 5) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot15
plot16 <- VlnPlot(CPos_Macrophage,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 5) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot16
merged_plot2 <- grid.arrange(plot15,plot16,nrow=1)
merged_plot2
ggsave("C57BL6Mice_Aging_Testis_Cmklr1PosMacrophage_Polariztion.tif", plot = merged_plot2, width = 6, height = 4, dpi = 300)

##CouplingCmklr1+Macrophage
CPos_Macrophage <- subset(seurat_merge, ident=c("0"))
CPos_Macrophage <- NormalizeData(CPos_Macrophage)
CPos_Macrophage <- FindVariableFeatures(CPos_Macrophage, nfeatures = 50)
CPos_Macrophage <- ScaleData(CPos_Macrophage)
CPos_Macrophage <- RunPCA(CPos_Macrophage, npc=15)
CPos_Macrophage <- FindNeighbors(CPos_Macrophage, dims = 1:10)
CPos_Macrophage <- FindClusters(CPos_Macrophage, resolution = 0.2, verbose = FALSE)
CPos_Macrophage <- RunTSNE(CPos_Macrophage, dims = 1:10)
mycol<- brewer.pal(12, 'Set1')
mycol0<- colorRampPalette(mycol)(10)
DimPlot(CPos_Macrophage, reduction = "tsne", pt.size = 1, label = T) 
DimPlot(CPos_Macrophage, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 1)
plot17 <- DimPlot(CPos_Macrophage, reduction = "tsne", label = FALSE, repel = TRUE, pt.size = 2, alpha=1) + 
  scale_color_manual(values = mycol0) + 
  scale_fill_manual(values = mycol0)
plot17
ggsave("C57BL6Mice_Aging_Testis_Cmklr1PosMacrophage_Subcluster.tif", plot = plot17, width = 6, height = 4, dpi = 300)
plot18 <- FeaturePlot(CPos_Macrophage, features = c("H2-Aa","Mrc1"), cols =c("lightgrey","red"),pt.size = 2)        
plot18
ggsave("C57BL6Mice_Aging_Testis_Cmklr1PosMacrophage_annotation.tif", plot = plot18, width = 8, height = 4, dpi = 300)

##SubclusteringCell
MHCIIHi_Macrophage <- subset(CPos_Macrophage, ident=c("2"))
CD206Hi_Macrophage <- subset(CPos_Macrophage, ident=c("0"))
Neg_Macrophage <- subset(CPos_Macrophage, ident=c("1"))

##OXPHOS
my_comparisons <- list( c("SDY", "SDO"))
plot19 <- VlnPlot(MHCIIHi_Macrophage,"Ndufc2", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot19
plot20 <- VlnPlot(CD206Hi_Macrophage,"Ndufc2", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot20
merged_plot3 <- grid.arrange(plot19,plot20,nrow=1)
merged_plot3
ggsave("C57BL6Mice_Aging_Testis_Cmklr1PosMacrophage_OXPHOS.tif", plot = merged_plot3, width = 4, height = 4, dpi = 300)

##OXPHOS
plot21 <- VlnPlot(MHCIIHi_Macrophage,"Cpt1b", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot21
plot22 <- VlnPlot(CD206Hi_Macrophage,"Cpt1b", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot22
merged_plot4 <- grid.arrange(plot21,plot22,nrow=1)
merged_plot4
ggsave("C57BL6Mice_Aging_Testis_Cmklr1PosMacrophage_FAO.tif", plot = merged_plot4, width = 4, height = 4, dpi = 300)

##Polarization
plot23 <- VlnPlot(MHCIIHi_Macrophage,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot23
plot24 <- VlnPlot(CD206Hi_Macrophage,"Mrc1", group.by = "orig.ident", pt.size = 0, 
                  combine = TRUE) + NoLegend() + geom_boxplot(fill="white", alpha = 0.9, width=0.1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  ylim(0, 6) + scale_fill_manual(values=c("#c2d0ea","#0070b8")) +
  theme(legend.position = 'none') + scale_x_discrete(limits = c("SDY", "SDO")) 
plot24
merged_plot5 <- grid.arrange(plot23,plot24,nrow=1)
merged_plot5
ggsave("C57BL6Mice_Aging_Testis_Cmklr1PosMacrophage_M2.tif", plot = merged_plot5, width = 4, height = 4, dpi = 300)