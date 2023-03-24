setwd("/work/lzy/project/onco_fetal/")
source("../utils/utils_plot.R")
source("../utils/utils_data_processing.R")
source("../utils/utils_color.R")
library(Seurat)
# library(SeuratObject)
library(sva)
library(Polychrome)
library(dplyr)
library(dendextend)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)

# >Color panel----
Cluster_levels <- c("B","CD8+ T","CD4+ T","Tregs","NK","Mast","DC","pDC","SPP1+ TAM2","MT1G+ TAM3","MYH11+ CAF","Other Mononuclear","Other Endothelium","Other CAF","FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Hepatocyte","Hepatocyte_P7_2","Hepatocyte_P15")
OF_Annotation_color_panel <- c("OF_CAF"="#EAA944","OF_Endo"="#69B4CE","OF_TAM"="#F3746C","OF_CAF_TAM"="#EF8F58","OF_CAF_Endo"="#A9AE89","OF_TAM_Endo"="#AE949D","OF_Endo_CAF_TAM"="#C29B80")
hcc01_color_panel <- 
  c("C1"="#ea776f","C2"="#eec98e","C3"="#ad5290","C4"="#77b077","C5"="#506ec3",
    "C6"="#f4bfb4","C7"="#92c2d2","C8"="#448732","C9"="#966fa7","C10"="#62bafd",
    "C11"="#e2b1cd")
hcc03_color_panel <- 
  c("C1"="#f4bfb4","C2"="#eec98e","C3"="#77b077","C4"="#966fa7","C5"="#ad5290",
    "C6"="#448732","C7"="#62bafd","C8"="#ea776f","C9"="#506ec3","C10"="#92c2d2",
    "C11"="#e2b1cd","C12"="#44979c")

# >Load single cell data----
HCC_seu <- readRDS("./02.processed_data/HCC_seu.rds")
# fetal_seu <- readRDS("./02.processed_data/fetal_scRNA_moreClusters_annotated.rds")
# Fib_seu <- readRDS("./02.processed_data/Fib_fetal.rds")

# >Generate Seurat object----
# source("../utils/utils_Gem2Seu.R")
# hcc01_gem <- data.table::fread("./00.data/BGI/hcc01.gem")
# Gem2Seu(data = hcc01_gem, binsize = 10, sample = "hcc01", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc01_gem, binsize = 25, sample = "hcc01", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc01_gem, binsize = 50, sample = "hcc01", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc01_gem, binsize = 100, sample = "hcc01", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc01_gem, binsize = 150, sample = "hcc01", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc01_gem, binsize = 200, sample = "hcc01", out.dir = "./02.processed_data/BGI");gc()
# hcc03_gem <- data.table::fread("./00.data/BGI/hcc03.gem")
# Gem2Seu(data = hcc03_gem, binsize = 10, sample = "hcc03", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc03_gem, binsize = 25, sample = "hcc03", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc03_gem, binsize = 50, sample = "hcc03", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc03_gem, binsize = 100, sample = "hcc03", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc03_gem, binsize = 150, sample = "hcc03", out.dir = "./02.processed_data/BGI");gc()
# Gem2Seu(data = hcc03_gem, binsize = 200, sample = "hcc03", out.dir = "./02.processed_data/BGI");gc()

# >Seurat pipeline----
# Reference
HCC_seu <- HCC_seu %>% subset(Sub_ClusterNew != "Bi-Potent")
HCC_seu <- SCTransform(HCC_seu, ncells = 10000, verbose = FALSE)
# hcc01_bin100
hcc01_bin100 <- readRDS("./02.processed_data/BGI/hcc01_bin100_seurat.rds")
hcc01_bin100@images$slice1@coordinates$row <- trunc(hcc01_bin100@images$slice1@coordinates$row/100) + 1
hcc01_bin100@images$slice1@coordinates$col <- trunc(hcc01_bin100@images$slice1@coordinates$col/100) + 1
hcc01_bin100@images$slice1@coordinates$imagerow <- trunc(hcc01_bin100@images$slice1@coordinates$imagerow/100) + 1
hcc01_bin100@images$slice1@coordinates$imagecol <- trunc(hcc01_bin100@images$slice1@coordinates$imagecol/100) + 1
hcc01_bin100@images$slice1@image <- hcc01_bin100@images$slice1@image[1:max(hcc01_bin100@images$slice1@coordinates$imagerow),1:max(hcc01_bin100@images$slice1@coordinates$imagecol)]
hcc01_bin100 <- SCTransform(hcc01_bin100, assay = "Spatial",  return.only.var.genes = FALSE, verbose = FALSE)
hcc01_bin100 <- NormalizeData(hcc01_bin100, assay = "Spatial", verbose = FALSE)
hcc01_bin100 <- ScaleData(hcc01_bin100, assay = "Spatial", verbose = FALSE)
hcc01_bin100 <- FindVariableFeatures(hcc01_bin100, assay = "Spatial", nfeatures = 2000)
hcc01_bin100 <- RunPCA(hcc01_bin100, assay = "SCT", verbose = FALSE)
hcc01_bin100 <- FindNeighbors(hcc01_bin100, reduction = "pca", dims = 1:30)
hcc01_bin100 <- FindClusters(hcc01_bin100, verbose = FALSE, resolution = .6)
hcc01_bin100 <- RunUMAP(hcc01_bin100, reduction = "pca", dims = 1:30, min.dist = .1, spread = 2)
hcc01_bin100@meta.data$UMAP_1 <- hcc01_bin100@reductions$umap@cell.embeddings[,1]
hcc01_bin100@meta.data$UMAP_2 <- hcc01_bin100@reductions$umap@cell.embeddings[,2]
hcc01_bin100$nCount <- log2(hcc01_bin100$nCount_Spatial + 1)
hcc01_bin100$nFeature <- log2(hcc01_bin100$nFeature_Spatial + 1)
hcc01_bin100_markers <- FindAllMarkers(hcc01_bin100)
write.csv(hcc01_bin100_markers, file = "./hcc01_bin100_markers.csv")
anchors <- FindTransferAnchors(reference = HCC_seu, query = hcc01_bin100, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = HCC_seu$Sub_ClusterNew, prediction.assay = TRUE, weight.reduction = hcc01_bin100[["pca"]], dims = 1:30)
hcc01_bin100[["predictions"]] <- predictions.assay
# hcc03_bin100
hcc03_bin100 <- readRDS("./02.processed_data/BGI/hcc03_bin100_seurat.rds")
hcc03_bin100@images$slice1@coordinates$row <- trunc(hcc03_bin100@images$slice1@coordinates$row/100) + 1
hcc03_bin100@images$slice1@coordinates$col <- trunc(hcc03_bin100@images$slice1@coordinates$col/100) + 1
hcc03_bin100@images$slice1@coordinates$imagerow <- trunc(hcc03_bin100@images$slice1@coordinates$imagerow/100) + 1
hcc03_bin100@images$slice1@coordinates$imagecol <- trunc(hcc03_bin100@images$slice1@coordinates$imagecol/100) + 1
hcc03_bin100@images$slice1@image <- hcc03_bin100@images$slice1@image[1:max(hcc03_bin100@images$slice1@coordinates$imagerow),1:max(hcc03_bin100@images$slice1@coordinates$imagecol)]
hcc03_bin100 <- SCTransform(hcc03_bin100, assay = "Spatial",  return.only.var.genes = FALSE, verbose = FALSE)
hcc03_bin100 <- NormalizeData(hcc03_bin100, assay = "Spatial", verbose = FALSE)
hcc03_bin100 <- FindVariableFeatures(hcc03_bin100, assay = "Spatial", nfeatures = 2000)
hcc03_bin100 <- ScaleData(hcc03_bin100, assay = "Spatial", verbose = FALSE)
hcc03_bin100 <- RunPCA(hcc03_bin100, assay = "SCT", verbose = FALSE)
hcc03_bin100 <- FindNeighbors(hcc03_bin100, reduction = "pca", dims = 1:30)
hcc03_bin100 <- FindClusters(hcc03_bin100, verbose = FALSE, resolution = .6)
hcc03_bin100 <- RunUMAP(hcc03_bin100, reduction = "pca", dims = 1:30)
hcc03_bin100@meta.data$UMAP_1 <- hcc03_bin100@reductions$umap@cell.embeddings[,1]
hcc03_bin100@meta.data$UMAP_2 <- hcc03_bin100@reductions$umap@cell.embeddings[,2]
hcc03_bin100$nCount <- log2(hcc03_bin100$nCount_Spatial + 1)
hcc03_bin100$nFeature <- log2(hcc03_bin100$nFeature_Spatial + 1)
hcc03_bin100_markers <- FindAllMarkers(hcc03_bin100)
anchors <- FindTransferAnchors(reference = HCC_seu, query = hcc03_bin100, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = HCC_seu$Sub_ClusterNew, prediction.assay = TRUE, weight.reduction = hcc03_bin100[["pca"]], dims = 1:30)
hcc03_bin100[["predictions"]] <- predictions.assay

# >Analyze----
HCC_seu <- readRDS("./02.processed_data/HCC_seu.rds")
hcc01_bin100 <- readRDS("./02.processed_data/BGI/hcc01_bin100_seurat_processed.rds")
hcc03_bin100 <- readRDS("./02.processed_data/BGI/hcc03_bin100_seurat_processed.rds")
hcc01_bin100@active.assay <- hcc03_bin100@active.assay <- "SCT"

# Basic data quality----
p1 <- SpatialFeaturePlot(hcc01_bin100, features = c("nFeature_Spatial"), pt.size.factor = 1.5, stroke = 0.2, min.cutoff = "q1")
p1.pdf <- ggAISpatial(p1) + ggtitle("HCC01")
p2 <- SpatialFeaturePlot(hcc03_bin100, features = c("nFeature_Spatial"), pt.size.factor = 1.5, stroke = 0.2, min.cutoff = "q1")
p2.pdf <- ggAISpatial(p2) + ggtitle("HCC03")
p <- plot_grid(p1.pdf, p2.pdf, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_nCount.pdf", width = 12, height = 5)

nGene.df <- data.frame(
  Sample = rep(c("HCC01","HCC03"), each = 6),
  BinSize = rep(c(10,25,50,100,150,200), 2),
  nGene = c(108.7663,556.363,1702.911,4418.332,6884.976,8881.547,
            94.8139,472.5612,1421.897,3470.576,5840.209,7643.424),
  nMID = c(145.7386,909.6103,3632.355,14436.48,32229.46,65940.78,
           147.0934,918.9376,3675.005,14573.38,32486.16,57319.09)
)
p1 <- ggplot(nGene.df, aes(x = BinSize, y = nGene)) +
  geom_point(aes(color = Sample)) +
  geom_line(aes(color = Sample)) +
  scale_y_log10() +
  theme_cowplot(font_size = 7)
p2 <- ggplot(nGene.df, aes(x = BinSize, y = nMID)) +
  geom_point(aes(color = Sample)) +
  geom_line(aes(color = Sample)) +
  scale_y_log10() +
  theme_cowplot(font_size = 7)
p <- p1 + p2
ggsave(p, file = "./04.figures/05.BGI_binsize_nCount.pdf", width = 5, height = 2)

# Gene expression of onco-fetal clusters----
hcc01_plot.list <- hcc03_plot.list <- list()
for(gene in c("FOLR2","CD68","PLVAP","PECAM1","POSTN","ACTA2","ENG","CD4","CTLA4","FOXP3","IL34","CSF1R","CSF1")){
  hcc01_plot.list[[gene]] <- SpatialFeaturePlot(hcc01_bin100, features = gene, pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
  hcc01_plot.list[[gene]] <- ggAISpatial(hcc01_plot.list[[gene]])
  hcc03_plot.list[[gene]] <- SpatialFeaturePlot(hcc03_bin100, features = gene, pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
  hcc03_plot.list[[gene]] <- ggAISpatial(hcc03_plot.list[[gene]])
}
p <- plot_grid(plotlist = hcc01_plot.list, nrow = 4)
ggsave(p, file = "./04.figures/05.BGI_hcc01_genes.pdf", width = 24, height = 20)
p <- plot_grid(plotlist = hcc03_plot.list, nrow = 4)
ggsave(p, file = "./04.figures/05.BGI_hcc03_genes.pdf", width = 24, height = 20)

OF_TAM_signature <- c("FOLR2","CD163","MRC1","TREM2")
OF_Endo_signature <- c("PLVAP","ACKR1","DLL4")
OF_CAF_signature <- c("POSTN","FAP","ENG","PTGDS")
hcc01_bin100@meta.data$OF_TAM_score <- colSums(hcc01_bin100@assays$SCT@scale.data[OF_TAM_signature,])
hcc01_bin100@meta.data$OF_Endo_score <- colSums(hcc01_bin100@assays$SCT@scale.data[OF_Endo_signature,])
hcc01_bin100@meta.data$OF_CAF_score <- colSums(hcc01_bin100@assays$SCT@scale.data[OF_CAF_signature,])
hcc03_bin100@meta.data$OF_TAM_score <- colSums(hcc03_bin100@assays$SCT@scale.data[OF_TAM_signature,])
hcc03_bin100@meta.data$OF_Endo_score <- colSums(hcc03_bin100@assays$SCT@scale.data[OF_Endo_signature,])
hcc03_bin100@meta.data$OF_CAF_score <- colSums(hcc03_bin100@assays$SCT@scale.data[OF_CAF_signature,])
hcc01_plot.list <- list(
  OF_TAM = ggAISpatial(SpatialFeaturePlot(hcc01_bin100, features = "OF_TAM_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')),
  OF_Endo = ggAISpatial(SpatialFeaturePlot(hcc01_bin100, features = "OF_Endo_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')),
  OF_CAF = ggAISpatial(SpatialFeaturePlot(hcc01_bin100, features = "OF_CAF_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99'))
)
p <- plot_grid(plotlist = hcc01_plot.list, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_hcc01_OF_score.pdf", width = 18, height = 5)
hcc03_plot.list <- list(
  OF_TAM = ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_TAM_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')),
  OF_Endo = ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_Endo_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')),
  OF_CAF = ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_CAF_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99'))
)
p <- plot_grid(plotlist = hcc03_plot.list, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_hcc03_OF_score.pdf", width = 18, height = 5)

# Spatial neighbors of OF clusters----
for(i in c("OF_Endo","OF_CAF","OF_TAM")){
  cutoff <- quantile(hcc01_bin100@meta.data[,paste0(i,"_score")], .975)
  high_cells <- hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data[,paste0(i,"_score")] > cutoff]
  high_neighbors <- FindNeighborSpots(hcc01_bin100, high_cells) 
  hcc01_bin100@meta.data[,paste0(i,"_high")] <- FALSE
  hcc01_bin100@meta.data[high_neighbors,paste0(i,"_high")] <- TRUE
  hcc01_bin100@meta.data[,paste0(i,"_high_self")] <- FALSE
  hcc01_bin100@meta.data[high_cells,paste0(i,"_high_self")] <- TRUE
}
for(i in c("OF_Endo","OF_CAF","OF_TAM")){
  cutoff <- quantile(hcc03_bin100@meta.data[,paste0(i,"_score")], .975)
  high_cells <- hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data[,paste0(i,"_score")] > cutoff]
  high_neighbors <- FindNeighborSpots(hcc03_bin100, high_cells) 
  hcc03_bin100@meta.data[,paste0(i,"_high")] <- FALSE
  hcc03_bin100@meta.data[high_neighbors,paste0(i,"_high")] <- TRUE
  hcc03_bin100@meta.data[,paste0(i,"_high_self")] <- FALSE
  hcc03_bin100@meta.data[high_cells,paste0(i,"_high_self")] <- TRUE
}

p1 <- ggplot(hcc01_bin100@meta.data, aes(x = OF_Endo_high, y = log2(OF_CAF_score - min(OF_CAF_score) + 1))) +
  geom_boxplot(aes(fill = OF_Endo_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "POSTN+ CAF score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "PLVAP+ Endo low","TRUE" = "PLVAP+ Endo high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p2 <- ggplot(hcc01_bin100@meta.data, aes(x = OF_TAM_high, y = log2(OF_CAF_score - min(OF_CAF_score) + 1))) +
  geom_boxplot(aes(fill = OF_TAM_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "POSTN+ CAF score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "FOLR2+ TAM1 low","TRUE" = "FOLR2+ TAM1 high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p <- p1 + p2
ggsave(p, file = "./04.figures/05.BGI_hcc01_OF_neighbor_enrich.pdf", width = 2.5, height = 3.5)

p1 <- ggplot(hcc03_bin100@meta.data, aes(x = OF_Endo_high, y = log2(OF_CAF_score - min(OF_CAF_score) + 1))) +
  geom_boxplot(aes(fill = OF_Endo_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "POSTN+ CAF score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "PLVAP+ Endo low","TRUE" = "PLVAP+ Endo high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p2 <- ggplot(hcc03_bin100@meta.data, aes(x = OF_TAM_high, y = log2(OF_CAF_score - min(OF_CAF_score) + 1))) +
  geom_boxplot(aes(fill = OF_TAM_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "POSTN+ CAF score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "FOLR2+ TAM1 low","TRUE" = "FOLR2+ TAM1 high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p <- p1 + p2
ggsave(p, file = "./04.figures/05.BGI_hcc03_OF_neighbor_enrich.pdf", width = 2.5, height = 3.5)

# Define OF high spots----
cells_highlight <- list(
  hcc01 = list(
    OF_TAM_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_TAM_high],
    OF_CAF_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_CAF_high],
    OF_Endo_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_Endo_high]),
  hcc03 = list(
    OF_TAM_high = hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_TAM_high],
    OF_CAF_high = hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_CAF_high],
    OF_Endo_high = hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_Endo_high])
)
cells_highlight$hcc01$OF_high <- intersect(intersect(cells_highlight$hcc01$OF_TAM_high, cells_highlight$hcc01$OF_CAF_high),cells_highlight$hcc01$OF_Endo_high)
hcc01_bin100@meta.data[,"OF_high"] <- FALSE
hcc01_bin100@meta.data[cells_highlight$hcc01$OF_high,"OF_high"] <- TRUE
cells_highlight$hcc03$OF_high <- intersect(intersect(cells_highlight$hcc03$OF_TAM_high, cells_highlight$hcc03$OF_CAF_high),cells_highlight$hcc03$OF_Endo_high)
hcc03_bin100@meta.data[,"OF_high"] <- FALSE
hcc03_bin100@meta.data[cells_highlight$hcc03$OF_high,"OF_high"] <- TRUE
plot.list <- list()
for(i in names(cells_highlight$hcc01)){
  plot.list[[i]] <- ggAISpatial(SpatialDimPlot(hcc01_bin100, cols.highlight = c("#F84141","#9FA0A3"),
                                               cells.highlight = cells_highlight$hcc01[[i]],
                                               pt.size.factor = 1.6, stroke = .1)) + ggtitle(i)
}
p <- plot_grid(plot.list$OF_TAM_high, plot.list$OF_CAF_high, plot.list$OF_Endo_high, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_hcc01_OF_cluster_high_spatial.pdf", width = 18, height = 5)
ggsave(plot.list$OF_high, file = "./04.figures/05.BGI_hcc01_OF_all_high_spatial.pdf", width = 6, height = 5)

plot.list <- list()
for(i in names(cells_highlight$hcc03)){
  plot.list[[i]] <- ggAISpatial(SpatialDimPlot(hcc03_bin100, cols.highlight = c("#F84141","#9FA0A3"),
                                               cells.highlight = cells_highlight$hcc03[[i]],
                                               pt.size.factor = 1.6, stroke = .1)) + ggtitle(i)
}
p <- plot_grid(plot.list$OF_TAM_high, plot.list$OF_CAF_high, plot.list$OF_Endo_high, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_hcc03_OF_cluster_high_spatial.pdf", width = 18, height = 5)
ggsave(plot.list$OF_high, file = "./04.figures/05.BGI_hcc03_OF_all_high_spatial.pdf", width = 6, height = 5)

# Calculate Neighborhood Index----
hcc03_bin100 <- CalNeighIndex(hcc03_bin100, hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_CAF_high_self], "OF_CAF")
hcc03_bin100 <- CalNeighIndex(hcc03_bin100, hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_TAM_high_self], "OF_TAM")
hcc03_bin100 <- CalNeighIndex(hcc03_bin100, hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_Endo_high_self], "OF_Endo")
hcc03_bin100 <- CalNeighIndex(hcc03_bin100, hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_high], "OF_High")
hcc03_bin100 <- CalProxiIndex(hcc03_bin100, hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_high], "OF_High")
hcc01_bin100 <- CalNeighIndex(hcc01_bin100, hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_CAF_high_self], "OF_CAF")
hcc01_bin100 <- CalNeighIndex(hcc01_bin100, hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_TAM_high_self], "OF_TAM")
hcc01_bin100 <- CalNeighIndex(hcc01_bin100, hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_Endo_high_self], "OF_Endo")
hcc01_bin100 <- CalNeighIndex(hcc01_bin100, hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_high], "OF_High")
hcc01_bin100 <- CalProxiIndex(hcc01_bin100, hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_high], "OF_High")
p1 <- ggAISpatial(SpatialFeaturePlot(hcc01_bin100, features = "OF_High_NeighIndex")) + ggtitle("HCC01")
p2 <- ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_High_NeighIndex")) + ggtitle("HCC03")
plot_list <- list(
  p1 = ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_CAF_NeighIndex")) + ggtitle("POSTN+ CAF"),
  p2 = ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_TAM_NeighIndex")) + ggtitle("FOLR2+ TAM1"),
  p3 = ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_Endo_NeighIndex")) + ggtitle("PLVAP+ Endo"),
  p4 = ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_High_NeighIndex")) + ggtitle("HCC03")
)
p <- plot_grid(p1, p2, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_OF_high_neighborhood_index.pdf", width = 11, height = 5)

p <- plot_grid(plot_list$p1, plot_list$p2, plot_list$p3, plot_list$p4, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_HCC03_OF_neighborhood_index.pdf", width = 21, height = 5)

# Ligands enriched in OF_high spots----
prioritized_ligands <- c("LTA","ADAM17","PGF","COL2A1","COL4A1","BMP2","ITGAM","TGFB3","VEGFA","BTLA")
hcc01_bin100@meta.data[,"Prioritized_ligands_score"] <- colSums(hcc01_bin100@assays$SCT@scale.data[prioritized_ligands,])
hcc03_bin100@meta.data[,"Prioritized_ligands_score"] <- colSums(hcc03_bin100@assays$SCT@scale.data[prioritized_ligands,])
p <- ggplot(hcc01_bin100@meta.data, aes(x = OF_high, y = log2(Prioritized_ligands_score - min(Prioritized_ligands_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "Prioritized ligands score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
ggsave(p, file = "./04.figures/05.BGI_hcc01_OF_prioritized_ligands_score.pdf", width = 1.5, height = 3.5)

p <- ggplot(hcc03_bin100@meta.data, aes(x = OF_high, y = log2(Prioritized_ligands_score - min(Prioritized_ligands_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "Prioritized ligands score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
ggsave(p, file = "./04.figures/05.BGI_hcc03_OF_prioritized_ligands_score.pdf", width = 1.5, height = 3.5)

# Zoom out OF-high regions----
# hcc01
hcc01_bin100@images$slice1@coordinates[,c("coord_col","coord_row")] <- stringr::str_split_fixed(row.names(hcc01_bin100@images$slice1@coordinates), pattern = ":|_", 3)[,c(2,3)]
hcc01_bin100@images$slice1@coordinates$coord_col <- as.numeric(hcc01_bin100@images$slice1@coordinates$coord_col)
hcc01_bin100@images$slice1@coordinates$coord_row <- as.numeric(hcc01_bin100@images$slice1@coordinates$coord_row)
hcc01_bin10 <- readRDS("./02.processed_data/BGI/hcc01_bin10_seurat.rds")
hcc01_bin10@images$slice1@coordinates[,c("coord_col","coord_row")] <- stringr::str_split_fixed(row.names(hcc01_bin10@images$slice1@coordinates), pattern = ":|_", 3)[,c(2,3)]
hcc01_bin10@images$slice1@coordinates$coord_col <- as.numeric(hcc01_bin10@images$slice1@coordinates$coord_col)
hcc01_bin10@images$slice1@coordinates$coord_row <- as.numeric(hcc01_bin10@images$slice1@coordinates$coord_row)
hcc01_bin10_OF_high_cells <- c()
for(cell in hcc01_OF_high){
  coord_x <- hcc01_bin100@images$slice1@coordinates[cell,"coord_row"]
  coord_y <- hcc01_bin100@images$slice1@coordinates[cell,"coord_col"]
  cells_used <- hcc01_bin10@images$slice1@coordinates %>% filter(coord_row>=coord_x, (coord_row<=coord_x+100), coord_col>=coord_y, (coord_col<=coord_y+100)) %>% row.names()
  hcc01_bin10_OF_high_cells <- c(hcc01_bin10_OF_high_cells, cells_used)
}
hcc01_bin10_OF_high_cells <- unique(hcc01_bin10_OF_high_cells)
hcc01_bin10@meta.data$OF_high <- hcc01_bin10@meta.data$cell %in% hcc01_bin10_OF_high_cells

OF_TAM_signature <- c("FOLR2","CD163","MRC1","TREM2")
OF_Endo_signature <- c("PLVAP","ACKR1","DLL4")
OF_CAF_signature <- c("POSTN","FAP","ENG","PTGDS")

cells_used1 <- hcc01_bin10@images$slice1@coordinates %>% filter(coord_row > 54800, coord_row < 55600, coord_col > 84950, coord_col < 85750) %>% row.names()
hcc01_bin10_subset1 <- subset(hcc01_bin10, subset = cell %in% cells_used1)
SpatialDimPlot(hcc01_bin10_subset1, group.by = "OF_high")
hcc01_bin10_subset1 <- SCTransform(hcc01_bin10_subset1, assay = "Spatial",  return.only.var.genes = FALSE, verbose = FALSE)
hcc01_bin10_subset1@meta.data$OF_CAF_score <- colSums(hcc01_bin10_subset1@assays$SCT@scale.data[intersect(OF_CAF_signature,row.names(hcc01_bin10_subset1)),])
hcc01_bin10_subset1@meta.data$OF_TAM_score <- colSums(hcc01_bin10_subset1@assays$SCT@scale.data[intersect(OF_TAM_signature,row.names(hcc01_bin10_subset1)),])
hcc01_bin10_subset1@meta.data$OF_Endo_score <- colSums(hcc01_bin10_subset1@assays$SCT@scale.data[intersect(OF_Endo_signature,row.names(hcc01_bin10_subset1)),])
for(i in c("OF_Endo","OF_CAF","OF_TAM")){
  mix.model <- mixtools::normalmixEM(hcc01_bin10_subset1@meta.data[,paste0(i,"_score")],k=2)
  cutoff <- mix.model$mu[2]
  high_cells <- hcc01_bin10_subset1@meta.data$cell[hcc01_bin10_subset1@meta.data[,paste0(i,"_score")] > cutoff]
  high_neighbors <- FindNeighborSpots(hcc01_bin10_subset1, high_cells)
  hcc01_bin10_subset1@meta.data[,paste0(i,"_high")] <- FALSE
  hcc01_bin10_subset1@meta.data[high_neighbors,paste0(i,"_high")] <- TRUE
}
hcc01_bin10_subset1@meta.data$OF_Annotation <- NA
hcc01_bin10_subset1@meta.data$OF_Annotation[hcc01_bin10_subset1@meta.data$OF_CAF_high] <- "OF_CAF"
hcc01_bin10_subset1@meta.data$OF_Annotation[hcc01_bin10_subset1@meta.data$OF_Endo_high] <- "OF_Endo"
hcc01_bin10_subset1@meta.data$OF_Annotation[hcc01_bin10_subset1@meta.data$OF_TAM_high] <- "OF_TAM"
hcc01_bin10_subset1@meta.data$OF_Annotation[hcc01_bin10_subset1@meta.data$OF_CAF_high & hcc01_bin10_subset1@meta.data$OF_TAM_high] <- "OF_CAF_TAM"
hcc01_bin10_subset1@meta.data$OF_Annotation[hcc01_bin10_subset1@meta.data$OF_CAF_high & hcc01_bin10_subset1@meta.data$OF_Endo_high] <- "OF_CAF_Endo"
hcc01_bin10_subset1@meta.data$OF_Annotation[hcc01_bin10_subset1@meta.data$OF_TAM_high & hcc01_bin10_subset1@meta.data$OF_Endo_high] <- "OF_TAM_Endo"
hcc01_bin10_subset1@meta.data$OF_Annotation[hcc01_bin10_subset1@meta.data$OF_Endo_high & hcc01_bin10_subset1@meta.data$OF_CAF_high & hcc01_bin10_subset1@meta.data$OF_TAM_high] <- "OF_Endo_CAF_TAM"

p0 <- ggAISpatial(SpatialDimPlot(hcc01_bin10_subset1, group.by = "OF_high", pt.size.factor = 3.5, stroke = .1))
p1 <- ggAISpatial(
  SpatialDimPlot(hcc01_bin10_subset1, group.by = "OF_Annotation", pt.size.factor = 3.5, stroke = .1) +
    scale_fill_manual(values = OF_Annotation_color_panel, na.value = "lightgrey") +
    guides(fill = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1)))
p <- plot_grid(p0, p1)
ggsave(p, file = "./04.figures/05.BGI_hcc01_bin10_subset1_OF_annotation.pdf", width = 10, height = 4)

cells_used2 <- hcc01_bin10@images$slice1@coordinates %>% filter(coord_row > 44550, coord_row < 45150, coord_col > 88950, coord_col < 89550) %>% row.names()
hcc01_bin10_subset2 <- subset(hcc01_bin10, subset = cell %in% cells_used2)
hcc01_bin10_subset2 <- SCTransform(hcc01_bin10_subset2, assay = "Spatial",  return.only.var.genes = FALSE, verbose = FALSE)
hcc01_bin10_subset2@meta.data$OF_CAF_score <- colSums(hcc01_bin10_subset2@assays$SCT@scale.data[intersect(OF_CAF_signature, row.names(hcc01_bin10_subset2)),])
hcc01_bin10_subset2@meta.data$OF_TAM_score <- colSums(hcc01_bin10_subset2@assays$SCT@scale.data[intersect(OF_TAM_signature, row.names(hcc01_bin10_subset2)),])
hcc01_bin10_subset2@meta.data$OF_Endo_score <- colSums(hcc01_bin10_subset2@assays$SCT@scale.data[intersect(OF_Endo_signature, row.names(hcc01_bin10_subset2)),])
for(i in c("OF_Endo","OF_CAF","OF_TAM")){
  mix.model <- mixtools::normalmixEM(hcc01_bin10_subset2@meta.data[,paste0(i,"_score")],k=2)
  cutoff <- mix.model$mu[2]
  high_cells <- hcc01_bin10_subset2@meta.data$cell[hcc01_bin10_subset2@meta.data[,paste0(i,"_score")] > cutoff]
  high_neighbors <- FindNeighborSpots(hcc01_bin10_subset2, high_cells) 
  hcc01_bin10_subset2@meta.data[,paste0(i,"_high")] <- FALSE
  hcc01_bin10_subset2@meta.data[high_neighbors,paste0(i,"_high")] <- TRUE
}
hcc01_bin10_subset2@meta.data$OF_Annotation <- NA
hcc01_bin10_subset2@meta.data$OF_Annotation[hcc01_bin10_subset2@meta.data$OF_CAF_high] <- "OF_CAF"
hcc01_bin10_subset2@meta.data$OF_Annotation[hcc01_bin10_subset2@meta.data$OF_Endo_high] <- "OF_Endo"
hcc01_bin10_subset2@meta.data$OF_Annotation[hcc01_bin10_subset2@meta.data$OF_TAM_high] <- "OF_TAM"
hcc01_bin10_subset2@meta.data$OF_Annotation[hcc01_bin10_subset2@meta.data$OF_CAF_high & hcc01_bin10_subset2@meta.data$OF_TAM_high] <- "OF_CAF_TAM"
hcc01_bin10_subset2@meta.data$OF_Annotation[hcc01_bin10_subset2@meta.data$OF_CAF_high & hcc01_bin10_subset2@meta.data$OF_Endo_high] <- "OF_CAF_Endo"
hcc01_bin10_subset2@meta.data$OF_Annotation[hcc01_bin10_subset2@meta.data$OF_TAM_high & hcc01_bin10_subset2@meta.data$OF_Endo_high] <- "OF_TAM_Endo"
hcc01_bin10_subset2@meta.data$OF_Annotation[hcc01_bin10_subset2@meta.data$OF_Endo_high & hcc01_bin10_subset2@meta.data$OF_CAF_high & hcc01_bin10_subset2@meta.data$OF_TAM_high] <- "OF_Endo_CAF_TAM"

p0 <- ggAISpatial(SpatialDimPlot(hcc01_bin10_subset2, group.by = "OF_high", pt.size.factor = 4, stroke = .1))
p1 <- ggAISpatial(
  SpatialDimPlot(hcc01_bin10_subset2, group.by = "OF_Annotation", pt.size.factor = 4, stroke = .1) +
    scale_fill_manual(values = OF_Annotation_color_panel, na.value = "lightgrey") +
    guides(fill = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1)))
p <- plot_grid(p0, p1)
ggsave(p, file = "./04.figures/05.BGI_hcc01_bin10_subset2_OF_annotation.pdf", width = 10, height = 4)

# hcc03
hcc03_bin100@images$slice1@coordinates[,c("coord_col","coord_row")] <- stringr::str_split_fixed(row.names(hcc03_bin100@images$slice1@coordinates), pattern = ":|_", 3)[,c(2,3)]
hcc03_bin100@images$slice1@coordinates$coord_col <- as.numeric(hcc03_bin100@images$slice1@coordinates$coord_col)
hcc03_bin100@images$slice1@coordinates$coord_row <- as.numeric(hcc03_bin100@images$slice1@coordinates$coord_row)
hcc03_bin10 <- readRDS("./02.processed_data/BGI/hcc03_bin10_seurat.rds")
hcc03_bin10@images$slice1@coordinates[,c("coord_col","coord_row")] <- stringr::str_split_fixed(row.names(hcc03_bin10@images$slice1@coordinates), pattern = ":|_", 3)[,c(2,3)]
hcc03_bin10@images$slice1@coordinates$coord_col <- as.numeric(hcc03_bin10@images$slice1@coordinates$coord_col)
hcc03_bin10@images$slice1@coordinates$coord_row <- as.numeric(hcc03_bin10@images$slice1@coordinates$coord_row)
hcc03_bin10_OF_high_cells <- c()
for(cell in hcc03_OF_high){
  coord_x <- hcc03_bin100@images$slice1@coordinates[cell,"coord_row"]
  coord_y <- hcc03_bin100@images$slice1@coordinates[cell,"coord_col"]
  cells_used <- hcc03_bin10@images$slice1@coordinates %>% filter(coord_row>=coord_x, (coord_row<=coord_x+100), coord_col>=coord_y, (coord_col<=coord_y+100)) %>% row.names()
  hcc03_bin10_OF_high_cells <- c(hcc03_bin10_OF_high_cells, cells_used)
}
hcc03_bin10_OF_high_cells <- unique(hcc03_bin10_OF_high_cells)
hcc03_bin10@meta.data$OF_high <- hcc03_bin10@meta.data$cell %in% hcc03_bin10_OF_high_cells

OF_TAM_signature <- c("FOLR2","CD163","MRC1","TREM2")
OF_Endo_signature <- c("PLVAP","ACKR1","DLL4")
OF_CAF_signature <- c("POSTN","FAP","ENG","PTGDS")

cells_used1 <- hcc03_bin10@images$slice1@coordinates %>% filter(coord_row > 20000, coord_row < 21500, coord_col > 14550, coord_col < 16050) %>% row.names()
hcc03_bin10_subset1 <- subset(hcc03_bin10, subset = cell %in% cells_used1)
SpatialDimPlot(hcc03_bin10_subset1, group.by = "OF_high")
hcc03_bin10_subset1 <- SCTransform(hcc03_bin10_subset1, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
hcc03_bin10_subset1@meta.data$OF_CAF_score <- colSums(hcc03_bin10_subset1@assays$SCT@scale.data[intersect(OF_CAF_signature,row.names(hcc03_bin10_subset1)),])
hcc03_bin10_subset1@meta.data$OF_TAM_score <- colSums(hcc03_bin10_subset1@assays$SCT@scale.data[intersect(OF_TAM_signature,row.names(hcc03_bin10_subset1)),])
hcc03_bin10_subset1@meta.data$OF_Endo_score <- colSums(hcc03_bin10_subset1@assays$SCT@scale.data[intersect(OF_Endo_signature,row.names(hcc03_bin10_subset1)),])
for(i in c("OF_Endo","OF_CAF","OF_TAM")){
  mix.model <- mixtools::normalmixEM(hcc03_bin10_subset1@meta.data[,paste0(i,"_score")],k=2)
  cutoff <- mix.model$mu[2]
  high_cells <- hcc03_bin10_subset1@meta.data$cell[hcc03_bin10_subset1@meta.data[,paste0(i,"_score")] > cutoff]
  high_neighbors <- FindNeighborSpots(hcc03_bin10_subset1, high_cells) 
  hcc03_bin10_subset1@meta.data[,paste0(i,"_high")] <- FALSE
  hcc03_bin10_subset1@meta.data[high_neighbors,paste0(i,"_high")] <- TRUE
}
hcc03_bin10_subset1@meta.data$OF_Annotation <- NA
hcc03_bin10_subset1@meta.data$OF_Annotation[hcc03_bin10_subset1@meta.data$OF_CAF_high] <- "OF_CAF"
hcc03_bin10_subset1@meta.data$OF_Annotation[hcc03_bin10_subset1@meta.data$OF_Endo_high] <- "OF_Endo"
hcc03_bin10_subset1@meta.data$OF_Annotation[hcc03_bin10_subset1@meta.data$OF_TAM_high] <- "OF_TAM"
hcc03_bin10_subset1@meta.data$OF_Annotation[hcc03_bin10_subset1@meta.data$OF_CAF_high & hcc03_bin10_subset1@meta.data$OF_TAM_high] <- "OF_CAF_TAM"
hcc03_bin10_subset1@meta.data$OF_Annotation[hcc03_bin10_subset1@meta.data$OF_CAF_high & hcc03_bin10_subset1@meta.data$OF_Endo_high] <- "OF_CAF_Endo"
hcc03_bin10_subset1@meta.data$OF_Annotation[hcc03_bin10_subset1@meta.data$OF_TAM_high & hcc03_bin10_subset1@meta.data$OF_Endo_high] <- "OF_TAM_Endo"
hcc03_bin10_subset1@meta.data$OF_Annotation[hcc03_bin10_subset1@meta.data$OF_Endo_high & hcc03_bin10_subset1@meta.data$OF_CAF_high & hcc03_bin10_subset1@meta.data$OF_TAM_high] <- "OF_Endo_CAF_TAM"
OF_Annotation_color_panel <- c("OF_CAF"="#EAA944","OF_Endo"="#69B4CE","OF_TAM"="#F3746C","OF_CAF_TAM"="#EF8F58","OF_CAF_Endo"="#A9AE89","OF_TAM_Endo"="#AE949D","OF_Endo_CAF_TAM"="#C29B80")
p0 <- ggAISpatial(SpatialDimPlot(hcc03_bin10_subset1, group.by = "OF_high", pt.size.factor = 2, stroke = .1))
p1 <- ggAISpatial(
  SpatialDimPlot(hcc03_bin10_subset1, group.by = "OF_Annotation", pt.size.factor = 2, stroke = .1) +
    scale_fill_manual(values = OF_Annotation_color_panel, na.value = "lightgrey") +
    guides(fill = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1)))
p <- plot_grid(p0, p1)
ggsave(p, file = "./04.figures/05.BGI_hcc03_bin10_subset1_OF_annotation.pdf", width = 10, height = 4)

cells_used2 <- hcc03_bin10@images$slice1@coordinates %>% filter(coord_row > 15550, coord_row < 17650, coord_col > 4700, coord_col < 6700) %>% row.names()
hcc03_bin10_subset2 <- subset(hcc03_bin10, subset = cell %in% cells_used2)
hcc03_bin10_subset2 <- SCTransform(hcc03_bin10_subset2, assay = "Spatial",  return.only.var.genes = FALSE, verbose = FALSE)
hcc03_bin10_subset2@meta.data$OF_CAF_score <- colSums(hcc03_bin10_subset2@assays$SCT@scale.data[intersect(OF_CAF_signature, row.names(hcc03_bin10_subset2)),])
hcc03_bin10_subset2@meta.data$OF_TAM_score <- colSums(hcc03_bin10_subset2@assays$SCT@scale.data[intersect(OF_TAM_signature, row.names(hcc03_bin10_subset2)),])
hcc03_bin10_subset2@meta.data$OF_Endo_score <- colSums(hcc03_bin10_subset2@assays$SCT@scale.data[intersect(OF_Endo_signature, row.names(hcc03_bin10_subset2)),])
for(i in c("OF_Endo","OF_CAF","OF_TAM")){
  mix.model <- mixtools::normalmixEM(hcc03_bin10_subset2@meta.data[,paste0(i,"_score")],k=2)
  cutoff <- mix.model$mu[2]
  high_cells <- hcc03_bin10_subset2@meta.data$cell[hcc03_bin10_subset2@meta.data[,paste0(i,"_score")] > cutoff]
  high_neighbors <- FindNeighborSpots(hcc03_bin10_subset2, high_cells)
  hcc03_bin10_subset2@meta.data[,paste0(i,"_high")] <- FALSE
  hcc03_bin10_subset2@meta.data[high_neighbors,paste0(i,"_high")] <- TRUE
}
hcc03_bin10_subset2@meta.data$OF_Annotation <- NA
hcc03_bin10_subset2@meta.data$OF_Annotation[hcc03_bin10_subset2@meta.data$OF_CAF_high] <- "OF_CAF"
hcc03_bin10_subset2@meta.data$OF_Annotation[hcc03_bin10_subset2@meta.data$OF_Endo_high] <- "OF_Endo"
hcc03_bin10_subset2@meta.data$OF_Annotation[hcc03_bin10_subset2@meta.data$OF_TAM_high] <- "OF_TAM"
hcc03_bin10_subset2@meta.data$OF_Annotation[hcc03_bin10_subset2@meta.data$OF_CAF_high & hcc03_bin10_subset2@meta.data$OF_TAM_high] <- "OF_CAF_TAM"
hcc03_bin10_subset2@meta.data$OF_Annotation[hcc03_bin10_subset2@meta.data$OF_CAF_high & hcc03_bin10_subset2@meta.data$OF_Endo_high] <- "OF_CAF_Endo"
hcc03_bin10_subset2@meta.data$OF_Annotation[hcc03_bin10_subset2@meta.data$OF_TAM_high & hcc03_bin10_subset2@meta.data$OF_Endo_high] <- "OF_TAM_Endo"
hcc03_bin10_subset2@meta.data$OF_Annotation[hcc03_bin10_subset2@meta.data$OF_Endo_high & hcc03_bin10_subset2@meta.data$OF_CAF_high & hcc03_bin10_subset2@meta.data$OF_TAM_high] <- "OF_Endo_CAF_TAM"

p0 <- ggAISpatial(SpatialDimPlot(hcc03_bin10_subset2, group.by = "OF_high", pt.size.factor = 1.5, stroke = .1))
p1 <- ggAISpatial(
  SpatialDimPlot(hcc03_bin10_subset2, group.by = "OF_Annotation", pt.size.factor = 1.5, stroke = .1) +
    scale_fill_manual(values = OF_Annotation_color_panel, na.value = "lightgrey") +
    guides(fill = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1)))
p <- plot_grid(p0, p1)
ggsave(p, file = "./04.figures/05.BGI_hcc03_bin10_subset2_OF_annotation.pdf", width = 10, height = 4)

 # Cluster based analysis----
hcc01_bin100@meta.data$Cluster <- plyr::revalue(
  hcc01_bin100@meta.data$SCT_snn_res.0.6,
  replace=c("0"="C2","1"="C6","2"="C4","3"="C1","4"="C3",
            "5"="C10","6"="C4","7"="C1","8"="C7","9"="C8",
            "10"="C5","11"="C9","12"="C4","13"="C10","14"="C11"))
hcc01_bin100@meta.data$Cluster <- factor(hcc01_bin100@meta.data$Cluster, levels = paste0("C",1:11))

p1 <- ggplot(hcc01_bin100@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 0.5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(name = "", values = hcc01_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p1.pdf <- ggAIplot(p1)
ggsave(p1.pdf, file = "./04.figures/05.BGI_hcc01_umap.pdf", width = 8, height = 6.5)
p2 <- ggAISpatial(SpatialDimPlot(hcc01_bin100, pt.size.factor = 1.5, stroke = 0, group.by = "Cluster") + scale_fill_manual(values = alpha(hcc01_color_panel, .9)))
ggsave(p2, file = "./04.figures/05.BGI_hcc01_cluster_spatial.pdf", width = 6, height = 5)

features_used <- c("OF_CAF_score","OF_TAM_score","OF_Endo_score")
hcc01.plot.data <- melt(hcc01_bin100@meta.data[,c(features_used,"UMAP_1","UMAP_2")],id.vars = c("UMAP_1","UMAP_2"))
hcc01.plot.data$value <- pmin(hcc01.plot.data$value, quantile(hcc01.plot.data$value,.995))
hcc01.plot.data$value <- pmax(hcc01.plot.data$value, quantile(hcc01.plot.data$value,.005))
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p <- 
  ggplot(hcc01.plot.data[,c("value","variable","UMAP_1","UMAP_2")] %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 1) +
  facet_wrap(~variable, nrow = 1) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("Exp", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/05.BGI_hcc01_OF_score_UMAP.pdf", width = 16, height = 5.75)

hcc03_bin100@meta.data$Cluster <- plyr::revalue(
  hcc03_bin100@meta.data$SCT_snn_res.0.6,
  replace=c("0"="C1","1"="C2","2"="C7","3"="C3","4"="C4",
            "5"="C8","6"="C7","7"="C12","8"="C6","9"="C5",
            "10"="C9","11"="C8","12"="C2","13"="C2","14"="C10",
            "15"="C11"))
hcc03_bin100@meta.data$Cluster <- factor(hcc03_bin100@meta.data$Cluster, levels = paste0("C",1:12))

p1 <- ggplot(hcc03_bin100@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 0.5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(name = "", values = hcc03_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p1.pdf <- ggAIplot(p1)
ggsave(p1.pdf, file = "./04.figures/05.BGI_hcc03_umap.pdf", width = 8, height = 6.5)
p2 <- ggAISpatial(SpatialDimPlot(hcc03_bin100, pt.size.factor = 1.5, stroke = 0, group.by = "Cluster") + scale_fill_manual(values = alpha(hcc03_color_panel,0.9)))
ggsave(p2, file = "./04.figures/05.BGI_hcc03_cluster_spatial.pdf", width = 6, height = 5)

features_used <- c("OF_CAF_score","OF_TAM_score","OF_Endo_score")
hcc03.plot.data <- melt(hcc03_bin100@meta.data[,c(features_used,"UMAP_1","UMAP_2")],id.vars = c("UMAP_1","UMAP_2"))
hcc03.plot.data$value <- pmin(hcc03.plot.data$value, quantile(hcc03.plot.data$value,.995))
hcc03.plot.data$value <- pmax(hcc03.plot.data$value, quantile(hcc03.plot.data$value,.005))
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p <- 
  ggplot(hcc03.plot.data[,c("value","variable","UMAP_1","UMAP_2")] %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 1) +
  facet_wrap(~variable, nrow = 1) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("Exp", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/05.BGI_hcc03_OF_score_UMAP.pdf", width = 16, height = 5.75)

# OF high related tumor genes expressed in spots----
degenes <- read.csv("./03.results/DEGenes/HCC_bulk_Recurrence_OF_high_vs_OF_low.csv", row.names = 1)
recur_OFhigh_genes <- intersect(degenes %>% arrange(desc(logFC)) %>% head(20) %>% pull(Symbol) %>% as.character(), row.names(hcc01_bin100@assays$SCT@scale.data))
hcc01_bin100@meta.data[,"recur_OFhigh_score"] <- colSums(hcc01_bin100@assays$SCT@scale.data[recur_OFhigh_genes,])
cutoff <- quantile(hcc01_bin100@meta.data[,"recur_OFhigh_score"], .975)
high_cells <- hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data[,"recur_OFhigh_score"] > cutoff]
hcc01_bin100@meta.data[,"recur_OFhigh"] <- FALSE
hcc01_bin100@meta.data[high_cells,"recur_OFhigh"] <- TRUE
hcc01_bin100 <- CalNeighIndex(hcc01_bin100, hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$recur_OFhigh], "recur_OFhigh")
hcc03_bin100@meta.data[,"recur_OFhigh_score"] <- colSums(hcc03_bin100@assays$SCT@scale.data[recur_OFhigh_genes,])
cutoff <- quantile(hcc03_bin100@meta.data[,"recur_OFhigh_score"], .975)
high_cells <- hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data[,"recur_OFhigh_score"] > cutoff]
hcc03_bin100@meta.data[,"recur_OFhigh"] <- FALSE
hcc03_bin100@meta.data[high_cells,"recur_OFhigh"] <- TRUE
hcc03_bin100 <- CalNeighIndex(hcc03_bin100, hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$recur_OFhigh], "recur_OFhigh")

hcc01_bin100@meta.data$OF_score <- hcc01_bin100@meta.data$OF_CAF_NeighIndex + hcc01_bin100@meta.data$OF_TAM_NeighIndex + hcc01_bin100@meta.data$OF_Endo_NeighIndex
p1 <- ggAISpatial(SpatialFeaturePlot(hcc01_bin100, features = "OF_score", pt.size.factor = 1.5, stroke = 0, max.cutoff = "q95") + viridis::scale_fill_viridis())
ggsave(p1, file = "./04.figures/05.BGI_hcc01_OF_NeighIndex_score.pdf", width = 6, height = 5)
p2 <- ggAISpatial(SpatialDimPlot(hcc01_bin100, group.by = "recur_OFhigh", pt.size.factor = 1.5, stroke = 0) + scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#F84141")))
ggsave(p2, file = "./04.figures/05.BGI_hcc01_recur_OFhigh_score.pdf", width = 6, height = 5)

hcc03_bin100@meta.data$OF_score <- hcc03_bin100@meta.data$OF_CAF_NeighIndex + hcc03_bin100@meta.data$OF_TAM_NeighIndex + hcc03_bin100@meta.data$OF_Endo_NeighIndex
p1 <- ggAISpatial(SpatialFeaturePlot(hcc03_bin100, features = "OF_score", pt.size.factor = 1.5, stroke = 0, max.cutoff = "q95") + viridis::scale_fill_viridis())
ggsave(p1, file = "./04.figures/05.BGI_hcc03_OF_NeighIndex_score.pdf", width = 6, height = 5)
p2 <- ggAISpatial(SpatialDimPlot(hcc03_bin100, group.by = "recur_OFhigh", pt.size.factor = 1.5, stroke = 0) + scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#F84141")))
ggsave(p2, file = "./04.figures/05.BGI_hcc03_recur_OFhigh_score.pdf", width = 6, height = 5)

hcc01_recur_dist_OF <- 
  DistFeature(seu = hcc01_bin100, 
              cells = hcc01_bin100@meta.data %>% filter(recur_OFhigh) %>% row.names(),
              feature = "OF_score",
              max.dist = 20)
hcc01_recur_dist_OF_random_list <- list()
for(i in 1:10){
  hcc01_recur_dist_OF_random_list[[i]] <- 
    DistFeature(seu = hcc01_bin100, 
                cells = sample(row.names(hcc01_bin100@meta.data), sum(hcc01_bin100@meta.data$recur_OFhigh)),
                feature = "OF_score",
                max.dist = 20)
}

hcc03_recur_dist_OF <- 
  DistFeature(seu = hcc03_bin100, 
              cells = hcc03_bin100@meta.data %>% filter(recur_OFhigh) %>% row.names(),
              feature = "OF_score",
              max.dist = 20)
hcc03_recur_dist_OF_random_list <- list()
for(i in 1:10){
  hcc03_recur_dist_OF_random_list[[i]] <- 
    DistFeature(seu = hcc03_bin100, 
                cells = sample(row.names(hcc03_bin100@meta.data), sum(hcc03_bin100@meta.data$recur_OFhigh)),
                feature = "OF_score",
                max.dist = 20)
}

save(hcc01_recur_dist_OF, hcc01_recur_dist_OF_random_list, hcc03_recur_dist_OF, hcc03_recur_dist_OF_random_list, file = "./03.results/STOmics_hcc_recur_dist_statistic.rda")

hcc01_dist_OF_random <- hcc03_dist_OF_random <- c()
for(i in 1:10){
  hcc01_dist_OF_random <- rbind(hcc01_dist_OF_random, colMeans(hcc01_recur_dist_OF_random_list[[i]]$dist.feature))
  hcc03_dist_OF_random <- rbind(hcc03_dist_OF_random, colMeans(hcc03_recur_dist_OF_random_list[[i]]$dist.feature))
}

p <- hcc03_dist_OF_random %>% melt() %>% 
  ggplot(aes(x = factor(Var2), y = value)) + 
  geom_boxplot(color = "#9FA0A3", outlier.stroke = 0) + 
  geom_point(data = data.frame(value = colMeans(hcc03_recur_dist_OF$dist.feature),
                               Vars = factor(1:20)),
             aes(x = Vars, y = value), color = "#F84141") +
  labs(x = "Distance", y = "Proportions of Onco-fetal High Bins") +
  theme_cowplot(font_size = 7, line_size = .25)
ggsave(p, file = "./04.figures/05.BGI_hcc03_recur_OFhigh_dist_score.pdf", width = 4, height = 3)

# TLS enriched in OF_high spots----
TLS_chemokine_genes <- intersect(c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13"), row.names(hcc01_bin100@assays$SCT@scale.data))
TLS_TFH_genes <- intersect(c("CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1"), row.names(hcc01_bin100@assays$SCT@scale.data))
TLS_TB_genes <- intersect(unique(c("CD4","CCR5","CXCR3","CSF2","IGSF6","IL2RA","CD38","CD40","CD5","MS4A1","SDC1","GFI1","IL1R1","IL1R2","IL10","CCL20","IRF4","TRAF6","STAT5A")), row.names(hcc01_bin100@assays$SCT@scale.data))
hcc01_bin100@meta.data[,"TLS_chemokine_score"] <- colSums(hcc01_bin100@assays$SCT@scale.data[TLS_chemokine_genes,])
hcc01_bin100@meta.data[,"TLS_TFH_score"] <- colSums(hcc01_bin100@assays$SCT@scale.data[TLS_TFH_genes,])
hcc01_bin100@meta.data[,"TLS_TB_score"] <- colSums(hcc01_bin100@assays$SCT@scale.data[TLS_TB_genes,])
hcc03_bin100@meta.data[,"TLS_chemokine_score"] <- colSums(hcc03_bin100@assays$SCT@scale.data[TLS_chemokine_genes,])
hcc03_bin100@meta.data[,"TLS_TFH_score"] <- colSums(hcc03_bin100@assays$SCT@scale.data[TLS_TFH_genes,])
hcc03_bin100@meta.data[,"TLS_TB_score"] <- colSums(hcc03_bin100@assays$SCT@scale.data[TLS_TB_genes,])
p1 <- ggplot(hcc01_bin100@meta.data, aes(x = OF_high, y = log2(TLS_chemokine_score - min(TLS_chemokine_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "TLS Signature (Chemokine)") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p2 <- ggplot(hcc01_bin100@meta.data, aes(x = OF_high, y = log2(TLS_TFH_score - min(TLS_TFH_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "TLS Signature (TFH)") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p3 <- ggplot(hcc01_bin100@meta.data, aes(x = OF_high, y = log2(TLS_TB_score - min(TLS_TB_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "TLS Signature (Th1&B cell)") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p <- p1 + p2 + p3
ggsave(p, file = "./04.figures/05.BGI_hcc01_OF_TLS_score.pdf", width = 3.5, height = 3.5)

p1 <- ggplot(hcc03_bin100@meta.data, aes(x = OF_high, y = log2(TLS_chemokine_score - min(TLS_chemokine_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "TLS Signature (Chemokine)") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p2 <- ggplot(hcc03_bin100@meta.data, aes(x = OF_high, y = log2(TLS_TFH_score - min(TLS_TFH_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "TLS Signature (TFH)") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p3 <- ggplot(hcc03_bin100@meta.data, aes(x = OF_high, y = log2(TLS_TB_score - min(TLS_TB_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "TLS Signature (Th1&B cell)") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p <- p1 + p2 + p3
ggsave(p, file = "./04.figures/05.BGI_hcc03_OF_TLS_score.pdf", width = 3.5, height = 3.5)

SpatialFeaturePlot(hcc01_bin100, features = "TLS_chemokine_score", max.cutoff = "q99") + theme(legend.position = "right", legend.title = element_blank()) + ggtitle("Chemokine")
SpatialFeaturePlot(hcc01_bin100, features = "TLS_TFH_score", max.cutoff = "q99") + theme(legend.position = "right", legend.title = element_blank()) + ggtitle("TFH")
SpatialFeaturePlot(hcc01_bin100, features = "TLS_TB_score", max.cutoff = "q99") + theme(legend.position = "right", legend.title = element_blank()) + ggtitle("TB")
SpatialFeaturePlot(hcc03_bin100, features = "TLS_chemokine_score", max.cutoff = "q99") + theme(legend.position = "right", legend.title = element_blank()) + ggtitle("Chemokine")
SpatialFeaturePlot(hcc03_bin100, features = "TLS_TFH_score", max.cutoff = "q99") + theme(legend.position = "right", legend.title = element_blank()) + ggtitle("TFH")
SpatialFeaturePlot(hcc03_bin100, features = "TLS_TB_score", max.cutoff = "q99") + theme(legend.position = "right", legend.title = element_blank()) + ggtitle("TB")


# B cells enriched in OF high spots----
HCC_genes_exp <- aggregate(t(HCC_seu@assays$RNA@data), list(Cluster = HCC_seu$Sub_ClusterNew), mean)
row.names(HCC_genes_exp) <- HCC_genes_exp$Cluster
HCC_genes_exp$Cluster <- c()
HCC_genes_exp <- t(HCC_genes_exp)
HCC_genes_prop <- aggregate(t(HCC_seu@assays$RNA@data), list(Cluster = HCC_seu$Sub_ClusterNew), function(x){sum(x > 0) / length(x)})
HCC_genes_prop$Cluster <- c()
HCC_genes_prop <- t(HCC_genes_prop)
OF_genes <- unique(c(
  row.names(HCC_genes_exp)[HCC_genes_prop[,"FOLR2+ TAM1"] > 0.3 & HCC_genes_exp[,"FOLR2+ TAM1"] > 0.75],
  row.names(HCC_genes_exp)[HCC_genes_prop[,"PLVAP+ EC"] > 0.3 & HCC_genes_exp[,"PLVAP+ EC"] > 0.75],
  row.names(HCC_genes_exp)[HCC_genes_prop[,"POSTN+ CAF"] > 0.3 & HCC_genes_exp[,"POSTN+ CAF"] > 0.75]
))

hcc01_OF_high <- hcc01_bin100@meta.data %>% filter(OF_high == TRUE) %>% pull(cell)
hcc01_degenes <- LIMMA(
  hcc01_bin100@assays$SCT@data,
  factor(hcc01_bin100$cell %in% hcc01_OF_high)
)
hcc03_OF_high <- hcc03_bin100@meta.data %>% filter(OF_high == TRUE) %>% pull(cell)
hcc03_degenes <- LIMMA(
  hcc03_bin100@assays$SCT@data,
  factor(hcc03_bin100$cell %in% hcc03_OF_high)
)

OF_high_upregulated <- intersect(row.names(HCC_seu), intersect(
  hcc01_degenes %>% filter(Symbol %ni% OF_genes, Grp == "C_TRUE", adj.P.Val < 0.05) %>% arrange(desc(abs(logFC))) %>% pull(Symbol) %>% head(50),
  hcc03_degenes %>% filter(Symbol %ni% OF_genes, Grp == "C_TRUE", adj.P.Val < 0.05) %>% arrange(desc(abs(logFC))) %>% pull(Symbol) %>% head(50)
))
cells_used <- HCC_seu@meta.data %>% filter(Sub_ClusterNew != "Bi-Potent") %>% pull(CellName)
p <- BoxHeatmap(as.matrix(HCC_seu@assays$RNA@data[OF_high_upregulated,cells_used]),factor(HCC_seu@meta.data[cells_used,"Global_Cluster"], levels = c("B cell","T cell","ILC","Mast","Mononuclear","Fibroblast","Endothelium","Hepatocyte")))
pdf("./04.figures/05.BGI_OF_high_upregulated.pdf", width = 8, height = 2.5)
grid.draw(p)
dev.off()

p <- data.frame(Bcell = hcc01_bin100@assays$predictions@data["B",], OF_high = hcc01_bin100$OF_high) %>%
  ggplot(aes(x = OF_high, y = Bcell)) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .5, outlier.size = 1.2, outlier.stroke = .25) +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  labs(y = "Predicted proportions of B cells", x = "") +
  theme_cowplot(font_size = 8) +
  RotatedAxis()
ggsave(p, file = "./04.figures/05.BGI_OF_high_Bcell_enriched.pdf", width = 1.75, height = 3.5)

hcc01_OF_high_neighbor <- FindNeighborSpots(hcc01_bin100, hcc01_OF_high)
hcc03_OF_high_neighbor <- FindNeighborSpots(hcc03_bin100, hcc03_OF_high)
hcc01_plot.list <- hcc03_plot.list <- list()
for(gene in c("CD79A","ITM2C","JCHAIN")){
  hcc01_plot.list[[gene]] <- SpatialFeaturePlot(hcc01_bin100 %>% subset(cell %in% hcc01_OF_high_neighbor), features = gene, min.cutoff = "q1", max.cutoff = 'q99')
  hcc01_plot.list[[gene]] <- ggAISpatial(hcc01_plot.list[[gene]])
  hcc03_plot.list[[gene]] <- SpatialFeaturePlot(hcc03_bin100 %>% subset(cell %in% hcc03_OF_high_neighbor), features = gene, min.cutoff = "q1", max.cutoff = 'q99')
  hcc03_plot.list[[gene]] <- ggAISpatial(hcc03_plot.list[[gene]])
}
p <- plot_grid(plotlist = hcc01_plot.list, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_hcc01_OF_high_Bcell_markers.pdf", width = 16, height = 4)
p <- plot_grid(plotlist = hcc03_plot.list, nrow = 1)
ggsave(p, file = "./04.figures/05.BGI_hcc03_OF_high_Bcell_markers.pdf", width = 16, height = 4)

# >Save rds files----
saveRDS(hcc01_bin100, file = "./02.processed_data/BGI/hcc01_bin100_seurat_processed.rds")
saveRDS(hcc03_bin100, file = "./02.processed_data/BGI/hcc03_bin100_seurat_processed.rds")

# RNAFISH data----
RNAFISH <- data.frame(
  Tissue = rep(rep(c("B017","Fetal_16wk"), each = 6), 2),
  Level = rep(rep(c("low","high"), each = 3), 4),
  Group = rep(c("FOLR2","PLVAP"), each = 12),
  Exp = c(0.005,0.206,0.065, 0.159,0.170,0.053,
          0.047,0.013,0.143, 0.249,0.459,0.340,
          0.000,0.014,0.038, 0.008,0.044,0.061,
          0.003,0.000,0.000, 0.016,0.021,0.063)
)
RNAFISH$Group2 <- factor(paste0(RNAFISH$Tissue,"_",RNAFISH$Level),levels = c("B017_low","B017_high","Fetal_16wk_low","Fetal_16wk_high"))
p1 <- ggplot(RNAFISH %>% filter(Group == "FOLR2"), aes(x = Group2, y = Exp)) +
  geom_boxplot(aes(fill = Level)) +
  scale_fill_manual(values = c("low" = "#9FA0A3", "high" = "#F84141")) +
  labs(x = "", y = "FOLR2 Signature") +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
p2 <- ggplot(RNAFISH %>% filter(Group == "PLVAP"), aes(x = Group2, y = Exp)) +
  geom_boxplot(aes(fill = Level)) +
  scale_fill_manual(values = c("low" = "#9FA0A3", "high" = "#F84141")) +
  labs(x = "", y = "PLVAP Signature") +
  theme_cowplot(font_size = 8) +
  RotatedAxis()
p <- p1 + p2
ggsave(p, file = "./04.figures/05.RNAFISH_statistics.pdf", width = 3.5, height = 4)
