setwd("/Volumes/ZiyiLi/PostDoc_Project/Onco-fetal Project/")
source("./onco-fetal_utils.R")
library(Seurat)
library(MERINGUE)
library(Polychrome)
library(dplyr)
library(dendextend)
library(ggplot2)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)

# Load data----
HCC_seu <- readRDS("./Onco-fetal CAF Identification/HCC_Ankur_release.rds")
hcc01_bin100 <- readRDS("./HCC01 and HCC03 Stereo-seq/hcc01_bin100.rds")
hcc01_bin100 <- ScaleData(hcc01_bin100, features = row.names(hcc01_bin100))
hcc03_bin100 <- readRDS("./HCC01 and HCC03 Stereo-seq/hcc03_bin100.rds")
hcc03_bin100 <- ScaleData(hcc03_bin100, features = row.names(hcc03_bin100))

# Fig.4a----
hcc03_plot.list <- list()
for(gene in c("FOLR2","PLVAP","POSTN")){
  hcc03_plot.list[[gene]] <- SpatialFeaturePlot(hcc03_bin100, features = gene, pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
}

# Fig.4b----
OF_TAM_signature <- c("FOLR2","CD163","MRC1","TREM2")
OF_Endo_signature <- c("PLVAP","ACKR1","DLL4")
OF_CAF_signature <- c("POSTN","FAP","ENG","PTGDS")
VSMC_signature <- c("MYH11","PLN","RERGL")
panCAF_signature <- c("TAGLN","MYL9","ACTA2","TPM2")
hcc01_bin100@meta.data$OF_TAM_score <- colSums(hcc01_bin100@assays$SCT@scale.data[OF_TAM_signature,])
hcc01_bin100@meta.data$OF_Endo_score <- colSums(hcc01_bin100@assays$SCT@scale.data[OF_Endo_signature,])
hcc01_bin100@meta.data$OF_CAF_score <- colSums(hcc01_bin100@assays$SCT@scale.data[OF_CAF_signature,])
hcc01_bin100@meta.data$VSMC_score <- colSums(hcc01_bin100@assays$SCT@scale.data[VSMC_signature,])
hcc01_bin100@meta.data$panCAF_score <- colSums(hcc01_bin100@assays$SCT@scale.data[panCAF_signature,])
hcc03_bin100@meta.data$OF_TAM_score <- colSums(hcc03_bin100@assays$SCT@scale.data[OF_TAM_signature,])
hcc03_bin100@meta.data$OF_Endo_score <- colSums(hcc03_bin100@assays$SCT@scale.data[OF_Endo_signature,])
hcc03_bin100@meta.data$OF_CAF_score <- colSums(hcc03_bin100@assays$SCT@scale.data[OF_CAF_signature,])
hcc03_bin100@meta.data$VSMC_score <- colSums(hcc03_bin100@assays$SCT@scale.data[VSMC_signature,])
hcc03_bin100@meta.data$panCAF_score <- colSums(hcc03_bin100@assays$SCT@scale.data[panCAF_signature,])
for(i in c("OF_Endo","OF_CAF","OF_TAM","VSMC","panCAF")){
  cutoff <- quantile(hcc01_bin100@meta.data[,paste0(i,"_score")], .975)
  high_cells <- hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data[,paste0(i,"_score")] > cutoff]
  high_neighbors <- FindNeighborSpots(hcc01_bin100, high_cells)
  hcc01_bin100@meta.data[,paste0(i,"_high")] <- FALSE
  hcc01_bin100@meta.data[high_neighbors,paste0(i,"_high")] <- TRUE
}
for(i in c("OF_Endo","OF_CAF","OF_TAM","VSMC","panCAF")){
  cutoff <- quantile(hcc03_bin100@meta.data[,paste0(i,"_score")], .975)
  high_cells <- hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data[,paste0(i,"_score")] > cutoff]
  high_neighbors <- FindNeighborSpots(hcc03_bin100, high_cells) 
  hcc03_bin100@meta.data[,paste0(i,"_high")] <- FALSE
  hcc03_bin100@meta.data[high_neighbors,paste0(i,"_high")] <- TRUE
}
cells_highlight <- list(
  hcc01 = list(
    OF_TAM_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_TAM_high],
    OF_CAF_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_CAF_high],
    OF_Endo_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$OF_Endo_high],
    VSMC_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$VSMC_high],
    panCAF_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$panCAF_high]),
  hcc03 = list(
    OF_TAM_high = hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_TAM_high],
    OF_CAF_high = hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_CAF_high],
    OF_Endo_high = hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$OF_Endo_high],
    VSMC_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$VSMC_high],
    panCAF_high = hcc01_bin100@meta.data$cell[hcc01_bin100@meta.data$panCAF_high])
)
cells_highlight$hcc01$OF_high <- intersect(intersect(cells_highlight$hcc01$OF_TAM_high, cells_highlight$hcc01$OF_CAF_high),cells_highlight$hcc01$OF_Endo_high)
hcc01_bin100@meta.data[,"OF_high"] <- FALSE
hcc01_bin100@meta.data[cells_highlight$hcc01$OF_high,"OF_high"] <- TRUE
cells_highlight$hcc03$OF_high <- intersect(intersect(cells_highlight$hcc03$OF_TAM_high, cells_highlight$hcc03$OF_CAF_high),cells_highlight$hcc03$OF_Endo_high)
hcc03_bin100@meta.data[,"OF_high"] <- FALSE
hcc03_bin100@meta.data[cells_highlight$hcc03$OF_high,"OF_high"] <- TRUE
plot.list <- list()
for(i in names(cells_highlight$hcc03)){
  plot.list[[i]] <- SpatialDimPlot(hcc03_bin100, cols.highlight = c("#F84141","#9FA0A3"),
                                               cells.highlight = cells_highlight$hcc03[[i]],
                                               pt.size.factor = 1.6, stroke = .1) + ggtitle(i)
}

# Fig.4c----
ggplot(hcc03_bin100@meta.data, aes(x = OF_Endo_high, y = log2(OF_CAF_score - min(OF_CAF_score) + 1))) +
  geom_boxplot(aes(fill = OF_Endo_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "POSTN+ CAF score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "PLVAP+ Endo low","TRUE" = "PLVAP+ Endo high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
ggplot(hcc03_bin100@meta.data, aes(x = OF_TAM_high, y = log2(OF_CAF_score - min(OF_CAF_score) + 1))) +
  geom_boxplot(aes(fill = OF_TAM_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "POSTN+ CAF score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "FOLR2+ TAM1 low","TRUE" = "FOLR2+ TAM1 high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()

cells_used <- hcc03_bin100@meta.data %>% filter(hcc03_bin100@meta.data$OF_High_NeighIndex > 0) %>% pull(cell)
posw <- hcc03_bin100@images$slice1@coordinates[cells_used,c("row","col")]
weight <- getSpatialNeighbors(posw, filterDist = 3)
gexp1 <- hcc03_bin100@meta.data[cells_used,"OF_CAF_score"] %>% setNames(cells_used)
gexp2 <- hcc03_bin100@meta.data[cells_used,"OF_Endo_score"] %>% setNames(cells_used)
gexp3 <- hcc03_bin100@meta.data[cells_used,"OF_TAM_score"] %>% setNames(cells_used)
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, ncores = 12)
spatialCrossCor(gexp1, gexp3, weight)
spatialCrossCorTest(gexp1, gexp3, weight, ncores = 12)

# Fig.4d----
hcc01_bin10 <- readRDS("./HCC01 and HCC03 Stereo-seq/hcc01_bin10.rds")
hcc01_OF_high <- hcc01_bin100@meta.data %>% filter(OF_high == TRUE) %>% pull(cell)
hcc01_bin100@images$slice1@coordinates[,c("coord_col","coord_row")] <- stringr::str_split_fixed(row.names(hcc01_bin100@images$slice1@coordinates), pattern = ":|_", 3)[,c(2,3)]
hcc01_bin100@images$slice1@coordinates$coord_col <- as.numeric(hcc01_bin100@images$slice1@coordinates$coord_col)
hcc01_bin100@images$slice1@coordinates$coord_row <- as.numeric(hcc01_bin100@images$slice1@coordinates$coord_row)
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
SpatialDimPlot(hcc01_bin10_subset1, group.by = "OF_high")
SpatialDimPlot(hcc01_bin10_subset1, group.by = "OF_Annotation", pt.size.factor = 2.5, stroke = .1) +
  guides(fill = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

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
SpatialDimPlot(hcc01_bin10_subset2, group.by = "OF_high", pt.size.factor = 4, stroke = .1)
SpatialDimPlot(hcc01_bin10_subset2, group.by = "OF_Annotation", pt.size.factor = 4, stroke = .1) +
  guides(fill = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

hcc03_bin10 <- readRDS("./HCC01 and HCC03 Stereo-seq/hcc03_bin10.rds")
hcc03_bin100@images$slice1@coordinates[,c("coord_col","coord_row")] <- stringr::str_split_fixed(row.names(hcc03_bin100@images$slice1@coordinates), pattern = ":|_", 3)[,c(2,3)]
hcc03_bin100@images$slice1@coordinates$coord_col <- as.numeric(hcc03_bin100@images$slice1@coordinates$coord_col)
hcc03_bin100@images$slice1@coordinates$coord_row <- as.numeric(hcc03_bin100@images$slice1@coordinates$coord_row)
hcc03_bin10@images$slice1@coordinates[,c("coord_col","coord_row")] <- stringr::str_split_fixed(row.names(hcc03_bin10@images$slice1@coordinates), pattern = ":|_", 3)[,c(2,3)]
hcc03_bin10@images$slice1@coordinates$coord_col <- as.numeric(hcc03_bin10@images$slice1@coordinates$coord_col)
hcc03_bin10@images$slice1@coordinates$coord_row <- as.numeric(hcc03_bin10@images$slice1@coordinates$coord_row)
hcc03_OF_high <- hcc03_bin100@meta.data %>% filter(OF_high == TRUE) %>% pull(cell)
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
SpatialDimPlot(hcc03_bin10_subset1, group.by = "OF_high", pt.size.factor = 2, stroke = .1)
SpatialDimPlot(hcc03_bin10_subset1, group.by = "OF_Annotation", pt.size.factor = 2, stroke = .1) +
    guides(fill = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

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
SpatialDimPlot(hcc03_bin10_subset2, group.by = "OF_high", pt.size.factor = 1.5, stroke = .1)
SpatialDimPlot(hcc03_bin10_subset2, group.by = "OF_Annotation", pt.size.factor = 1.5, stroke = .1) +
    guides(fill = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

# Fig.4f----
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
ggplot(RNAFISH %>% filter(Group == "FOLR2"), aes(x = Group2, y = Exp)) +
  geom_boxplot(aes(fill = Level)) +
  scale_fill_manual(values = c("low" = "#9FA0A3", "high" = "#F84141")) +
  labs(x = "", y = "FOLR2 Signature") +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
ggplot(RNAFISH %>% filter(Group == "PLVAP"), aes(x = Group2, y = Exp)) +
  geom_boxplot(aes(fill = Level)) +
  scale_fill_manual(values = c("low" = "#9FA0A3", "high" = "#F84141")) +
  labs(x = "", y = "PLVAP Signature") +
  theme_cowplot(font_size = 8) +
  RotatedAxis()

# Extended Data Fig.5a-b----
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
p1 + p2

# Extended Data Fig.5c----
SpatialFeaturePlot(hcc01_bin100, features = c("nFeature_Spatial"), pt.size.factor = 1.5, stroke = 0.2, min.cutoff = "q1")
SpatialFeaturePlot(hcc03_bin100, features = c("nFeature_Spatial"), pt.size.factor = 1.5, stroke = 0.2, min.cutoff = "q1")

# Extended Data Fig.5d----
hcc01_plot.list <- list()
for(gene in c("FOLR2","PLVAP","POSTN")){
  hcc01_plot.list[[gene]] <- SpatialFeaturePlot(hcc01_bin100, features = gene, pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
}

# Extended Data Fig.5e----
plot.list <- list()
for(i in names(cells_highlight$hcc01)){
  plot.list[[i]] <- ggAISpatial(SpatialDimPlot(hcc01_bin100, cols.highlight = c("#F84141","#9FA0A3"),
                                               cells.highlight = cells_highlight$hcc01[[i]],
                                               pt.size.factor = 1.6, stroke = .1)) + ggtitle(i)
}

# Extended Data Fig.5f----
hcc01_plot.list <- hcc03_plot.list <- list()
for(gene in c("CD68","CSF1R","PECAM1","ACTA2","ENG","CD4","FOXP3")){
  hcc01_plot.list[[gene]] <- SpatialFeaturePlot(hcc01_bin100, features = gene, pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
  hcc03_plot.list[[gene]] <- SpatialFeaturePlot(hcc03_bin100, features = gene, pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
}

# Extended Data Fig.5g----
ggplot(hcc01_bin100@meta.data, aes(x = OF_Endo_high, y = log2(VSMC_score - min(VSMC_score) + 1))) +
  geom_boxplot(aes(fill = OF_Endo_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "VSMC score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "PLVAP+ Endo low","TRUE" = "PLVAP+ Endo high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
ggplot(hcc01_bin100@meta.data, aes(x = OF_TAM_high, y = log2(VSMC_score - min(VSMC_score) + 1))) +
  geom_boxplot(aes(fill = OF_TAM_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "VSMC score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "FOLR2+ TAM1 low","TRUE" = "FOLR2+ TAM1 high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()

# Extended Data Fig.5h----
DotPlot(HCC_seu %>% subset(Global_Cluster == "Fibroblast"), features = c(OF_CAF_signature, VSMC_signature, panCAF_signature)) + 
  theme_cowplot(font_size = 7) +
  RotatedAxis()

# Extended Data Fig.5i----
hcc03_plot.list <- list(
  OF_TAM = SpatialFeaturePlot(hcc03_bin100, features = "OF_TAM_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99'),
  OF_Endo = SpatialFeaturePlot(hcc03_bin100, features = "OF_Endo_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99'),
  OF_CAF = SpatialFeaturePlot(hcc03_bin100, features = "OF_CAF_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99'),
  VSMC = SpatialFeaturePlot(hcc03_bin100, features = "VSMC_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99'),
  panCAF = SpatialFeaturePlot(hcc03_bin100, features = "panCAF_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
)

# Extended Data Fig.5j----
ggplot(hcc03_bin100@meta.data, aes(x = OF_Endo_high, y = log2(VSMC_score - min(VSMC_score) + 1))) +
  geom_boxplot(aes(fill = OF_Endo_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "VSMC score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "PLVAP+ Endo low","TRUE" = "PLVAP+ Endo high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
ggplot(hcc03_bin100@meta.data, aes(x = OF_TAM_high, y = log2(VSMC_score - min(VSMC_score) + 1))) +
  geom_boxplot(aes(fill = OF_TAM_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "VSMC score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "FOLR2+ TAM1 low","TRUE" = "FOLR2+ TAM1 high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()

ggplot(hcc03_bin100@meta.data, aes(x = OF_Endo_high, y = log2(panCAF_score - min(panCAF_score) + 1))) +
  geom_boxplot(aes(fill = OF_Endo_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "panCAF score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "PLVAP+ Endo low","TRUE" = "PLVAP+ Endo high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
ggplot(hcc03_bin100@meta.data, aes(x = OF_TAM_high, y = log2(panCAF_score - min(panCAF_score) + 1))) +
  geom_boxplot(aes(fill = OF_TAM_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "panCAF score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "FOLR2+ TAM1 low","TRUE" = "FOLR2+ TAM1 high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()

cells_used <- hcc03_bin100@meta.data %>% filter(hcc03_bin100@meta.data$OF_High_NeighIndex > 0) %>% pull(cell)
posw <- hcc03_bin100@images$slice1@coordinates[cells_used,c("row","col")]
weight <- getSpatialNeighbors(posw, filterDist = 3)
gexp2 <- hcc03_bin100@meta.data[cells_used,"OF_Endo_score"] %>% setNames(cells_used)
gexp3 <- hcc03_bin100@meta.data[cells_used,"OF_TAM_score"] %>% setNames(cells_used)
gexp4 <- hcc03_bin100@meta.data[cells_used,"VSMC_score"] %>% setNames(cells_used)
gexp5 <- hcc03_bin100@meta.data[cells_used,"panCAF_score"] %>% setNames(cells_used)
spatialCrossCor(gexp4, gexp2, weight)
spatialCrossCorTest(gexp4, gexp2, weight, ncores = 12)
spatialCrossCor(gexp4, gexp3, weight)
spatialCrossCorTest(gexp4, gexp3, weight, ncores = 12)
spatialCrossCor(gexp5, gexp2, weight)
spatialCrossCorTest(gexp5, gexp2, weight, ncores = 12)
spatialCrossCor(gexp5, gexp3, weight)
spatialCrossCorTest(gexp5, gexp3, weight, ncores = 12)

# Extended Data Fig.6a----
core_ligands <- c("LTA","PGF","BMP2","ITGAM","VEGFA","ADAM17","COL2A1","COL4A1","TGFB3","BTLA")
hcc03_bin100@meta.data$Core_Ligand_score <- colSums(hcc03_bin100@assays$SCT@scale.data[core_ligands,])
SpatialFeaturePlot(hcc03_bin100, features = "Core_Ligand_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')

# Extended Data Fig.6b----
prioritized_ligands <- c("LTA","ADAM17","PGF","COL2A1","COL4A1","BMP2","ITGAM","TGFB3","VEGFA","BTLA")
hcc01_bin100@meta.data[,"Prioritized_ligands_score"] <- colSums(hcc01_bin100@assays$SCT@scale.data[prioritized_ligands,])
hcc03_bin100@meta.data[,"Prioritized_ligands_score"] <- colSums(hcc03_bin100@assays$SCT@scale.data[prioritized_ligands,])
ggplot(hcc01_bin100@meta.data, aes(x = OF_high, y = log2(Prioritized_ligands_score - min(Prioritized_ligands_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "Prioritized ligands score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()
ggplot(hcc03_bin100@meta.data, aes(x = OF_high, y = log2(Prioritized_ligands_score - min(Prioritized_ligands_score) + 1))) +
  geom_boxplot(aes(fill = OF_high), alpha = .8, outlier.alpha = .8, outlier.size = .5) +
  labs(x = "", y = "Prioritized ligands score") +
  ggsignif::geom_signif(comparisons = list(c("TRUE","FALSE")), textsize = 3) +
  scale_fill_manual(values = c("TRUE" = "#F84141", "FALSE" = "#9FA0A3")) +
  scale_x_discrete(labels = c("FALSE" = "Onco-fetal low","TRUE" = "Onco-fetal high")) +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "none") +
  RotatedAxis()

# Extended Data Fig.6c----
hcc03_core_ligand_plot.list <- list()
for(gene in c("PGF","BMP2","ITGAM","TGFB3")){
  hcc03_core_ligand_plot.list[[gene]] <- SpatialFeaturePlot(hcc03_bin100, features = gene, pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
}

plot.list <- list()
core_ligands <- c("PGF","BMP2","ITGAM","TGFB3")
for(gene in core_ligands){
  plot.list[[gene]] <- 
    data.frame(exp = as.vector(hcc03_bin100@assays$SCT@scale.data[gene,]),
               OF_niche = hcc03_bin100@meta.data$OF_high) %>% 
    filter(exp > 0) %>%
    ggplot(aes(x = exp, color = OF_niche)) +
    geom_density(aes(fill = OF_niche), alpha = .2) +
    ggtitle(gene) +
    scale_color_manual(values = c("TRUE" = "#F84141","FALSE" = "#9FA0A3")) +
    scale_fill_manual(values = c("TRUE" = "#F84141","FALSE" = "#9FA0A3")) +
    theme_cowplot(font_size = 7) +
    theme(legend.position = "none")
}

cells_used <- hcc03_bin100@meta.data %>% filter(hcc03_bin100@meta.data$OF_High_NeighIndex > 0) %>% pull(cell)
posw <- hcc03_bin100@images$slice1@coordinates[cells_used,c("row","col")]
weight <- getSpatialNeighbors(posw, filterDist = 3)

gexp1 <- hcc03_bin100@assays$SCT@data["PGF",cells_used]
gexp2 <- rowSums(hcc03_bin100@meta.data[cells_used,c("OF_TAM_score","OF_Endo_score","OF_CAF_score")])
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, plot=TRUE)

gexp1 <- hcc03_bin100@assays$SCT@data["TGFB3",cells_used]
gexp2 <- rowSums(hcc03_bin100@meta.data[cells_used,c("OF_TAM_score","OF_Endo_score","OF_CAF_score")])
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, plot=TRUE)

gexp1 <- hcc03_bin100@assays$SCT@data["BMP2",cells_used]
gexp2 <- rowSums(hcc03_bin100@meta.data[cells_used,c("OF_TAM_score","OF_Endo_score","OF_CAF_score")])
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, plot=TRUE)

gexp1 <- hcc03_bin100@assays$SCT@data["ITGAM",cells_used]
gexp2 <- rowSums(hcc03_bin100@meta.data[cells_used,c("OF_TAM_score","OF_Endo_score","OF_CAF_score")])
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, plot=TRUE)

# Extended Data Fig.6d-g----
color_panel <- c("CXCL12" = "red3", "CXCR4" = "green3", "POSTN" = "yellow", "FOLR2" = "blue3", "DLL4" = "red3",  "NOTCH3" = "green3",  "PLVAP" = "blue3")
plot.list <- list(hcc01 = list(), hcc03 = list())
for(gene in c("CXCL12","CXCR4","POSTN","FOLR2","DLL4","NOTCH3","PLVAP")){
  plot.list[["hcc01"]][[gene]] <- SpatialFeaturePlot(hcc01_bin100, features = gene, pt.size.factor = 1.5, stroke = 0, min.cutoff = "q1", max.cutoff = 'q99') + scale_fill_gradient(low = "white", high = color_panel[gene])
  plot.list[["hcc03"]][[gene]] <- SpatialFeaturePlot(hcc03_bin100, features = gene, pt.size.factor = 1.5, stroke = 0, min.cutoff = "q1", max.cutoff = 'q99') + scale_fill_gradient(low = "white", high = color_panel[gene])
}

cells_used <- hcc01_bin100@meta.data %>% filter(hcc01_bin100@meta.data$OF_High_NeighIndex > 0) %>% pull(cell)
posw <- hcc01_bin100@images$slice1@coordinates[cells_used,c("row","col")]
weight <- getSpatialNeighbors(posw, filterDist = 3)
gexp1 <- hcc01_bin100@assays$SCT@data["CXCL12",cells_used]
gexp2 <- hcc01_bin100@assays$SCT@data["CXCR4",cells_used]
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, plot=TRUE)

gexp1 <- hcc01_bin100@assays$SCT@data["DLL4",cells_used]
gexp2 <- hcc01_bin100@assays$SCT@data["NOTCH3",cells_used]
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, plot=TRUE)

cells_used <- hcc03_bin100@meta.data %>% filter(hcc03_bin100@meta.data$OF_High_NeighIndex > 0) %>% pull(cell)
posw <- hcc03_bin100@images$slice1@coordinates[cells_used,c("row","col")]
weight <- getSpatialNeighbors(posw, filterDist = 3)
gexp1 <- hcc03_bin100@assays$SCT@data["CXCL12",cells_used]
gexp2 <- hcc03_bin100@assays$SCT@data["CXCR4",cells_used]
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, plot=TRUE)

gexp1 <- hcc03_bin100@assays$SCT@data["DLL4",cells_used]
gexp2 <- hcc03_bin100@assays$SCT@data["NOTCH3",cells_used]
spatialCrossCor(gexp1, gexp2, weight)
spatialCrossCorTest(gexp1, gexp2, weight, plot=TRUE)

# Supplementary Fig.6a----
hcc01_bin100@meta.data$Cluster <- plyr::revalue(
  hcc01_bin100@active.ident,
  replace=c("0"="C2","1"="C6","2"="C4","3"="C1","4"="C3",
            "5"="C10","6"="C4","7"="C1","8"="C7","9"="C8",
            "10"="C5","11"="C9","12"="C4","13"="C10","14"="C11"))
hcc01_bin100@meta.data$Cluster <- factor(hcc01_bin100@meta.data$Cluster, levels = paste0("C",1:11))
hcc01_bin100@meta.data[,c("UMAP_1","UMAP_2")] <- hcc01_bin100@reductions$umap@cell.embeddings[,c(1,2)]
ggplot(hcc01_bin100@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 0.5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
SpatialDimPlot(hcc01_bin100, pt.size.factor = 1.5, stroke = 0, group.by = "Cluster")

hcc03_bin100@meta.data$Cluster <- plyr::revalue(
  hcc03_bin100@meta.data$SCT_snn_res.0.6,
  replace=c("0"="C1","1"="C2","2"="C7","3"="C3","4"="C4",
            "5"="C8","6"="C7","7"="C12","8"="C6","9"="C5",
            "10"="C9","11"="C8","12"="C2","13"="C2","14"="C10",
            "15"="C11"))
hcc03_bin100@meta.data$Cluster <- factor(hcc03_bin100@meta.data$Cluster, levels = paste0("C",1:12))
hcc03_bin100@meta.data[,c("UMAP_1","UMAP_2")] <- hcc03_bin100@reductions$umap@cell.embeddings[,c(1,2)]
ggplot(hcc03_bin100@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 0.5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
SpatialDimPlot(hcc03_bin100, pt.size.factor = 1.5, stroke = 0, group.by = "Cluster")

# Supplementary Fig.6b----
features_used <- c("OF_CAF_score","OF_TAM_score","OF_Endo_score")
hcc01.plot.data <- melt(hcc01_bin100@meta.data[,c(features_used,"UMAP_1","UMAP_2")],id.vars = c("UMAP_1","UMAP_2"))
hcc01.plot.data$value <- pmin(hcc01.plot.data$value, quantile(hcc01.plot.data$value,.995))
hcc01.plot.data$value <- pmax(hcc01.plot.data$value, quantile(hcc01.plot.data$value,.005))
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
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

features_used <- c("OF_CAF_score","OF_TAM_score","OF_Endo_score")
hcc03.plot.data <- melt(hcc03_bin100@meta.data[,c(features_used,"UMAP_1","UMAP_2")],id.vars = c("UMAP_1","UMAP_2"))
hcc03.plot.data$value <- pmin(hcc03.plot.data$value, quantile(hcc03.plot.data$value,.995))
hcc03.plot.data$value <- pmax(hcc03.plot.data$value, quantile(hcc03.plot.data$value,.005))
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
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