setwd("/work/lzy/project/onco_fetal/")
source("../utils/utils_plot.R")
source("../utils/utils_data_processing.R")
source("../utils/utils_color.R")
library(Seurat)
library(Polychrome)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(cowplot)
library(tidyverse)
library(ComplexHeatmap)

# >Load dataset----
Fib_fetal <- readRDS("./02.processed_data/Fib_fetal.rds")
Fib_fetal_color_panel <- c(
  "Fib" = "#F69459", "ACTA2+ Fib" = "#84B8D7", "COL1A1+ Fib" = "#8DC594",
  "BGN+ Fib" = "#68AB9F", "APOC3+ Fib" = "#B294C7", "ITM2C+ Fib" = "#815AA8",
  "Mesenchymal" = "#E93A3B", "PLEK+ Fib" = "#F58584"
)

# >Cluster and annotation----
# >>UMAP plot----
p <- ggplot(Fib_fetal@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 4, alpha = .8) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(name = "", values = Fib_fetal_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/02.Fib_fetal_UMAP.pdf", width = 8, height = 6.5)

# >>Marker genes heatmap----
Fib_fetal <- subset(Fib_fetal, subset = Sub_Cluster != "Mesenchymal")
Fib_fetal_markers <- FindAllMarkers(Fib_fetal)
write.csv(Fib_fetal_markers, file = "./Results/DEGenes/Fib_fetal_all_markers.csv")
genes_used <- Fib_fetal_markers %>% group_by(cluster) %>% filter(pct.2 < 0.8, pct.1 > 0.5) %>% slice_max(n = 25, order_by = avg_logFC) %>% pull(gene) %>% unique()
genes_used <- genes_used[-grep("^MT",genes_used)]
genes_used <- c(genes_used, "ACTA2", "PLEK")
markers_mean_exp <- 
  aggregate(
    as.matrix(t(Fib_fetal@assays$RNA@data[genes_used,])),
    list(Cluster = Fib_fetal$Sub_Cluster),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Cluster
markers_mean_exp$Cluster <- c()
markers_plot_matrix <- apply(markers_mean_exp, 2, zscore) %>% t()
color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), max(markers_plot_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
max.avg <- apply(markers_plot_matrix, 1, which.max)
gene_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(gene = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), cluster = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
  if(sum(max.avg == i) == 1){
    temp <- data.frame(gene = names(max.avg)[max.avg == i], cluster = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
}
markers_plot_matrix_quantile <- quantile(markers_plot_matrix, c(0.01, 0.98))
markers_plot_matrix <- pmax(markers_plot_matrix, markers_plot_matrix_quantile[1])
markers_plot_matrix <- pmin(markers_plot_matrix, markers_plot_matrix_quantile[2])
genes_labeled <- c("ACTA2","HBM","HBG1","COL1A1","COL1A2","COL3A1","JUNB","FOS","MEG3","BGN","HSPB1","HSPA1A","APOA1","ALB","APOC3","APOA2","HBB","HLA-C","CD99","KRT18","FCGRT","RAMP1","CTSC","ITM2C","HLA-B","NR1H4","DCN","IFITM3","ECM1","COLEC10","VIM","CALD1","PLEK","S100A4","FABP5","SCL25A5","RPL36")
p <- Heatmap(markers_plot_matrix[gene_order$gene,],
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 6),
             col = color_used,
             name = "Exp") +
  rowAnnotation(link = anno_mark(
    at = which(gene_order$gene %in% genes_labeled),
    labels = gene_order$gene[which(gene_order$gene %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
pdf("./04.figures/02.Fib_fetal_markers_heatmap.pdf", width = 3.5, height = 9.5)
draw(p)
dev.off()

# >>Expression of marker genes----
genes_used <- c("ACTA2","COL1A1","BGN","APOC3","ITM2C","CD99","HLA-C","S100A10","PLEK")
Fib_fetal.cells.used <- Fib_fetal@meta.data %>% filter(Sub_Cluster != "Mesenchymal") %>% pull(CellID)
Fib_fetal.plot.data <- cbind(t(as.matrix(Fib_fetal@assays$RNA@data[genes_used,Fib_fetal.cells.used])), Fib_fetal@meta.data[Fib_fetal.cells.used, c("UMAP_1","UMAP_2")])
Fib_fetal.plot.data <- melt(Fib_fetal.plot.data, id.vars = c("UMAP_1","UMAP_2"))
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
p <- 
  ggplot(Fib_fetal.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 3) +
  facet_wrap(~variable, nrow = 3) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/02.Fib_fetal_markers_UMAP.pdf", width = 16, height = 15.5)

# >>Proportion----
# Cluster by patient
p <- 
  ggplot(Fib_fetal@meta.data, aes(x = Sub_Cluster, fill = PatientID)) +
  geom_bar(position = "fill", alpha = 0.8, color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 16) + 
  scale_fill_manual(values = c54[25:28]) +
  scale_x_discrete(limits = sort(unique(Fib_fetal@meta.data$Sub_Cluster))) +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(override.aes = list(size = 6), ncol = 1)) 
ggsave(p, file = "./04.figures/02.Fib_fetal_patient_proportions.pdf", width = 5, height = 6)
# Ro/e
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:6])
p <- 
  ROIE(table(Fib_fetal@meta.data[,c("Sub_Cluster","PatientID")])) %>% melt() %>% filter(Var1 != "Mesenchymal") %>% mutate(value = pmin(value, 3)) %>% mutate(Var1 = factor(as.character(Var1), levels = c("PLEK+ Fib","Fib","BGN+ Fib","ACTA2+ Fib","COL1A1+ Fib","APOC3+ Fib","ITM2C+ Fib"))) %>% 
  ggplot(data = ., aes(Var2, forcats::fct_rev(Var1), fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = myPalette(100)) +
  labs(x = "", y = "") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.line = element_blank(), 
        axis.ticks = element_blank()) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 12))
ggsave(p, file = "./04.figures/02.Fib_fetal_RoIe.pdf", width = 4, height = 4)

# >Trajectory----
# >>Monocle3----
library(monocle3)
library(SeuratWrappers)
library(magrittr)
Fib_fetal_UMAP.cds <- as.cell_data_set(Fib_fetal)
Fib_fetal_UMAP.cds <- cluster_cells(Fib_fetal_UMAP.cds, k = 30, partition_qval = 0.05, resolution = 0.01, reduction_method = "UMAP")
colData(Fib_fetal_UMAP.cds)$Cluster <- colData(Fib_fetal_UMAP.cds)$Sub_Cluster
Fib_fetal_UMAP.cds@clusters$UMAP$clusters <- as.character(Fib_fetal_UMAP.cds@colData$Cluster)
names(Fib_fetal_UMAP.cds@clusters$UMAP$clusters) <- Fib_fetal_UMAP.cds@colData$index
plot_grid(
  plot_cells(
    Fib_fetal_UMAP.cds,
    color_cells_by = "Cluster",
    show_trajectory_graph = FALSE),
  plot_cells(
    Fib_fetal_UMAP.cds, 
    color_cells_by = "partition", 
    show_trajectory_graph = FALSE)
)

Fib_fetal_UMAP.cds <- learn_graph(Fib_fetal_UMAP.cds, learn_graph_control = list(minimal_branch_len = 3))
Fib_fetal_UMAP.cds <- order_cells(Fib_fetal_UMAP.cds)
save(Fib_fetal_UMAP.cds, file = "./03.results/Trajectory/Fib_fetal_Monocle3.rda")

p <- ggplot(Fib_fetal@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 4, alpha = .8) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(name = "", values = Fib_fetal_color_panel) +
  theme(legend.position = "none")
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/02.Fib_fetal_monocle3_UMAP.pdf", width = 4, height = 4)

p1 <- plot_cells(
  Fib_fetal_UMAP.cds, 
  color_cells_by = "Sub_Cluster", 
  cell_size = 0,
  trajectory_graph_segment_size = 1,
  label_groups_by_cluster = FALSE, 
  label_cell_groups = FALSE,
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  scale_color_manual(values = Fib_fetal_color_panel) +
  labs(x = "UMAP1", y = "UMAP2")
p1 <- plot_grid(p1 + theme(legend.position = "none"))
ggsave(p1, file = "./04.figures/02.Fib_fetal_monocle3_UMAP_trajectory.pdf", width = 4, height = 4)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
p2 <- plot_cells(
  Fib_fetal_UMAP.cds, 
  color_cells_by = "pseudotime",
  cell_size = 2,
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  scale_color_gradientn(colours = myPalette(100)) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1, "cm"))
p2.pdf <- ggAIplot(p2 + theme(legend.position = "none"))
legend <- get_legend(p2)
p <- plot_grid(p2.pdf, legend, nrow = 1, rel_widths = c(4,1))
ggsave(p, file = "./04.figures/02.Fib_fetal_monocle3_UMAP_pseudotime.pdf", width = 5, height = 4)

Fib_fetal$Monocle3_pseudotime <- Fib_fetal_UMAP.cds@principal_graph_aux@listData$UMAP$pseudotime

p <- Fib_fetal@meta.data %>% filter(Sub_Cluster != "Mesenchymal") %>% mutate(Grp = plyr::revalue(Sub_Cluster, replace = c(
  "Fib" = "Grp1","ACTA2+ Fib" = "Grp3","PLEK+ Fib" = "Grp1","BGN+ Fib" = "Grp2","APOC3+ Fib" = "Grp4","ITM2C+ Fib" = "Grp4","COL1A1+ Fib" = "Grp3"
))) %>% 
  mutate(Grp = factor(Grp, levels = c("Grp1","Grp2","Grp3","Grp4"))) %>%
  mutate(Sub_Cluster = factor(as.character(Sub_Cluster),levels = rev(c("PLEK+ Fib","Fib","ACTA2+ Fib","BGN+ Fib","COL1A1+ Fib","APOC3+ Fib","ITM2C+ Fib")))) %>%
  ggplot(aes(x = Monocle3_pseudotime, y = Sub_Cluster)) +
  geom_boxplot(aes(color = Sub_Cluster)) +
  geom_point(aes(color = Sub_Cluster), position = position_jitter(h=0.2), size = .8) +
  facet_grid(Grp~., scales = "free", space = "free")  +
  scale_color_manual(name = "", values = Fib_fetal_color_panel) +
  theme_bw(base_size = 10) +
  ylab("") + xlab("Monocle3 pseudotime") +
  theme(legend.position = "none", 
        strip.text = element_blank(),
        strip.background = element_blank())
ggsave(p, file = "./04.figures/02.Fib_fetal_Monocle3_pseudotime_boxplot.pdf", width = 5, height = 4)

# >>SCENIC (server)----
library(dplyr)
library(SCENIC)
library(Seurat)
library(SCopeLoomR)
# expression matrix
cluster_used <- c("ACTA2+ Fib", "COL1A1+ Fib", "BGN+ Fib", "APOC3+ Fib", "ITM2C+ Fib", "PLEK+ Fib")
cells_used <- Fib_fetal@meta.data %>% filter(Sub_Cluster %in% cluster_used) %>% pull(CellID)
exprMat <- as.matrix(Fib_fetal@assays$RNA@data[,cells_used])
setwd("./03.results/SCENIC/Fib_fetal_tpm/")
dir.create("./int")
# cell information
cellInfo <- Fib_fetal@meta.data[cells_used,]
saveRDS(cellInfo, file="./int/cellInfo.Rds")
# initialize settings
org <- "hgnc"
dbDir <- "/work/lzy/tools/SCENIC/"
myDatasetTitle <- "SCENIC analysis of Fib fetal"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
scenicOptions@inputDatasetInfo$cellInfo <- "./int/cellInfo.Rds"
saveRDS(scenicOptions, file="./int/scenicOptions.Rds") 
scenicOptions <- readRDS("./int/scenicOptions.Rds")
# gene filter
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
# Correlation
runCorrelation(exprMat_filtered, scenicOptions)
# Genie3
runGenie3(exprMat_filtered, scenicOptions)
# build and score the GRN
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)
runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# >>>>Cluster-specific regulons heatmap----
library(SCENIC)
library(AUCell)
setwd("./03.results/SCENIC/Fib_fetal_tpm/")
scenicOptions <- readRDS("./int/scenicOptions.Rds")
cellInfo <- data.frame(CellType=Idents(Fib_fetal %>% subset(subset = Sub_Cluster %in% c("Fib", "ACTA2+ Fib", "COL1A1+ Fib", "BGN+ Fib", "APOC3+ Fib", "ITM2C+ Fib", "PLEK+ Fib"))))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

FindDEGenes(
  expression_matrix = as.matrix(regulonAUC@assays@data$AUC),
  groupid = cellInfo[colnames(regulonAUC@assays@data$AUC),"CellType"],
  out.prefix = "./Fib_fetal_SCENIC/Fib_fetal",
  cutoff = 1,
  logFC = 0
)

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
topRegulators_used <- topRegulators %>% filter(CellType %in% cluster_used) %>% group_by(CellType) %>% arrange(desc(RelativeActivity), .by_group = T) %>% top_n(15, wt = RelativeActivity) %>% pull(Regulon) %>% unique() %>% as.character()
Fib_fetal_SCENIC_degenes <- read.csv("./Fib_fetal_SCENIC/Fib_fetal_cutoff1_all_de_genes.csv", row.names = 1)
regulons_used <- Fib_fetal_SCENIC_degenes %>% filter(Group %in% cluster_used) %>% group_by(Group) %>% filter(Exp.Per.Out < 0.8, Exp.Per.In > 0.2) %>% arrange(desc(AUC), .by_group = T) %>% top_n(20, wt = AUC) %>% pull(Symbol) %>% unique()
regulons_used <- unique(c(as.character(regulons_used), topRegulators_used))

regulons_mean_exp <- 
  aggregate(
    as.matrix(t(regulonAUC@assays@data$AUC[regulons_used,cells_used])),
    list(Cluster = Fib_fetal@meta.data[cells_used,"Sub_Cluster"]),
    mean)
row.names(regulons_mean_exp) <- regulons_mean_exp$Cluster
regulons_mean_exp$Cluster <- c()
regulons_plot_matrix <- apply(regulons_mean_exp, 2, zscore) %>% t()
max.avg <- apply(regulons_plot_matrix, 1, which.max)
regulons_order <- c()
for(i in 1:ncol(regulons_plot_matrix)){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(gene = names(sort(regulons_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), cluster = colnames(regulons_plot_matrix)[i], stringsAsFactors = F)
    regulons_order <- rbind(regulons_order, temp)
  }
  if(sum(max.avg == i) == 1){
    temp <- data.frame(gene = names(max.avg)[max.avg == i], cluster = colnames(regulons_plot_matrix)[i], stringsAsFactors = F)
    regulons_order <- rbind(regulons_order, temp)
  }
}
regulons_plot_matrix_quantile <- quantile(regulons_plot_matrix, c(0.02, 0.98))
regulons_plot_matrix <- pmax(regulons_plot_matrix, regulons_plot_matrix_quantile[1])
regulons_plot_matrix <- pmin(regulons_plot_matrix, regulons_plot_matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(regulons_plot_matrix), max(regulons_plot_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
p <- Heatmap(regulons_plot_matrix[regulons_order$gene,],
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10),
             col = color_used,
             name = "Regulon Activity")
pdf("../../../04.figures/02.Fib_fetal_SCENIC_regulon_activity_heatmap.pdf", width = 5, height = 10)
draw(p)
dev.off()

# >>>>Regulon Specificity Score(RSS)----
library(scFunctions)
binary_regulon <- loadInt(scenicOptions, "aucell_binary_nonDupl")
binary_regulon <- binary_regulon[onlyNonDuplicatedExtended(rownames(binary_regulon)),]
rrs_df <- calculate_rrs(Fib_fetal@meta.data[cells_used,], binary_regulon, "Sub_Cluster")
rrs_df <- rrs_df %>% subset(!grepl("extended", regulon))

rrs_df_sub <- rrs_df %>% subset(cell_type == "ITM2C+ Fib") %>% 
  arrange(desc(RSS))
rrs_df_sub <- rrs_df_sub %>% mutate(rank = as.numeric(rownames(rrs_df_sub)))
rrs_ranking_plot <- 
  ggplot(rrs_df_sub, aes(rank, RSS, label = regulon)) + 
  geom_point(data = subset(rrs_df_sub,rank > 12),color = "#619CFF", size = 1) + 
  geom_point(data = subset(rrs_df_sub,rank <= 12), color = "#F87269", size = .8) + 
  labs(x = "Rank", y = "Regulon Specificity Score", title = "Top regulons of\nPOSTN+ CAF") +
  theme_cowplot(font_size = 12, line_size = 1)
rrs_ranking_plot <- 
  ggAIplot(rrs_ranking_plot, width = 3, height = 5) + 
  geom_text_repel(data = subset(rrs_df_sub, rank <= 12),
                  nudge_x = 35 - subset(rrs_df_sub, rank <= 12)$rank,
                  segment.size = 0.2,
                  segment.color = "#F87269", color = "#F87269",
                  direction = "y", hjust = 0)
ggsave(rrs_ranking_plot, file = "../../../Figures/01.CAF_HCC_SCENIC_POSTN_RSS.pdf", width = 3, height = 5)