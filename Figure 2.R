setwd("/Volumes/ZiyiLi/PostDoc_Project/Onco-fetal Project/")
source("./onco-fetal_utils.R")
library(Seurat)
library(Polychrome)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)

# Load data----
CAF_HCC <- readRDS("./Onco-fetal CAF Identification/CAF_HCC_release.rds")
Fib_fetal <- readRDS("./Onco-fetal CAF Identification/Fib_fetal_release.rds")
# The "Fib" cluster in Fib_fetal data was "IGFBP3+ Fib" cluster in manuscript
# The "Mesenchymal" cluster in Fib_fetal data was "Mesothelial cell" cluster in manuscript
Fib_integrated <- readRDS("./Onco-fetal CAF Identification/Fib_integrated_release.rds")

# Fig.2a----
ggplot(CAF_HCC@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

# Fig.2b----
genes_used <- c("MYH11","HSPA6","SDC2","FAP","POSTN","COL4A1","ENG","APOA2","MT1M","ABCA8")
CAF_HCC.cells.used <- CAF_HCC$CellID
CAF_HCC.plot.data <- cbind(t(as.matrix(CAF_HCC@assays$RNA@data[genes_used,CAF_HCC.cells.used])), CAF_HCC@meta.data[CAF_HCC.cells.used, c("UMAP_1","UMAP_2")])
CAF_HCC.plot.data <- melt(CAF_HCC.plot.data, id.vars = c("UMAP_1","UMAP_2"))
CAF_HCC.plot.data$value <- pmin(CAF_HCC.plot.data$value, 3.5)
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p <- 
  ggplot(CAF_HCC.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 1) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))

# Fig.2c----
ggplot(Fib_fetal@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 4, alpha = .8) +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

# Fig.2d----
genes_used <- c("IGFBP3","ACTA2","COL1A1","BGN","APOC3","ITM2C","PLEK","MSLN","CD99","POSTN")
Fib_fetal.cells.used <- Fib_fetal@meta.data %>% pull(CellID)
Fib_fetal.plot.data <- cbind(t(as.matrix(Fib_fetal@assays$RNA@data[genes_used,Fib_fetal.cells.used])), Fib_fetal@meta.data[Fib_fetal.cells.used, c("UMAP_1","UMAP_2")])
Fib_fetal.plot.data <- melt(Fib_fetal.plot.data, id.vars = c("UMAP_1","UMAP_2"))
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
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

# Fig.2e----
ggplot(Fib_integrated@meta.data %>% filter(Source == "Healthy"), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 3, alpha = .8) +
  theme_cowplot(font_size = 7) +
  lims(x = c(-7.5,6), y = c(-6,5.5)) +
  NoAxes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")
ggplot(Fib_integrated@meta.data %>% filter(Source == "Fetal"), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 3, alpha = .8) +
  theme_cowplot(font_size = 7) +
  lims(x = c(-7.5,6), y = c(-6,5.5)) +
  NoAxes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")
ggplot(Fib_integrated@meta.data %>% filter(Source == "Tumor"), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 3, alpha = .8) +
  theme_cowplot(font_size = 7) +
  lims(x = c(-7.5,6), y = c(-6,5.5)) +
  NoAxes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")
ggplot(Fib_integrated@meta.data %>% filter(Source == "Adj Normal"), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 3, alpha = .8) +
  theme_cowplot(font_size = 7) +
  lims(x = c(-7.5,6), y = c(-6,5.5)) + 
  NoAxes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")

# Fig.2f----
plot.data <- Fib_integrated@meta.data %>% filter(Integrated_cluster %in% c(3,4,7), Source != "Healthy") %>% droplevels()
ggplot(plot.data, aes(x = Integrated_cluster, fill = Sub_Cluster)) +
  geom_bar(position = "fill", color = NA) +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 1))
cross_tab <- table(plot.data[,c("Integrated_cluster","Sub_Cluster")])
cross_tab_prop <- cross_tab/rowSums(cross_tab)

# Fig.2g----
Fib_markers <- FindMarkers(Fib_fetal %>% subset(subset = Sub_Cluster != "Mesenchymal"), ident.1 = "ITM2C+ Fib", group.by = "Sub_Cluster", logfc.threshold = 0, min.pct = 0)
Fib_markers$gene <- row.names(Fib_markers)
CAF_markers <- FindMarkers(CAF_HCC %>% subset(subset = NTF == "Tumor"), ident.1 = "POSTN+ CAF", group.by = "Sub_Cluster", logfc.threshold = 0, min.pct = 0)
CAF_markers$gene <- row.names(CAF_markers)
genes_used <- intersect(Fib_markers$gene, CAF_markers$gene)
plot.df <- cbind(
  Fib_markers[genes_used,c("avg_log2FC","p_val_adj","pct.1","pct.2")],
  CAF_markers[genes_used,c("avg_log2FC","p_val_adj","pct.1","pct.2","gene")]
)
colnames(plot.df) <- c("logFC_fetal","p_fetal","pct_in_fetal","pct_our_fetal","logFC_tumor","p_tumor","pct_in_tumor","pct_out_tumor","Symbol")
plot.df$Sig <- "Not up-regulated"
plot.df[plot.df$logFC_fetal > 0 & plot.df$logFC_tumor > 0 & plot.df$p_fetal < 0.05 & plot.df$p_tumor < 0.05,"Sig"] <- "Both low up-regulated"
plot.df[plot.df$logFC_fetal > 0.5 & plot.df$logFC_tumor > 0.5 & plot.df$p_fetal < 0.05 & plot.df$p_tumor < 0.05,"Sig"] <- "Both high up-regulated"
plot.df[abs(plot.df$logFC_fetal) < 0.5 & plot.df$logFC_tumor > 0.5 & plot.df$p_fetal > 0.05 & plot.df$p_tumor < 0.05,"Sig"] <- "Tumor-specific up-regulated"
plot.df[plot.df$logFC_fetal < 0 & plot.df$p_fetal < 0.05 & plot.df$logFC_tumor < 0 & plot.df$p_tumor < 0.05,"Sig"] <- "Down-regulated"
plot.df$Label <- c()
labeled <- (plot.df$logFC_fetal > .5 & plot.df$logFC_tumor > .5)
labeled[plot.df$Symbol %in% c("PTGDS","LXN","LUM","CYP1B1","LTBP2","COL6A3","THBS2","CCDC80","TNFSF13B","GPNMB","PDGFRA","FBLN1","INHBA","CTGF","PTGIS","MMP2","TNFSF10","CTSS","VCAM1","IL32","THY1","ICAM1")] <- TRUE
plot.df[labeled,"Label"] <- plot.df[labeled,"Symbol"]
ggplot(plot.df %>% arrange(desc(Sig)), aes(x = logFC_fetal, y = logFC_tumor)) +
  geom_point(aes(color = Sig), size = 1.5, alpha = .8) +
  labs(x = "logFC of ITM2C+ Fib", y = "logFC of POSTN+ CAF") +
  scale_color_manual(name = "", 
                     values = c("Not up-regulated" = "lightgrey",
                                "Tumor-specific up-regulated" = "#58B44D",
                                "Both low up-regulated" = "orange",
                                "Both high up-regulated" = "#E51B1A",
                                "Down-regulated" = "#4480BF")) +
  xlim(-2,2) + ylim(-2,2) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = c(0,0.5), linetype = "dashed") + 
  geom_vline(xintercept = c(0,0.5), linetype = "dashed") +
  geom_text_repel(aes(label = Label)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 3)), keyheight = 1.5)

# Fig.2h----
# Analyzed in GEPIA2 (http://gepia2.cancer-pku.cn/#index)

# Extended Data Fig.2a----
ggplot(CAF_HCC@meta.data %>% filter(NTF == "Tumor"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  lims(x = c(-5,6.5), y = c(-7,5.5)) + 
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
ggplot(CAF_HCC@meta.data %>% filter(NTF == "Adj Normal"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  lims(x = c(-5,6.5), y = c(-7,5.5)) + 
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

# Extended Data Fig.2b----
CAF_HCC_markers <- FindAllMarkers(CAF_HCC)
genes_used <- CAF_HCC_markers %>% group_by(cluster) %>% filter(pct.2 < 0.8, pct.1 > 0.5) %>% slice_max(n = 8, order_by = avg_log2FC) %>% pull(gene) %>% unique()
genes_used <- unique(c(genes_used, "FAP", "SDC2", "MYH11", "CD36", "POSTN", "ENG", "PI16","CD34","HAS1","PDGFRA","MCAM"))
genes_used <- genes_used[!genes_used %in% c("CASQ2","PBXIP1","MRGPRF","GPRC5C","COX4I2","EPAS1","ZFAND2A","CACYBP","XIST","PTP4A3","SFRP4","PRSS23","POLR2L","IFI27","ZFP36","ADIRF","FTH1","MEG3","ITM2A","ADH1B","C7")]
gene.mean.matrix <- 
  aggregate(t(as.matrix(CAF_HCC@assays$RNA@data[genes_used,])),
            list(Cluster = CAF_HCC$Sub_Cluster), 
            mean)
gene.mean.zscore <- apply(gene.mean.matrix[,2:ncol(gene.mean.matrix)], 2, function(x){(x - mean(x)) / sd(x)})
row.names(gene.mean.zscore) <- gene.mean.matrix[,1]
gene.mean.zscore.df <- 
  data.frame(Group = gene.mean.matrix[,1], gene.mean.zscore, check.names = F) %>% 
  melt(id.vars = "Group", variable.name = "Gene", value.name = "Exp")
gene.mean.zscore <- t(gene.mean.zscore)
max.avg <- apply(gene.mean.zscore, 1, which.max)
gene_order <- c()
for(i in 1:ncol(gene.mean.zscore)){
  if(sum(max.avg == i) > 0){
    temp <- data.frame(Gene = names(sort(gene.mean.zscore[names(max.avg)[max.avg == i],i], decreasing = T)), Gene.Group = colnames(gene.mean.zscore)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
}
gene.per <- 
  aggregate(t(as.matrix(CAF_HCC@assays$RNA@data[genes_used,])),
            list(Group = CAF_HCC$Sub_Cluster),
            function(x){sum(x > 2) / length(x)}) %>% 
  melt(id.vars = "Group", variable.name = "Gene", value.name = "Per")
plot.data <- merge(merge(gene.mean.zscore.df, gene.per), gene_order)
plot.data$Gene <- as.character(plot.data$Gene)
plot.data$Gene.Group <- factor(plot.data$Gene.Group, levels = c("CAF","MYH11+ CAF","HSP+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF","Fib","MT1M+ Fib","ABCAB+ Fib"))
plot.data$Group <- factor(plot.data$Group, levels = rev(c("CAF","MYH11+ CAF","HSP+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF","Fib","MT1M+ Fib","ABCAB+ Fib")))
plot.data.order <- plot.data[order(plot.data$Gene.Group, plot.data$Exp),]
plot.data$Gene <- factor(plot.data$Gene, levels = unique(plot.data.order$Gene))
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
p <- ggplot(plot.data, aes(x = Group, y = Gene)) +
  geom_point(aes(size = Per, fill = Exp, color = Exp), shape = 21) +
  facet_grid(Gene.Group~., scales = "free", space = "free") +
  theme_minimal() + 
  scale_x_discrete(limits = rev(levels(plot.data$Group))) +
  theme(
    axis.ticks = element_line(size = .5),
    axis.ticks.length = unit(0.05,"in"),
    strip.text = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", angle = 90, 
                               hjust = 1, vjust = 0.5, size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "lightgrey",
      linetype = "dashed",
      size = 0.2
    )) + 
  labs(x = "", y = "") +
  scale_fill_gradientn("Exp",
                       colours = myColorPalette(100), 
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black", barwidth = 0.8)) +
  scale_color_gradientn("Exp",
                        colours = myColorPalette(100), 
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks.colour = "black", barwidth = 0.8)) +
  scale_size_continuous("Percentage",
                        breaks = seq(0, 0.8, 0.2), range = c(1,5))

# Extended Data Fig.2c----
# Data from HCC atlas
ggplot(HCC_atlas_fibroblast@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster == "Fibroblast") %>% arrange(desc(Sub_Cluster)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 4, alpha = .8) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Extended Data Fig.2d-f----
library(monocle3)
library(SeuratWrappers)
library(magrittr)
CAF_HCC_UMAP.cds <- as.cell_data_set(CAF_HCC)
CAF_HCC_UMAP.cds <- cluster_cells(CAF_HCC_UMAP.cds, k = 30, partition_qval = 0.05, resolution = 0.01, reduction_method = "UMAP")
colData(CAF_HCC_UMAP.cds)$Cluster <- colData(CAF_HCC_UMAP.cds)$Sub_Cluster
CAF_HCC_UMAP.cds@clusters$UMAP$clusters <- as.character(CAF_HCC_UMAP.cds@colData$Cluster)
names(CAF_HCC_UMAP.cds@clusters$UMAP$clusters) <- CAF_HCC_UMAP.cds@colData$index
CAF_HCC_UMAP.cds <- learn_graph(CAF_HCC_UMAP.cds, learn_graph_control = list(minimal_branch_len = 3))
CAF_HCC_UMAP.cds <- order_cells(CAF_HCC_UMAP.cds)
p1 <- plot_cells(
  CAF_HCC_UMAP.cds, 
  color_cells_by = "Sub_Cluster", 
  cell_size = 1,
  trajectory_graph_segment_size = 1,
  label_groups_by_cluster = FALSE, 
  label_cell_groups = FALSE,
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  labs(x = "UMAP1", y = "UMAP2")

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
p2 <- plot_cells(
  CAF_HCC_UMAP.cds, 
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

CAF_HCC <- AddMetaData(CAF_HCC, metadata = CAF_HCC_UMAP.cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "Monocle3_pseudotime")
plot.data <- CAF_HCC@meta.data %>% filter(Sub_Cluster != "APOA2+ CAF") %>% filter(!(Sub_Cluster == "CAF" & Monocle3_pseudotime > 8)) %>% mutate(Grp = plyr::revalue(Sub_Cluster, replace = c(
  "ABCAB+ Fib" = "Grp1", "MT1M+ Fib" = "Grp1",
  "Fib" = "Grp1", "CAF" = "Grp2", 
  "HSP+ CAF" = "Grp4", "MYH11+ CAF" = "Grp4",
  "SDC2+ CAF" = "Grp5", "POSTN+ CAF" = "Grp5"
))) %>% 
  mutate(Grp = factor(Grp, levels = c("Grp1","Grp2","Grp3","Grp4","Grp5"))) %>%
  mutate(Sub_Cluster = factor(as.character(Sub_Cluster),levels = c("MT1M+ Fib","Fib","ABCAB+ Fib","CAF","MYH11+ CAF","HSP+ CAF","POSTN+ CAF","SDC2+ CAF"))) %>%
  select(Monocle3_pseudotime, Sub_Cluster, Grp)
ggplot(plot.data, aes(x = Monocle3_pseudotime, y = Sub_Cluster)) +
  geom_boxplot(aes(color = Sub_Cluster)) +
  geom_point(aes(color = Sub_Cluster), position = position_jitter(h=0.2), size = .8) +
  facet_grid(Grp~., scales = "free", space = "free")  +
  theme_bw(base_size = 10) +
  ylab("") + xlab("Monocle3 pseudotime") +
  theme(legend.position = "none", 
        strip.text = element_blank(),
        strip.background = element_blank())

# Extended Data Fig.2g----
# SCENIC(R) results not provided
library(scFunctions)
binary_regulon <- loadInt(scenicOptions, "aucell_binary_nonDupl")
binary_regulon <- binary_regulon[onlyNonDuplicatedExtended(rownames(binary_regulon)),]
rrs_df <- calculate_rrs(CAF_HCC@meta.data, binary_regulon, "Sub_Cluster")
rrs_df <- rrs_df %>% subset(!grepl("extended", regulon))
rrs_df_sub1 <- rrs_df %>% subset(cell_type == "POSTN+ CAF") %>% 
  arrange(desc(RSS))
rrs_df_sub1 <- rrs_df_sub1 %>% mutate(rank = as.numeric(rownames(rrs_df_sub1)))
rrs_ranking_plot1 <- 
  ggplot(rrs_df_sub1, aes(rank, RSS, label = regulon)) + 
  geom_point(data = subset(rrs_df_sub1,rank > 12),color = "#619CFF", size = 1) + 
  geom_point(data = subset(rrs_df_sub1,rank <= 12), color = "#F87269", size = .8) + 
  labs(x = "Rank", y = "Regulon Specificity Score", title = "Top regulons of\nPOSTN+ CAF") +
  geom_text_repel(data = subset(rrs_df_sub1, rank <= 12),
                  nudge_x = 35 - subset(rrs_df_sub1, rank <= 12)$rank,
                  segment.size = 0.2,
                  segment.color = "#F87269", color = "#F87269",
                  direction = "y", hjust = 0) +
  theme_cowplot(font_size = 12, line_size = 1)
rrs_df_sub2 <- rrs_df %>% subset(cell_type == "MYH11+ CAF") %>% 
  arrange(desc(RSS))
rrs_df_sub2 <- rrs_df_sub2 %>% mutate(rank = as.numeric(rownames(rrs_df_sub2)))
rrs_ranking_plot2 <- 
  ggplot(rrs_df_sub2, aes(rank, RSS, label = regulon)) + 
  geom_point(data = subset(rrs_df_sub2,rank > 8),color = "#619CFF", size = 1) + 
  geom_point(data = subset(rrs_df_sub2,rank <= 8), color = "#F87269", size = .8) + 
  labs(x = "Rank", y = "Regulon Specificity Score", title = "Top regulons of\nMYH11+ CAF") +
  geom_text_repel(data = subset(rrs_df_sub, rank <= 8),
                  nudge_x = 20 - subset(rrs_df_sub, rank <= 8)$rank,
                  segment.size = 0.2,
                  segment.color = "#F87269", color = "#F87269",
                  direction = "y", hjust = 0) +
  theme_cowplot(font_size = 12, line_size = 1)

# Extended Data Fig.2h----
cells_used <- CAF_HCC@meta.data %>% filter(Sub_Cluster %in% c("POSTN+ CAF", "MYH11+ CAF")) %>% pull(CellID)
degenes <- LIMMA(
  CAF_HCC@assays$RNA@data[,cells_used],
  stringr::str_split_fixed(CAF_HCC@meta.data[cells_used,"Sub_Cluster"],"\\+",2)[,1]
)
degenes$Sig <- FALSE
degenes[degenes$adj.P.Val < 0.05 & abs(degenes$logFC) > 0.5,"Sig"] <- TRUE
genes_labeled <- c(
  degenes %>% filter(Sig == TRUE) %>% slice_max(n = 20, order_by = logFC) %>% pull(Symbol),
  degenes %>% filter(Sig == TRUE) %>% slice_min(n = 20, order_by = logFC) %>% pull(Symbol)
)
genes_labeled <- c(genes_labeled,"POSTN","ENG")
degenes[genes_labeled,"label"] <- genes_labeled
volcano_plot <- 
  ggplot(degenes, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = Binomial_color_panel) +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  theme_cowplot(font_size = 16) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1)) +
  geom_vline(xintercept = c(-0.5,0.5), color = "grey", linetype = "dashed", lwd = 1) +
  geom_text_repel(data = degenes, aes(label = label), size = 4)

# Extended Data Fig.3a----
Fib_fetal <- subset(Fib_fetal, subset = Sub_Cluster != "Mesenchymal")
Fib_fetal_markers <- FindAllMarkers(Fib_fetal)
genes_used <- Fib_fetal_markers %>% group_by(cluster) %>% filter(pct.2 < 0.8, pct.1 > 0.5) %>% slice_max(n = 25, order_by = avg_log2FC) %>% pull(gene) %>% unique()
genes_used <- genes_used[-grep("^MT",genes_used)]
genes_used <- c(genes_used, "ACTA2", "PLEK")
markers_mean_exp <- 
  aggregate(
    t(as.matrix(Fib_fetal@assays$RNA@data[genes_used,])),
    list(Cluster = Fib_fetal$Sub_Cluster),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Cluster
markers_mean_exp$Cluster <- c()
markers_plot_matrix <- apply(markers_mean_exp, 2, scale) %>% t()
colnames(markers_plot_matrix) <- row.names(markers_mean_exp)
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

# Extended Data Fig.3b-c----
library(monocle3)
library(SeuratWrappers)
library(magrittr)
Fib_fetal_UMAP.cds <- as.cell_data_set(Fib_fetal)
Fib_fetal_UMAP.cds <- cluster_cells(Fib_fetal_UMAP.cds, k = 30, partition_qval = 0.05, resolution = 0.01, reduction_method = "UMAP")
colData(Fib_fetal_UMAP.cds)$Cluster <- colData(Fib_fetal_UMAP.cds)$Sub_Cluster
Fib_fetal_UMAP.cds@clusters$UMAP$clusters <- as.character(Fib_fetal_UMAP.cds@colData$Cluster)
names(Fib_fetal_UMAP.cds@clusters$UMAP$clusters) <- Fib_fetal_UMAP.cds@colData$index
Fib_fetal_UMAP.cds <- learn_graph(Fib_fetal_UMAP.cds, learn_graph_control = list(minimal_branch_len = 3))
Fib_fetal_UMAP.cds <- order_cells(Fib_fetal_UMAP.cds)
plot_cells(
  Fib_fetal_UMAP.cds, 
  color_cells_by = "Sub_Cluster", 
  cell_size = 1,
  trajectory_graph_segment_size = 1,
  label_groups_by_cluster = FALSE, 
  label_cell_groups = FALSE,
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  labs(x = "UMAP1", y = "UMAP2")

plot_cells(
  Fib_fetal_UMAP.cds, 
  color_cells_by = "pseudotime",
  cell_size = 1,
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  scale_color_gradientn(colours = myPalette(100)) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1, "cm"))

Fib_fetal$Monocle3_pseudotime <- Fib_fetal_UMAP.cds@principal_graph_aux@listData$UMAP$pseudotime
plot.data <- Fib_fetal@meta.data %>% filter(Sub_Cluster != "Mesenchymal") %>% mutate(Grp = plyr::revalue(Sub_Cluster, replace = c(
  "Fib" = "Grp1","ACTA2+ Fib" = "Grp3","PLEK+ Fib" = "Grp2","BGN+ Fib" = "Grp1","APOC3+ Fib" = "Grp4","ITM2C+ Fib" = "Grp4","COL1A1+ Fib" = "Grp3"
))) %>% 
  mutate(Grp = factor(Grp, levels = c("Grp1","Grp2","Grp3","Grp4"))) %>%
  mutate(Sub_Cluster = factor(as.character(Sub_Cluster),levels = rev(c("Fib","BGN+ Fib","ACTA2+ Fib","PLEK+ Fib","COL1A1+ Fib","APOC3+ Fib","ITM2C+ Fib"))))
ggplot(plot.data, aes(x = Monocle3_pseudotime, y = Sub_Cluster)) +
  geom_boxplot(aes(color = Sub_Cluster)) +
  geom_point(aes(color = Sub_Cluster), position = position_jitter(h=0.2), size = .8) +
  facet_grid(Grp~., scales = "free", space = "free")  +
  scale_color_manual(name = "", values = Fib_fetal_color_panel) +
  theme_bw(base_size = 10) +
  ylab("") + xlab("Monocle3 pseudotime") +
  theme(legend.position = "none", 
        strip.text = element_blank(),
        strip.background = element_blank())

# Extended Data Fig.3d----
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:6])
plot.data <- ROIE(table(Fib_fetal@meta.data[,c("Sub_Cluster","PatientID")])) %>% melt() %>% filter(Var1 != "Mesenchymal") %>% mutate(value = pmin(value, 3)) %>% mutate(Var1 = factor(as.character(Var1), levels = c("Fib","BGN+ Fib","PLEK+ Fib","ACTA2+ Fib","COL1A1+ Fib","APOC3+ Fib","ITM2C+ Fib")))
ggplot(plot.data, aes(Var2, forcats::fct_rev(Var1), fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = myPalette(100)) +
  labs(x = "", y = "") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.line = element_blank(), 
        axis.ticks = element_blank()) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 12))

# Extended Data Fig.3e----
ggplot(Fib_integrated@meta.data, aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 2, alpha = .8) +
  theme_cowplot(font_size = 7) +
  labs(x = "UMAP1", y = "UMAP2") +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

# Extended Data Fig.3f----
plot.data <- Fib_integrated@meta.data
plot.data$Sub_Cluster <- as.character(plot.data$Sub_Cluster)
plot.data$Sub_Cluster[!(plot.data$Sub_Cluster %in% c("POSTN+ CAF","ITM2C+ Fib","PLEK+ Fib","Fib(Tumor)","Fib(Fetal)"))] <- "A0"
ggplot(plot.data %>% arrange(Sub_Cluster), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 2, alpha = .6) +
  theme_cowplot(font_size = 7) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_color_manual(name = "", values = c("A0" = "lightgrey", "Fib(Tumor)" = "#68AB9F", "Fib(Fetal)" = "#F69459", "PLEK+ Fib" = "#F58584", "ITM2C+ Fib" = "#815AA8", "POSTN+ CAF" = "#693D99")) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

library(monocle3)
library(SeuratWrappers)
library(magrittr)
Fib_integrated_UMAP.cds <- as.cell_data_set(Fib_integrated, assay = "RNA")
Fib_integrated_UMAP.cds <- cluster_cells(Fib_integrated_UMAP.cds)
Fib_integrated_UMAP.cds <- learn_graph(Fib_integrated_UMAP.cds, learn_graph_control = list(minimal_branch_len = 5))
root_cells <- data.frame(colData(Fib_integrated_UMAP.cds)) %>% filter(Sub_Cluster == "HSC") %>% row.names()
Fib_integrated_UMAP.cds <- order_cells(Fib_integrated_UMAP.cds, root_cells = root_cells)
plot_cells(
  Fib_integrated_UMAP.cds, 
  color_cells_by = "Sub_Cluster", 
  cell_size = 0,
  trajectory_graph_segment_size = 1,
  label_groups_by_cluster = FALSE, 
  label_cell_groups = FALSE,
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  labs(x = "UMAP1", y = "UMAP2")

# Extended Data Fig.3g----
ggplot(Fib_integrated@meta.data, aes(x = Integrated_cluster, fill = Source)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 6), ncol = 1)) +
  scale_x_discrete(labels=paste0(
    names(table(Fib_integrated$Integrated_cluster))," (",
    table(Fib_integrated$Integrated_cluster),")")) +
  coord_flip()

# Extended Data Fig.3h----
library(scibetR)
CAF_HCC_cells_used <- CAF_HCC@meta.data %>% pull(CellID)
Fib_fetal_cells_used <- Fib_fetal@meta.data %>% pull(CellID)
scibet_test_set <- data.frame(
  t(as.matrix(CAF_HCC@assays$RNA@counts[,CAF_HCC_cells_used])), 
  label = CAF_HCC@meta.data[CAF_HCC_cells_used,"Sub_Cluster"]
)
scibet_train_set <- data.frame(
  t(as.matrix(Fib_fetal@assays$RNA@counts[,Fib_fetal_cells_used])), 
  label = Fib_fetal@meta.data[Fib_fetal_cells_used,"Sub_Cluster"]
)
scibet_prd <- SciBet_R(scibet_train_set, scibet_test_set[,-ncol(scibet_test_set)])
confusion_matrix <- table(as.character(scibet_test_set$label),scibet_prd)
confusion_matrix <- confusion_matrix / rowSums(confusion_matrix)
confusion_matrix <- as.data.frame(confusion_matrix)
confusion_matrix <- dcast(confusion_matrix, Var1~scibet_prd, value.var = "Freq")
row.names(confusion_matrix) <- confusion_matrix$Var1
confusion_matrix$Var1 <- c()
confusion_matrix <- as.matrix(confusion_matrix)
test_order <- c("CAF","APOA2+ CAF","HSP+ CAF","MYH11+ CAF","SDC2+ CAF","POSTN+ CAF","Fib","MT1M+ Fib","ABCAB+ Fib") 
train_order <- c("Fib","ACTA2+ Fib","COL1A1+ Fib","BGN+ Fib","APOC3+ Fib","ITM2C+ Fib","PLEK+ Fib","Mesenchymal")
color_used <- 
  circlize::colorRamp2(c(seq(0,0.5,length.out = 50),0.7), 
                       c(colorRampPalette(rev(brewer.pal(6,'Blues')))(100)[51:100],"red2"), space = "RGB")
p <- Heatmap(
  confusion_matrix[test_order, train_order],
  name = 'Cell Proportion',
  column_title = "Fetal liver (train data)", row_title = "HCC (test data)",
  col = color_used,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_side = "left", row_dend_side = "right",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_title_gp = gpar(fontsize = 16), column_title_gp = gpar(fontsize = 16),
  row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 16),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 16), 
                              labels_gp = gpar(fontsize = 16),
                              legend_height = unit(3, "cm"))
)

# Extended Data Fig.3i----
library(dendextend)
library(scibetR)
library(plyr)
cluster <- c(
  paste0("Tumor_",as.character(CAF_HCC@meta.data[CAF_HCC_cells_used,"Sub_Cluster"])), 
  paste0("Fetal_",as.character(Fib_fetal@meta.data[Fib_fetal_cells_used,"Sub_Cluster"]))
)
CAF_HCC <- ScaleData(CAF_HCC, features = row.names(CAF_HCC))
Fib_fetal <- ScaleData(Fib_fetal, features = row.names(Fib_fetal))
train_genes_used <- SelectGene_R(scibet_train_set, k = 1000)
test_genes_used <- SelectGene_R(scibet_test_set, k = 1000)
genes_used <- intersect(train_genes_used, test_genes_used)
genes_used <- genes_used[!grepl("[.]", genes_used)]
expression_used <- cbind(CAF_HCC@assays$RNA@scale.data[genes_used,CAF_HCC_cells_used], Fib_fetal@assays$RNA@scale.data[genes_used,Fib_fetal_cells_used])
avg_expression <- aggregate(t(expression_used), list(cluster), mean)
row.names(avg_expression) <- avg_expression$Group.1
avg_expression$Group.1 <- c()
M <- (1- cor(t(avg_expression),method="pearson"))/2
par(mfcol=c(1,1))
hc <- hclust(as.dist(M),method="complete")
plot(hc)

Fib_fetal@meta.data$tempCluster <- plyr::revalue(
  Fib_fetal@meta.data$Sub_Cluster,
  replace = c("Fib" = "C3",
              "ACTA2+ Fib" = "C1",
              "COL1A1+ Fib" = "C4",
              "BGN+ Fib" = "C4",
              "APOC3+ Fib" = "C2",
              "ITM2C+ Fib" = "C2",
              "Mesenchymal" = "C1",
              "PLEK+ Fib" = "C1"))
Fib_fetal <- SetIdent(Fib_fetal, value = Fib_fetal@meta.data$tempCluster)
Fib_fetal.markers <- FindAllMarkers(Fib_fetal)
CAF_HCC@meta.data$tempCluster <- plyr::revalue(
  CAF_HCC@meta.data$Sub_Cluster,
  replace = c("HSP+ CAF" = "C1",
              "MYH11+ CAF" = "C1",
              "CAF" = "C1",
              "POSTN+ CAF" = "C2",
              "SDC2+ CAF" = "C2",
              "APOA2+ CAF" = "C3",
              "Fib" = "C4",
              "MT1M+ Fib" = "C4",
              "ABCAB+ Fib" = "C4"))
CAF_HCC <- SetIdent(CAF_HCC, value = CAF_HCC@meta.data$tempCluster)
CAF_HCC.markers <- FindAllMarkers(CAF_HCC)
genes_used <- c(
  intersect(Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C1", avg_log2FC > 0) %>% top_n(n = 250, wt = avg_log2FC) %>% pull(gene),
            CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C1", avg_log2FC > 0) %>% top_n(n = 250, wt = avg_log2FC) %>% pull(gene)),
  intersect(Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C2", avg_log2FC > 0) %>% top_n(n = 250, wt = avg_log2FC) %>% pull(gene),
            CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C2", avg_log2FC > 0) %>% top_n(n = 250, wt = avg_log2FC) %>% pull(gene)),
  intersect(Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C3", avg_log2FC > 0) %>% top_n(n = 250, wt = avg_log2FC) %>% pull(gene),
            CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C3", avg_log2FC > 0) %>% top_n(n = 250, wt = avg_log2FC) %>% pull(gene)),
  intersect(Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C4", avg_log2FC > 0) %>% top_n(n = 250, wt = avg_log2FC) %>% pull(gene),
            CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C4", avg_log2FC > 0) %>% top_n(n = 250, wt = avg_log2FC) %>% pull(gene)),
  intersect(Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C1", avg_log2FC < 0) %>% top_n(n = 250, wt = -avg_log2FC) %>% pull(gene),
            CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C1", avg_log2FC < 0) %>% top_n(n = 250, wt = -avg_log2FC) %>% pull(gene)),
  intersect(Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C2", avg_log2FC < 0) %>% top_n(n = 250, wt = -avg_log2FC) %>% pull(gene),
            CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C2", avg_log2FC < 0) %>% top_n(n = 250, wt = -avg_log2FC) %>% pull(gene)),
  intersect(Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C3", avg_log2FC < 0) %>% top_n(n = 250, wt = -avg_log2FC) %>% pull(gene),
            CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C3", avg_log2FC < 0) %>% top_n(n = 250, wt = -avg_log2FC) %>% pull(gene)),
  intersect(Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C4", avg_log2FC < 0) %>% top_n(n = 200, wt = -avg_log2FC) %>% pull(gene),
            CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene), cluster == "C4", avg_log2FC < 0) %>% top_n(n = 250, wt = -avg_log2FC) %>% pull(gene))
)
genes_used <- unique(genes_used)
CAF_HCC_exp <- aggregate(t(CAF_HCC@assays$RNA@scale.data[genes_used,]), list(Cluster = CAF_HCC$Sub_Cluster), mean)
row.names(CAF_HCC_exp) <- CAF_HCC_exp$Cluster
CAF_HCC_exp$Cluster <- c()
CAF_HCC_exp <- t(CAF_HCC_exp)
colnames(CAF_HCC_exp) <- paste0("Tumor_",colnames(CAF_HCC_exp))
Fib_fetal_exp <- aggregate(t(Fib_fetal@assays$RNA@scale.data[genes_used,]), list(Cluster = Fib_fetal$Sub_Cluster), mean)
row.names(Fib_fetal_exp) <- Fib_fetal_exp$Cluster
Fib_fetal_exp$Cluster <- c()
Fib_fetal_exp <- t(Fib_fetal_exp)
colnames(Fib_fetal_exp) <- paste0("Fetal_",colnames(Fib_fetal_exp))
CAF_all_exp_agg <- cbind(Fib_fetal_exp, CAF_HCC_exp)
CAF_all_exp_agg <- apply(CAF_all_exp_agg, 2, scale)
CAF_all_exp_agg_quantile <- quantile(CAF_all_exp_agg,c(0.01,0.98))
CAF_all_exp_agg <- pmax(CAF_all_exp_agg, CAF_all_exp_agg_quantile[1])
CAF_all_exp_agg <- pmin(CAF_all_exp_agg, CAF_all_exp_agg_quantile[2])
CAF_all_exp_agg <- CAF_all_exp_agg[,hc$labels[hc$order]]
row.names(CAF_all_exp_agg) <- genes_used
max.avg <- apply(CAF_all_exp_agg, 1, which.max)
gene_order <- c()
for(i in ncol(CAF_all_exp_agg):1){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(gene = names(sort(CAF_all_exp_agg[names(max.avg)[max.avg == i],i], decreasing = F)), cluster = colnames(CAF_all_exp_agg)[i], stringsAsFactors = F)
  }
  if(sum(max.avg == i) == 1){
    temp <- data.frame(gene = names(max.avg)[max.avg == i], cluster = colnames(CAF_all_exp_agg)[i], stringsAsFactors = F)
  }
  gene_order <- rbind(gene_order, temp)
}
CAF_all_exp_agg <- CAF_all_exp_agg[gene_order$gene,]
CAF_all_exp_agg <- CAF_all_exp_agg[,hc$labels]
dend <- hc
dend <- dendextend::color_branches(dend, k = 4)
color_used <- circlize::colorRamp2(seq(min(CAF_all_exp_agg, na.rm = TRUE), max(CAF_all_exp_agg, na.rm = TRUE), length = 11), rev(RColorBrewer::brewer.pal(11,"RdBu")))
genes_labeled <- c("COL1A2","COL1A1","ATF4","FOS","NFKNIZ","CD63","ITM2B","IER2","LMNA","NR4A1","CCNL1","JUND","ATF3","MYL6","ACTB","B2M","ALB","MYL12A","COLEC11","CCL21","CTSC","ITM2C","TCF4","COL6A1","COL5A1","CXCL12","APOC3","TGFBI","ITGA1","MYH10","PLSCR4","MARCKS","PLD3","PRNP","LAMB1","ASPN","POSTN","COL6A3","CD81","COL14A1","IGFBP4","IFITM3","COL4A1","APOA2","TSPAN3","CHMP2A","CD151","HSP90AA1","HSPA1A","VDAC2","PPP1CA","HMGB1","PGK1","PLP2","BCAM","IGFBP5","CD9")
p <- Heatmap(CAF_all_exp_agg,
             cluster_rows = F,col = color_used,
             column_dend_side = "bottom",
             row_names_gp = gpar(fontsize = 6),
             cluster_columns = dend,
             heatmap_legend_param = list(
               at = c(-1.5, -0.75, 0, 0.75, 1.5, 2.25),
               labels = c("<-1.5", "-0.75", "0", "0.75", "1.5", ">2.25"),
               title = "Normalized expression",
               legend_height = unit(4, "cm"),
               title_position = "lefttop-rot")) +
  rowAnnotation(link = anno_mark(
    at = which(gene_order$gene %in% genes_labeled),
    labels = gene_order$gene[which(gene_order$gene %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

# Extended Data Fig.3j-k----
library(richR)
hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
high_up_gene_name <- plot.df %>% filter(p_fetal < 0.05 & p_tumor < 0.05 & (logFC_tumor > 0.05 | logFC_fetal > 0.05)) %>% pull(Symbol)
high_up_go <- richGO(high_up_gene_name,godata = hsago,ontology ="BP")
high_up_go@result %>% filter(Padj < 0.05) %>% pull(Term)
# load("./Onco-fetal CAF Identification/Co_up_reg_pathway.rda")
# The results might change based on the annotation used
GO_to_show <- c(1,5,7,12,17,40,44,61,65,70,119,138,140,152,166,180,181)
pathway_plot.df <- data.frame(
  pathway = high_up_go@result$Term[GO_to_show],
  pvalue = high_up_go@result$Padj[GO_to_show],
  generatio = (high_up_go@result$Significant/high_up_go@result$Annotated)[GO_to_show],
  database = "GO (BP)"
)
pathway_plot.df$logP <- -log10(pathway_plot.df$pvalue)
pathway_plot.df$logP <- pmin(pathway_plot.df$logP,10)
ggplot(pathway_plot.df, aes(x = reorder(pathway,generatio), y = generatio)) +
  geom_bar(aes(fill = logP), stat = "identity") +
  scale_fill_gradientn(name = "-log10 (P-value)", limits = c(0,10), colors = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)) +
  theme_bw(base_size = 12) +
  labs(y = "Gene Ratio", x = "") +
  coord_flip() +
  theme(strip.text.x = element_blank(),strip.background = element_rect(fill=NA, color=NA))

tumor_up_gene_name <- plot.df %>% filter(p_fetal > 0.05 & p_tumor < 0.05 & logFC_tumor > 0.5 & abs(logFC_fetal) < 0.5) %>% pull(Symbol)
tumor_up_go <- richGO(tumor_up_gene_name,godata = hsago,ontology ="BP")
tumor_up_go@result %>% filter(Padj < 0.05) %>% pull(Term)
# load("./Onco-fetal CAF Identification/Tumor_up_reg_pathway.rda")
# The results might change based on the annotation used
GO_to_show <- c(1,3,7,12,15,20,24,39,62,78,109,120,124,127,129,171)
pathway_plot.df <- data.frame(
  pathway = tumor_up_go@result$Term[GO_to_show],
  pvalue = tumor_up_go@result$Padj[GO_to_show],
  generatio = (tumor_up_go@result$Significant/tumor_up_go@result$Annotated)[GO_to_show],
  database = "GO (BP)"
)
pathway_plot.df$logP <- -log10(pathway_plot.df$pvalue)
pathway_plot.df$logP <- pmin(pathway_plot.df$logP,10)
ggplot(pathway_plot.df, aes(x = reorder(pathway,generatio), y = generatio)) +
  geom_bar(aes(fill = logP), stat = "identity") +
  scale_fill_gradientn(name = "-log10 (P-value)", limits = c(0,10), colors = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)) +
  theme_bw(base_size = 12) +
  labs(y = "Gene Ratio", x = "") +
  coord_flip() +
  theme(strip.text.x = element_blank(),strip.background = element_rect(fill=NA, color=NA))

# Supplementary Fig.4----
DefaultAssay(Fib_integrated) <- "RNA"
fib_integrated_markers <- list()
for(cluster in 0:11){
  fib_integrated_markers[[cluster+1]] <- FindConservedMarkers(Fib_integrated,ident.1 = cluster, grouping.var = "Source")
}
genes_used <- c()
for(i in 1:12){
  genes_used <- c(genes_used, row.names(fib_integrated_markers[[i]])[1:25])
}
genes_used <- unique(genes_used[-grep("^MT|[.]",genes_used)])
plot.list <- list()
for(tissue in c("Adj Normal","Fetal","Healthy","Tumor")){
  cells_used <- Fib_integrated@meta.data %>% filter(Source == tissue) %>% row.names()
  num.clusters <- table(Fib_integrated@meta.data[cells_used,"Integrated_cluster"])
  cluster.used <- names(num.clusters[num.clusters > 10])
  cells_used <- Fib_integrated@meta.data %>% filter(Source == tissue, Integrated_cluster %in% cluster.used) %>% row.names()
  markers_mean_exp <- 
    aggregate(
      t(as.matrix(Fib_integrated@assays$RNA@data[genes_used,cells_used])),
      list(Cluster = Fib_integrated@meta.data[cells_used,"Integrated_cluster"]),
      mean)
  row.names(markers_mean_exp) <- markers_mean_exp$Cluster
  markers_mean_exp$Cluster <- c()
  markers_plot_matrix <- apply(markers_mean_exp, 2, scale) %>% t()
  markers_plot_matrix <- markers_plot_matrix[!is.na(rowSums(markers_plot_matrix)),]
  colnames(markers_plot_matrix) <- row.names(markers_mean_exp)
  gene_order <- c()
  max.avg <- apply(markers_plot_matrix, 1, which.max)
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
  markers_plot_matrix_quantile <- quantile(markers_plot_matrix, c(0.01, 0.99))
  markers_plot_matrix <- pmax(markers_plot_matrix, markers_plot_matrix_quantile[1])
  markers_plot_matrix <- pmin(markers_plot_matrix, markers_plot_matrix_quantile[2])
  color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), max(markers_plot_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
  genes_labeled <- c("BCAM","PDGFRB","ALB","COLEC11","ACTB","APOC3","B2M","COL3A1","EGR1","CTSD","FOSB","IGFBP7","DCN","TIMP1","TPM2","MYL9","DSTN","SPARC","MYH11","IFITM2","IRF1","RGS5","HSPA1A","JUN","HSP90AB1","TAGLN","COL6A3","COL5A1","IGFBP3","COL1A2","PPIC","FOS","COL1A1","ITM2C","CTSC","NR1H4","CD9","VIM","TCF4","LUM","ITM2A","APOA1","APOA2","COL6A2","DNAJA1","IGFBP5","RPS3","RPL13","S100A4","S100A6","ENO1","TUBB4B","IDH2","S100A8","CXCL8","TFF2","CCL2","FCN3","CLDN5","KDR","IL33","FLT1","CD93","CTSL","TGFBR2","VAMP5")
  plot.list[[tissue]] <- Heatmap(markers_plot_matrix[gene_order$gene,],
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
}
