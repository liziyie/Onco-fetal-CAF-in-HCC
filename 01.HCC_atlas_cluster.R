setwd("/work/lzy/project/onco_fetal2/")
source("../utils/utils_plot.R")
source("../utils/utils_data_processing.R")
source("../utils/utils_color.R")
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# >Load dataset----
HCC_atlas <- readRDS("./02.processed_data/RPCA_integration/HCC_atlas.combined.rds")

DimPlot(HCC_atlas, label = T, group.by = "Anno_R1")
DimPlot(HCC_atlas, group.by = "Global_Cluster", label = T, split.by = "orig.ident", ncol = 3)

HCC_atlas@meta.data$Anno_R1 <- plyr::revalue(
  HCC_atlas@meta.data$integrated_snn_res.0.2,
  replace = c("0" = "T/NK",
              "1" = "Mononuclear",
              "2" = "Hepatocyte",
              "3" = "Hepatocyte",
              "4" = "Endothelium",
              "5" = "B cells",
              "6" = "Hepatocyte",
              "7" = "Fibroblast",
              "8" = "Cycling cells",
              "9" = "Hepatocyte",
              "10" = "Hepatocyte",
              "11" = "Hepatocyte")
)
DefaultAssay(HCC_atlas) <- "RNA"

HCC_atlas_hepatocyte <- HCC_atlas %>% subset(Anno_R1 == "Hepatocyte")

# **Round 2----
# >Fibroblast clustering----
# HCC_atlas_fibroblast <- HCC_atlas %>% subset(Anno_R1 == "Fibroblast")
# saveRDS(HCC_atlas_fibroblast, file = "./02.processed_data/RPCA_integration/HCC_atlas_fibroblast.rds")
HCC_atlas_fibroblast <- readRDS("./02.processed_data/RPCA_integration/HCC_atlas_fibroblast.rds")

HCC_atlas_fibroblast.list <- SplitObject(HCC_atlas_fibroblast, split.by = "orig.ident")
HCC_atlas_fibroblast.list  <- lapply(X = HCC_atlas_fibroblast.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_fibroblast.list)
HCC_atlas_fibroblast.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_fibroblast.list, anchor.features = features, k.filter = 50)
HCC_atlas_fibroblast <- IntegrateData(anchorset = HCC_atlas_fibroblast.anchors)
DefaultAssay(HCC_atlas_fibroblast) <- "integrated"
HCC_atlas_fibroblast <- ScaleData(HCC_atlas_fibroblast, verbose = FALSE)
HCC_atlas_fibroblast <- RunPCA(HCC_atlas_fibroblast, npcs = 30, verbose = FALSE)
HCC_atlas_fibroblast <- RunUMAP(HCC_atlas_fibroblast, reduction = "pca", dims = 1:10, min.dist = .1, spread = 1)
HCC_atlas_fibroblast <- FindNeighbors(HCC_atlas_fibroblast, reduction = "pca", dims = 1:10)
HCC_atlas_fibroblast <- FindClusters(HCC_atlas_fibroblast, resolution = 0.4)
saveRDS(HCC_atlas_fibroblast, file = "./02.processed_data/CCA_integration/HCC_atlas_fibroblast.rds")

HCC_atlas_fibroblast <- HCC_atlas_fibroblast %>% subset(seurat_clusters %in% c(0,1,2,3,4,5,6,7,9)) %>% subset(orig.ident %ni% c("HCC_Fan","HCC_Zhang"))
HCC_atlas_fibroblast.list <- SplitObject(HCC_atlas_fibroblast, split.by = "orig.ident")
HCC_atlas_fibroblast.list  <- lapply(X = HCC_atlas_fibroblast.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_fibroblast.list)
HCC_atlas_fibroblast.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_fibroblast.list, anchor.features = features)
HCC_atlas_fibroblast <- IntegrateData(anchorset = HCC_atlas_fibroblast.anchors)
DefaultAssay(HCC_atlas_fibroblast) <- "integrated"
HCC_atlas_fibroblast <- ScaleData(HCC_atlas_fibroblast, verbose = FALSE)
HCC_atlas_fibroblast <- RunPCA(HCC_atlas_fibroblast, npcs = 30, verbose = FALSE)
HCC_atlas_fibroblast <- RunUMAP(HCC_atlas_fibroblast, reduction = "pca", dims = 1:10)
HCC_atlas_fibroblast <- FindNeighbors(HCC_atlas_fibroblast, reduction = "pca", dims = 1:10, force.recalc = TRUE)
HCC_atlas_fibroblast <- FindClusters(HCC_atlas_fibroblast, resolution = 0.4)
DefaultAssay(HCC_atlas_fibroblast) <- "RNA"
fibroblast_markers <- FindAllMarkers(HCC_atlas_fibroblast)

HCC_atlas_fibroblast <- HCC_atlas_fibroblast %>% subset(seurat_clusters %ni% c(8,9))
HCC_atlas_fibroblast <- RunUMAP(HCC_atlas_fibroblast, reduction = "pca", dims = 1:10)
HCC_atlas_fibroblast@meta.data[,c("UMAP_1","UMAP_2")] <- HCC_atlas_fibroblast@reductions$umap@cell.embeddings[,c(1,2)]
HCC_atlas_fibroblast@meta.data$integrated_snn_res.0.2 <- HCC_atlas_fibroblast@meta.data$integrated_snn_res.0.4 <- HCC_atlas_fibroblast@meta.data$Anno_R1 <- c()

fibroblast_markers <- FindAllMarkers(HCC_atlas_fibroblast)
HCC_atlas_fibroblast@meta.data$Release_Cluster <- plyr::mapvalues(
  HCC_atlas_fibroblast@meta.data$seurat_clusters,
  from = c(0:7,10),
  to = c("F1","F7","F2","F3","F4","F5","F8","F9","F6")
)
HCC_atlas_fibroblast@meta.data$Release_Cluster <- factor(as.character(HCC_atlas_fibroblast@meta.data$Release_Cluster), levels = paste0("F",1:9))
HCC_atlas_fibroblast@meta.data[HCC_atlas_fibroblast@meta.data$orig.ident == 'HCC_Filliol' & HCC_atlas_fibroblast@meta.data$Tissue == "Adj Normal","Tissue"] <- "Cirrhosis"

CAF_color_panel <- c(
  "F1" = "#e8dae3", "F2" = "#de5178", "F3" = "#e4a9d1",
  "F4" = "#9c5cad", "F5" = "#72477d", "F6" = "#a9d8c7",
  "F7" = "#6ebb6e", "F8" = "#daf2ae", "F9" = "#33583a")
Tissue_color_panel <- c(
  "Adj Normal" = "#3f9b71", "Tumor" = "#d02c26", "Lymphnode" = "#a600d1", "Cirrhosis" = "#ec8c00"
)
Dataset_color_panel <- c(
  "HCC_Ankur" = "#caacd7", "HCC_Lu" = "#efc96d", "HCC_Ma" = "#ef9d9b", "HCC_Filliol" = "#4b69b4"
)
CAF_HCC_color_panel <- c(
  "CAF" = "#62A3C7", "Fib" = "#68AB9F", "SDC2+ CAF" = "#CAB2D6",
  "MYH11+ CAF" = "#F69459", "HSP+ CAF" = "#1C76B3", "POSTN+ CAF" = "#693D99", 
  "APOA2+ CAF" = "#FB9A99", "MT1M+ Fib" = "#B2DE89", "ABCAB+ Fib" = "#4DAC47"
)

p <- ggplot(HCC_atlas_fibroblast@meta.data %>% arrange(desc(Release_Cluster)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = 2.5, alpha = .8) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = CAF_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_CAF_cluster_UMAP.pdf", width = 3, height = 2.5)

p <- ggplot(HCC_atlas_fibroblast@meta.data %>% arrange(Tissue), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Tissue), size = 2, alpha = 1) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Tissue_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_CAF_tissue_UMAP.pdf", width = 3, height = 2.5)

human_gene_plot.df <- data.frame(
  PI16 = HCC_atlas_fibroblast@assays$RNA@data["PI16",],
  CD34 = HCC_atlas_fibroblast@assays$RNA@data["CD34",],
  DPT = HCC_atlas_fibroblast@assays$RNA@data["DPT",],
  Group = HCC_atlas_fibroblast@meta.data$Release_Cluster) %>% melt
p <- 
  ggplot(human_gene_plot.df, aes(x = Group, y = value)) +
  geom_violin(aes(color = Group, fill = Group), scale = "width") +
  scale_fill_manual(values = CAF_color_panel) +
  scale_color_manual(values = CAF_color_panel) +
  labs(y = "log2(TPM)") +
  theme_cowplot(font_size = 7) +
  facet_grid(variable ~ .) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "null")
ggsave(p, file = "./04.figures/R2_01_CAF_UNI_marker_violin.pdf", width = 2.25, height = 1.5)

p <- ggplot(HCC_atlas_fibroblast@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster == "Fibroblast") %>% arrange(desc(Sub_Cluster)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 4, alpha = .8) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_CAF_Ankur_Sub_Cluster_UMAP.pdf", width = 3.25, height = 2.5)

p <- 
  ggplot(HCC_atlas_fibroblast@meta.data, aes(x = Release_Cluster, fill = orig.ident)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  scale_fill_manual(values = Dataset_color_panel) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))
ggsave(p, file = "./04.figures/R2_01_CAF_dataset_bar.pdf", width = 2.5, height = 4)

p <- 
  ggplot(HCC_atlas_fibroblast@meta.data, aes(x = Tissue, fill = Release_Cluster)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  scale_fill_manual(values = CAF_color_panel) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))
ggsave(p, file = "./04.figures/R2_01_CAF_Tissue_bar.pdf", width = 2.5, height = 4)

HCC_atlas_fibroblast <- SetIdent(HCC_atlas_fibroblast, value = HCC_atlas_fibroblast@meta.data$Release_Cluster)
fibroblast_markers <- FindAllMarkers(HCC_atlas_fibroblast)
write.csv(fibroblast_markers, file = "./03.results/DEGenes/HCC_atlas/HCC_atlas_fibroblast_markers.csv")

genes_used <- fibroblast_markers %>% group_by(cluster) %>% filter(pct.2 < 0.8, pct.1 > 0.5) %>% slice_max(n = 8, order_by = avg_logFC) %>% pull(gene) %>% unique()
genes_used <- c(genes_used, "ABCA8")
gene.mean.matrix <- 
  aggregate(t(as.matrix(HCC_atlas_fibroblast@assays$RNA@data[genes_used,])),
            list(Cluster = HCC_atlas_fibroblast$Release_Cluster), 
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
  if(sum(max.avg == i) > 1){
    temp <- data.frame(Gene = names(sort(gene.mean.zscore[names(max.avg)[max.avg == i],i], decreasing = T)), Gene.Group = colnames(gene.mean.zscore)[i], stringsAsFactors = F)
  }else if(sum(max.avg == i) == 1){
    temp <- data.frame(Gene = names(max.avg)[max.avg == i], Gene.Group = colnames(gene.mean.zscore)[i], stringsAsFactors = F)
  }
  gene_order <- rbind(gene_order, temp)
}
gene.per <- 
  aggregate(t(as.matrix(HCC_atlas_fibroblast@assays$RNA@data[genes_used,])),
            list(Group = HCC_atlas_fibroblast@meta.data$Release_Cluster),
            function(x){sum(x > 2) / length(x)}) %>% 
  melt(id.vars = "Group", variable.name = "Gene", value.name = "Per")
plot.data <- merge(merge(gene.mean.zscore.df, gene.per), gene_order)
plot.data$Gene <- as.character(plot.data$Gene)
plot.data$Gene.Group <- factor(plot.data$Gene.Group, levels = paste0("F",9:1))
plot.data$Group <- factor(plot.data$Group, levels = paste0("F",9:1))
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
    axis.text.y = element_text(color = "black", size = 7),
    axis.text.x = element_text(color = "black", angle = 90, 
                               hjust = 1, vjust = 0.5, size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
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
                        breaks = seq(0, 0.8, 0.2), range = c(1,4))
ggsave(p, file = "./04.figures/R2_01_CAF_markers_bubble_heatmap.pdf", width = 3, height = 10.5)

universe_entrezid <- bitr(row.names(HCC_atlas_fibroblast@assays$RNA@counts), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
universe_entrezid <- unique(universe_entrezid$ENTREZID)
gene_list <- list()
for(cluster_used in as.character(unique(fibroblast_markers$cluster))){
  genes_symbol <- fibroblast_markers %>% filter(cluster == cluster_used, p_val_adj < 0.05, avg_logFC > 0) %>% top_n(n = 100, wt = avg_logFC) %>% pull(gene)
  gene.df <- bitr(genes_symbol, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  gene_list[[cluster_used]] <- unique(gene.df$ENTREZID)
}
cluster_go <- compareCluster(gene_list, fun = "enrichGO",
                             OrgDb = org.Hs.eg.db, ont = "ALL", 
                             universe = universe_entrezid,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = TRUE, pool = TRUE)
description_used <- c(
  "DNA-binding transcription activator activity, RNA polymerase II-specific",
  "transcription factor complex",
  "response to temperature stimulus",
  "muscle system process",
  "muscle contraction",
  "muscle tissue development",
  "myofibril assembly",
  "respiratory chain complex",
  "ATP synthesis coupled electron transport",
  "oxidative phosphorylation",
  "focal adhesion",
  "triglyceride catabolic process",
  "triglyceride metabolic process",
  "neutral lipid catabolic process",
  "lipid localization",
  "lipid transport",
  "blood microparticle",
  "vesicle lumen",
  "platelet degranulation",
  "high-density lipoprotein particle",
  "collagen-containing extracellular matrix",
  "extracellular matrix organization",
  "integrin binding",
  "cell adhesion molecule binding",
  "complex of collagen trimers",
  "heparin binding",
  "cell-matrix adhesion",
  "response to transforming growth factor beta",
  "regulation of ERK1 and ERK2 cascade",
  "response to unfolded protein",
  "protein localization to endoplasmic reticulum",
  "protein folding",
  "ribosomal subunit",
  "MHC class II protein complex",
  "MHC protein complex",
  "complement activation, lectin pathway",
  "growth factor binding"
)
cluster_go_plot <- cluster_go@compareClusterResult %>% 
  filter(Description %in% description_used)
dim_x <- unique(as.character(cluster_go_plot$Cluster))
dim_y <- description_used
map_x <- setNames(seq_along(dim_x), dim_x)
map_y <- setNames(seq_along(dim_y), dim_y)
cluster_go_plot_mat <- sparseMatrix(
  i=map_x[as.character(cluster_go_plot$Cluster)], 
  j=map_y[as.character(cluster_go_plot$Description)], 
  x=-log10(cluster_go_plot$p.adjust), 
  dims=c(length(dim_x), length(dim_y)), 
  dimnames=list(dim_x, dim_y)
)
cluster_go_plot_mat <- as.matrix(cluster_go_plot_mat)
cluster_go_plot_mat <- cluster_go_plot_mat[paste0("F",1:9),]
cluster_go_plot_mat_quantile <- quantile(cluster_go_plot_mat, c(0.05, 0.92))
cluster_go_plot_mat <- pmax(cluster_go_plot_mat, cluster_go_plot_mat_quantile[1])
cluster_go_plot_mat <- pmin(cluster_go_plot_mat, cluster_go_plot_mat_quantile[2])
color_used <- circlize::colorRamp2(seq(min(cluster_go_plot_mat), max(cluster_go_plot_mat), length = 9), rev(RColorBrewer::brewer.pal(11,"RdYlBu")[2:10]))
p <- Heatmap(cluster_go_plot_mat %>% t(),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "-log10(adjusted P-value)")
pdf("./04.figures/R2_01_CAF_pathway_heatmap.pdf", width = 5.5, height = 4)
draw(p)
dev.off()

cells_used <- HCC_atlas_fibroblast@meta.data %>% filter(Release_Cluster %in% c("F1","F2","F3","F7","F8")) %>% pull(CellName)
degenes <- LIMMA(
  HCC_atlas_fibroblast@assays$RNA@data[,cells_used],
  c("VSMC","CAF")[as.factor(HCC_atlas_fibroblast@meta.data[cells_used,"Release_Cluster"] %ni% c("F1","F2","F3"))]
)
degenes$Sig <- FALSE
degenes[degenes$adj.P.Val < 0.05 & abs(degenes$logFC) > 0.5,"Sig"] <- TRUE
genes_labeled <- c(
  degenes %>% filter(Sig == TRUE) %>% slice_max(n = 20, order_by = logFC) %>% pull(Symbol),
  degenes %>% filter(Sig == TRUE) %>% slice_min(n = 20, order_by = logFC) %>% pull(Symbol)
)
genes_labeled <- c(genes_labeled,"POSTN")
degenes[genes_labeled,"label"] <- genes_labeled
degenes[degenes$adj.P.Val == 0,"adj.P.Val"] <- 1e-319 + 1e-320*rnorm(sum(degenes$adj.P.Val == 0))
volcano_plot <- 
  ggplot(degenes, aes(x = -logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = Binomial_color_panel) +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
legend <- get_legend(volcano_plot)
volcano_plot <- 
  ggAIplot(volcano_plot + theme(legend.position = "none")) +
  geom_vline(xintercept = c(-0.5,0.5), color = "grey", linetype = "dashed", lwd = 1) +
  geom_text_repel(data = degenes, aes(label = label), size = 2.5)
volcano_plot <- plot_grid(volcano_plot, legend, rel_widths = c(4,1))
ggsave(volcano_plot, file = "./04.figures/R2_01_VSMC_vs_CAF_DEGenes_volcano.pdf", width = 5, height = 4)

universe_entrezid <- bitr(unique(degenes$Symbol), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
universe_entrezid <- unique(universe_entrezid$ENTREZID)
gene_list <- list()
for(cluster_used in as.character(unique(degenes$Grp))){
  genes_symbol <- degenes %>% filter(Grp == cluster_used, adj.P.Val < 0.05, abs(logFC) > 0.5) %>% top_n(n = 100, wt = abs(logFC)) %>% pull(Symbol)
  gene.df <- bitr(genes_symbol, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  gene_list[[cluster_used]] <- unique(gene.df$ENTREZID)
}
cluster_go <- compareCluster(gene_list, fun = "enrichGO",
                             OrgDb = org.Hs.eg.db, ont = "ALL", 
                             universe = universe_entrezid,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = TRUE, pool = TRUE)
description_used <- c(
  "cytochrome-c oxidase activity",
  "oxidative phosphorylation",
  "muscle system process",
  "muscle contraction",
  "vascular smooth muscle cell proliferation",
  "tissue remodeling",
  "extracellular matrix organization",
  "growth factor binding",
  "cell-matrix adhesion",
  "response to transforming growth factor beta",
  "axonogenesis",
  "leukocyte migration"
)
cluster_go_plot <- cluster_go@compareClusterResult %>% 
  filter(Description %in% description_used)
p1 <- cluster_go_plot %>% filter(Cluster == "VSMC") %>% 
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description,p.adjust))) +
  geom_bar(stat = "identity", aes(fill = Cluster)) +
  theme_cowplot(font_size = 7) +
  theme(axis.title.y = element_blank())
p2 <- cluster_go_plot %>% filter(Cluster == "CAF") %>% 
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description,-p.adjust))) +
  geom_bar(stat = "identity", aes(fill = Cluster)) +
  theme_cowplot(font_size = 7) +
  theme(axis.title.y = element_blank())
p <- p2 + p1 + plot_layout(ncol = 1)
ggsave(p, file = "./04.figures/R2_01_VSMC_vs_CAF_GO_bar.pdf", width = 5, height = 3)

genes_used <- c("MYH11","PDGFRA")
plot.data <- cbind(t(as.matrix(HCC_atlas_fibroblast@assays$RNA@data[genes_used,])), HCC_atlas_fibroblast@meta.data[, c("UMAP_1","UMAP_2")])
plot.data <- melt(plot.data, id.vars = c("UMAP_1","UMAP_2"))
plot.data$value <- pmin(plot.data$value, 3.5)
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p <- 
  ggplot(plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 2) +
  facet_wrap(~variable, nrow = 1) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/R2_01_CAF_MYH11_PDGFRA_exp_UMAP.pdf", width = 6.5, height = 3)

saveRDS(HCC_atlas_fibroblast, file = "./02.processed_data/HCC_atlas_fibroblast_release.rds")

# myHSC vs cyHSC----
HCC_atlas_fibroblast <- readRDS("./02.processed_data/HCC_atlas_fibroblast_release.rds")
myHSC_signature <- c("SPP1","MGP","TIMP1","SERPINE2","PTN","COL15A1","COL1A1","ELN","LPL","LTBP2","DPT","COL1A2","LGALS1","SERPINA3","MFAP4","SPARCL1","HTRA1","FBN2","TNC","TIMP3","PLAT","COL6A3","IGF1","LOXL1","NCAM1","GAS6","COL8A1","COL5A2","LOX","COL4A5","CXCL14","MFAP2","COL3A1","FSTL1","MDK","FAM180A","PDGFRL","AEBP1","COL6A2","PLA1A","CRISPLD2","IGFBP4","MMP2")
cyHSC_signature <- c("COLEC11","ANGPTL6","ECM1","COLEC10","RELN","MASP1","CXCL12","SOD3","COL14A1","GDF2","DCN","IGFBP3","WFDC1","TGFBI","RSPO3","CCL2","HGF","PRELP","BMP5","TNFRSF11B","EFEMP1","CXCL1","LGALS3BP","B2M","BMP10","CFB","ADM","IFNAR2","ISG15","TINAGL1","LAMA1","PF4","NGF","APOE","RARRES2")
HCC_atlas_fibroblast <- AddModuleScore(HCC_atlas_fibroblast, features = list(
  myHSC = myHSC_signature,
  cyHSC = cyHSC_signature
))

human_gene_plot.df <- data.frame(
  myHSC = HCC_atlas_fibroblast@meta.data$Cluster1,
  cyHSC = HCC_atlas_fibroblast@meta.data$Cluster2,
  Group = HCC_atlas_fibroblast@meta.data$Release_Cluster) %>% melt
p <- 
  ggplot(human_gene_plot.df, aes(x = Group, y = value)) +
  geom_violin(aes(color = Group, fill = Group), scale = "width") +
  scale_fill_manual(values = CAF_color_panel) +
  scale_color_manual(values = CAF_color_panel) +
  labs(y = "log2(TPM)") +
  theme_cowplot(font_size = 7) +
  facet_grid(variable ~ .) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "null")
ggsave(p, file = "./04.figures/R2_01_CAF_myHSC_cyHSC_violin.pdf", width = 2.25, height = 1.5)

# CAF cluster DEGenes----
HCC_atlas_fibroblast@meta.data$Anno <- plyr::revalue(
  HCC_atlas_fibroblast@meta.data$Release_Cluster,
  replace = c("F1"="VSMC","F2"="VSMC","F3"="VSMC","F4"="LA","F5"="LA","F6"="UNI","F7"="EM","F8"="EM","F9"="IM")
)
HCC_atlas_fibroblast <- SetIdent(HCC_atlas_fibroblast, value = HCC_atlas_fibroblast@meta.data$Anno)
HCC_atlas_fibroblast_markers <- FindAllMarkers(HCC_atlas_fibroblast)
write.csv(HCC_atlas_fibroblast_markers, file = "./03.results/DEGenes/HCC_atlas/HCC_atlas_fibroblast_anno_markers.csv")

# >Endothelium clustering----
#HCC_atlas_endothelium <- HCC_atlas %>% subset(Anno_R1 == "Endothelium")
#saveRDS(HCC_atlas_endothelium, file = "./02.processed_data/RPCA_integration/HCC_atlas_endothelium.rds")
HCC_atlas_endothelium <- readRDS("./02.processed_data/RPCA_integration/HCC_atlas_endothelium.rds")

HCC_atlas_endothelium.list <- SplitObject(HCC_atlas_endothelium, split.by = "orig.ident")
HCC_atlas_endothelium.list  <- lapply(X = HCC_atlas_endothelium.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_endothelium.list)
HCC_atlas_endothelium.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_endothelium.list, anchor.features = features, k.filter = 50)
HCC_atlas_endothelium <- IntegrateData(anchorset = HCC_atlas_endothelium.anchors)
DefaultAssay(HCC_atlas_endothelium) <- "integrated"
HCC_atlas_endothelium <- ScaleData(HCC_atlas_endothelium, verbose = FALSE)
HCC_atlas_endothelium <- RunPCA(HCC_atlas_endothelium, npcs = 30, verbose = FALSE)
HCC_atlas_endothelium <- RunUMAP(HCC_atlas_endothelium, reduction = "pca", dims = 1:10, min.dist = .1, spread = 1)
HCC_atlas_endothelium <- FindNeighbors(HCC_atlas_endothelium, reduction = "pca", dims = 1:10)
HCC_atlas_endothelium <- FindClusters(HCC_atlas_endothelium, resolution = 0.4)
saveRDS(HCC_atlas_endothelium, file = "./02.processed_data/CCA_integration/HCC_atlas_endothelium.rds")

DimPlot(HCC_atlas_endothelium, label = T)
DimPlot(HCC_atlas_endothelium, label = T, split.by = "orig.ident", ncol = 3)
DimPlot(HCC_atlas_endothelium %>% subset(orig.ident == "HCC_Ankur" & Global_Cluster != "Endothelium"), label = T, group.by = "Sub_Cluster")
FeaturePlot(HCC_atlas_endothelium, features = c("PECAM1"))
DimPlot(HCC_atlas_endothelium, group.by = "Tissue")

HCC_atlas_endothelium <- HCC_atlas_endothelium %>% subset(seurat_clusters %in% c(0,1,2,3,4,5,6,7,9)) %>% subset(orig.ident %ni% c("HCC_Fan","HCC_Zhang"))
HCC_atlas_endothelium.list <- SplitObject(HCC_atlas_endothelium, split.by = "orig.ident")
HCC_atlas_endothelium.list  <- lapply(X = HCC_atlas_endothelium.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_endothelium.list)
HCC_atlas_endothelium.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_endothelium.list, anchor.features = features)
HCC_atlas_endothelium <- IntegrateData(anchorset = HCC_atlas_endothelium.anchors)
DefaultAssay(HCC_atlas_endothelium) <- "integrated"
HCC_atlas_endothelium <- ScaleData(HCC_atlas_endothelium, verbose = FALSE)
HCC_atlas_endothelium <- RunPCA(HCC_atlas_endothelium, npcs = 30, verbose = FALSE)
HCC_atlas_endothelium <- RunUMAP(HCC_atlas_endothelium, reduction = "pca", dims = 1:10)
HCC_atlas_endothelium <- FindNeighbors(HCC_atlas_endothelium, reduction = "pca", dims = 1:10, force.recalc = TRUE)
HCC_atlas_endothelium <- FindClusters(HCC_atlas_endothelium, resolution = 0.4)
DefaultAssay(HCC_atlas_endothelium) <- "RNA"
endothelium_markers <- FindAllMarkers(HCC_atlas_endothelium)

HCC_atlas_endothelium <- RunUMAP(HCC_atlas_endothelium, reduction = "pca", dims = 1:8, min.dist = .5, seed.use = 417)
HCC_atlas_endothelium@meta.data[,c("UMAP_1","UMAP_2")] <- HCC_atlas_endothelium@reductions$umap@cell.embeddings[,c(1,2)]
HCC_atlas_endothelium@meta.data$integrated_snn_res.0.2 <- HCC_atlas_endothelium@meta.data$integrated_snn_res.0.4 <- HCC_atlas_endothelium@meta.data$Anno_R1 <- HCC_atlas_endothelium@meta.data$Anno_R2 <- c()
HCC_atlas_endothelium@meta.data$Release_Cluster <- plyr::mapvalues(
  HCC_atlas_endothelium@meta.data$seurat_clusters,
  from = c(0:10),
  to = c("E8","E6","E7","E1","E10","E2","E3","E5","E4","E9","E11")
)
HCC_atlas_endothelium@meta.data$Release_Cluster <- factor(as.character(HCC_atlas_endothelium@meta.data$Release_Cluster), levels = paste0("E",1:11))
HCC_atlas_endothelium@meta.data[HCC_atlas_endothelium@meta.data$orig.ident == 'HCC_Filliol' & HCC_atlas_endothelium@meta.data$Tissue == "Adj Normal","Tissue"] <- "Cirrhosis"

Endo_color_panel <- c(
  "E1" = "#fbe390", "E2" = "#c19595", "E3" = "#fdb069",
  "E4" = "#bf98f0", "E5" = "#c4926d", "E6" = "#666e58",
  "E7" = "#dca48a", "E8" = "#aa4f4c", "E9" = "#ffc9d4",
  "E10" = "#6d5154", "E11" = "#88c4f6")
Tissue_color_panel <- c(
  "Adj Normal" = "#3f9b71", "Tumor" = "#d02c26", "Lymphnode" = "#a600d1", "Cirrhosis" = "#ec8c00"
)
Dataset_color_panel <- c(
  "HCC_Ankur" = "#caacd7", "HCC_Lu" = "#efc96d", "HCC_Ma" = "#ef9d9b", "HCC_Filliol" = "#4b69b4"
)
Endo_HCC_color_panel <- c(
  "IGFBP3+ EC" = "#fbe390", "CD9+ EC" = "#fdb069", "PLPP3+ EC" = "#aa4f4c",
  "PLVAP+ EC_3" = "#c4926d", "PLVAP+ EC_4" = "#666e58", "PLVAP+ EC_9" = "#dca48a",
  "CD320+ EC" = "#6d5154", "TFF3+ EC" = "#88c4f6")

p <- ggplot(HCC_atlas_endothelium@meta.data %>% arrange(Release_Cluster), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = 1.5, alpha = .9) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Endo_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_Endothelium_cluster_UMAP.pdf", width = 3, height = 2.5)

p <- ggplot(HCC_atlas_endothelium@meta.data %>% arrange(Tissue), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Tissue), size = 1.5, alpha = .5) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Tissue_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_Endothelium_tissue_UMAP.pdf", width = 3, height = 2.5)

p <- ggplot(HCC_atlas_endothelium@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster == "Endothelium") %>% arrange(desc(Sub_Cluster)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 2, alpha = .8) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Endo_HCC_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_Endothelium_Ankur_Sub_Cluster_UMAP.pdf", width = 3.25, height = 2.5)

p <- 
  ggplot(HCC_atlas_endothelium@meta.data, aes(x = orig.ident, fill = Release_Cluster)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  scale_fill_manual(values = Endo_color_panel) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))
ggsave(p, file = "./04.figures/R2_01_Endothelium_dataset_bar.pdf", width = 2.5, height = 4)

genes_used <- c("PLVAP","HLA-DRA","ACKR1","IGFBP3","CD9","CD320","TFF3","PLPP3","NRP1","CRHBP","FCN3")
endo.plot.data <- cbind(t(as.matrix(HCC_atlas_endothelium@assays$RNA@data[genes_used,])), HCC_atlas_endothelium@meta.data[, c("UMAP_1","UMAP_2")])
endo.plot.data <- melt(endo.plot.data, id.vars = c("UMAP_1","UMAP_2"))
endo.plot.data$value <- pmin(endo.plot.data$value, 5)
myColorPalette <- colorRampPalette(c("grey",brewer.pal(9, "YlOrRd")))
p <- 
  ggplot(endo.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 1.5) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/R2_01_Endothelium_markers_UMAP.pdf", width = 29, height = 10.5)

saveRDS(HCC_atlas_endothelium, file = "./02.processed_data/HCC_atlas_endothelium_release.rds")

# >Mononuclear clustering----
#HCC_atlas_mononuclear <- HCC_atlas %>% subset(Anno_R1 == "Mononuclear")
#saveRDS(HCC_atlas_mononuclear, file = "./02.processed_data/RPCA_integration/HCC_atlas_mononuclear.rds")
HCC_atlas_mononuclear <- readRDS("./02.processed_data/RPCA_integration/HCC_atlas_mononuclear.rds")

HCC_atlas_mononuclear.list <- SplitObject(HCC_atlas_mononuclear, split.by = "orig.ident")
HCC_atlas_mononuclear.list  <- lapply(X = HCC_atlas_mononuclear.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_mononuclear.list)
HCC_atlas_mononuclear.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_mononuclear.list, anchor.features = features)
HCC_atlas_mononuclear <- IntegrateData(anchorset = HCC_atlas_mononuclear.anchors)
DefaultAssay(HCC_atlas_mononuclear) <- "integrated"
HCC_atlas_mononuclear <- ScaleData(HCC_atlas_mononuclear, verbose = FALSE)
HCC_atlas_mononuclear <- RunPCA(HCC_atlas_mononuclear, npcs = 30, verbose = FALSE)
HCC_atlas_mononuclear <- RunUMAP(HCC_atlas_mononuclear, reduction = "pca", dims = 1:10)
HCC_atlas_mononuclear <- FindNeighbors(HCC_atlas_mononuclear, reduction = "pca", dims = 1:10)
HCC_atlas_mononuclear <- FindClusters(HCC_atlas_mononuclear, resolution = 0.4)
saveRDS(HCC_atlas_mononuclear, file = "./02.processed_data/CCA_integration/HCC_atlas_mononuclear.rds")

HCC_atlas_mononuclear <- HCC_atlas_mononuclear %>% subset(seurat_clusters %ni% c(8,9,11,13))
HCC_atlas_mononuclear.list <- SplitObject(HCC_atlas_mononuclear, split.by = "orig.ident")
HCC_atlas_mononuclear.list  <- lapply(X = HCC_atlas_mononuclear.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_mononuclear.list)
HCC_atlas_mononuclear.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_mononuclear.list, anchor.features = features)
HCC_atlas_mononuclear <- IntegrateData(anchorset = HCC_atlas_mononuclear.anchors)
DefaultAssay(HCC_atlas_mononuclear) <- "integrated"
HCC_atlas_mononuclear <- ScaleData(HCC_atlas_mononuclear, verbose = FALSE)
HCC_atlas_mononuclear <- RunPCA(HCC_atlas_mononuclear, npcs = 30, verbose = FALSE)
HCC_atlas_mononuclear <- RunUMAP(HCC_atlas_mononuclear, reduction = "pca", dims = 1:10)
HCC_atlas_mononuclear <- FindNeighbors(HCC_atlas_mononuclear, reduction = "pca", dims = 1:10, force.recalc = TRUE)
HCC_atlas_mononuclear <- FindClusters(HCC_atlas_mononuclear, resolution = 0.4)
DefaultAssay(HCC_atlas_mononuclear) <- "RNA"
mononuclear_markers <- FindAllMarkers(HCC_atlas_mononuclear)

HCC_atlas_mononuclear <- HCC_atlas_mononuclear %>% subset(seurat_clusters %ni% c(10,11))
HCC_atlas_mononuclear <- RunUMAP(HCC_atlas_mononuclear, reduction = "pca", dims = 1:10, seed.use = 903)
HCC_atlas_mononuclear@meta.data[,c("UMAP_1","UMAP_2")] <- HCC_atlas_mononuclear@reductions$umap@cell.embeddings[,c(1,2)]
HCC_atlas_mononuclear@meta.data$integrated_snn_res.0.2 <- HCC_atlas_mononuclear@meta.data$integrated_snn_res.0.4 <- HCC_atlas_mononuclear@meta.data$Anno_R1 <- c()

# >TNK clustering----
# HCC_atlas_TNK <- HCC_atlas %>% subset(Anno_R1 == "T/NK")
# saveRDS(HCC_atlas_TNK, file = "./02.processed_data/RPCA_integration/HCC_atlas_TNK.rds")
HCC_atlas_TNK <- readRDS("./02.processed_data/RPCA_integration/HCC_atlas_TNK.rds")

HCC_atlas_TNK.list <- SplitObject(HCC_atlas_TNK, split.by = "orig.ident")
HCC_atlas_TNK.list  <- lapply(X = HCC_atlas_TNK.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_TNK.list)
HCC_atlas_TNK.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_TNK.list, anchor.features = features)
HCC_atlas_TNK <- IntegrateData(anchorset = HCC_atlas_TNK.anchors)
DefaultAssay(HCC_atlas_TNK) <- "integrated"
HCC_atlas_TNK <- ScaleData(HCC_atlas_TNK, verbose = FALSE)
HCC_atlas_TNK <- RunPCA(HCC_atlas_TNK, npcs = 30, verbose = FALSE)
HCC_atlas_TNK <- RunUMAP(HCC_atlas_TNK, reduction = "pca", dims = 1:10, min.dist = .1, spread = 1)
HCC_atlas_TNK <- FindNeighbors(HCC_atlas_TNK, reduction = "pca", dims = 1:10)
HCC_atlas_TNK <- FindClusters(HCC_atlas_TNK, resolution = 0.4)
saveRDS(HCC_atlas_TNK, file = "./02.processed_data/CCA_integration/HCC_atlas_TNK.rds")

DimPlot(HCC_atlas_TNK, label = T)
DimPlot(HCC_atlas_TNK, label = T, split.by = "orig.ident", ncol = 3)
DimPlot(HCC_atlas_TNK, group.by = "Global_Cluster", split.by = "orig.ident", ncol = 3, label = T)
DimPlot(HCC_atlas_TNK %>% subset(orig.ident == "HCC_Zhang" & Global_Cluster %in% c("Lymphoid-T","Lymphoid-NK","Lymphoid-T-NK-cycling","Myeloid-Mast","Lymphoid-B","Lymphoid-B-Plasma")),  group.by = "Sub_Cluster", label = T) + scale_color_manual(values = c54)
DimPlot(HCC_atlas_TNK, group.by = "Tissue")

DefaultAssay(HCC_atlas_TNK) <- "RNA"
TNK_markers <- FindAllMarkers(HCC_atlas_TNK, max.cells.per.ident = 2000)

TNK_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC) -> top10
DotPlot(HCC_atlas_TNK, features = unique(top10$gene), dot.scale = 6) + 
  labs(x = "", y = "") + 
  theme_cowplot(font_size = 12) +
  NoLegend() + RotatedAxis()
DotPlot(HCC_atlas_fibroblast, features = fibroblast_markers %>% filter(cluster == 5) %>% top_n(n = 20, wt = avg_logFC) %>% pull(gene), dot.scale = 6) + 
  labs(x = "", y = "") + 
  theme_cowplot(font_size = 12) +
  NoLegend() + RotatedAxis()
DotPlot(HCC_atlas, features = fibroblast_markers %>% filter(cluster == 7) %>% top_n(n = 20, wt = avg_logFC) %>% pull(gene), dot.scale = 6) + 
  labs(x = "", y = "") + 
  theme_cowplot(font_size = 12) +
  NoLegend() + RotatedAxis()

# >Bcell clustering----
# HCC_atlas_Bcell <- HCC_atlas %>% subset(Anno_R1 == "B cells")
# saveRDS(HCC_atlas_Bcell, file = "./02.processed_data/RPCA_integration/HCC_atlas_Bcell.rds")
HCC_atlas_Bcell <- readRDS("./02.processed_data/RPCA_integration/HCC_atlas_Bcell.rds")

HCC_atlas_Bcell.list <- SplitObject(HCC_atlas_Bcell, split.by = "orig.ident")
HCC_atlas_Bcell.list  <- lapply(X = HCC_atlas_Bcell.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_Bcell.list)
HCC_atlas_Bcell.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_Bcell.list, anchor.features = features)
HCC_atlas_Bcell <- IntegrateData(anchorset = HCC_atlas_Bcell.anchors)
DefaultAssay(HCC_atlas_Bcell) <- "integrated"
HCC_atlas_Bcell <- ScaleData(HCC_atlas_Bcell, verbose = FALSE)
HCC_atlas_Bcell <- RunPCA(HCC_atlas_Bcell, npcs = 30, verbose = FALSE)
HCC_atlas_Bcell <- RunUMAP(HCC_atlas_Bcell, reduction = "pca", dims = 1:10, min.dist = .1, spread = 1)
HCC_atlas_Bcell <- FindNeighbors(HCC_atlas_Bcell, reduction = "pca", dims = 1:10)
HCC_atlas_Bcell <- FindClusters(HCC_atlas_Bcell, resolution = 0.4)
saveRDS(HCC_atlas_Bcell, file = "./02.processed_data/CCA_integration/HCC_atlas_Bcell.rds")

# **Round 3----
# >HCC myeloid----
HCC_atlas_myeloid <- merge(HCC_atlas_mononuclear, list(HCC_atlas_TNK %>% subset(seurat_clusters == 12), HCC_atlas_Bcell %>% subset(seurat_clusters %in% c(6,7))))
HCC_atlas_myeloid.list <- SplitObject(HCC_atlas_myeloid, split.by = "orig.ident")
HCC_atlas_myeloid.list  <- lapply(X = HCC_atlas_myeloid.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_myeloid.list)
HCC_atlas_myeloid.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_myeloid.list, anchor.features = features)
HCC_atlas_myeloid <- IntegrateData(anchorset = HCC_atlas_myeloid.anchors)
DefaultAssay(HCC_atlas_myeloid) <- "integrated"
HCC_atlas_myeloid <- ScaleData(HCC_atlas_myeloid, verbose = FALSE)
HCC_atlas_myeloid <- RunPCA(HCC_atlas_myeloid, npcs = 30, verbose = FALSE)
HCC_atlas_myeloid <- RunUMAP(HCC_atlas_myeloid, reduction = "pca", dims = 1:10)
HCC_atlas_myeloid <- FindNeighbors(HCC_atlas_myeloid, reduction = "pca", dims = 1:10, force.recalc = TRUE)
HCC_atlas_myeloid <- FindClusters(HCC_atlas_myeloid, resolution = 0.4)
DefaultAssay(HCC_atlas_myeloid) <- "RNA"

HCC_atlas_myeloid <- HCC_atlas_myeloid %>% subset(seurat_clusters != 12)
HCC_atlas_myeloid.list <- SplitObject(HCC_atlas_myeloid, split.by = "orig.ident")
HCC_atlas_myeloid.list  <- lapply(X = HCC_atlas_myeloid.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_myeloid.list)
HCC_atlas_myeloid.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_myeloid.list, anchor.features = features)
HCC_atlas_myeloid <- IntegrateData(anchorset = HCC_atlas_myeloid.anchors)
DefaultAssay(HCC_atlas_myeloid) <- "integrated"
HCC_atlas_myeloid <- ScaleData(HCC_atlas_myeloid, verbose = FALSE)
HCC_atlas_myeloid <- RunPCA(HCC_atlas_myeloid, npcs = 30, verbose = FALSE)
HCC_atlas_myeloid <- RunUMAP(HCC_atlas_myeloid, reduction = "pca", dims = 1:10)
HCC_atlas_myeloid <- FindNeighbors(HCC_atlas_myeloid, reduction = "pca", dims = 1:10, force.recalc = TRUE)
HCC_atlas_myeloid <- FindClusters(HCC_atlas_myeloid, resolution = 0.6)

DimPlot(HCC_atlas_myeloid, label = T) +
  DimPlot(HCC_atlas_myeloid, label = T, group.by = "Release_Cluster")

HCC_atlas_myeloid@meta.data$Release_Cluster <- plyr::mapvalues(
  HCC_atlas_myeloid@meta.data$seurat_clusters,
  from = c(0:16),
  to = c("M4","M9","M2","M3","M5","M1","M8","M6","M3","M3","M10","M1","M7","M12","M13","M2","M11")
)
HCC_atlas_myeloid@meta.data$Release_Cluster <- factor(as.character(HCC_atlas_myeloid@meta.data$Release_Cluster), levels = paste0("M",1:13))
HCC_atlas_myeloid@meta.data[HCC_atlas_myeloid@meta.data$orig.ident == 'HCC_Filliol' & HCC_atlas_myeloid@meta.data$Tissue == "Adj Normal","Tissue"] <- "Cirrhosis"
HCC_atlas_myeloid@meta.data[HCC_atlas_myeloid@meta.data$Tissue %in% c("Tumor_core","Tumor_edge"),"Tissue"] <- "Tumor"
HCC_atlas_myeloid@meta.data[,c("UMAP_1","UMAP_2")] <- HCC_atlas_myeloid@reductions$umap@cell.embeddings[,c(1,2)]
HCC_atlas_myeloid@meta.data$integrated_snn_res.0.2 <- HCC_atlas_myeloid@meta.data$integrated_snn_res.0.4 <- HCC_atlas_myeloid@meta.data$integrated_snn_res.0.6 <- HCC_atlas_myeloid@meta.data$Anno_R1 <- HCC_atlas_myeloid@meta.data$nCount_integrated <- HCC_atlas_myeloid@meta.data$nFeature_integrated <- c()

DefaultAssay(HCC_atlas_myeloid) <- "RNA"
myeloid_markers <- FindAllMarkers(HCC_atlas_myeloid)

Myeloid_color_panel <- c(
  "M1" = "#415f96", "M2" = "#f2a7b8", "M3" = "#fb942a", "M4" = "#c8ebe8", "M5" = "#516dd3",
  "M6" = "#f5f179", "M7" = "#fbc160", "M8" = "#5cbac8", "M9" = "#d945a0", "M10" = "#97c6df",
  "M11" = "#f97d3c", "M12" = "#8456c0", "M13" = "#dbc5d5"
)
Tissue_color_panel <- c(
  "Adj Normal" = "#3f9b71", "Tumor" = "#d02c26", "Lymphnode" = "#a600d1", "Cirrhosis" = "#ec8c00", "Blood" = "#b5f3ff", "Ascites" = "#f49ddf"
)
Dataset_color_panel <- c(
  "HCC_Ankur" = "#caacd7", "HCC_Lu" = "#efc96d", "HCC_Ma" = "#ef9d9b", "HCC_Filliol" = "#4b69b4", "HCC_Fan" = "#b0e861", "HCC_Zhang" = "#f19a97"
)
Mono_HCC_color_panel <- c(
  "DC1" = "#97c6df", "DC2" = "#d945a0", "pDC" = "#8456c0", 
  "FOLR2+ TAM1_1" = "#c8ebe8", "FOLR2+ TAM1_7" = "#c8ebe8",
  "Mo-derived cells" = "#f2a7b8", "Monocyte" = "#415f96", 
  "MT1G+ TAM3" = "#fbc160", "SPP1+ TAM2" = "#fb942a",
  "Mast" = "#dbc5d5"
)

p <- ggplot(HCC_atlas_myeloid@meta.data %>% arrange(Release_Cluster), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = .8, alpha = .9) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Myeloid_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_myeloid_cluster_UMAP.pdf", width = 3, height = 2.5)

p <- ggplot(HCC_atlas_myeloid@meta.data %>% arrange(Tissue), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Tissue), size = .8, alpha = .5) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Tissue_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_myeloid_tissue_UMAP.pdf", width = 3, height = 2.5)

p <- ggplot(HCC_atlas_myeloid@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster %in% c("Mononuclear","Mast")) %>% arrange(desc(Sub_Cluster)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 1.5, alpha = .8) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Mono_HCC_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_myeloid_Ankur_Sub_Cluster_UMAP.pdf", width = 3.25, height = 2.5)

p <- 
  ggplot(HCC_atlas_myeloid@meta.data, aes(x = orig.ident, fill = Release_Cluster)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  scale_fill_manual(values = Myeloid_color_panel) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))
ggsave(p, file = "./04.figures/R2_01_myeloid_dataset_bar.pdf", width = 2.5, height = 4)

genes_used <- c("FOLR2","FCN1","SPP1","CXCL9","MARCO","MT1G","APOA2","FCER1A","CLEC9A","LAMP3","LILRA4","CPA3")
myeloid.plot.data <- cbind(t(as.matrix(HCC_atlas_myeloid@assays$RNA@data[genes_used,])), HCC_atlas_myeloid@meta.data[, c("UMAP_1","UMAP_2")])
myeloid.plot.data <- melt(myeloid.plot.data, id.vars = c("UMAP_1","UMAP_2"))
myeloid.plot.data$value <- pmin(myeloid.plot.data$value, 5)
myColorPalette <- colorRampPalette(c("grey",brewer.pal(9, "YlOrRd")))
p <- 
  ggplot(myeloid.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 1) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/R2_01_myeloid_markers_UMAP.pdf", width = 29, height = 10.5)

saveRDS(HCC_atlas_myeloid, file = "./02.processed_data/HCC_atlas_myeloid_release.rds")

# >HCC TNK----
DefaultAssay(HCC_atlas_TNK) <- "RNA"
HCC_atlas_TNK <- HCC_atlas_TNK %>% subset(seurat_clusters %ni% c(9,10,12))
HCC_atlas_TNK.list <- SplitObject(HCC_atlas_TNK, split.by = "orig.ident")
HCC_atlas_TNK.list  <- lapply(X = HCC_atlas_TNK.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_TNK.list)
HCC_atlas_TNK.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_TNK.list, anchor.features = features)
HCC_atlas_TNK <- IntegrateData(anchorset = HCC_atlas_TNK.anchors)
DefaultAssay(HCC_atlas_TNK) <- "integrated"
HCC_atlas_TNK <- ScaleData(HCC_atlas_TNK, verbose = FALSE)
HCC_atlas_TNK <- RunPCA(HCC_atlas_TNK, npcs = 30, verbose = FALSE)
HCC_atlas_TNK <- RunUMAP(HCC_atlas_TNK, reduction = "pca", dims = 1:10)
HCC_atlas_TNK <- FindNeighbors(HCC_atlas_TNK, reduction = "pca", dims = 1:10, force.recalc = TRUE)
HCC_atlas_TNK <- FindClusters(HCC_atlas_TNK, resolution = 0.4)

HCC_atlas_TNK <- HCC_atlas_TNK %>% subset(seurat_clusters != 3)
DefaultAssay(HCC_atlas_TNK) <- "integrated"
HCC_atlas_TNK <- ScaleData(HCC_atlas_TNK, verbose = FALSE)
HCC_atlas_TNK <- RunPCA(HCC_atlas_TNK, npcs = 30, verbose = FALSE)
HCC_atlas_TNK <- RunUMAP(HCC_atlas_TNK, reduction = "pca", dims = 1:15, seed.use = 903)
HCC_atlas_TNK <- FindNeighbors(HCC_atlas_TNK, reduction = "pca", dims = 1:15, force.recalc = TRUE)
HCC_atlas_TNK <- FindClusters(HCC_atlas_TNK, resolution = 0.4)
DefaultAssay(HCC_atlas_TNK) <- "RNA"
TNK_markers <- FindAllMarkers(HCC_atlas_TNK, max.cells.per.ident = 2000)

HCC_atlas_TNK@meta.data$Release_Cluster <- plyr::mapvalues(
  HCC_atlas_TNK@meta.data$seurat_clusters,
  from = c(0:10),
  to = c("T1","T4","T9","T6","T8","T2","T3","T10","T7","T11","T5")
)
HCC_atlas_TNK@meta.data$Release_Cluster <- factor(as.character(HCC_atlas_TNK@meta.data$Release_Cluster), levels = paste0("T",1:11))
HCC_atlas_TNK@meta.data[HCC_atlas_TNK@meta.data$orig.ident == 'HCC_Filliol' & HCC_atlas_TNK@meta.data$Tissue == "Adj Normal","Tissue"] <- "Cirrhosis"
HCC_atlas_TNK@meta.data[HCC_atlas_TNK@meta.data$Tissue %in% c("Tumor_core","Tumor_edge"),"Tissue"] <- "Tumor"
HCC_atlas_TNK@meta.data[,c("UMAP_1","UMAP_2")] <- HCC_atlas_TNK@reductions$umap@cell.embeddings[,c(1,2)]
HCC_atlas_TNK@meta.data$integrated_snn_res.0.2 <- HCC_atlas_TNK@meta.data$integrated_snn_res.0.4 <- HCC_atlas_TNK@meta.data$integrated_snn_res.0.6 <- HCC_atlas_TNK@meta.data$Anno_R1 <- c()

TNK_color_panel <- c(
  "T1" = "#f2a7b8", "T2" = "#415f96", "T3" = "#fb942a", "T4" = "#c8ebe8", "T5" = "#516dd3",
  "T6" = "#97c6df", "T7" = "#f5f179", "T8" = "#8457c0", "T9" = "#5cbac8", "T10" = "#fbc160",
  "T11" = "#d945a0"
)
Tissue_color_panel <- c(
  "Adj Normal" = "#3f9b71", "Tumor" = "#d02c26", "Lymphnode" = "#a600d1", "Cirrhosis" = "#ec8c00", "Blood" = "#b5f3ff", "Ascites" = "#f49ddf"
)
Dataset_color_panel <- c(
  "HCC_Ankur" = "#caacd7", "HCC_Lu" = "#efc96d", "HCC_Ma" = "#ef9d9b", "HCC_Filliol" = "#4b69b4", "HCC_Fan" = "#b0e861", "HCC_Zhang" = "#f19a97"
)

p <- ggplot(HCC_atlas_TNK@meta.data %>% arrange(Release_Cluster), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = .1, alpha = .8) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = TNK_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_TNK_cluster_UMAP.pdf", width = 3, height = 2.5)

saveRDS(HCC_atlas_TNK, file = "./02.processed_data/HCC_atlas_TNK_release.rds")

# >HCC B cell----
HCC_atlas_Bcell <- merge(HCC_atlas_Bcell %>% subset(seurat_clusters %in% c(0,4,5,2,3,8,9)), HCC_atlas_TNK %>% subset(seurat_clusters == 10))
HCC_atlas_Bcell.list <- SplitObject(HCC_atlas_Bcell, split.by = "orig.ident")
HCC_atlas_Bcell.list <- lapply(X = HCC_atlas_Bcell.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_Bcell.list)
HCC_atlas_Bcell.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_Bcell.list, anchor.features = features)
HCC_atlas_Bcell <- IntegrateData(anchorset = HCC_atlas_Bcell.anchors)
DefaultAssay(HCC_atlas_Bcell) <- "integrated"
HCC_atlas_Bcell <- ScaleData(HCC_atlas_Bcell, verbose = FALSE)
HCC_atlas_Bcell <- RunPCA(HCC_atlas_Bcell, npcs = 30, verbose = FALSE)
HCC_atlas_Bcell <- RunUMAP(HCC_atlas_Bcell, reduction = "pca", dims = 1:10, min.dist = .1, spread = 1)
HCC_atlas_Bcell <- FindNeighbors(HCC_atlas_Bcell, reduction = "pca", dims = 1:10)
HCC_atlas_Bcell <- FindClusters(HCC_atlas_Bcell, resolution = 0.4)

HCC_atlas_Bcell <- HCC_atlas_Bcell %>% subset(seurat_clusters %in% c(0,1,3,4,7,8))
DefaultAssay(HCC_atlas_Bcell) <- "RNA"
HCC_atlas_Bcell.list <- SplitObject(HCC_atlas_Bcell, split.by = "orig.ident")
HCC_atlas_Bcell.list <- lapply(X = HCC_atlas_Bcell.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_Bcell.list)
HCC_atlas_Bcell.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_Bcell.list, anchor.features = features)
HCC_atlas_Bcell <- IntegrateData(anchorset = HCC_atlas_Bcell.anchors)
DefaultAssay(HCC_atlas_Bcell) <- "integrated"
HCC_atlas_Bcell <- ScaleData(HCC_atlas_Bcell, verbose = FALSE)
HCC_atlas_Bcell <- RunPCA(HCC_atlas_Bcell, npcs = 30, verbose = FALSE)
HCC_atlas_Bcell <- RunUMAP(HCC_atlas_Bcell, reduction = "pca", dims = 1:10, min.dist = .1, spread = 1)
HCC_atlas_Bcell <- FindNeighbors(HCC_atlas_Bcell, reduction = "pca", dims = 1:10)
HCC_atlas_Bcell <- FindClusters(HCC_atlas_Bcell, resolution = 0.4)

HCC_atlas_Bcell <- HCC_atlas_Bcell %>% subset(seurat_clusters %ni% c(8,9,10,11))
DefaultAssay(HCC_atlas_Bcell) <- "RNA"
HCC_atlas_Bcell.list <- SplitObject(HCC_atlas_Bcell, split.by = "orig.ident")
HCC_atlas_Bcell.list <- lapply(X = HCC_atlas_Bcell.list , FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = HCC_atlas_Bcell.list)
HCC_atlas_Bcell.anchors <- FindIntegrationAnchors(object.list = HCC_atlas_Bcell.list, anchor.features = features)
HCC_atlas_Bcell <- IntegrateData(anchorset = HCC_atlas_Bcell.anchors)
DefaultAssay(HCC_atlas_Bcell) <- "integrated"
HCC_atlas_Bcell <- ScaleData(HCC_atlas_Bcell, verbose = FALSE)
HCC_atlas_Bcell <- RunPCA(HCC_atlas_Bcell, npcs = 30, verbose = FALSE)
HCC_atlas_Bcell <- RunUMAP(HCC_atlas_Bcell, reduction = "pca", dims = 1:10, spread = 3)
HCC_atlas_Bcell <- FindNeighbors(HCC_atlas_Bcell, reduction = "pca", dims = 1:10)
HCC_atlas_Bcell <- FindClusters(HCC_atlas_Bcell, resolution = 0.4)
DefaultAssay(HCC_atlas_Bcell) <- "RNA"

HCC_atlas_Bcell@meta.data$Release_Cluster <- plyr::mapvalues(
  HCC_atlas_Bcell@meta.data$seurat_clusters,
  from = c(0:7),
  to = c("B1","B6","B5","B7","B4","B8","B2","B3")
)
HCC_atlas_Bcell@meta.data$Release_Cluster <- factor(as.character(HCC_atlas_Bcell@meta.data$Release_Cluster), levels = paste0("B",1:8))
HCC_atlas_Bcell@meta.data[HCC_atlas_Bcell@meta.data$orig.ident == 'HCC_Filliol' & HCC_atlas_Bcell@meta.data$Tissue == "Adj Normal","Tissue"] <- "Cirrhosis"
HCC_atlas_Bcell@meta.data[HCC_atlas_Bcell@meta.data$Tissue %in% c("Tumor_core","Tumor_edge"),"Tissue"] <- "Tumor"
HCC_atlas_Bcell@meta.data[,c("UMAP_1","UMAP_2")] <- HCC_atlas_Bcell@reductions$umap@cell.embeddings[,c(1,2)]
HCC_atlas_Bcell@meta.data$integrated_snn_res.0.2 <- HCC_atlas_Bcell@meta.data$integrated_snn_res.0.4 <- HCC_atlas_Bcell@meta.data$Anno_R1 <- HCC_atlas_Bcell@meta.data$nCount_integrated <- HCC_atlas_Bcell@meta.data$nFeature_integrated <- c()

Bcell_color_panel <- c(
  "B1" = "#fbe390", "B2" = "#c19595", "B3" = "#fdb069",
  "B4" = "#bf98f0", "B5" = "#c4926d", "B6" = "#666e58",
  "B7" = "#dca48a", "B8" = "#aa4f4c"
)

p <- ggplot(HCC_atlas_Bcell@meta.data %>% arrange(Release_Cluster), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = 1.5, alpha = .9) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Bcell_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_Bcell_cluster_UMAP.pdf", width = 3, height = 2.5)

saveRDS(HCC_atlas_Bcell, file = "./02.processed_data/HCC_atlas_Bcell_release.rds")

# >All cells umap----
HCC_atlas <- readRDS("./02.processed_data/RPCA_integration/HCC_atlas.combined.rds")
HCC_atlas_hepatocyte <- HCC_atlas %>% subset(integrated_snn_res.0.2 %in% c(2,3,6,9,10,11))
HCC_atlas_hepatocyte@meta.data$integrated_snn_res.0.2 <- HCC_atlas_hepatocyte@meta.data$seurat_clusters <- HCC_atlas_hepatocyte@meta.data$Anno_R1 <- c()
HCC_atlas_hepatocyte@meta.data$Release_Global_Cluster <- "Hepatocyte"
HCC_atlas_TNK <- readRDS("./02.processed_data/HCC_atlas_TNK_release.rds")
HCC_atlas_TNK@meta.data$Release_Global_Cluster <- "T/NK"
HCC_atlas_Bcell <- readRDS("./02.processed_data/HCC_atlas_Bcell_release.rds")
HCC_atlas_Bcell@meta.data$Release_Global_Cluster <- "B"
HCC_atlas_myeloid <- readRDS("./02.processed_data/HCC_atlas_myeloid_release.rds")
HCC_atlas_myeloid@meta.data$Release_Global_Cluster <- "Myeloid"
HCC_atlas_fibroblast <- readRDS("./02.processed_data/HCC_atlas_fibroblast_release.rds")
HCC_atlas_fibroblast@meta.data$Release_Global_Cluster <- "Fibroblast"
HCC_atlas_endothelium <- readRDS("./02.processed_data/HCC_atlas_endothelium_release.rds")
HCC_atlas_endothelium@meta.data$Release_Global_Cluster <- "Endothelium"
cells_used <- c(colnames(HCC_atlas_hepatocyte), colnames(HCC_atlas_TNK), colnames(HCC_atlas_Bcell), colnames(HCC_atlas_myeloid), colnames(HCC_atlas_fibroblast), colnames(HCC_atlas_endothelium))
HCC_atlas <- HCC_atlas[,cells_used]
HCC_atlas@meta.data$Release_Global_Cluster <- c(
  rep("Hepatocyte",ncol(HCC_atlas_hepatocyte)),
  rep("T/NK",ncol(HCC_atlas_TNK)),
  rep("B",ncol(HCC_atlas_Bcell)),
  rep("Myeloid",ncol(HCC_atlas_myeloid)),
  rep("Fibroblast",ncol(HCC_atlas_fibroblast)),
  rep("Endothelium",ncol(HCC_atlas_endothelium))
)
DefaultAssay(HCC_atlas) <- "integrated"
HCC_atlas <- RunUMAP(HCC_atlas, reduction = "pca", dims = 1:10)

HCC_atlas@meta.data$integrated_snn_res.0.2 <- c()
HCC_atlas@meta.data$Release_Cluster <- c(
  rep("Hepatocyte",ncol(HCC_atlas_hepatocyte)),
  as.character(HCC_atlas_TNK@meta.data$Release_Cluster),
  as.character(HCC_atlas_Bcell@meta.data$Release_Cluster),
  as.character(HCC_atlas_myeloid@meta.data$Release_Cluster),
  as.character(HCC_atlas_fibroblast@meta.data$Release_Cluster),
  as.character(HCC_atlas_endothelium@meta.data$Release_Cluster)
)
HCC_atlas@meta.data[,c("Global_UMAP_1","Global_UMAP_2")] <- HCC_atlas@reductions$umap@cell.embeddings[,c(1,2)]
HCC_atlas@meta.data[,c("UMAP_1","UMAP_2")] <- rbind(
  data.frame(UMAP_1 = rep("NA",ncol(HCC_atlas_hepatocyte)),
             UMAP_2 = rep("NA",ncol(HCC_atlas_hepatocyte)), 
             stringsAsFactors = F),
  HCC_atlas_TNK@meta.data[,c("UMAP_1","UMAP_2")],
  HCC_atlas_Bcell@meta.data[,c("UMAP_1","UMAP_2")],
  HCC_atlas_myeloid@meta.data[,c("UMAP_1","UMAP_2")],
  HCC_atlas_fibroblast@meta.data[,c("UMAP_1","UMAP_2")],
  HCC_atlas_endothelium@meta.data[,c("UMAP_1","UMAP_2")]
)

Global_color_panel <- c(
  "T/NK" = "#415f96", "B" = "#f2a7b8", "Myeloid" = "#fb942a", 
  "Endothelium" = "#5cbac8", "Fibroblast" = "#8456c0", "Hepatocyte" = "#fbc160"
)
Dataset_color_panel <- c(
  "HCC_Ankur" = "#caacd7", "HCC_Lu" = "#efc96d", "HCC_Ma" = "#ef9d9b", "HCC_Filliol" = "#4b69b4", "HCC_Fan" = "#97c6df", "HCC_Zhang" = "#f97d3c"
)

p <- ggplot(HCC_atlas@meta.data %>% arrange(Release_Global_Cluster), aes(x = Global_UMAP_1, y = Global_UMAP_2)) +
  geom_point(aes(color = Release_Global_Cluster), size = .1, alpha = .9) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Global_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_Allcell_cluster_UMAP.pdf", width = 3, height = 2.5)

p <- ggplot(HCC_atlas@meta.data, aes(x = Global_UMAP_1, y = Global_UMAP_2)) +
  geom_point(aes(color = orig.ident), size = .1, alpha = .9) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = Dataset_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/R2_01_Allcell_dataset_UMAP.pdf", width = 3, height = 2.5)

p <- 
  ggplot(HCC_atlas@meta.data, aes(x = orig.ident, fill = Release_Global_Cluster)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  scale_fill_manual(values = Global_color_panel) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))
ggsave(p, file = "./04.figures/R2_01_Allcell_dataset_bar.pdf", width = 2, height = 2)

saveRDS(HCC_atlas, file = "./02.processed_data/HCC_atlas_release.rds")

HCC_atlas_metadata <- HCC_atlas@meta.data
write.csv(HCC_atlas_metadata, file = "./02.processed_data/HCC_atlas_release_metadata.csv")

table(HCC_atlas_metadata$orig.ident, HCC_atlas_metadata$Release_Global_Cluster)
