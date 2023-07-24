setwd("/Volumes/ZiyiLi/PostDoc_Project/Onco-fetal Project/")
source("./onco-fetal_utils.R")
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggalluvial)
library(RColorBrewer)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load data----
HCC_atlas <- readRDS("./HCC atlas/HCC_atlas_all_release.rds")
HCC_atlas_Bcell <- readRDS("./HCC atlas/HCC_atlas_Bcell_release.rds")
HCC_atlas_TNK <- readRDS("./HCC atlas/HCC_atlas_TNK_release.rds")
HCC_atlas_endothelium <- readRDS("./HCC atlas/HCC_atlas_endothelium_release.rds")
HCC_atlas_fibroblast <- readRDS("./HCC atlas/HCC_atlas_fibroblast_release.rds")
HCC_atlas_myeloid <- readRDS("./HCC atlas/HCC_atlas_myeloid_release.rds")

# Fig.1b----
ggplot(HCC_atlas@meta.data %>% arrange(Release_Global_Cluster), aes(x = Global_UMAP_1, y = Global_UMAP_2)) +
  geom_point(aes(color = Release_Global_Cluster), size = .1, alpha = .9) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

ggplot(HCC_atlas_endothelium@meta.data %>% arrange(Release_Cluster), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = 1.5, alpha = .9) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

ggplot(HCC_atlas_fibroblast@meta.data %>% arrange(desc(Release_Cluster)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = 2.5, alpha = .8) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

ggplot(HCC_atlas_myeloid@meta.data %>% arrange(Release_Cluster), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = .8, alpha = .9) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Fig.1c----
ggplot(HCC_atlas_endothelium@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster == "Endothelium") %>% arrange(desc(Sub_Cluster)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 2, alpha = .8) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Fig.1d----
genes_used <- c("PLVAP","HLA-DRA","ACKR1","IGFBP3","CD9","CD320","TFF3","PLPP3","NRP1","CRHBP","FCN3")
endo.plot.data <- cbind(t(as.matrix(HCC_atlas_endothelium@assays$RNA@data[genes_used,])), HCC_atlas_endothelium@meta.data[, c("UMAP_1","UMAP_2")])
endo.plot.data <- melt(endo.plot.data, id.vars = c("UMAP_1","UMAP_2"))
endo.plot.data$value <- pmin(endo.plot.data$value, 5)
myColorPalette <- colorRampPalette(c("grey",brewer.pal(9, "YlOrRd")))
ggplot(endo.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 0.5) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))

# Fig.1e----
ggplot(HCC_atlas_myeloid@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster %in% c("Mononuclear","Mast")) %>% arrange(desc(Sub_Cluster)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 1.5, alpha = .8) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Fig.1f----
genes_used <- c("FOLR2","FCN1","SPP1","CXCL9","MARCO","MT1G","FCER1A","CLEC9A","LAMP3","LILRA4","CPA3")
myeloid.plot.data <- cbind(t(as.matrix(HCC_atlas_myeloid@assays$RNA@data[genes_used,])), HCC_atlas_myeloid@meta.data[, c("UMAP_1","UMAP_2")])
myeloid.plot.data <- melt(myeloid.plot.data, id.vars = c("UMAP_1","UMAP_2"))
myeloid.plot.data$value <- pmin(myeloid.plot.data$value, 5)
myColorPalette <- colorRampPalette(c("grey",brewer.pal(9, "YlOrRd")))
ggplot(myeloid.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 0.5) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))

# Fig.1g----
HCC_atlas_fibroblast <- SetIdent(HCC_atlas_fibroblast, value = HCC_atlas_fibroblast@meta.data$Release_Cluster)
fibroblast_markers <- FindAllMarkers(HCC_atlas_fibroblast)
universe_entrezid <- bitr(row.names(HCC_atlas_fibroblast@assays$RNA@counts), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
universe_entrezid <- unique(universe_entrezid$ENTREZID)
gene_list <- list()
for(cluster_used in as.character(unique(fibroblast_markers$cluster))){
  genes_symbol <- fibroblast_markers %>% filter(cluster == cluster_used, p_val_adj < 0.05, avg_log2FC > 0) %>% top_n(n = 100, wt = avg_log2FC) %>% pull(gene)
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

# Fig.1h----
genes_used <- c("MYH11","PDGFRA")
plot.data <- cbind(t(as.matrix(HCC_atlas_fibroblast@assays$RNA@data[genes_used,])), HCC_atlas_fibroblast@meta.data[, c("UMAP_1","UMAP_2")])
plot.data <- melt(plot.data, id.vars = c("UMAP_1","UMAP_2"))
plot.data$value <- pmin(plot.data$value, 3.5)
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
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

# Fig.1j----
cells_used <- HCC_atlas_fibroblast@meta.data %>% filter(Release_Cluster %in% c("F1","F2","F3","F7","F8")) %>% pull(CellName)
degenes <- LIMMA(
  HCC_atlas_fibroblast@assays$RNA@data[,cells_used],
  c("VSMC","CAF")[as.factor(HCC_atlas_fibroblast@meta.data[cells_used,"Release_Cluster"] %ni% c("F1","F2","F3"))]
)
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
p2 + p1 + plot_layout(ncol = 1)

# Fig.1k----
ggplot(HCC_atlas_fibroblast@meta.data, aes(x = Tissue, fill = Release_Cluster)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))

# Extended Data Fig.1a----
ggplot(HCC_atlas_fibroblast@meta.data, aes(x = Release_Cluster, fill = orig.ident)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))

# Extended Data Fig.1b----
HCC_atlas_fibroblast <- SetIdent(HCC_atlas_fibroblast, value = HCC_atlas_fibroblast@meta.data$Release_Cluster)
fibroblast_markers <- FindAllMarkers(HCC_atlas_fibroblast)
genes_used <- fibroblast_markers %>% group_by(cluster) %>% filter(pct.2 < 0.8, pct.1 > 0.5) %>% slice_max(n = 8, order_by = avg_log2FC) %>% pull(gene) %>% unique()
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
ggplot(plot.data, aes(x = Group, y = Gene)) +
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

# Extended Data Fig.1c----
human_gene_plot.df <- data.frame(
  PI16 = HCC_atlas_fibroblast@assays$RNA@data["PI16",],
  CD34 = HCC_atlas_fibroblast@assays$RNA@data["CD34",],
  DPT = HCC_atlas_fibroblast@assays$RNA@data["DPT",],
  Group = HCC_atlas_fibroblast@meta.data$Release_Cluster) %>% melt
ggplot(human_gene_plot.df, aes(x = Group, y = value)) +
  geom_violin(aes(color = Group, fill = Group), scale = "width") +
  labs(y = "log2(TPM)") +
  theme_cowplot(font_size = 7) +
  facet_grid(variable ~ .) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "null")

# Extended Data Fig.1d----
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
ggplot(degenes, aes(x = -logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = Binomial_color_panel) +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1)) +
  geom_vline(xintercept = c(-0.5,0.5), color = "grey", linetype = "dashed", lwd = 1) +
  geom_text_repel(data = degenes, aes(label = label), size = 2.5)

# Extended Data Fig.1e----
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
ggplot(human_gene_plot.df, aes(x = Group, y = value)) +
  geom_violin(aes(color = Group, fill = Group), scale = "width") +
  labs(y = "log2(TPM)") +
  theme_cowplot(font_size = 7) +
  facet_grid(variable ~ .) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "null")

# Extended Data Fig.1f----
plot.list <- list()
for(i in c("1L","1T","1N","2L","2T","2N")){
  plot.list[[i]] <- list()
  for(gene in c("MYH11","POSTN","ABCA8","RELN","DPT")){
    plot.list[[i]][[gene]] <- SpatialFeaturePlot(HCC_GJ[[i]], features = gene, pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
  }
}

# Supplementary Fig.1a----
HCC_atlas_metadata <- read.csv("./HCC atlas/HCC_atlas_metadata_batch_effect.csv", row.names = 1)
plot.list1 <- plot.list2 <- list()
for(algorithm in c("raw","fastmnn","harmony","liger","RPCA","scanorama")){
  plot.list1[[algorithm]] <- ggplot(HCC_atlas_metadata %>% arrange(Release_Global_Cluster), aes_string(x = paste0(algorithm,"_UMAP_1"), y = paste0(algorithm,"_UMAP_2"))) +
    geom_point(aes(color = Release_Global_Cluster), size = .1, alpha = .9) +
    theme_cowplot(font_size = 7) +
    guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
  
  plot.list2[[algorithm]] <- ggplot(HCC_atlas_metadata, aes_string(x = paste0(algorithm,"_UMAP_1"), y = paste0(algorithm,"_UMAP_2"))) +
    geom_point(aes(color = orig.ident), size = .1, alpha = .9) +
    theme_cowplot(font_size = 7) +
    guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
}

# Supplementary Fig.1b----
HCC_atlas <- SetIdent(HCC_atlas, value = HCC_atlas@meta.data$Release_Global_Cluster)
Global_marker_genes <- FindAllMarkers(HCC_atlas, max.cells.per.ident = 5000)
genes_for_plot1 <- Global_marker_genes %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 30) %>% pull(gene) %>% unique()
genes_for_plot2 <- c(
  "CD3D","CD3E","CD8A","CD8B","LEF1","CD4","IL7R","TIGIT","FOXP3","GZMB","KLRF1","NKG7","CD160",
  "CD79A","CD19","IGHG1","JCHAIN",
  "CPA3","CD74","CD163","LYZ","C1QA","C1QB","S100A8","S100A9","CD1C","XCR1","HLA-DPA1","HLA-DQA1","HLA-DRA",
  "TAGLN","MYL9","ACTA2","MYH11","COL3A1","COL4A1","COL1A1",
  "PLVAP","FCN3","PLPP1","PECAM1","CD9","TFF3","CAV1","VWF","SPARC",
  "ALB","APOC3","APOA2","FABP1","FGA","FGB","ORM2","ORM1")
genes_for_plot <- unique(c(genes_for_plot1, genes_for_plot2))
markers_mean_exp <- 
  aggregate(
    as.matrix(t(HCC_atlas@assays$RNA@data[genes_for_plot,])),
    list(Cluster = HCC_atlas@meta.data$Release_Global_Cluster),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Cluster
markers_mean_exp$Cluster <- c()
markers_plot_matrix <- apply(markers_mean_exp, 2, zscore) %>% t()
markers_plot_matrix <- markers_plot_matrix[,c("Hepatocyte","Fibroblast","Endothelium","Myeloid","B","T/NK")]
color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), max(markers_plot_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
max.avg <- apply(markers_plot_matrix, 1, which.max)
gene_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) > 0){
    temp <- data.frame(gene = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), cluster = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
}
genes_labeled <- genes_for_plot2
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

# Supplementary Fig.1c-d----
batch_entropy <- list("raw" = list(), "fastmnn" = list(), "harmony" = list(),
                      "liger" = list(), "RPCA" = list(), "scanorama" = list())
cells_used <- list()
for(i in 1:10){
  cat("Round ",i,"\n")
  cells_used[[i]] <- sample(HCC_atlas_metadata$CellName, 50000)
  for(algorithm in c("raw","fastmnn","harmony","liger","RPCA","scanorama")){
    cat("Algorithm:",algorithm,"\n")
    cat("Calculate batch entropy.\n")
    batch_entropy[[algorithm]][[i]] <- BatchEntropy(
      input.data = HCC_atlas_metadata[cells_used[[i]],c(paste0(algorithm,"_UMAP_1"),paste0(algorithm,"_UMAP_2"))],
      group.id = list(Dataset = HCC_atlas_metadata[cells_used[[i]],"orig.ident"],
                      Cluster = HCC_atlas_metadata[cells_used[[i]],"Release_Global_Cluster"]),
      k.used = 30,
      dimension.used = "raw")
  }
}
batch_entropy_df <- data.frame(algorithm = NA, rep = NA, batch = NA, median = NA, mean = NA)
for(algorithm in c("raw","fastmnn","harmony","liger","RPCA","scanorama")){
  for(rep in 1:10){
    batch_entropy_df <- rbind(batch_entropy_df, c(algorithm, rep, "Dataset", median(batch_entropy[[algorithm]][[rep]][["Dataset"]]), mean(batch_entropy[[algorithm]][[rep]][["Dataset"]])))
    batch_entropy_df <- rbind(batch_entropy_df, c(algorithm, rep, "Cluster", median(batch_entropy[[algorithm]][[rep]][["Cluster"]]), mean(batch_entropy[[algorithm]][[rep]][["Cluster"]])))
  }
}
batch_entropy_df <- batch_entropy_df[-1,]
batch_entropy_df$median <- as.numeric(batch_entropy_df$median)
batch_entropy_df$mean <- as.numeric(batch_entropy_df$mean)
ggplot(batch_entropy_df %>% filter(batch == "Dataset"), aes(x = reorder(algorithm,mean), y = mean)) +
  geom_boxplot(lwd = .2, outlier.size = .1) +
  labs(x = "", y = "Mean of Dataset-based Entropy") +
  theme_cowplot(font_size = 7) +
  RotatedAxis()
ggplot(batch_entropy_df %>% filter(batch == "Cluster"), aes(x = reorder(algorithm,mean), y = mean)) +
  geom_boxplot(lwd = .2, outlier.size = .1) +
  labs(x = "", y = "Mean of Cluster-based Entropy") +
  theme_cowplot(font_size = 7) +
  RotatedAxis()

# Supplementary Fig.1e----
ggplot(HCC_atlas@meta.data, aes(x = orig.ident, fill = Release_Global_Cluster)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))

# Supplementary Fig.1f----
ggplot(HCC_atlas_TNK@meta.data %>% arrange(Release_Cluster), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = .1, alpha = .8) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Supplementary Fig.1g----
ggplot(HCC_atlas_Bcell@meta.data %>% arrange(Release_Cluster), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Release_Cluster), size = 1.5, alpha = .9) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Supplementary Fig.2a----
ggplot(HCC_atlas_endothelium@meta.data %>% arrange(Tissue), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Tissue), size = 1) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Supplementary Fig.2b----
# see Fig.1d

# Supplementary Fig.2c----
ggplot(HCC_atlas_endothelium@meta.data, aes(x = orig.ident, fill = Release_Cluster)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))

# Supplementary Fig.2d----
plot.data <- HCC_atlas_endothelium@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster == "Endothelium") %>% dplyr::select(Release_Cluster, Sub_Cluster) %>% table() %>% as.data.frame() %>% filter(Freq > 0) %>% ggplot(aes(axis1 = Release_Cluster, axis2 = Sub_Cluster, y = Freq)) +
  geom_alluvium(aes(fill = Sub_Cluster), curve_type = "sigmoid", width = 1/8) +
  geom_stratum(width = 1/8) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),
                   expand = c(0.15, 0.05)) +
  theme_void() +
  theme(legend.position = "none")

# Supplementary Fig.2e----
cells_used <- HCC_atlas_endothelium@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster == "Endothelium") %>% pull(CellName)
Confusion_heatmap_new(HCC_atlas_endothelium@meta.data[cells_used,"Sub_Cluster"], HCC_atlas_endothelium@meta.data[cells_used,"Release_Cluster"])
cross_table <- as.matrix(table(HCC_atlas_endothelium@meta.data[cells_used,"Sub_Cluster"], HCC_atlas_endothelium@meta.data[cells_used,"Release_Cluster"]))
cross_table_per <- cross_table/rowSums(cross_table)

# Supplementary Fig.2f----
ggplot(HCC_atlas_myeloid@meta.data %>% arrange(Tissue), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Tissue), size = .5) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Supplementary Fig.2g----
# see Fig.1f

# Supplementary Fig.2h----
ggplot(HCC_atlas_myeloid@meta.data, aes(x = orig.ident, fill = Release_Cluster)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))

# Supplementary Fig.2i----
HCC_atlas_myeloid@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster == "Mononuclear") %>% dplyr::select(Release_Cluster, Sub_Cluster) %>% table() %>% as.data.frame() %>% filter(Freq > 0) %>%
  ggplot(aes(axis1 = Release_Cluster, axis2 = Sub_Cluster, y = Freq)) +
  geom_alluvium(aes(fill = Sub_Cluster), curve_type = "sigmoid", width = 1/8) +
  geom_stratum(width = 1/8) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),
                   expand = c(0.15, 0.05)) +
  theme_void() +
  theme(legend.position = "none")

# Supplementary Fig.2j----
cells_used <- HCC_atlas_myeloid@meta.data %>% filter(orig.ident == "HCC_Ankur", Global_Cluster == "Mononuclear") %>% pull(CellName)
Confusion_heatmap_new(HCC_atlas_myeloid@meta.data[cells_used,"Sub_Cluster"], HCC_atlas_myeloid@meta.data[cells_used,"Release_Cluster"])