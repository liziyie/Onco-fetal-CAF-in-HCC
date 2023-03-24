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
library(ComplexHeatmap)

# >Load dataset----
CAF_HCC <- readRDS("./02.processed_data/CAF_HCC.rds")
CAF_HCC_color_panel <- c(
  "CAF" = "#62A3C7", "Fib" = "#68AB9F", "SDC2+ CAF" = "#CAB2D6",
  "MYH11+ CAF" = "#F69459", "HSP+ CAF" = "#1C76B3", "POSTN+ CAF" = "#693D99", 
  "APOA2+ CAF" = "#FB9A99", "MT1M+ Fib" = "#B2DE89", "ABCAB+ Fib" = "#4DAC47"
)

# Plot global cluster in all HCC cells----
HCC_seu <- readRDS("./02.processed_data/HCC_seu.rds")
Global_Cluster_color_panel <- c("Lymphocyte" = "#69B4CE","B cell" = "#69B4CE", "T cell" = "#ECAFCF", "ILC" = "#A0D7C9", "Mast" = "#BB4A94", "Mononuclear" = "#F3746C", "Hepatocyte" = "#CAA57D", "Fibroblast" = "#C35338", "Doublet" = "#C35338", "Endothelium" = "#EAA944")
p <- ggplot(HCC_seu@meta.data %>% filter(Global_Cluster != "Doublet"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Global_Cluster), size = 0.5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(name = "", values = Global_Cluster_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/01.HCC_seu_UMAP.pdf", width = 8, height = 6.5)

HCC_seu@active.ident <- factor(HCC_seu$Global_Cluster)
genes_used <- c("PTPRC","PECAM1","ALB","PDGFRB","TAGLN","TPSAB1","EMCN","KRT18","RGS5","ACTA2")
HCC_seu.cells.used <- HCC_seu@meta.data %>% filter(Global_Cluster != "Doublet") %>% pull(CellName)
HCC_seu.plot.data <- cbind(t(as.matrix(HCC_seu@assays$RNA@data[genes_used,HCC_seu.cells.used])), HCC_seu@meta.data[HCC_seu.cells.used, c("UMAP_1","UMAP_2")])
HCC_seu.plot.data <- melt(HCC_seu.plot.data, id.vars = c("UMAP_1","UMAP_2"))
HCC_seu.plot.data$value <- pmin(HCC_seu.plot.data$value, 3.5)
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p <- 
  ggplot(HCC_seu.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = .5) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/01.HCC_seu_markers_UMAP.pdf", width = 27, height = 10.5)

# >Cluster and annotation----
# >>UMAP plot----
p <- ggplot(CAF_HCC@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/01.CAF_HCC_UMAP.pdf", width = 8, height = 6.5)

p1 <- ggplot(CAF_HCC@meta.data %>% filter(NTF == "Tumor"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  lims(x = c(-5,6.5), y = c(-7,5.5)) + 
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p1.pdf <- ggAIplot(p1 + theme(legend.position = "none"))
legend <- get_legend(p1)
p2 <- ggplot(CAF_HCC@meta.data %>% filter(NTF == "Adj Normal"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 5, alpha = .8) +
  theme_cowplot(font_size = 12) +
  lims(x = c(-5,6.5), y = c(-7,5.5)) + 
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p2.pdf <- ggAIplot(p2 + theme(legend.position = "none"))
p <- plot_grid(p1.pdf, p2.pdf, legend, rel_widths = c(5,5,1), nrow = 1)
ggsave(p, file = "./04.figures/01.CAF_HCC_UMAP_Tissue.pdf", width = 11, height = 5)

# >>Marker genes heatmap----
CAF_HCC_markers <- FindAllMarkers(CAF_HCC)
write.csv(CAF_HCC_markers, file = "./Results/DEGenes/CAF_HCC_all_markers.csv")
genes_used <- CAF_HCC_markers %>% group_by(cluster) %>% filter(pct.2 < 0.8, pct.1 > 0.5) %>% slice_max(n = 20, order_by = avg_logFC) %>% pull(gene) %>% unique()
genes_used <- unique(c(genes_used, "FAP", "PDPN", "SDC2", "MYH11", "CD36", "POSTN","MYH11"))
markers_mean_exp <- 
  aggregate(
    as.matrix(t(CAF_HCC@assays$RNA@data[genes_used,])),
    list(Cluster = CAF_HCC$Sub_Cluster),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Cluster
markers_mean_exp$Cluster <- c()
markers_plot_matrix <- apply(markers_mean_exp, 2, zscore) %>% t()
color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), max(markers_plot_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
max.avg <- apply(markers_plot_matrix, 1, which.max)
gene_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) > 0){
    temp <- data.frame(gene = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), cluster = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
}
markers_plot_matrix_quantile <- quantile(markers_plot_matrix, c(0.01, 0.98))
markers_plot_matrix <- pmax(markers_plot_matrix, markers_plot_matrix_quantile[1])
markers_plot_matrix <- pmin(markers_plot_matrix, markers_plot_matrix_quantile[2])
genes_labeled <- c("COX4I2","MEF2C","BCAM","MYH11","COX6C","CD151","S100A4","HSPA6","BAG3","DNAJB1","HSPH1","HSPA1B","HSPA1A","FOS","HSP90AB1","LHFP","LMNA","HES1","RHOB","PTGDS","CCDC80","COL1A1","POSTN","COL3A1","VCAN","S100A10","ASPN","FAP","COL1A2","FN1","CSF3","CCL2","IGFBP6","DCN","HLA-A","SDC2","HINT1","APOA2","APOC3","IGKC","ALB","CD36","MT1A","MT1M","MYC","IRF1","BTG2","APOD","ABCA8","ITM2A","JUNB","JUN")
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
pdf("./04.figures/01.CAF_HCC_markers_heatmap.pdf", width = 4.5, height = 9.5)
draw(p)
dev.off()

# >>Marker genes bubble heatmap----
CAF_HCC_markers <- read.csv("./03.results/DEGenes/CAF_HCC_all_markers.csv", row.names = 1, stringsAsFactors = F)
genes_used <- CAF_HCC_markers %>% group_by(cluster) %>% filter(pct.2 < 0.8, pct.1 > 0.5) %>% slice_max(n = 8, order_by = avg_logFC) %>% pull(gene) %>% unique()
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
ggsave(p, file = "./04.figures/01.CAF_HCC_markers_bubble_heatmap.pdf", width = 4, height = 12)

# iCAF and apCAF signal----
iCAF_markers <- c("PDGFRA","CXCL12","CFD","DPT")
iCAF_markers <- intersect(iCAF_markers, row.names(CAF_HCC))
myCAF_markers <- c("ACTA2","TAGLN","MYH11","TPM1")
myCAF_markers <- intersect(myCAF_markers, row.names(CAF_HCC))
apCAF_markers <- c("PTGIS","HLA-DQB1","CD74")
apCAF_markers <- intersect(apCAF_markers, row.names(CAF_HCC))

CAF_HCC@meta.data$iCAF_score <- colMeans(CAF_HCC@assays$RNA@data[iCAF_markers,])
CAF_HCC@meta.data$myCAF_score <- colMeans(CAF_HCC@assays$RNA@data[myCAF_markers,])
CAF_HCC@meta.data$apCAF_score <- colMeans(CAF_HCC@assays$RNA@data[apCAF_markers,])
p <- 
  CAF_HCC@meta.data[,c("Sub_Cluster","iCAF_score","myCAF_score","apCAF_score")] %>% filter(Sub_Cluster %in% c("CAF","MYH11+ CAF","HSP+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF")) %>% melt(.,id.vars = "Sub_Cluster") %>% 
  ggplot(., aes(x = Sub_Cluster, y = value)) +
  geom_boxplot(aes(color = Sub_Cluster, fill = Sub_Cluster), alpha = .4, outlier.size = -1) +
  facet_wrap(~variable, scales = "free_y", nrow = 3) +
  scale_color_manual(values = CAF_HCC_color_panel) +
  scale_fill_manual(values = CAF_HCC_color_panel) +
  theme_cowplot(font_size = 12) +
  ggsignif::geom_signif(comparisons = list(
    c("MYH11+ CAF","HSP+ CAF"),
    c("MYH11+ CAF","POSTN+ CAF"),
    c("SDC2+ CAF","POSTN+ CAF")
  ), step = .1) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(p, file = "./04.figures/01.iCAF_apCAF_score.pdf", width = 2.5, height = 8)

p <- ggplot(CAF_HCC@meta.data, aes(x = iCAF_score, y = apCAF_score)) +
  geom_jitter(aes(color = Sub_Cluster), size = 5, alpha = .8, width = .01, height = .02) +
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p) +
  geom_vline(xintercept = .3, linetype = "dashed") + 
  geom_hline(yintercept = .3, linetype = "dashed")
ggsave(p.pdf, file = "./04.figures/01.iCAF_apCAF_score_scatter.pdf", width = 5.5, height = 4)

# >>MYH11 vs POSTN----
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
genes_labeled <- c(genes_labeled,"POSTN")
degenes[genes_labeled,"label"] <- genes_labeled
volcano_plot <- 
  ggplot(degenes, aes(x = -logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = Binomial_color_panel) +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  theme_cowplot(font_size = 16) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
legend <- get_legend(volcano_plot)
volcano_plot <- 
  ggAIplot(volcano_plot + theme(legend.position = "none")) +
  geom_vline(xintercept = c(-0.5,0.5), color = "grey", linetype = "dashed", lwd = 1) +
  geom_text_repel(data = degenes, aes(label = label), size = 4)
volcano_plot <- plot_grid(volcano_plot, legend, rel_widths = c(4,1))
ggsave(volcano_plot, file = "./04.figures/01.myCAF_vs_apCAF_DEGenes_volcano.pdf", width = 7.75, height = 6)

# >>Proliferation ability----
pro_genes <- intersect(c("FEN1","FOXM1","MCM2","MCM3","MCM4","MCM5","MCM6","MKI67","PCNA","PLK1","TOP2A","TYMS","CCNB1","CCNE1"), row.names(CAF_HCC@assays$RNA@data))
CAF_HCC$Pro_Score <- colMeans(CAF_HCC@assays$RNA@data[pro_genes,])
p <- ggplot(CAF_HCC@meta.data %>% filter(!(Sub_Cluster %in% c("Fib","MT1M+ Fib","ABCAB+ Fib")), Pro_Score < 0.5), aes(x = Sub_Cluster, y = Pro_Score)) +
  geom_boxplot(aes(color = Sub_Cluster)) +
  scale_color_manual(values = CAF_HCC_color_panel) +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ggsignif::geom_signif(comparisons = list(c("POSTN+ CAF","MYH11+ CAF"),
                                           c("POSTN+ CAF","SDC2+ CAF")),
                        step_increase = .1)
ggsave(p, file = "./04.figures/01.CAF_HCC_proliferation_boxplot.pdf", height = 5, width = 3)

# >>Expression of marker genes----
genes_used <- c("MYH11","HSPA6","SDC2","FAP","POSTN","COL4A1","ENG","APOA2","MT1M","ABCA8")
CAF_HCC.cells.used <- CAF_HCC$CellID
CAF_HCC.plot.data <- cbind(t(as.matrix(CAF_HCC@assays$RNA@data[genes_used,CAF_HCC.cells.used])), CAF_HCC@meta.data[CAF_HCC.cells.used, c("UMAP_1","UMAP_2")])
CAF_HCC.plot.data <- melt(CAF_HCC.plot.data, id.vars = c("UMAP_1","UMAP_2"))
CAF_HCC.plot.data$value <- pmin(CAF_HCC.plot.data$value, 3.5)
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p <- 
  ggplot(CAF_HCC.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 3) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/01.CAF_HCC_markers_UMAP.pdf", width = 25, height = 10.5)

# >>GSVA pathway enrichment----
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(clusterProfiler)
data(c2BroadSets)
exp <- as.matrix(CAF_HCC@assays$RNA@data)
exp_entrezid <- bitr(row.names(CAF_HCC@assays$RNA@data), fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
exp <- as.matrix(exp[exp_entrezid$SYMBOL,])
row.names(exp) <- exp_entrezid$ENTREZID
gsva.matrix <- gsva(exp, c2BroadSets)
gsva.matrix.mean <- aggregate(t(gsva.matrix), list(as.character(CAF_HCC@meta.data[colnames(gsva.matrix),"Sub_Cluster"])), mean)
row.names(gsva.matrix.mean) <- gsva.matrix.mean$Group.1
gsva.matrix.mean$Group.1 <- c()
gsva.matrix.mean <- t(gsva.matrix.mean)

FindDEGenes(expression_matrix = gsva.matrix,
            groupid = as.character(CAF_HCC@meta.data[colnames(gsva.matrix),"Sub_Cluster"]),
            out.prefix = "./03.results/DEGenes/GSVA/CAF_HCC_GSVA",
            cutoff = 1,
            logFC = 0)
GSVA_deg <- read.csv("./03.results/DEGenes/GSVA/CAF_HCC_GSVA_cutoff1_all_de_genes_raw.csv", stringsAsFactors = F)
GSVA_deg$dataset <- stringr::str_split_fixed(GSVA_deg$Symbol,"_",2)[,1]
pathway_used <- GSVA_deg %>% filter(dataset %in% c("KEGG","REACTOME")) %>% group_by(Group) %>% slice_max(order_by = AUC, n = 35) %>% pull(Symbol) %>% unique()

gsva.matrix.mean.filter <- gsva.matrix.mean[pathway_used,]
mat_for_plot <- t(apply(gsva.matrix.mean.filter,1,scale))
mat_quan <- quantile(gsva.matrix.mean.filter, c(.001,.999))
mat_for_plot <- pmax(mat_for_plot, mat_quan[1])
mat_for_plot <- pmin(mat_for_plot, mat_quan[2])
mat_for_plot[mat_for_plot < 0.1 & mat_for_plot > -0.1] <- 0
colnames(mat_for_plot) <- colnames(gsva.matrix.mean.filter)
mat_for_plot <- mat_for_plot[,c("CAF","MYH11+ CAF","HSP+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF","Fib","MT1M+ Fib","ABCAB+ Fib")]

max.avg <- apply(mat_for_plot, 1, which.max)
pathway_order <- c()
for(i in 1:ncol(mat_for_plot)){
  if(sum(max.avg == i) > 0){
    temp <- data.frame(Pathway = names(sort(mat_for_plot[names(max.avg)[max.avg == i],i], decreasing = T)), cluster = colnames(mat_for_plot)[i], stringsAsFactors = F)
    pathway_order <- rbind(pathway_order, temp)
  }
}
mat_for_plot <- mat_for_plot[pathway_order$Pathway,]
color_used = circlize::colorRamp2(seq(min(mat_for_plot), max(mat_for_plot), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))

pathway_labeled <- unique(c(
  "REACTOME_MUSCLE_CONTRACTION","REACTOME_SMOOTH_MUSCLE_CONTRACTION",
  "KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION",
  "REACTOME_SIGNALING_BY_NOTCH",
  "REACTOME_GLUCOSE_REGULATION_OF_INSULIN_SECRETION",
  "REACTOME_INTEGRATION_OF_ENERGY_METABOLISM",
  "REACTOME_SIGNALING_BY_WNT", "REACTOME_STABILIZATION_OF_P53",
  "REACTOME_APOPTOSIS", "KEGG_MAPK_SIGNALING_PATHWAY", 
  "REACTOME_SIGNALLING_BY_NGF",
  "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "REACTOME_NCAM1_INTERACTIONS", "REACTOME_SIGNALING_BY_PDGF",
  "REACTOME_PLATELET_ADHESION_TO_EXPOSED_COLLAGEN", 
  "KEGG_FOCAL_ADHESION", "REACTOME_PLATELET_ACTIVATION",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES", "KEGG_P53_SIGNALING_PATHWAY",
  "REACTOME_GPCR_LIGAND_BINDING", "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "KEGG_NITROGEN_METABOLISM",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "REACTOME_COMPLEMENT_CASCADE",
  "KEGG_PPAR_SIGNALING_PATHWAY","KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY", "REACTOME_METABOLISM_OF_AMINO_ACIDS",
  "KEGG_PURINE_METABOLISM", "KEGG_GLUTATHIONE_METABOLISM",
  "KEGG_NITROGEN_METABOLISM","KEGG_PROGESTERONE_MEDIATED_OOCYTE_MATURATION",
  "REACTOME_ADENYLATE_CYCLASE_ACTIVATING_PATHWAY","KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION",
  "KEGG_CARDIAC_MUSCLE_CONTRACTION","KEGG_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_GLUCOSE_REGULATION_OF_INSULIN_SECRETION","REACTOME_SIGNALING_BY_NOTCH",
  "KEGG_MAPK_SIGNALING_PATHWAY","BIOCARTA_CCR5_PATHWAY","KEGG_REGULATION_OF_ACTIN_CYTOSKELETON",
  "REACTOME_SEMAPHORIN_INTERACTIONS","KEGG_N_GLYCAN_BIOSYNTHESIS","REACTOME_NCAM1_INTERACTIONS",
  "REACTOME_SIGNALING_BY_PDGF","KEGG_ECM_RECEPTOR_INTERACTION",
  "REACTOME_PLATELET_ADHESION_TO_EXPOSED_COLLAGEN","REACTOME_SIGNALING_BY_WNT",
  "BIOCARTA_P53HYPOXIA_PATHWAY","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "REACTOME_REGULATION_OF_BETA_CELL_DEVELOPMENT","KEGG_PPAR_SIGNALING_PATHWAY",
  "REACTOME_LIPOPROTEIN_METABOLISM","REACTOME_HDL_MEDIATED_LIPID_TRANSPORT",
  "REACTOME_METABOLISM_OF_LIPIDS_AND_LIPOPROTEINS","BIOCARTA_THELPER_PATHWAY",
  "REACTOME_STEROID_METABOLISM","BIOCARTA_DC_PATHWAY","BIOCARTA_IL5_PATHWAY","BIOCARTA_CTLA4_PATHWAY",
  "REACTOME_DOWNSTREAM_EVENTS_IN_GPCR_SIGNALING","KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_INOSITOL_PHOSPHATE_METABOLISM","KEGG_ERBB_SIGNALING_PATHWAY","KEGG_P53_SIGNALING_PATHWAY",
  "REACTOME_INNATE_IMMUNITY_SIGNALING","REACTOME_CLASSICAL_ANTIBODY_MEDIATED_COMPLEMENT_ACTIVATION",
  "BIOCARTA_41BB_PATHWAY","BIOCARTA_TOLL_PATHWAY","KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "BIOCARTA_IL12_PATHWAY","KEGG_JAK_STAT_SIGNALING_PATHWAY"
))
ha <-  HeatmapAnnotation(
  Cluster = factor(colnames(mat_for_plot),levels = c("CAF","MYH11+ CAF","HSP+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF","Fib","MT1M+ Fib","ABCAB+ Fib")),
  col = list(Cluster = CAF_HCC_color_panel))
p <- Heatmap(mat_for_plot, 
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = F,
             show_column_names = F,
             column_names_gp = gpar(fontsize = 5),
             col = color_used,
             top_annotation = ha,
             name = "GSVA score") +
  rowAnnotation(link = anno_mark(
    at = which(pathway_order$Pathway %in% pathway_labeled),
    labels =  pathway_order$Pathway[which(pathway_order$Pathway %in% pathway_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
pdf("./04.figures/01.CAF_HCC_GSVA_heatmap2.pdf", width = 10, height = 7)
draw(p)
dev.off()

save(gsva.matrix, file = "./02.processed_data/CAF_HCC_GSVA.rda")

# >>Progeny----
library(progeny)
CAF_HCC <- progeny(CAF_HCC, scale=FALSE, organism="Human", top=500, perm=1, 
                   return_assay = TRUE)
CAF_HCC <- ScaleData(CAF_HCC, assay = "progeny") 
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(CAF_HCC, slot = "scale.data", 
                               assay = "progeny"))) %>%
  tibble::rownames_to_column("Cell") %>%
  tidyr::gather(Pathway, Activity, -Cell)
progeny_scores_df$Sub_Cluster <- CAF_HCC@meta.data[progeny_scores_df$Cell,"Sub_Cluster"]
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, Sub_Cluster) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  tidyr::spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
paletteLength <- 100
myColor <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))(paletteLength)
progenyBreaks <- c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
pheatmap::pheatmap(t(summarized_progeny_scores_df[,-1]),
                   fontsize=7, fontsize_row = 7,
                   color=myColor, breaks = progenyBreaks,
                   main = "PROGENy", angle_col = 45,
                   treeheight_col = 0,  border_color = NA,
                   filename = "./04.figures/01.CAF_HCC_PROGENy_heatmap.pdf",
                   width = 3.5, height = 3)

# >>Proportion----
# Cluster by tissue
p <- 
  ggplot(CAF_HCC@meta.data, aes(x = Sub_Cluster, fill = NTF)) +
  geom_bar(position = "fill", alpha = 0.7, color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 16) + 
  scale_fill_manual(values = Tissue_color_panel) +
  scale_x_discrete(limits = sort(unique(CAF_HCC@meta.data$Sub_Cluster))) +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(override.aes = list(size = 6), nrow = 2)) 
ggsave(p, file = "./04.figures/01.CAF_HCC_tissue_proportions.pdf", width = 6, height = 6)
# Cluster by patient
p <- 
  ggplot(CAF_HCC@meta.data, aes(x = Sub_Cluster, fill = PatientID)) +
  geom_bar(position = "fill", alpha = 0.7, color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 16) + 
  scale_fill_manual(values = c54) +
  scale_x_discrete(limits = sort(unique(CAF_HCC@meta.data$Sub_Cluster))) +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(override.aes = list(size = 6), ncol = 2)) 
ggsave(p, file = "./04.figures/01.CAF_HCC_patient_proportions.pdf", width = 6, height = 6)

# >Trajectory analysis----
# >>Monocle3----
library(monocle3)
library(SeuratWrappers)
library(magrittr)
CAF_HCC_UMAP.cds <- as.cell_data_set(CAF_HCC)
CAF_HCC_UMAP.cds <- cluster_cells(CAF_HCC_UMAP.cds, k = 30, partition_qval = 0.05, resolution = 0.01, reduction_method = "UMAP")
colData(CAF_HCC_UMAP.cds)$Cluster <- colData(CAF_HCC_UMAP.cds)$Sub_Cluster
CAF_HCC_UMAP.cds@clusters$UMAP$clusters <- as.character(CAF_HCC_UMAP.cds@colData$Cluster)
names(CAF_HCC_UMAP.cds@clusters$UMAP$clusters) <- CAF_HCC_UMAP.cds@colData$index
plot_grid(
  plot_cells(
    CAF_HCC_UMAP.cds,
    color_cells_by = "Cluster",
    show_trajectory_graph = FALSE),
  plot_cells(
    CAF_HCC_UMAP.cds, 
    color_cells_by = "partition", 
    show_trajectory_graph = FALSE)
)

CAF_HCC_UMAP.cds <- learn_graph(CAF_HCC_UMAP.cds, learn_graph_control = list(minimal_branch_len = 3))
CAF_HCC_UMAP.cds <- order_cells(CAF_HCC_UMAP.cds)
save(CAF_HCC_UMAP.cds, file = "./Results/Trajectory/CAF_HCC_Monocle3.rda")

p <- ggplot(CAF_HCC@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 4, alpha = .8) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  theme(legend.position = "none")
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/01.CAF_HCC_monocle3_UMAP.pdf", width = 4, height = 4)

p1 <- plot_cells(
  CAF_HCC_UMAP.cds, 
  color_cells_by = "Sub_Cluster", 
  cell_size = 0,
  trajectory_graph_segment_size = 1,
  label_groups_by_cluster = FALSE, 
  label_cell_groups = FALSE,
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  scale_color_manual(values = CAF_HCC_color_panel) +
  labs(x = "UMAP1", y = "UMAP2")
p1 <- plot_grid(p1 + theme(legend.position = "none"))
ggsave(p1, file = "./04.figures/01.CAF_HCC_monocle3_UMAP_trajectory.pdf", width = 4, height = 4)

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
p2.pdf <- ggAIplot(p2 + theme(legend.position = "none"))
legend <- get_legend(p2)
p <- plot_grid(p2.pdf, legend, nrow = 1, rel_widths = c(4,1))
ggsave(p, file = "./04.figures/01.CAF_HCC_monocle3_UMAP_pseudotime.pdf", width = 5, height = 4)

CAF_HCC <- AddMetaData(CAF_HCC, metadata = CAF_HCC_UMAP.cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "Monocle3_pseudotime")

p <- CAF_HCC@meta.data %>% filter(Sub_Cluster != "APOA2+ CAF") %>% filter(!(Sub_Cluster == "CAF" & Monocle3_pseudotime > 8)) %>% mutate(Grp = plyr::revalue(Sub_Cluster, replace = c(
  "ABCAB+ Fib" = "Grp1", "MT1M+ Fib" = "Grp1",
  "Fib" = "Grp1", "CAF" = "Grp2", 
  "HSP+ CAF" = "Grp4", "MYH11+ CAF" = "Grp4",
  "SDC2+ CAF" = "Grp5", "POSTN+ CAF" = "Grp5"
))) %>% 
  mutate(Grp = factor(Grp, levels = c("Grp1","Grp2","Grp3","Grp4","Grp5"))) %>%
  mutate(Sub_Cluster = factor(as.character(Sub_Cluster),levels = c("MT1M+ Fib","Fib","ABCAB+ Fib","CAF","MYH11+ CAF","HSP+ CAF","POSTN+ CAF","SDC2+ CAF"))) %>%
  ggplot(aes(x = Monocle3_pseudotime, y = Sub_Cluster)) +
  geom_boxplot(aes(color = Sub_Cluster)) +
  geom_point(aes(color = Sub_Cluster), position = position_jitter(h=0.2), size = .8) +
  facet_grid(Grp~., scales = "free", space = "free")  +
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  theme_bw(base_size = 10) +
  ylab("") + xlab("Monocle3 pseudotime") +
  theme(legend.position = "none", 
        strip.text = element_blank(),
        strip.background = element_blank())
ggsave(p, file = "./04.figures/01.CAF_HCC_Monocle3_pseudotime_boxplot.pdf", width = 5, height = 4)

# >>Diffusion map----
library(destiny)
library(scater)
library(car)
library(rgl)
CAF_HCC_markers <- read.csv("./Results/DEGenes/CAF_HCC_all_markers.csv", head = T, row.names = 1)
genes_used <- CAF_HCC_markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05, abs(avg_logFC) > 0.8) %>% pull(gene) %>% unique()

set.seed(123)
CAF_HCC.ct <- data.frame(t(as.matrix(CAF_HCC@assays$RNA@counts[genes_used,])))
CAF_HCC.ct <- cbind(CellName = row.names(CAF_HCC.ct), CAF_HCC.ct)
CAF_HCC.es <- as.ExpressionSet(CAF_HCC.ct)
CAF_HCC.es$Cluster <- CAF_HCC@meta.data$Sub_Cluster
CAF_HCC.dm <- DiffusionMap(CAF_HCC.es, rotate = T)
CAF_HCC.dpt <- DPT(CAF_HCC.dm)
save(CAF_HCC.dm, CAF_HCC.dpt, file = "./Results/Trajectory/CAF_HCC_Diffusion_map.rda")

CAF_HCC$DC1 <- CAF_HCC.dm$DC1;CAF_HCC$DC2 <- CAF_HCC.dm$DC2;CAF_HCC$DC3 <- CAF_HCC.dm$DC3
CAF_HCC$DPT <- CAF_HCC.dpt$DPT2

p <- ggplot(CAF_HCC@meta.data, aes(x = DC1, y = DC2)) +
  geom_point(aes(color = Sub_Cluster), size = 3) +
  theme_cowplot(font_size = 12) +
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/01.CAF_HCC_Diffusion_map.pdf", width = 8, height = 6.5)

p <- ggplot(CAF_HCC@meta.data, aes(x = DC1, y = DC2)) +
  geom_point(aes(color = DPT), size = 3) +
  theme_cowplot(font_size = 12) +
  scale_color_gradientn(colours = myPalette(100))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/01.CAF_HCC_Diffusion_pseudotime.pdf", width = 7, height = 6.5)

p <- CAF_HCC@meta.data %>% mutate(Grp = plyr::revalue(Sub_Cluster, replace = c(
  "ABCAB+ Fib" = "Grp1", "MT1M+ Fib" = "Grp1",
  "Fib" = "Grp2", "CAF" = "Grp2", "APOA2+ CAF" = "Grp2",
  "HSP+ CAF" = "Grp3", "MYH11+ CAF" = "Grp4",
  "SDC2+ CAF" = "Grp5", "POSTN+ CAF" = "Grp5"
))) %>% 
  mutate(Grp = factor(Grp, levels = c("Grp1","Grp2","Grp3","Grp4","Grp5"))) %>%
  mutate(Sub_Cluster = factor(as.character(Sub_Cluster),levels = rev(levels(Sub_Cluster)))) %>%
  ggplot(aes(x = DPT, y = Sub_Cluster)) +
  geom_boxplot(aes(color = Sub_Cluster)) +
  geom_point(aes(color = Sub_Cluster), position = position_jitter(h=0.2), size = .8) +
  facet_grid(Grp~., scales = "free", space = "free")  +
  scale_color_manual(name = "", values = CAF_HCC_color_panel) +
  theme_bw(base_size = 10) +
  ylab("") + xlab("diffusion pseudotime (dpt)") +
  theme(legend.position = "none", 
        strip.text = element_blank(),
        strip.background = element_blank())
ggsave(p, file = "./04.figures/01.CAF_HCC_Diffusion_pseudotime_boxplot.pdf", width = 5, height = 4)

# DPT DEGenes (SDC2+CAF to POSTN+ CAF)----
cells_used <- CAF_HCC@meta.data %>% filter(Sub_Cluster %in% c("SDC2+ CAF","POSTN+ CAF")) %>% pull(CellID)
genes_used <- row.names(CAF_HCC)[rowSums(CAF_HCC[,cells_used]) > 20]
genes_used <- genes_used[-grep("^RP|^MT-",genes_used)]
CAF_HCC_subset <- CAF_HCC[genes_used, cells_used]
CAF_HCC_subset$dpt_order = rank(CAF_HCC_subset$DPT)
CAF_HCC_subset$dpt_order_perc = CAF_HCC_subset$dpt_order / max(CAF_HCC_subset$dpt_order)
p = ggplot(CAF_HCC_subset@meta.data, aes(x=dpt_order_perc, fill=Sub_Cluster)) +
  geom_density(size=0.2, adjust=1.5, n=200, trim=F, color="black", alpha=0.75) +
  ggpubr::theme_classic2() +
  scale_fill_manual(values=CAF_HCC_color_panel)
ggsave(p, file = "./04.figures/01.CAF_HCC_Diffusion_pseudotime_density.pdf", width = 6, height = 2)

gam.pval <- c()
for(gene in genes_used){
  dat = data.frame(exp=CAF_HCC_subset@assays$RNA@data[gene,], time=CAF_HCC_subset$dpt_order_perc)
  fit = mgcv::gam(exp ~ time, data=dat)
  summ = summary(fit)
  coef = unname(coef(fit)[2])
  p = summ[4][[1]][2]
  df = data.frame(geneSymbol=gene, coef=coef, pval=p)
  gam.pval <- rbind(gam.pval, df)
}
colnames(gam.pval) = c("geneSymbol", "coef", "pval")
gam.pval = gam.pval[order(gam.pval$pval,decreasing=F),]
gam.pval$adj.p = p.adjust(gam.pval$pval, method = "BH")
saveRDS(gam.pval, file="./Results/DEGenes/CAF_HCC_DPT_DEGenes.rds")
gam.pval.subset = gam.pval[(gam.pval$adj.p<0.01 & abs(gam.pval$coef)>.9),]
gam.pval.subset = arrange(gam.pval.subset, adj.p, desc(abs(coef)))
plot.gene = as.character(gam.pval.subset$geneSymbol)
plot.gene = plot.gene[plot.gene != "MALAT1"]
highlight.genes = c("SPARCL1","ADIRF","HSPA1A","ANXA5","FKBP1A","TMEM14C","ANXA1","FAP","LTBP2","COL10A1","POSTN","COL5A1","CTSK","LAMP2","HSPA6","S100A13","CTSA","MDH1","PLD3","LAP3","CTLA","INHBA","LY6E","SCP2","IGFBP6","VCAN","S100A10","IL32")
annot=rowAnnotation(gene=anno_mark(at=match(highlight.genes, plot.gene), labels=highlight.genes))
colSet <- list()
colSet[['dpt_order_perc']] = circlize::colorRamp2(c(0,max(CAF_HCC_subset$dpt_order_perc)),c("white","darkblue"))
colSet[['cluster.name']] = structure(CAF_HCC_color_panel, names=names(CAF_HCC_color_panel))
set.seed(1)
sce <- as.SingleCellExperiment(CAF_HCC_subset)
sce$meta.cluster = as.character(sce$Sub_Cluster)
sce$cluster.name = factor(as.character(sce$Sub_Cluster))
hm.obj = heatmap_sm(sce[plot.gene,], assay.name="logcounts",
                    colSet=colSet,
                    out.prefix="./04.figures/01.CAF_HCC_DPT_heatmap.pdf",
                    columns=c("dpt_order_perc"), 
                    ncell.downsample=500,
                    use_raster=T, raster_quality=2,
                    columns.order="dpt_order_perc",
                    show_column_names=F, show_row_names=F,
                    row_names_side="right",
                    row_gap = unit(1, "mm"), column_gap = unit(0, "mm"),
                    border=F, ann.bar.height=0.5, palette.name="RdYlBu", z.hi = 2.5,
                    pdf.width=10.5, pdf.height=9, do.scale=F,
                    do.clustering.row=T, do.clustering.col=F,
                    dend.row=T,  right_annotation=annot,
                    clustering.distance="cosine",
                    clustering.method="ward.D2",
                    row_km=5, row_km_repeats=100,
)
genes_clu = row_order(hm.obj)
if(is.null(names(genes_clu))){
  names(genes_clu) = 1:length(genes_clu)
}
for (i in names(genes_clu)){
  genes_clu[[i]] = plot.gene[genes_clu[[i]]]
}

# DPT DEGenes (MYH11+CAF vs POSTN+ CAF)----
cells_used_MYH11 <- CAF_HCC@meta.data %>% filter(Sub_Cluster %in% c("CAF","MYH11+ CAF")) %>% pull(CellID)
genes_used_MYH11 <- row.names(CAF_HCC)[rowSums(CAF_HCC[,cells_used_MYH11]) > 20]
genes_used_MYH11 <- genes_used_MYH11[-grep("^RP|^MT-",genes_used_MYH11)]
CAF_HCC_subset_MYH11 <- CAF_HCC[genes_used_MYH11, cells_used_MYH11]
CAF_HCC_subset_MYH11$dpt_order = rank(CAF_HCC_subset_MYH11$DPT)
CAF_HCC_subset_MYH11$dpt_order_perc = CAF_HCC_subset_MYH11$dpt_order / max(CAF_HCC_subset_MYH11$dpt_order)
gam.pval_MYH11 <- c()
for(gene in genes_used_MYH11){
  dat = data.frame(exp=CAF_HCC_subset_MYH11@assays$RNA@data[gene,], time=CAF_HCC_subset_MYH11$dpt_order_perc)
  fit = mgcv::gam(exp ~ time, data=dat)
  summ = summary(fit)
  coef = unname(coef(fit)[2])
  p = summ[4][[1]][2]
  df = data.frame(geneSymbol=gene, coef=coef, pval=p)
  gam.pval_MYH11 <- rbind(gam.pval_MYH11, df)
}
colnames(gam.pval_MYH11) = c("geneSymbol", "coef", "pval")
gam.pval_MYH11 = gam.pval_MYH11[order(gam.pval_MYH11$pval,decreasing=F),]
gam.pval_MYH11$adj.p = p.adjust(gam.pval_MYH11$pval, method = "BH")

cells_used_POSTN <- CAF_HCC@meta.data %>% filter(Sub_Cluster %in% c("CAF","POSTN+ CAF")) %>% pull(CellID)
genes_used_POSTN <- row.names(CAF_HCC)[rowSums(CAF_HCC[,cells_used_POSTN]) > 20]
genes_used_POSTN <- genes_used_POSTN[-grep("^RP|^MT-",genes_used_POSTN)]
CAF_HCC_subset_POSTN <- CAF_HCC[genes_used_POSTN, cells_used_POSTN]
CAF_HCC_subset_POSTN$dpt_order = rank(CAF_HCC_subset_POSTN$DPT)
CAF_HCC_subset_POSTN$dpt_order_perc = CAF_HCC_subset_POSTN$dpt_order / max(CAF_HCC_subset_POSTN$dpt_order)
gam.pval_POSTN <- c()
for(gene in genes_used_POSTN){
  dat = data.frame(exp=CAF_HCC_subset_POSTN@assays$RNA@data[gene,], time= CAF_HCC_subset_POSTN$dpt_order_perc)
  fit = mgcv::gam(exp ~ time, data=dat)
  summ = summary(fit)
  coef = unname(coef(fit)[2])
  p = summ[4][[1]][2]
  df = data.frame(geneSymbol=gene, coef=coef, pval=p)
  gam.pval_POSTN <- rbind(gam.pval_POSTN, df)
}
colnames(gam.pval_POSTN) = c("geneSymbol", "coef", "pval")
gam.pval_POSTN = gam.pval_POSTN[order(gam.pval_POSTN$pval,decreasing=F),]
gam.pval_POSTN$adj.p = p.adjust(gam.pval_POSTN$pval, method = "BH")

gam.pval.subset_MYH11 = gam.pval_MYH11 %>% filter(adj.p < 0.01 & abs(coef) > 1.2) %>% arrange(adj.p, desc(abs(coef)))
gam.pval.subset_POSTN = gam.pval_POSTN %>% filter(adj.p < 0.01 & abs(coef) > 1.2) %>% arrange(adj.p, desc(abs(coef)))
plot.gene = unique(c(as.character(gam.pval.subset_POSTN$geneSymbol),as.character(gam.pval.subset_MYH11$geneSymbol)))
plot.gene = plot.gene[!plot.gene %in% c("SPARCL1","TPM2","MALAT1","PTMA","TAGLN")]
highlight.genes = c("TMEM176B","LY6E","POSTN","FAP","PRSS23","S100A13","INHBA","FSTL1","VCAN","S100A10","PLTP","COL6A3","LUM","THY1","DCN","ANXA1","CD99","BST2","COL1A2","C1R","COL1A1","IL32","BCAM","CCND1","SLC25A4","MYH11","PGAM1","TUBB","SLC25A3","CD59","HSPA5","CTSD","CCT3","ILK","VDAC2")
annot=rowAnnotation(gene=anno_mark(at=match(highlight.genes, plot.gene), labels=highlight.genes))

ncells_sampled <- sum(CAF_HCC$Sub_Cluster == "POSTN+ CAF")
cells_used <- c(
  CAF_HCC@meta.data %>% filter(Sub_Cluster %in% c("CAF")) %>% slice_min(n = 36, order_by = DPT) %>% pull(CellID),
  CAF_HCC@meta.data %>% filter(Sub_Cluster %in% c("MYH11+ CAF")) %>% slice_max(n = ncells_sampled, order_by = DPT) %>% pull(CellID),
  CAF_HCC@meta.data %>% filter(Sub_Cluster %in% c("POSTN+ CAF")) %>% slice_max(n = ncells_sampled, order_by = DPT) %>% pull(CellID)
)
CAF_HCC_subset <- CAF_HCC[plot.gene, cells_used]
MYH11_cells <- CAF_HCC_subset@meta.data %>% filter(Sub_Cluster == "MYH11+ CAF") %>% pull(CellID)
CAF_HCC_subset@meta.data[MYH11_cells, "DPT"] <- -CAF_HCC_subset@meta.data[MYH11_cells, "DPT"]
CAF_HCC_subset$dpt_order = rank(CAF_HCC_subset$DPT)
CAF_HCC_subset$dpt_order_perc = CAF_HCC_subset$dpt_order / max(CAF_HCC_subset$dpt_order)
sce <- as.SingleCellExperiment(CAF_HCC_subset)
sce$meta.cluster = as.character(sce$Sub_Cluster)
sce$cluster.name = factor(as.character(sce$Sub_Cluster))
colSet <- list()
colSet[['dpt_order_perc']] = circlize::colorRamp2(c(0,0.5,1),c("darkblue","white","darkblue"))
colSet[['cluster.name']] = structure(CAF_HCC_color_panel, names=names(CAF_HCC_color_panel))
set.seed(1)
hm.obj = heatmap_sm(sce, assay.name="logcounts",
                    colSet=colSet,
                    out.prefix="./04.figures/01.CAF_HCC_DPT_heatmap2.pdf",
                    columns=c("dpt_order_perc"), 
                    ncell.downsample=500, n_smooth = 20,
                    use_raster=T, raster_quality=2,
                    columns.order="dpt_order_perc",
                    show_column_names=F, show_row_names=T,
                    row_names_side="right", 
                    row_gap = unit(2, "mm"), column_gap = unit(0, "mm"),
                    border=F, ann.bar.height=0.5, palette.name="RdYlBu",
                    z.hi = 3,
                    pdf.width=10.5, pdf.height=10, do.scale=F,
                    do.clustering.row=T, do.clustering.col=F,
                    dend.row=T, right_annotation=annot,
                    clustering.distance="cosine",
                    clustering.method="ward.D2",
                    row_km=5, row_km_repeats=100)

# >>scVelo (server)----
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
ldat <- ReadVelocity(file = "./Data/livermerged_2.loom")
CAF_HCC_vel <- as.Seurat(x = ldat)
CAF_HCC_vel$index <- stringr::str_split_fixed(row.names(CAF_HCC_vel@meta.data),":",2)[,2]
CAF_HCC_vel$barcode <- stringr::str_split_fixed(row.names(CAF_HCC_vel@meta.data),":",2)[,1]
index_used <- paste0(stringr::str_split_fixed(row.names(CAF_HCC@meta.data),"_|-",3)[,2],"x")
CAF_HCC_vel <- subset(CAF_HCC_vel, subset = )
CAF_HCC_vel[["RNA"]] <- CAF_HCC_vel[["spliced"]]
CAF_HCC_vel <- SCTransform(CAF_HCC_vel)
CAF_HCC_vel <- RunPCA(CAF_HCC_vel)
CAF_HCC_vel <- RunUMAP(CAF_HCC_vel, dims = 1:20)
CAF_HCC_vel <- FindNeighbors(CAF_HCC_vel, dims = 1:20)
CAF_HCC_vel <- FindClusters(CAF_HCC_vel)
DefaultAssay(CAF_HCC_vel) <- "RNA"
SaveH5Seurat(CAF_HCC_vel, filename = "mouseCAF_HCC_vel.h5Seurat")
Convert("mouseCAF_HCC_vel.h5Seurat", dest = "h5ad")

# >>SCENIC (server)----
library(dplyr)
library(SCENIC)
library(Seurat)
library(SCopeLoomR)
# Calculation
exprMat <- as.matrix(CAF_HCC@assays$RNA@counts)
setwd("/data2/lzy/oncofetal/Analysis/Results/SCENIC/CAF_HCC/")
dir.create("./int")
cellInfo <- CAF_HCC@meta.data
saveRDS(cellInfo, file="./int/cellInfo.Rds")
org <- "hgnc"
dbDir <- "/data2/lzy/SCENIC/cisTarget_databases"
myDatasetTitle <- "SCENIC analysis of CAF HCC"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
scenicOptions@inputDatasetInfo$cellInfo <- "./int/cellInfo.Rds"
saveRDS(scenicOptions, file="./int/scenicOptions.Rds") 
scenicOptions <- readRDS("./int/scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions)
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
library(SCENIC, lib.loc = "D:/Softwares/R-4.0.2/library2/")
library(AUCell)
setwd("./Results/SCENIC/CAF_HCC_tpm/")
scenicOptions <- readRDS("./int/scenicOptions.Rds")
cellInfo <- data.frame(CellType=Idents(CAF_HCC))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
cluster_used <- c("CAF","MYH11+ CAF","HSP+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF")
cells_used <- CAF_HCC@meta.data %>% filter(Sub_Cluster %in% cluster_used) %>% pull(CellID)

FindDEGenes(
  expression_matrix = as.matrix(regulonAUC@assays@data$AUC),
  groupid = CAF_HCC$Sub_Cluster,
  out.prefix = "./CAF_HCC_SCENIC/CAF_HCC",
  cutoff = 1,
  logFC = 0
)

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
topRegulators_used <- topRegulators %>% filter(CellType %in% cluster_used) %>% group_by(CellType) %>% arrange(desc(RelativeActivity), .by_group = T) %>% top_n(15, wt = RelativeActivity) %>% pull(Regulon) %>% unique() %>% as.character()
CAF_HCC_SCENIC_degenes <- read.csv("./CAF_HCC_SCENIC/CAF_HCC_cutoff1_all_de_genes.csv", row.names = 1)
regulons_used <- CAF_HCC_SCENIC_degenes %>% filter(Group %in% cluster_used) %>% group_by(Group) %>% filter(Exp.Per.Out < 0.8, Exp.Per.In > 0.2) %>% arrange(desc(AUC), .by_group = T) %>% top_n(20, wt = AUC) %>% pull(Symbol) %>% unique()
regulons_used <- unique(c(regulons_used, topRegulators_used))

regulons_mean_exp <- 
  aggregate(
    as.matrix(t(regulonAUC@assays@data$AUC[regulons_used,cells_used])),
    list(Cluster = CAF_HCC@meta.data[cells_used,"Sub_Cluster"]),
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
pdf("../../../04.figures/01.CAF_HCC_SCENIC_regulon_activity_heatmap.pdf", width = 5, height = 10)
draw(p)
dev.off()

# >>>>Regulon Specificity Score(RSS)----
library(scFunctions)
binary_regulon <- loadInt(scenicOptions, "aucell_binary_nonDupl")
binary_regulon <- binary_regulon[onlyNonDuplicatedExtended(rownames(binary_regulon)),]
rrs_df <- calculate_rrs(CAF_HCC@meta.data, binary_regulon, "Sub_Cluster")
rrs_df <- rrs_df %>% subset(!grepl("extended", regulon))

rrs_df_sub <- rrs_df %>% subset(cell_type == "POSTN+ CAF") %>% 
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
ggsave(rrs_ranking_plot, file = "../../../04.figures/01.CAF_HCC_SCENIC_POSTN_RSS.pdf", width = 3, height = 5)

rrs_df_sub <- rrs_df %>% subset(cell_type == "MYH11+ CAF") %>% 
  arrange(desc(RSS))
rrs_df_sub <- rrs_df_sub %>% mutate(rank = as.numeric(rownames(rrs_df_sub)))
rrs_ranking_plot2 <- 
  ggplot(rrs_df_sub, aes(rank, RSS, label = regulon)) + 
  geom_point(data = subset(rrs_df_sub,rank > 8),color = "#619CFF", size = 1) + 
  geom_point(data = subset(rrs_df_sub,rank <= 8), color = "#F87269", size = .8) + 
  labs(x = "Rank", y = "Regulon Specificity Score", title = "Top regulons of\nMYH11+ CAF") +
  theme_cowplot(font_size = 12, line_size = 1)
rrs_ranking_plot2 <- 
  ggAIplot(rrs_ranking_plot2, width = 3, height = 5) + 
  geom_text_repel(data = subset(rrs_df_sub, rank <= 8),
                  nudge_x = 20 - subset(rrs_df_sub, rank <= 8)$rank,
                  segment.size = 0.2,
                  segment.color = "#F87269", color = "#F87269",
                  direction = "y", hjust = 0)
ggsave(rrs_ranking_plot2, file = "../../../04.figures/01.CAF_HCC_SCENIC_MYH11_RSS.pdf", width = 3, height = 5)

# >>>>Regulon activity in UMAP----
regulons_used <- c("MSC (16g)", "STAT2 (388g)", "STAT1 (200g)","RUNX1 (13g)", "BHLHE40 (340g)", "FOXO3 (49g)")
plot_mat <- regulonAUC@assays@data$AUC[regulons_used,]
plot_mat <- t(apply(plot_mat, 1, zscore))
plot_mat_quantile <- quantile(as.matrix(plot_mat), c(0.02, 0.98))
plot_mat <- pmax(plot_mat, plot_mat_quantile[1])
plot_mat <- pmin(plot_mat, plot_mat_quantile[2])
plot_df <- data.frame(
  t(plot_mat), 
  UMAP1 = CAF_HCC$UMAP_1, UMAP2 = CAF_HCC$UMAP_2, 
  check.names = F
)
myColorPalette <- colorRampPalette(brewer.pal(9, "OrRd"))
p <- plot_df %>% melt(id.vars = c("UMAP1","UMAP2")) %>% arrange(desc(value)) %>%
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = value), size = 2.5) +
  facet_wrap(~variable, nrow = 2) +
  scale_colour_gradientn("Regulon\nActivity", colors = myColorPalette(100))+
  theme_cowplot(font_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p <- ggAIplot.grid(p, facet.by = "variable")
ggsave(p, file = "../../../04.figures/01.CAF_HCC_SCENIC_regulon_activity_UMAP.pdf", width = 10, height = 7)

genes_used <- c("MSC","STAT2","STAT1","RUNX1","BHLHE40","FOXO3")
CAF_HCC.cells.used <- CAF_HCC$CellID
CAF_HCC.plot.data <- cbind(t(as.matrix(CAF_HCC@assays$RNA@data[genes_used,CAF_HCC.cells.used])), CAF_HCC@meta.data[CAF_HCC.cells.used, c("UMAP_1","UMAP_2")])
CAF_HCC.plot.data <- melt(CAF_HCC.plot.data, id.vars = c("UMAP_1","UMAP_2"))
plot_mat_quantile <- quantile(as.matrix(CAF_HCC.plot.data$value), c(0.02, 0.98))
CAF_HCC.plot.data$value <- pmax(CAF_HCC.plot.data$value, plot_mat_quantile[1])
CAF_HCC.plot.data$value <- pmin(CAF_HCC.plot.data$value, plot_mat_quantile[2])
myColorPalette <- colorRampPalette(brewer.pal(9, "OrRd"))
p <- 
  ggplot(CAF_HCC.plot.data %>% arrange(value), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = value), alpha = 0.8, size = 3) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("log2\n(TPM+1)", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "../../../04.figures/01.CAF_HCC_SCENIC_TF_UMAP.pdf", width = 10, height = 7)

# COMET (server)----
cells_used <- CAF_HCC@meta.data %>% filter(NTF == "Tumor") %>% pull(CellID)
write.table(CAF_HCC@assays$RNA@data[,cells_used], file = "./Data/CAF_HCC_COMET_exp.txt", quote = FALSE, sep = "\t")

vis_data <- CAF_HCC@meta.data[cells_used,c("CellID","UMAP_1","UMAP_2")]
vis_data$CellID <- gsub("-",".",vis_data$CellID)
write.table(vis_data, file = "./Data/CAF_HCC_COMET_vis.txt", quote = FALSE, sep = "\t", col.names = F, row.names = F)

cluster_data <- CAF_HCC@meta.data[cells_used,c("CellID","Sub_Cluster")]
cluster_data$CellID <- gsub("-",".",cluster_data$CellID)
cluster_data$Sub_Cluster <- as.numeric(factor(cluster_data$Sub_Cluster))
write.table(cluster_data, file = "./Data/CAF_HCC_COMET_cluster.txt", quote = FALSE, sep = "\t", col.names = F, row.names = F)

# Projection CAF clusters to pan-CAF study----
c3_genes <- intersect(c("CTHRC1","COL1A1","HSPA6","FN1","COL3A1","POSTN","COL1A2","GJB2","COL11A1","HSPA1B","COL8A1","CRYAB","HLA-B","CST1","CTSK"),row.names(CAF_HCC))
c8_genes <- intersect(c("C7","CFD","PRSS1","GSN","ADH1B","CLPS","APOD","PRSS2","PNLIP","CCL2","IGFBP5","COLEC11","A2M","CELA3B","MGP","CTRB2","SEPP1","CPB1","CPA1","RPS4Y1","FMO2"),row.names(CAF_HCC))
CAF_HCC@meta.data$c3_signature <- colMeans(CAF_HCC@assays$RNA@data[c3_genes[1:10],])
CAF_HCC@meta.data$c8_signature <- colMeans(CAF_HCC@assays$RNA@data[c8_genes[1:10],])
plot.matrix <- CAF_HCC@meta.data[,c("UMAP_1","UMAP_2","c3_signature","c8_signature")]
plot.matrix$c3_signature <- pmin(plot.matrix$c3_signature, quantile(plot.matrix$c3_signature, .99))
plot.matrix$c3_signature <- pmax(plot.matrix$c3_signature, quantile(plot.matrix$c3_signature, .01))
plot.matrix$c8_signature <- pmin(plot.matrix$c8_signature, quantile(plot.matrix$c3_signature, .99))
plot.matrix$c8_signature <- pmax(plot.matrix$c8_signature, quantile(plot.matrix$c3_signature, .01))
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p1 <- ggplot(plot.matrix %>% arrange(c3_signature), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = c3_signature), size = 3) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("CAF (c3) score") +
  scale_colour_gradientn("Exp", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = .2, barheight = 12))
p1.pdf <- ggAIplot(p1)
p2 <- ggplot(plot.matrix %>% arrange(c8_signature), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = c8_signature), size = 3) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("NAT (c8) score") +
  scale_colour_gradientn("Exp", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = .2, barheight = 12))
p2.pdf <- ggAIplot(p2)
p <- p1.pdf + p2.pdf
ggsave(p, file = "./04.figures/01.CAF_HCC_projection_of_c3&8_signature.pdf", width = 6.5, height = 3)

# Save data----
saveRDS(CAF_HCC, file = "./Data/CAF_HCC.rds")
