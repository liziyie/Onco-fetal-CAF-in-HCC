setwd("/work/lzy/project/onco_fetal/")
source("../utils/utils_plot.R")
source("../utils/utils_data_processing.R")
source("../utils/utils_color.R")
library(Seurat)
library(Polychrome)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)

# >Load dataset----
CAF_HCC <- readRDS("./02.processed_data/CAF_HCC.rds")
CAF_HCC$Sub_Cluster <- as.character(CAF_HCC$Sub_Cluster)
CAF_HCC$Sub_Cluster[CAF_HCC$Sub_Cluster == "Fib"] <- "Fib(Tumor)"
CAF_HCC$Sub_Cluster <- factor(CAF_HCC$Sub_Cluster, levels = c( "CAF","HSP+ CAF","MYH11+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF","Fib(Tumor)","MT1M+ Fib","ABCAB+ Fib"))
CAF_HCC_color_panel <- c(
  "CAF" = "#62A3C7", "Fib(Tumor)" = "#68AB9F", "SDC2+ CAF" = "#CAB2D6",
  "MYH11+ CAF" = "#F69459", "HSP+ CAF" = "#1C76B3", "POSTN+ CAF" = "#693D99", 
  "APOA2+ CAF" = "#FB9A99", "MT1M+ Fib" = "#B2DE89", "ABCAB+ Fib" = "#4DAC47"
)

Fib_fetal <- readRDS("./02.processed_data/Fib_fetal.rds")
Fib_fetal$Sub_Cluster <- as.character(Fib_fetal$Sub_Cluster)
Fib_fetal$Sub_Cluster[Fib_fetal$Sub_Cluster == "Fib"] <- "Fib(Fetal)"
Fib_fetal$Sub_Cluster[Fib_fetal$Sub_Cluster == "Mesenchymal"] <- "Mesenchymal(Fetal)"
Fib_fetal$Sub_Cluster <- factor(Fib_fetal$Sub_Cluster, levels = c("Fib(Fetal)","ACTA2+ Fib","COL1A1+ Fib","BGN+ Fib","APOC3+ Fib","ITM2C+ Fib","Mesenchymal(Fetal)","PLEK+ Fib"))
Fib_fetal_color_panel <- c(
  "Fib(Fetal)" = "#F69459", "ACTA2+ Fib" = "#84B8D7", "COL1A1+ Fib" = "#8DC594",
  "BGN+ Fib" = "#68AB9F", "APOC3+ Fib" = "#B294C7", "ITM2C+ Fib" = "#815AA8",
  "Mesenchymal(Fetal)" = "#E93A3B", "PLEK+ Fib" = "#F58584"
)

Healthy_color_panel <- c("VSMC" = "#F69459", "HSC" = "#62A3C7", "SAMes" = "#815AA8", "Mesenchymal(Healthy)" = "#E93A3B")

Fib_integrated <- readRDS("./02.processed_data/CAF_integration_onco-fetal.rds")
Fib_integrated@meta.data$Sub_Cluster <- factor(Fib_integrated$Sub_Cluster, levels = c(names(CAF_HCC_color_panel), names(Fib_fetal_color_panel),"VSMC","HSC","SAMes"))

# >Data integration----
# >>UMAP plot----
integrated_color_panel <- c(
  "#66CEF6",
  "#EAA944","#A0D7C9","#D6D4EB","#74517B","#F8F4A8",
  "#B9A96B","#EF5276","#EACF68","#F3746C","#7A8DBB",
  "#69B4CE")
names(integrated_color_panel) <- 0:11
p <- ggplot(Fib_integrated@meta.data, aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 2, alpha = .8) +
  theme_cowplot(font_size = 7) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_color_manual(name = "", values = integrated_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/03.Integrated_Cluster_UMAP.pdf", width = 5, height = 4.5)

p <- ggplot(Fib_integrated@meta.data, aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Source), size = 2, alpha = .8) +
  theme_cowplot(font_size = 7) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_color_manual(name = "", values = Tissue_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/03.Integrated_Source_UMAP.pdf", width = 5, height = 4.5)

p1 <- ggplot(Fib_integrated@meta.data %>% filter(Source == "Adj Normal"), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 3, alpha = .8) +
  theme_cowplot(font_size = 7) +
  lims(x = c(-7.5,6), y = c(-6,5.5)) + 
  scale_color_manual(name = "", values = integrated_color_panel) +
  theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")
p1.pdf <- ggAIplot(p1)
p2 <- ggplot(Fib_integrated@meta.data %>% filter(Source == "Fetal"), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 3, alpha = .8) +
  theme_cowplot(font_size = 7) +
  lims(x = c(-7.5,6), y = c(-6,5.5)) +
  scale_color_manual(name = "", values = integrated_color_panel) +
  theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")
p2.pdf <- ggAIplot(p2)
p3 <- ggplot(Fib_integrated@meta.data %>% filter(Source == "Healthy"), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 3, alpha = .8) +
  theme_cowplot(font_size = 7) +
  lims(x = c(-7.5,6), y = c(-6,5.5)) +
  scale_color_manual(name = "", values = integrated_color_panel) +
  theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")
p3.pdf <- ggAIplot(p3)
p4 <- ggplot(Fib_integrated@meta.data %>% filter(Source == "Tumor"), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Integrated_cluster), size = 3, alpha = .8) +
  theme_cowplot(font_size = 7) +
  lims(x = c(-7.5,6), y = c(-6,5.5)) +
  scale_color_manual(name = "", values = integrated_color_panel) +
  theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")
p4.pdf <- ggAIplot(p4)
p <- plot_grid(p1.pdf, p2.pdf, p3.pdf, p4.pdf, nrow = 2)
ggsave(p, file = "./04.figures/03.Integrated_UMAP_Source.pdf", width = 8, height = 8)

# >>Cell Proportions----
p <- 
  ggplot(Fib_integrated@meta.data, aes(x = Integrated_cluster, fill = Source)) +
  geom_bar(position = "fill", color = "white") +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  scale_fill_manual(values = Tissue_color_panel) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 6), ncol = 1)) +
  scale_x_discrete(labels=paste0(
    names(table(Fib_integrated$Integrated_cluster))," (",
    table(Fib_integrated$Integrated_cluster),")")) +
  coord_flip()
ggsave(p, file = "./04.figures/03.Integrated_Source_proportions.pdf", width = 5, height = 2.5)

p <- 
  ggplot(Fib_integrated@meta.data %>% filter(Integrated_cluster %in% c(3,4,7), Source != "Healthy"), aes(x = Integrated_cluster, fill = Sub_Cluster)) +
  geom_bar(position = "fill", color = NA) +
  labs(x = "", y = "Percentage") + theme_cowplot(font_size = 7) + 
  scale_fill_manual(values = c(CAF_HCC_color_panel, Fib_fetal_color_panel)) +
  theme(legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 1))
ggsave(p, file = "./04.figures/03.Integrated_Sub_Cluster_proportions.pdf", width = 2, height = 3)

per.df <- Fib_integrated@meta.data %>% filter(Integrated_cluster == 4) %>% dplyr::select(Source) %>% table() %>% data.frame()
colnames(per.df)[1] <- "Source"
per.df$Freq <- per.df$Freq / sum(per.df$Freq) * 100
p <- 
  ggplot(per.df, aes(x = Source, y = Freq, fill = Source)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=round(Freq,2)), vjust=0, color="black", size=5)+
  scale_x_discrete(limits=c("Healthy","Tumor","Fetal")) +
  xlab("") + ylab("Frequency") +
  scale_fill_manual(name = "", values=Tissue_color_panel)+
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(p, file = "./04.figures/03.Integrated_Cluster4_Source_proportions.pdf", width = 3, height = 6)

# >>Highlight fibroblast clusters----
plot.data <- Fib_integrated@meta.data
plot.data$Sub_Cluster <- as.character(plot.data$Sub_Cluster)
plot.data$Sub_Cluster[!(plot.data$Sub_Cluster %in% c("POSTN+ CAF","ITM2C+ Fib","PLEK+ Fib","Fib(Tumor)","Fib(Fetal)"))] <- "A0"
p <- ggplot(plot.data %>% arrange(Sub_Cluster), aes(x = integrated_UMAP_1, y = integrated_UMAP_2)) +
  geom_point(aes(color = Sub_Cluster), size = 2, alpha = .6) +
  theme_cowplot(font_size = 7) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_color_manual(name = "", values = c(CAF_HCC_color_panel, Fib_fetal_color_panel, Healthy_color_panel, "A0" = "lightgrey")) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/03.Integrated_Cluster_Highlight_UMAP.pdf", width = 5.5, height = 4.5)

library(monocle3)
library(SeuratWrappers)
library(magrittr)
Fib_integrated_UMAP.cds <- as.cell_data_set(Fib_integrated, assay = "RNA")
Fib_integrated_UMAP.cds <- cluster_cells(Fib_integrated_UMAP.cds)

Fib_integrated_UMAP.cds <- learn_graph(Fib_integrated_UMAP.cds, learn_graph_control = list(minimal_branch_len = 5))
root_cells <- data.frame(colData(Fib_integrated_UMAP.cds)) %>% filter(Sub_Cluster == "HSC") %>% row.names()
Fib_integrated_UMAP.cds <- order_cells(Fib_integrated_UMAP.cds, root_cells = root_cells)
save(Fib_integrated_UMAP.cds, file = "./03.results/Trajectory/Fib_integrated_Monocle3.rda")

p1 <- plot_cells(
  Fib_integrated_UMAP.cds, 
  color_cells_by = "Sub_Cluster", 
  cell_size = 0,
  trajectory_graph_segment_size = 1,
  label_groups_by_cluster = FALSE, 
  label_cell_groups = FALSE,
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  labs(x = "UMAP1", y = "UMAP2")
p1 <- plot_grid(p1 + theme(legend.position = "none"))
ggsave(p1, file = "./04.figures/03.Fib_integrated_monocle3_UMAP_trajectory.pdf", width = 4, height = 4)

# >>Similarity----
library(scibetR)
CAF_HCC_cells_used <- CAF_HCC@meta.data %>% pull(CellID)
Fib_fetal_cells_used <- Fib_fetal@meta.data %>% pull(CellID)
scibet_test_set <- data.frame(
  t(CAF_HCC@assays$RNA@counts[,CAF_HCC_cells_used]), 
  label = CAF_HCC@meta.data[CAF_HCC_cells_used,"Sub_Cluster"]
)
scibet_train_set <- data.frame(
  t(Fib_fetal@assays$RNA@counts[,Fib_fetal_cells_used]), 
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
test_order <- c("CAF","APOA2+ CAF","HSP+ CAF","MYH11+ CAF","SDC2+ CAF","POSTN+ CAF","Fib(Tumor)","MT1M+ Fib","ABCAB+ Fib") 
train_order <- c("Fib(Fetal)","ACTA2+ Fib","COL1A1+ Fib","BGN+ Fib","APOC3+ Fib","ITM2C+ Fib","PLEK+ Fib","Mesenchymal")

color_used <- 
  colorRamp2(c(seq(0,0.5,length.out = 50),0.7), 
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
pdf("./Figures/03.HCC_to_fetal_similarity.pdf", width = 6, height = 4)
p
dev.off()

# >>Combined heatmap----
# dendrogram
library(dendextend)
library(plyr)
cluster <- c(
  paste0("Tumor_",as.character(CAF_HCC@meta.data[CAF_HCC_cells_used,"Sub_Cluster"])), 
  paste0("Fetal_",as.character(Fib_fetal@meta.data[Fib_fetal_cells_used,"Sub_Cluster"]))
)
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

Fib_fetal.markers <- read.csv("./Results/DEGenes/Fib_fetal_all_markers.csv", row.names = 1)
Fib_fetal_genes <- Fib_fetal.markers %>% filter(!grepl("^RP|^MT",gene)) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) %>% filter(p_val < 0.05, avg_logFC > 0)
CAF_HCC.markers <- read.csv("./Results/DEGenes/CAF_HCC_all_markers.csv", row.names = 1)
CAF_HCC_genes <- CAF_HCC.markers %>% filter(!grepl("^RP|^MT",gene)) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) %>% filter(p_val < 0.05, avg_logFC > 0)
features_used <- intersect(CAF_HCC_genes$gene, Fib_fetal_genes$gene)
features_used <- unique(c(features_used, plot.df %>% filter(Sig == "Both high up-regulated") %>% pull(Symbol))) #high_up_gene_name defined in next panel
CAF_HCC_exp <- aggregate(t(CAF_HCC@assays$RNA@scale.data[features_used,]), list(Cluster = CAF_HCC$Sub_Cluster), mean)
row.names(CAF_HCC_exp) <- CAF_HCC_exp$Cluster
CAF_HCC_exp$Cluster <- c()
CAF_HCC_exp <- t(CAF_HCC_exp)
colnames(CAF_HCC_exp) <- paste0("Tumor_",colnames(CAF_HCC_exp))
Fib_fetal_exp <- aggregate(t(Fib_fetal@assays$RNA@scale.data[features_used,]), list(Cluster = Fib_fetal$Sub_Cluster), mean)
row.names(Fib_fetal_exp) <- Fib_fetal_exp$Cluster
Fib_fetal_exp$Cluster <- c()
Fib_fetal_exp <- t(Fib_fetal_exp)
colnames(Fib_fetal_exp) <- paste0("Fetal_",colnames(Fib_fetal_exp))
CAF_all_exp_agg <- cbind(Fib_fetal_exp, CAF_HCC_exp)
CAF_all_exp_agg <- apply(CAF_all_exp_agg, 2, zscore)
CAF_all_exp_agg_quantile <- quantile(CAF_all_exp_agg,c(0.01,0.98))
CAF_all_exp_agg <- pmax(CAF_all_exp_agg, CAF_all_exp_agg_quantile[1])
CAF_all_exp_agg <- pmin(CAF_all_exp_agg, CAF_all_exp_agg_quantile[2])
CAF_all_exp_agg <- CAF_all_exp_agg[,hc$labels[hc$order]]
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
genes_labeled <- c("COL3A1","PRELP","ICAM1","BGN","GPC3","ITM2A","CCNL1","SOCS3","HLA-E","ACTB","B2M","FABP1","ALB","APOA1","CCL2","CXCL14","CCPG1","CCL21","FCGRT","ECM1","ITM2C","COL5A1","THY1","TSPAN4","POSTN","VCAN","HLA-B","CD81","IFITM3","FN1","COL4A1","ACTA2","POLR2L","HSPA6","BAG3","S100A4","S100A6","PLAC9","LY6E")
p <- Heatmap(CAF_all_exp_agg,
             cluster_rows = F,col = color_used,
             column_dend_side = "bottom",
             row_names_gp = gpar(fontsize = 6),
             cluster_columns = dend,
             heatmap_legend_param = list(
               at = c(-2, -1, 0, 1, 2, 2.3),
               labels = c("<-2", "-1", "0", "1", "2", ">2.3"),
               title = "Normalized expression",
               legend_height = unit(4, "cm"),
               title_position = "lefttop-rot")) +
  rowAnnotation(link = anno_mark(
    at = which(gene_order$gene %in% genes_labeled),
    labels = gene_order$gene[which(gene_order$gene %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
pdf("./Figures/03.Integrated_heatmap.pdf", width = 5, height = 12)
draw(p)
dev.off()

# >>Upregulated genes----
# Differentially expressed genes
Fib_markers <- FindMarkers(Fib_fetal %>% subset(subset = Sub_Cluster != "Mesenchymal"), ident.1 = "ITM2C+ Fib", group.by = "Sub_Cluster", logfc.threshold = 0, min.pct = 0)
Fib_markers$gene <- row.names(Fib_markers)
CAF_markers <- FindMarkers(CAF_HCC %>% subset(subset = NTF == "Tumor"), ident.1 = "POSTN+ CAF", group.by = "Sub_Cluster", logfc.threshold = 0, min.pct = 0)
CAF_markers$gene <- row.names(CAF_markers)
genes_used <- intersect(Fib_markers$gene, CAF_markers$gene)
plot.df <- cbind(
  Fib_markers[genes_used,c("avg_logFC","p_val_adj","pct.1","pct.2")],
  CAF_markers[genes_used,c("avg_logFC","p_val_adj","pct.1","pct.2","gene")]
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
p <- ggplot(plot.df %>% arrange(desc(Sig)), aes(x = logFC_fetal, y = logFC_tumor)) +
  geom_point(aes(color = Sig), size = 2.5, alpha = .8) +
  labs(x = "logFC of ITM2C+ Fib", y = "logFC of POSTN+ CAF") +
  scale_color_manual(name = "", 
                     values = c("Not up-regulated" = "lightgrey",
                                "Tumor-specific up-regulated" = "#58B44D",
                                "Both low up-regulated" = "orange",
                                "Both high up-regulated" = "#E51B1A",
                                "Down-regulated" = "#4480BF")) +
  xlim(-2,2) + ylim(-2,2)
p <- ggAIplot(p + theme(legend.position = "none")) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = c(0,0.5), linetype = "dashed") + 
  geom_vline(xintercept = c(0,0.5), linetype = "dashed") +
  geom_text_repel(aes(label = Label)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 3)), keyheight = 1.5)
ggsave(p, file = "./Figures/03.Shared_up_regulated_genes.pdf", width = 6, height = 6.5)

write.csv(plot.df, file = "./Results/DEGenes/Co_up_reg_oncofetal_CAF_genes.csv")

# Pathway enrichment
library(richR)
hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
high_up_gene_name <- plot.df %>% filter(p_fetal < 0.05 & p_tumor < 0.05 & (logFC_tumor > 0.05 | logFC_fetal > 0.05)) %>% pull(Symbol)
high_up_go <- richGO(high_up_gene_name,godata = hsago,ontology ="BP")
high_up_go@result %>% filter(Padj < 0.05) %>% pull(Term)
save(high_up_go, file = "./Results/DEGenes/Co_up_reg_pathway.rda")
GO_to_show <- c(1,5,7,12,17,40,44,61,65,70,119,138,140,152,166,180,181)
pathway_plot.df <- data.frame(
    pathway = high_up_go@result$Term[GO_to_show],
    pvalue = high_up_go@result$Padj[GO_to_show],
    generatio = (high_up_go@result$Significant/high_up_go@result$Annotated)[GO_to_show],
    database = "GO (BP)"
)
p <- ggplot(pathway_plot.df, aes(x = reorder(pathway,-pvalue), y = -log10(pvalue))) +
  geom_bar(aes(fill = generatio), stat = "identity") +
  scale_fill_gradientn(name = "Gene Ratio", limits = c(0,0.25), colors = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)) +
  theme_bw(base_size = 12) +
  labs(y = "-log10 (P-Value)", x = "") +
  coord_flip() +
  theme(strip.text.x = element_blank(),strip.background = element_rect(fill=NA, color=NA))
ggsave(p, file = "./Figures/03.High_up_reg_genes_pathway_enrichment.pdf", width = 8, height = 4)

tumor_up_gene_name <- plot.df %>% filter(p_fetal > 0.05 & p_tumor < 0.05 & logFC_tumor > 0.5 & abs(logFC_fetal) < 0.5) %>% pull(Symbol)
tumor_up_go <- richGO(tumor_up_gene_name,godata = hsago,ontology ="BP")
tumor_up_go@result %>% filter(Padj < 0.05) %>% pull(Term)
save(tumor_up_go, file = "./Results/DEGenes/Tumor_up_reg_pathway.rda")
GO_to_show <- c(1,3,7,12,15,20,24,39,62,78,109,120,124,127,129,171)
pathway_plot.df <- data.frame(
  pathway = tumor_up_go@result$Term[GO_to_show],
  pvalue = tumor_up_go@result$Padj[GO_to_show],
  generatio = (tumor_up_go@result$Significant/tumor_up_go@result$Annotated)[GO_to_show],
  database = "GO (BP)"
)
p <- ggplot(pathway_plot.df, aes(x = reorder(pathway,-pvalue), y = -log10(pvalue))) +
  geom_bar(aes(fill = generatio), stat = "identity") +
  scale_fill_gradientn(name = "Gene Ratio", limits = c(0,0.11), colors = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)) +
  theme_bw(base_size = 12) +
  labs(y = "-log10 (P-Value)", x = "") +
  coord_flip() +
  theme(strip.text.x = element_blank(),strip.background = element_rect(fill=NA, color=NA))
ggsave(p, file = "./Figures/03.tumor_up_reg_genes_pathway_enrichment.pdf", width = 8, height = 4)
