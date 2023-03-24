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
Tissue_color_panel <- c("Tumor" = "#ED780F", "Normal" = "#2CAFB9", "Fetal" = "#B289BD")
Sample_color_panel <- c("B015" = "#56B4E9", "B017" = "#E69F00", "14w" = "#F0CE39", "16w" = "#A16FAD")
Tissue_box_color_panel <- c("Tumor" = "#F84141","Normal" = "#9FA0A3")
Cluster_levels <- c("B","T cells","Mast","DC","pDC","SPP1+ TAM2","MT1G+ TAM3","Other Mononuclear","Other Endothelium","Other CAF","FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Hepatocyte","Hepatocyte_P7_2","Hepatocyte_P15")
Cluster_color_panel <- c("B"="#F8BFAF","T cells"="#63B472","Mast"="#DEEAB1","DC"="#588198","pDC"="#D5E7F7",
                         "SPP1+ TAM2"="#D6D4EB","MT1G+ TAM3"="#EACF68",
                         "Other Mononuclear"="#F8F4A8","Other Endothelium"="#C9BDB2","Other CAF"="#7A8DBB",
                         "FOLR2+ TAM1"="#F3746C","PLVAP+ EC"="#69B4CE","POSTN+ CAF"="#EAA944",
                         "Hepatocyte"="#A0D7C9","Hepatocyte_P7_2"="#74517B","Hepatocyte_P15"="#EF5276")

# >>Nanostring data----
HCC_seu <- readRDS("./02.processed_data/HCC_seu.rds")
# >Load nanostring wta matrix----
load("./02.processed_data/Nanostring/nanostring.rda")
# ns_exp <- read.csv("./02.processed_data/Nanostring/nanostring_wta.csv", row.names = 1)
# ns_exp <- log2(ns_exp + 1)
# ns_sample <- read.csv("./02.processed_data/Nanostring/nanostring_sample.csv", row.names = 1, stringsAsFactors = F)
# ns_tumor_exp <- ns_exp[,ns_sample$Tissue %in% c("Tumor", "Normal")]
# ns_fetal_exp <- ns_exp[,ns_sample$Tissue == "Fetal"]
# ns_tumor_sample <- ns_sample[colnames(ns_tumor_exp),]
# ns_fetal_sample <- ns_sample[colnames(ns_fetal_exp),]
# ns_tumor_exp_norm <- ComBat(as.matrix(ns_tumor_exp), as.numeric(as.factor(ns_tumor_sample$Sample)), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
# ns_fetal_exp_norm <- ComBat(as.matrix(ns_fetal_exp), as.numeric(as.factor(ns_fetal_sample$Sample)), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
# ns_exp_norm <- ComBat(as.matrix(ns_exp), as.numeric(as.factor(ns_sample$Sample)), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
# ns_tumor_sample$Group <- factor(paste0(ns_tumor_sample$Sample,"_",ns_tumor_sample$Tissue),levels = c("B015_Normal","B017_Normal","B015_Tumor","B017_Tumor"))
# save(ns_exp_norm, ns_tumor_exp_norm, ns_fetal_exp_norm, ns_sample, ns_tumor_sample, ns_fetal_sample, file = "./02.processed_data/Nanostring/nanostring.rda")

# >Load CIBERSORTx deconvolution results----
tumor_prop <- read.table("./02.processed_data/Nanostring/seperated_cluster/Tumor/CIBERSORTx_Results.txt", sep = "\t", head = T, check.names = F, row.names = 1)
row.names(tumor_prop) <- ns_tumor_sample[ns_tumor_sample$Tissue == "Tumor", "NewID"]
tumor_prop2 <- tumor_prop[,c("CD4+ T","CD8+ T")]
for(cluster in unique(HCC_seu@meta.data$Sub_ClusterNew)){
  sub_cluster <- intersect(colnames(tumor_prop), c(unique(HCC_seu@meta.data %>% filter(Sub_ClusterNew == cluster) %>% pull(Sub_Cluster)), cluster))
  if(length(sub_cluster) == 1){
    tumor_prop2[,cluster] <- tumor_prop[,sub_cluster]
  }else{
    tumor_prop2[,cluster] <- rowSums(tumor_prop[,sub_cluster])
  }
}
tumor_prop2[,"T cells"] <- tumor_prop2[,"CD4+ T"] + tumor_prop2[,"CD8+ T"] + tumor_prop2[,"Tregs"]
tumor_prop2$`CD4+ T` <- tumor_prop2$`CD8+ T` <- tumor_prop2$Tregs <- c()
tumor_prop2[,"All Hepatocyte"] <- tumor_prop2[,"Hepatocyte_P15"] + tumor_prop2[,"Hepatocyte_P7_2"] + tumor_prop2[,"Hepatocyte"]
tumor_prop2 <- tumor_prop2[,Cluster_levels]

# >Analyze tumor and normal samples first----
# PCA analysis----
ns_tumor_sd <- sort(apply(ns_tumor_exp_norm,1,sd),decreasing = T)
features_used <- names(ns_tumor_sd)[1:1000]
pca_analysis <- prcomp(t(ns_tumor_exp_norm[features_used,]))
pca_analysis <- summary(pca_analysis)
pca_plot <- data.frame(
  pca_analysis$x[,1:2],
  Tissue = ns_tumor_sample$Tissue,
  Sample = ns_tumor_sample$Sample)
p <- ggplot(pca_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Tissue, shape = Sample), size = 3, alpha = .8) +
  scale_color_manual(values = Tissue_color_panel) +
  xlab(sprintf("PC1: %s%% variance",round(pca_analysis$importance[2,1],2) * 100)) + 
  ylab(sprintf("PC2: %s%% variance",round(pca_analysis$importance[2,2],2) * 100)) +
  theme_cowplot()
ggsave(p, file = "./04.figures/07.ns_tumor_pca.pdf", width = 6, height = 5)

# Compare onco-fetal cluster signatures----
OF_TAM_signature <- c("FOLR2","CD163","MRC1","TREM2")
OF_Endo_signature <- c("PLVAP","ACKR1","DLL4")
OF_CAF_signature <- c("POSTN","FAP","ENG","PTGDS")
ns_tumor_sample$OF_TAM_score <- colMeans(ns_tumor_exp_norm[OF_TAM_signature,])
ns_tumor_sample$OF_Endo_score <- colMeans(ns_tumor_exp_norm[OF_Endo_signature,])
ns_tumor_sample$OF_CAF_score <- colMeans(ns_tumor_exp_norm[OF_CAF_signature,])
p <- ns_tumor_sample %>% select(Tissue, Group, OF_TAM_score, OF_Endo_score, OF_CAF_score) %>% reshape2::melt() %>%
  ggplot(aes(x = Group, y = value)) +
  geom_boxplot(aes(fill = Tissue), alpha = .8) +
  scale_fill_manual(values = Tissue_box_color_panel) +
  facet_wrap(~variable, scales = "free") +
  labs(y = "Expression", x = "") +
  theme_cowplot() +
  ggsignif::geom_signif(comparisons = list(c("B015_Normal","B017_Normal"),
                                           c("B015_Tumor","B017_Tumor"),
                                           c("B015_Normal","B015_Tumor"),
                                           c("B017_Normal","B017_Tumor")),
                        step_increase = .1, tip_length = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(p, file = "./04.figures/07.ns_patient_score_comparison.pdf", width = 6, height = 5)

p <- data.frame(Sample = ns_tumor_sample[row.names(tumor_prop),"Sample"],
                tumor_prop2[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")],
                check.names = F) %>% reshape2::melt() %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_boxplot(aes(color = Sample), alpha = .8) +
  scale_color_manual(values = Sample_color_panel) +
  facet_wrap(~variable, scales = "free") +
  labs(y = "Predicted proportions", x = "") +
  theme_cowplot() +
  ggsignif::geom_signif(comparisons = list(c("B015","B017")),
                        step_increase = .1, tip_length = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(p, file = "./04.figures/07.ns_patient_proportion_comparison.pdf", width = 6, height = 5)

p <- data.frame(Sample = ns_tumor_sample[row.names(tumor_prop),"Sample"],
                OF_score = rowSums(tumor_prop2[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")]),
                check.names = F) %>%
  ggplot(aes(x = Sample, y = OF_score)) +
  geom_boxplot(aes(color = Sample), alpha = .8) +
  scale_color_manual(values = Sample_color_panel) +
  labs(y = "Onco-fetal Score", x = "") +
  theme_cowplot() +
  ggsignif::geom_signif(comparisons = list(c("B015","B017")),
                        step_increase = .1, tip_length = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(p, file = "./04.figures/07.ns_patient_OF_score_comparison.pdf", width = 3, height = 5)

# Correlation between onco-fetal clusters----
ns_tumor_sample_temp <- ns_tumor_sample %>% filter(Tissue == "Tumor")
ns_tumor_sample_temp <- cbind(ns_tumor_sample_temp, tumor_prop2[row.names(ns_tumor_sample_temp),])
EC_prop_cutoff <- quantile(ns_tumor_sample_temp$`PLVAP+ EC`,c(0.5,0.5))
ns_tumor_sample_temp[ns_tumor_sample_temp$`PLVAP+ EC` > EC_prop_cutoff[2],"Group2"] <- "PLVAP+ EC high"
ns_tumor_sample_temp[ns_tumor_sample_temp$`PLVAP+ EC` <= EC_prop_cutoff[1],"Group2"] <- "PLVAP+ EC low"
ns_tumor_sample_temp$Group2 <- factor(ns_tumor_sample_temp$Group2, levels = c("PLVAP+ EC low","PLVAP+ EC high",NA))
TAM_prop_cutoff <- quantile(ns_tumor_sample_temp$`FOLR2+ TAM1`,c(0.5,0.5))
ns_tumor_sample_temp[ns_tumor_sample_temp$`FOLR2+ TAM1` > TAM_prop_cutoff[2],"Group3"] <- "FOLR2+ TAM1 high"
ns_tumor_sample_temp[ns_tumor_sample_temp$`FOLR2+ TAM1` <= TAM_prop_cutoff[1],"Group3"] <- "FOLR2+ TAM1 low"
ns_tumor_sample_temp$Group3 <- factor(ns_tumor_sample_temp$Group3, levels = c("FOLR2+ TAM1 low","FOLR2+ TAM1 high",NA))
p1 <- ggplot(ns_tumor_sample_temp %>% filter(!is.na(Group2)), aes(x = Group2, y = `POSTN+ CAF`)) +
  geom_boxplot(aes(fill = Group2), alpha = .8) +
  theme_cowplot() +
  labs(y = "Predicted proportions of\n POSTN+ CAF", x = "") +
  scale_fill_manual(values = c("PLVAP+ EC high" = "#F84141", "PLVAP+ EC low" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("PLVAP+ EC high","PLVAP+ EC low"))) +
  RotatedAxis() +
  theme(legend.position = "none")
p2 <- ggplot(ns_tumor_sample_temp %>% filter(!is.na(Group3)), aes(x = Group3, y = `POSTN+ CAF`)) +
  geom_boxplot(aes(fill = Group3), alpha = .8) +
  theme_cowplot() +
  labs(y = "Predicted proportions of\n POSTN+ CAF", x = "") +
  scale_fill_manual(values = c("FOLR2+ TAM1 high" = "#F84141", "FOLR2+ TAM1 low" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("FOLR2+ TAM1 high","FOLR2+ TAM1 low"))) +
  RotatedAxis() +
  theme(legend.position = "none")
p <- p1 + p2
ggsave(p, file = "./04.figures/07.ns_tumor_POSTN+CAF_proportion_ECorTAM1highlow_boxplot.pdf", width = 4, height = 5)

# Cluster tumor samples according to proportions----
hc <- hclust(dist(tumor_prop2, method = "euclidean"))
dend.hc <- as.dendrogram(hc)
pdf("./04.figures/07.ns_tumor_dendogram.pdf", width = 6, height = 3)
dend.hc %>% set("labels_col", Sample_color_panel[ns_tumor_sample[labels(dend.hc),"Sample"]]) %>% plot(center = T)
dend.hc %>% rect.dendrogram(k=3, border=8, lty=5, lwd=2) 
dev.off()

plot.df <- tumor_prop2[,colSums(tumor_prop2)!=0] %>% mutate(Sample = row.names(.)) %>% melt() %>% mutate(Sample = factor(Sample, levels = labels(dend.hc)), variable = factor(variable, levels = Cluster_levels))
p <- ggplot(plot.df, aes(x = Sample, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity") +
  labs(y = "Predicted proportions", x = "") +
  scale_fill_manual(name = "", values = Cluster_color_panel) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(p, file = "./04.figures/07.ns_tumor_proportions_barplot.pdf", width = 10, height = 6)

# DEGenes heatmap----
ns_tumor_sample$Cluster <- c()
ns_tumor_sample[names(cutree(hc,3)),"Cluster"] <- cutree(hc,3)
sample_used <- ns_tumor_sample %>% filter(Cluster %in% c(1,2)) %>% pull(NewID)
degenes_12 <- LIMMA(ns_tumor_exp_norm[,sample_used],
                    ns_tumor_sample[sample_used,"Cluster"])
sample_used <- ns_tumor_sample %>% filter(Cluster %in% c(1,3)) %>% pull(NewID)
degenes_13 <- LIMMA(ns_tumor_exp_norm[,sample_used],
                    ns_tumor_sample[sample_used,"Cluster"])
sample_used <- ns_tumor_sample %>% filter(Cluster %in% c(2,3)) %>% pull(NewID)
degenes_23 <- LIMMA(ns_tumor_exp_norm[,sample_used],
                    ns_tumor_sample[sample_used,"Cluster"])
degenes <- unique(c(
  degenes_12 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 80) %>% pull(Symbol),
  degenes_12 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 80) %>% pull(Symbol),
  degenes_13 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes_13 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes_23 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes_23 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 25) %>% pull(Symbol),
  OF_TAM_signature,OF_Endo_signature,OF_CAF_signature
))
genes_labeled <- unique(c("AGXT","SLC38A4","APOA2","SLC7A2","F13B","ALB","CNTN1","A2M","GNMT","IGHG1","COL1A1","S100A11","LYZ","TMSB10","PLTP","BGN","PTPRC","CYBB","IL7R","THBS2","HLA-DRB1","S100A6","CCL2","C1QA","SPARC","HLA-A","ITGB3","ENG","CTSS","FABP1","SPP1","POSTN",OF_TAM_signature,OF_Endo_signature,OF_CAF_signature))
sample_used <- ns_tumor_sample %>% filter(Tissue == "Tumor") %>% pull(NewID)
mat_for_plot <- ns_tumor_exp_norm[degenes,sample_used]
mat_for_plot <- t(apply(mat_for_plot, 1, scale))
colnames(mat_for_plot) <- sample_used
mat_for_plot <- pmin(mat_for_plot, quantile(mat_for_plot,.975))
mat_for_plot <- pmax(mat_for_plot, quantile(mat_for_plot,.025))
color_used <- circlize::colorRamp2(seq(min(mat_for_plot), max(mat_for_plot), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
ha <-  HeatmapAnnotation(
  Sample = ns_tumor_sample[sample_used,"Sample"],
  col = list(Tissue = Tissue_color_panel,
             Sample = Sample_color_panel))
p <- Heatmap(mat_for_plot,
             show_column_names = F,
             show_row_names = F,
             cluster_rows = T,
             cluster_columns = dend.hc,
             column_split = 3,
             show_row_dend = F,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             top_annotation = ha,
             name = "Exp") +
  rowAnnotation(link = anno_mark(
    at = which(degenes %in% genes_labeled),
    labels = degenes[which(degenes %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
pdf("./04.figures/07.ns_tumor_heatmap.pdf", width = 10, height = 8)
draw(p)
dev.off()

# Proportions of onco-fetal clusters between cluster1 and cluster2----
ns_tumor_sample_temp <- ns_tumor_sample %>% filter(!is.na(Cluster))
ns_tumor_sample_temp$Cluster <- paste0("C",ns_tumor_sample_temp$Cluster)
ns_tumor_sample_temp <- cbind(ns_tumor_sample_temp, tumor_prop2[row.names(ns_tumor_sample_temp), ])
ns_tumor_sample_temp$OF_score <- ns_tumor_sample_temp$`FOLR2+ TAM1` + ns_tumor_sample_temp$`PLVAP+ EC` + ns_tumor_sample_temp$`POSTN+ CAF`
p <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = OF_score)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Onco-fetal Score", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
ggsave(p, file = "./04.figures/07.ns_tumor_onco-fetal_score_C12_boxplot.pdf", width = 1.25, height = 2)
p1 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = `FOLR2+ TAM1`)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Predicted proportions of\n FOLR2+ TAM1", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p2 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = `PLVAP+ EC`)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Predicted proportions of\n PLVAP+ EC", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p3 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = `POSTN+ CAF`)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Predicted proportions of\n POSTN+ CAF", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p4 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = `B`)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Predicted proportions of\n B Cells", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p5 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = `T cells`)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Predicted proportions of\n T Cells", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p6 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = `All Hepatocyte`)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Predicted proportions of\n Hepatocytes", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)
ggsave(p, file = "./04.figures/07.ns_tumor_onco-fetal_proportion_C12_boxplot.pdf", width = 2.5, height = 6)

# DEGenes between cluster 1 and cluster 2----
degenes_12$Sig <- FALSE
degenes_12$label <- NA
degenes_12[degenes_12$adj.P.Val < 0.05 & abs(degenes_12$logFC) > 0.5,"Sig"] <- TRUE
genes_labeled <- c(
  degenes_12 %>% filter(Sig == TRUE) %>% slice_max(n = 25, order_by = logFC) %>% pull(Symbol),
  degenes_12 %>% filter(Sig == TRUE) %>% slice_min(n = 25, order_by = logFC) %>% pull(Symbol)
)
degenes_12[genes_labeled,"label"] <- genes_labeled
write.csv(degenes_12, file = "./03.results/DEGenes/Nanostring_C1_C2_degenes.csv")
volcano_plot <- 
  ggplot(degenes_12, aes(x = -logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = Binomial_color_panel) +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  theme_cowplot(font_size = 16) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
legend <- get_legend(volcano_plot)
volcano_plot <- 
  ggAIplot(volcano_plot + theme(legend.position = "none")) +
  geom_vline(xintercept = c(-0.5,0.5), color = "grey", linetype = "dashed", lwd = 1) +
  geom_text_repel(data = degenes_12, aes(label = label), size = 4)
volcano_plot <- plot_grid(volcano_plot, legend, rel_widths = c(4,1))
ggsave(volcano_plot, file = "./04.figures/07.ns_tumor_degenes12_volcano.pdf", width = 7.75, height = 6)

library(AUCell)
load("./03.results/05_Bulk/HCC_tumor_cells_ranking.RData")
C2_AUC <- AUCell_calcAUC(geneSets = list(G1 = degenes_12 %>% filter(Sig == TRUE) %>% slice_min(n = 50, order_by = logFC) %>% pull(Symbol)), rankings = tumor_cells_rankings, aucMaxRank = nrow(tumor_cells_rankings)*0.05)
C1_AUC <- AUCell_calcAUC(geneSets = list(G1 = degenes_12 %>% filter(Sig == TRUE) %>% slice_max(n = 50, order_by = logFC) %>% pull(Symbol)), rankings = tumor_cells_rankings, aucMaxRank = nrow(tumor_cells_rankings)*0.05)
C2_AUC_matrix <- getAUC(C2_AUC)
C2_AUC_matrix <- pmin(C2_AUC_matrix, quantile(C2_AUC_matrix, .999))
C2_AUC_matrix <- pmax(C2_AUC_matrix, quantile(C2_AUC_matrix, .001))
HCC_seu@meta.data[colnames(C2_AUC_matrix),"C2_AUC"] <- C2_AUC_matrix[1,]
C1_AUC_matrix <- getAUC(C1_AUC)
C1_AUC_matrix <- pmin(C1_AUC_matrix, quantile(C1_AUC_matrix, .990))
C1_AUC_matrix <- pmax(C1_AUC_matrix, quantile(C1_AUC_matrix, .001))
HCC_seu@meta.data[colnames(C1_AUC_matrix),"C1_AUC"] <- C1_AUC_matrix[1,]
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p1 <- ggplot(HCC_seu@meta.data %>% filter(NTF == "Tumor") %>% arrange(C2_AUC), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = C2_AUC), size = .1) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("C2 AUC score", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = .2, barheight = 12))
p1.pdf <- ggAIplot(p1)
p2 <- ggplot(HCC_seu@meta.data %>% filter(NTF == "Tumor") %>% arrange(C1_AUC), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = C1_AUC), size = .1) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("C1 AUC score", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = .2, barheight = 12))
p2.pdf <- ggAIplot(p2)
p <- plot_grid(p1.pdf, p2.pdf)
ggsave(p, file = "./04.figures/07.ns_tumor_degenes12_AUC.pdf", width = 8, height = 3.5)

library(gprofiler2)
query = list(
  "C2" = degenes_12 %>% filter(Sig == TRUE) %>% slice_min(n = 50, order_by = logFC) %>% pull(Symbol),
  "C1" = degenes_12 %>% filter(Sig == TRUE) %>% slice_max(n = 50, order_by = logFC) %>% pull(Symbol)
)
gostres <- gost(
  query = query,
  organism = "mmusculus", ordered_query = FALSE,
  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
  measure_underrepresentation = FALSE, evcodes = FALSE,
  user_threshold = 0.05, correction_method = "g_SCS",
  domain_scope = "annotated", custom_bg = NULL,
  numeric_ns = "", sources = NULL, as_short_link = FALSE)
p1 <- gostres$result %>% filter(query == "C1") %>% arrange(desc(precision)) %>% head(30) %>% 
  ggplot(aes(x = -log10(p_value), y = reorder(term_name, - p_value))) +
  geom_bar(stat = "identity", aes(fill = source)) +
  labs(y = "") + ggtitle("C1") +
  scale_fill_nejm() +
  theme_cowplot()
p2 <- gostres$result %>% filter(query == "C2") %>% arrange(desc(precision)) %>% head(30) %>% 
  ggplot(aes(x = -log10(p_value), y = reorder(term_name, - p_value))) +
  geom_bar(stat = "identity", aes(fill = source)) +
  labs(y = "") + ggtitle("C2") +
  scale_fill_nejm() +
  theme_cowplot()
p <- plot_grid(p1, p2, nrow = 1)
ggsave(p, file = "./04.figures/07.ns_tumor_degenes12_pathway.pdf", width = 18, height = 6)

# Enrichement of predicted ligands in C1 and C2----
prioritized_ligands <- c("LTA","ADAM17","PGF","COL2A1","COL4A1","BMP2","ITGAM","TGFB3","VEGFA","BTLA")
ns_tumor_sample_temp$prioritized_ligands_score <- colMeans(ns_tumor_exp_norm[prioritized_ligands,row.names(ns_tumor_sample_temp)])
p <- ggplot(ns_tumor_sample_temp %>% filter(Cluster != "C3"), aes(x = Cluster, y = prioritized_ligands_score)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Shared Prioritized Ligands", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
ggsave(p, file = "./04.figures/07.ns_tumor_prioritized_ligands_C12_boxplot.pdf", width = 1, height = 3)

# Enrichement of onco-fetal low/high tumor specific genes in C1 and C2----
of_high_tumor_genes <- c("CXCL5","AGR2","MUC6","CXCL1","KRT7","TEX15","CXCL6","DUOX2","KRT19","MMP7","SLC6A19","DMBT1","SYT13","MISP","MUC5B","CA9","KRT80","CLDN4","WNK2","CTNND2")
of_low_tumor_genes <- c("PCK1","CYP2A7","FXYD1","CYP2A6","NR1I3","F9","CYP2E1","SLC25A47","CYP3A4","ADH1B","BHMT","THRSP","CYP8B1","CYP17A1","LECT2","CYP3A43","SLCO1B1")
ns_tumor_sample_temp$of_high_tumor_score <- colMeans(ns_tumor_exp_norm[of_high_tumor_genes,row.names(ns_tumor_sample_temp)])
ns_tumor_sample_temp$of_low_tumor_score <- colMeans(ns_tumor_exp_norm[of_low_tumor_genes,row.names(ns_tumor_sample_temp)])
p1 <- ggplot(ns_tumor_sample_temp, aes(x = Cluster, y = of_high_tumor_score)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Tumor Upregulated Genes", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2"),c("C1","C3"),c("C2","C3")), 
                        step_increase = .1, textsize = 3) +
  theme(legend.position = "none")
p2 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster != "C3") , aes(x = Cluster, y = of_low_tumor_score)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Tumor Upregulated Genes", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p <- p1 + p2 + plot_layout(widths = c(3,2))
ggsave(p, file = "./04.figures/07.ns_tumor_of_related_tumor_boxplot.pdf", width = 2.5, height = 3)

# >Analyze tumor, normal and fetal samples together----
# PCA analysis----
ns_tumor_sd <- sort(apply(ns_tumor_exp_norm,1,sd),decreasing = T)
ns_fetal_sd <- sort(apply(ns_fetal_exp_norm,1,sd),decreasing = T)
core_features <- unique(c(
  names(ns_tumor_sd)[1:500],
  names(ns_fetal_sd)[1:500]
))
pca_analysis <- prcomp(t(ns_exp_norm[core_features,]))
pca_analysis <- summary(pca_analysis)

samples_used <- ns_sample %>% pull(NewID)
hc <- hclust(dist(t(ns_exp_norm[core_features,samples_used]),method = "euclidean"))
dend.hc <- as.dendrogram(hc)
pdf("./04.figures/07.ns_all_dendogram.pdf", width = 14, height = 3)
dend.hc %>% set("labels_col", Tissue_color_panel[ns_sample[labels(dend.hc),"Tissue"]]) %>% plot()
dend.hc %>% rect.dendrogram(k=6, border=8, lty=5, lwd=2) 
dev.off()

cluster <- c(rep("C1",24),rep("C2",10),rep("C3_1",9),rep("C3_2",24),rep("C4",2),rep("C5",3))
ns_sample$Cluster <- c()
ns_sample[labels(dend.hc),"Cluster"] <- cluster

pca_plot <- data.frame(
  pca_analysis$x[,1:2],
  Tissue = ns_sample$Tissue,
  Sample = ns_sample$Sample,
  Cluster = ns_sample$Cluster,
  NewID = ns_sample$NewID)
p1 <- ggplot(pca_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Tissue, shape = Sample), size = 3, alpha = .8) +
  scale_color_manual(values = Tissue_color_panel) +
  xlab(sprintf("PC1: %s%% variance",round(pca_analysis$importance[2,1],2) * 100)) + 
  ylab(sprintf("PC2: %s%% variance",round(pca_analysis$importance[2,2],2) * 100)) +
  theme_cowplot()
p2 <- ggplot(pca_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Cluster), size = 3, alpha = .8) +
  scale_color_manual(values = c54[16:21]) +
  xlab(sprintf("PC1: %s%% variance",round(pca_analysis$importance[2,1],2) * 100)) + 
  ylab(sprintf("PC2: %s%% variance",round(pca_analysis$importance[2,2],2) * 100)) +
  theme_cowplot()
p <- p1 + p2
ggsave(p, file = "./04.figures/07.ns_all_pca.pdf", width = 12, height = 5)

# DEGenes heatmap----
sample_used <- ns_sample %>% filter(Cluster %in% c("C1","C3_2")) %>% pull(NewID)
degenes1 <- LIMMA(ns_exp_norm[,sample_used],
                  ns_sample[sample_used,"Cluster"])
sample_used <- ns_sample %>% filter(Cluster %in% c("C2","C3_1")) %>% pull(NewID)
degenes2 <- LIMMA(ns_exp_norm[,sample_used],
                  ns_sample[sample_used,"Cluster"])
sample_used <- ns_sample %>% filter(Cluster %in% c("C2","C3_2")) %>% pull(NewID)
degenes3 <- LIMMA(ns_exp_norm[,sample_used],
                  ns_sample[sample_used,"Cluster"])
sample_used <- ns_sample %>% filter(Cluster %in% c("C3_1","C3_2")) %>% pull(NewID)
degenes4 <- LIMMA(ns_exp_norm[,sample_used],
                  ns_sample[sample_used,"Cluster"])
sample_used <- ns_sample %>% filter(Cluster %in% c("C3_2","C4")) %>% pull(NewID)
degenes5 <- LIMMA(ns_exp_norm[,sample_used],
                  ns_sample[sample_used,"Cluster"])
sample_used <- ns_sample %>% filter(Cluster %in% c("C4","C5")) %>% pull(NewID)
degenes6 <- LIMMA(ns_exp_norm[,sample_used],
                  ns_sample[sample_used,"Cluster"])
degenes <- unique(c(
  degenes1 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes1 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes2 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes2 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes3 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes3 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes4 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 50) %>% pull(Symbol),
  degenes4 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 50) %>% pull(Symbol),
  degenes5 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes5 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes6 %>% filter(adj.P.Val < 0.05) %>% slice_max(order_by = logFC, n = 25) %>% pull(Symbol),
  degenes6 %>% filter(adj.P.Val < 0.05) %>% slice_min(order_by = logFC, n = 25) %>% pull(Symbol),
  OF_TAM_signature,OF_Endo_signature,OF_CAF_signature
))
genes_labeled <- unique(c("IGHA1","JCHAIN","MYH11","ADH4","CCL14","COL4A1","COL1A1","FGL1","TIMP1","IGF2","LY6E","FOS","LYZ","MRC1","SPP1","IGHG4","VIM","IL7R","A2M","ID1","NR4A1","THBS1","EGR1","CLEC1B","APOA2","SCD","CXCL1","CXCL6","FN1","VTN",OF_TAM_signature,OF_Endo_signature,OF_CAF_signature))
sample_used <- ns_sample %>% pull(NewID)
mat_for_plot <- ns_exp_norm[degenes,sample_used]
mat_for_plot <- t(apply(mat_for_plot, 1, scale))
colnames(mat_for_plot) <- sample_used
mat_for_plot <- pmin(mat_for_plot, quantile(mat_for_plot,.975))
mat_for_plot <- pmax(mat_for_plot, quantile(mat_for_plot,.025))
color_used <- circlize::colorRamp2(seq(min(mat_for_plot), max(mat_for_plot), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
ha <-  HeatmapAnnotation(
  Tissue = ns_sample[sample_used,"Tissue"],
  Sample = ns_sample[sample_used,"Sample"],
  col = list(Tissue = Tissue_color_panel,
             Sample = Sample_color_panel))
p <- Heatmap(mat_for_plot,
             show_column_names = F,
             show_row_names = F,
             cluster_rows = T,
             cluster_columns = dend.hc,
             column_split = 5,
             show_row_dend = F,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             top_annotation = ha,
             name = "Exp") +
  rowAnnotation(link = anno_mark(
    at = which(degenes %in% genes_labeled),
    labels = degenes[which(degenes %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
pdf("./04.figures/07.ns_all_heatmap.pdf", width = 10, height = 8)
draw(p)
dev.off()

# Proportions----
ns_sample[ns_sample$Tissue == "Fetal","Tissue"]
p <- ggplot(ns_sample, aes(x = Cluster)) +
  geom_bar(aes(fill = Tissue), color = "white", position = "fill") +
  scale_fill_manual(values = Tissue_color_panel) +
  labs(y = "Frequency (%)") +
  theme_cowplot() + 
  RotatedAxis()
ggsave(p, file = "./04.figures/07.ns_all_cluster_proportions.pdf", width = 4, height = 5)
