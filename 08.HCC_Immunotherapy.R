setwd("/work/lzy/project/onco_fetal2/")
source("../utils/utils_plot.R")
source("../utils/utils_color.R")
source("../utils/utils_data_processing.R")
library(Seurat)
library(Polychrome)
library(dplyr)
library(dendextend)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ComplexHeatmap)

# >Color panel----
Cluster_levels <- c("B","NK","T cells","Mast","DC","pDC","SPP1+ TAM2","MT1G+ TAM3","MYH11+ CAF","Other Mononuclear","Other Endothelium","Other CAF","FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","All Hepatocyte")
Cluster_color_panel <- c("B"="#F8BFAF","NK"="#EF5276","T cells"="#63B472","Mast"="#DEEAB1","DC"="#588198","pDC"="#D5E7F7",
                         "SPP1+ TAM2"="#D6D4EB","MT1G+ TAM3"="#EACF68","MYH11+ CAF"="#74517B",
                         "Other Mononuclear"="#F8F4A8","Other Endothelium"="#C9BDB2","Other CAF"="#7A8DBB",
                         "FOLR2+ TAM1"="#F3746C","PLVAP+ EC"="#69B4CE","POSTN+ CAF"="#EAA944",
                         "All Hepatocyte"="#A0D7C9")
Cluster_color_panel_sc <- c("B"="#F8BFAF","NK"="#EF5276","CD4+ T"="#63B472","CD8+ T"="#76CFA0","Tregs"="#3BB374",
                            "Mast"="#DEEAB1","DC"="#588198","pDC"="#D5E7F7",
                            "SPP1+ TAM2"="#D6D4EB","MT1G+ TAM3"="#EACF68","MYH11+ CAF"="#74517B",
                            "Other Mononuclear"="#F8F4A8","Other Endothelium"="#C9BDB2","Other CAF"="#7A8DBB",
                            "FOLR2+ TAM1"="#F3746C","PLVAP+ EC"="#69B4CE","POSTN+ CAF"="#EAA944",
                            "Hepatocyte"="#A0D7C9","Hepatocyte_P15"="#8CB3DC","Hepatocyte_P7_2"="#1279DF",
                            "Bi-Potent"="grey")
Tissue_box_color_panel <- c("T" = "#F84141","N" = "#9FA0A3")

# >Load scRNA-seq data----
HCC_seu <- readRDS("../onco_fetal/02.processed_data/HCC_seu.rds")

# >Load bulk RNA-seq data----
ABRS_exp <- read.table("./00.raw_data/HCC_Immunotherapy/genes.TPM_counts.txt", sep = "\t", head = T, row.names = 1, check.names = F)
ABRS_exp <- log2(ABRS_exp + 1)

ABRS_metadata <- read.csv("./00.raw_data/HCC_Immunotherapy/anon_clinical_data_20210913.csv", check.names = F, stringsAsFactors = F)
row.names(ABRS_metadata) <- ABRS_metadata$anon_sampleId
ABRS_metadata <- ABRS_metadata[colnames(ABRS_exp),]
ABRS_metadata$Response[ABRS_metadata$`Confirmed Response_IRF` %in% c("CR","PR")] <- "R"
ABRS_metadata$Response[ABRS_metadata$`Confirmed Response_IRF` %in% c("PD","SD")] <- "NR"
ABRS_metadata[,"PFS censoring\n(1=cens,0=evt)\nIRF"] <- abs(ABRS_metadata[,"PFS censoring\n(1=cens,0=evt)\nIRF"] - 1)

ABRS_gender <- read.csv("./00.raw_data/HCC_Immunotherapy/sample_gender.csv", head = F, stringsAsFactors = F)
ABRS_metadata <- merge(ABRS_metadata, ABRS_gender[,c("V1","V2","V4")], by.x = "anon_sampleId", by.y = "V1")
ABRS_metadata[ABRS_metadata$`Treatment Group` == "A" & ABRS_metadata$V2 == "Pre-treatment","Group"] <- "IMbrave150_AB"
ABRS_metadata[ABRS_metadata$`Treatment Group` == "A" & ABRS_metadata$V2 == "SCRN","Group"] <- "GO30140_A"
ABRS_metadata[ABRS_metadata$`Treatment Group` == "B" & ABRS_metadata$V2 == "Pre-treatment","Group"] <- "IMbrave150_S"
ABRS_metadata[ABRS_metadata$`Treatment Group` == "F1" & ABRS_metadata$V2 == "SCRN","Group"] <- "GO30140_F1"
ABRS_metadata[ABRS_metadata$`Treatment Group` == "F2" & ABRS_metadata$V2 == "SCRN","Group"] <- "GO30140_F2"
colnames(ABRS_metadata)[colnames(ABRS_metadata) == "V4"] <- "Gender"
row.names(ABRS_metadata) <- ABRS_metadata$anon_sampleId

ABRS_cibersort <- read.table("./02.processed_data/CIBERSORTx/Output/CIBERSORTx_Results.txt", sep = "\t", head = T, check.names = F)
row.names(ABRS_cibersort) <- ABRS_cibersort$Mixture
ABRS_cibersort$Mixture <- c()
ABRS_cibersort[,"T cells"] <- ABRS_cibersort[,"CD4+ T"] + ABRS_cibersort[,"CD8+ T"] + ABRS_cibersort[,"Tregs"]
ABRS_cibersort[,"All Hepatocyte"] <- ABRS_cibersort[,"Hepatocyte_P15"] + ABRS_cibersort[,"Hepatocyte_P7_2"] + ABRS_cibersort[,"Hepatocyte"]
ABRS_cibersort$`P-value` <- ABRS_cibersort$Correlation <- ABRS_cibersort$RMSE <- c()
colnames(ABRS_cibersort)[colnames(ABRS_cibersort) == "Other Endotehlium"] <- "Other Endothelium"
ABRS_metadata <- cbind(ABRS_metadata, ABRS_cibersort[row.names(ABRS_metadata),])
ABRS_metadata$OF_score <- ABRS_metadata$`PLVAP+ EC` + ABRS_metadata$`POSTN+ CAF` + ABRS_metadata$`FOLR2+ TAM1`

ABRS_metadata <- ABRS_metadata %>% filter(!is.na(Group))
ABRS_exp <- ABRS_exp[,row.names(ABRS_metadata)]

# >Onco-fetal signature score----
OF_related_genes <- read.csv("/work/lzy/project/onco_fetal/03.results/DEGenes/HCC_bulk_Recurrence_OF_high_vs_OF_low.csv", row.names = 1, stringsAsFactors = F)
OF_TAM_signature <- c("FOLR2","CD163","MRC1","TREM2")
OF_Endo_signature <- c("PLVAP","ACKR1","DLL4")
OF_CAF_signature <- c("POSTN","FAP","ENG","PTGDS")

ABRS_genes <- c("CXCR2P1","ICOS","TIMD4","CTLA4","PAX5","KLRC3","FCRL3","AIM2","GBP5","CCL4")
ABRS_metadata$ABRS_score <- colMeans(ABRS_exp[ABRS_genes,])
ABRS_metadata$OF_TAM_score <- colMeans(ABRS_exp[OF_TAM_signature,])
ABRS_metadata$OF_Endo_score <- colMeans(ABRS_exp[OF_Endo_signature,])
ABRS_metadata$OF_CAF_score <- colMeans(ABRS_exp[OF_CAF_signature,])
ABRS_metadata$OFhiR_score <- colMeans(ABRS_exp[intersect(row.names(ABRS_exp),OF_related_genes %>% slice_max(logFC, n = 20) %>% pull(Symbol)),])
ABRS_metadata$OFlowR_score <- colMeans(ABRS_exp[intersect(row.names(ABRS_exp),OF_related_genes %>% slice_min(logFC, n = 20) %>% pull(Symbol)),])

Treg_genes <- c("FOXP3", "CTLA4", "CCR8", "BATF", "TNFRSF4", "TNFRSF18", "IKZF2", "IL2RA")
ABRS_metadata$Treg_score <- colMeans(ABRS_exp[Treg_genes,])

Tex_genes <- c("ENTPD1","LAYN","PDCD1","HAVCR2","TIGIT","TNFRSF9")
ABRS_metadata$Tex_score <- colMeans(ABRS_exp[Tex_genes,])

# >ABRS and Treg Signature score in scRNA-seq----
HCC_seu@meta.data$ABRS_score <- colMeans(HCC_seu@assays$RNA@data[intersect(ABRS_genes,row.names(HCC_seu)),])
HCC_seu@meta.data$Treg_score <- colMeans(HCC_seu@assays$RNA@data[intersect(Treg_genes,row.names(HCC_seu)),])
HCC_seu@meta.data$temp_Cluster <- as.character(HCC_seu@meta.data$Global_Cluster)
HCC_seu@meta.data$temp_Cluster[HCC_seu@meta.data$Sub_Cluster == "CD8+ T"] <- "CD8+ T"
HCC_seu@meta.data$temp_Cluster[HCC_seu@meta.data$Sub_Cluster == "CD4+ T"] <- "CD4+ T"
HCC_seu@meta.data$temp_Cluster[HCC_seu@meta.data$Sub_Cluster == "Tregs"] <- "Tregs"
HCC_seu@meta.data$temp_Cluster <- factor(HCC_seu@meta.data$temp_Cluster, levels = c("CD8+ T","CD4+ T","Tregs","B cell","ILC","Mononuclear","Mast","Endothelium","Fibroblast","Hepatocyte","Double"))

plot.matrix <- aggregate(t(HCC_seu@assays$RNA@data[intersect(ABRS_genes,row.names(HCC_seu)),]),list(HCC_seu@meta.data$temp_Cluster),mean)
row.names(plot.matrix) <- plot.matrix$Group.1
plot.matrix$Group.1 <- c()
plot.matrix <- plot.matrix[1:10,]
plot.matrix <- apply(plot.matrix, 2, zscore)
plot.matrix_quantile <- quantile(plot.matrix, c(0.05, 0.95))
plot.matrix <- pmax(plot.matrix, plot.matrix_quantile[1])
plot.matrix <- pmin(plot.matrix, plot.matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(plot.matrix), max(plot.matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
p <- Heatmap(plot.matrix,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp")
pdf("./04.figures/R2_02_ABRS_genes_heatmap.pdf", width = 3, height = 2.5)
draw(p)
dev.off()

plot.matrix <- aggregate(t(HCC_seu@assays$RNA@data[intersect(Treg_genes,row.names(HCC_seu)),]),list(HCC_seu@meta.data$temp_Cluster),mean)
row.names(plot.matrix) <- plot.matrix$Group.1
plot.matrix$Group.1 <- c()
plot.matrix <- plot.matrix[1:10,]
plot.matrix <- apply(plot.matrix, 2, zscore)
plot.matrix_quantile <- quantile(plot.matrix, c(0.05, 0.95))
plot.matrix <- pmax(plot.matrix, plot.matrix_quantile[1])
plot.matrix <- pmin(plot.matrix, plot.matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(plot.matrix), max(plot.matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
p <- Heatmap(plot.matrix,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp")
pdf("./04.figures/R2_02_Treg_genes_heatmap.pdf", width = 3, height = 2.5)
draw(p)
dev.off()

p <- 
  ggplot(HCC_seu@meta.data %>% filter(temp_Cluster != "Double"), aes(x = temp_Cluster, y = ABRS_score)) +
  ggbeeswarm::geom_quasirandom(col = "grey", cex = 1, alpha = 0.5) +
  labs(y = "ABRS Score") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, vjust = .5),
        strip.text = element_blank(),
        axis.title.x = element_blank()) 
p.pdf  <- 
  ggAIplot(p, width = 10, height = 10) +
  geom_boxplot(outlier.size = -1, alpha = 0)
ggsave(p.pdf, file = "./04.figures/R2_02_ABRS_genes_boxplot.pdf", width = 2.5, height = 2)

p <- 
  ggplot(HCC_seu@meta.data %>% filter(temp_Cluster != "Double"), aes(x = temp_Cluster, y = Treg_score)) +
  ggbeeswarm::geom_quasirandom(col = "grey", cex = 1, alpha = 0.5) +
  labs(y = "Treg Score") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, vjust = .5),
        strip.text = element_blank(),
        axis.title.x = element_blank()) 
p.pdf  <- 
  ggAIplot(p, width = 10, height = 10) +
  geom_boxplot(outlier.size = -1, alpha = 0)
ggsave(p.pdf, file = "./04.figures/R2_02_Treg_genes_boxplot.pdf", width = 2.5, height = 2)

boxplot_list <- list()

# Response vs non-response group----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
boxplot_list[[1]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Response, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("R","NR")))

p0.1 <- ggplot(ABRS_metadata[sample_used,], aes(x = Response, y = ABRS_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 26) +
  ggsignif::geom_signif(comparisons = list(c("R","NR")))

p0.3 <- ggplot(ABRS_metadata[sample_used,] %>% mutate(`Confirmed Response_IRF` = factor(`Confirmed Response_IRF`, levels = c("CR","PR","SD","PD"))), aes(x = `Confirmed Response_IRF`, y = ABRS_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 26) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("CR","SD"),
                                           c("CR","PD"),
                                           c("PR","SD"),
                                           c("PR","PD")),
                        step_increase = .1)

p0.5 <- ggplot(ABRS_metadata[sample_used,] %>% mutate(`Confirmed Response_IRF` = factor(`Confirmed Response_IRF`, levels = c("CR","PR","SD","PD"))), aes(x = `Confirmed Response_IRF`, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 26) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("CR","SD"),
                                           c("CR","PD"),
                                           c("PR","SD"),
                                           c("PR","PD")),
                        step_increase = .1)

p0.6 <- ggplot(ABRS_metadata[sample_used,] %>% mutate(`Confirmed Response_IRF` = factor(`Confirmed Response_IRF`, levels = c("CR","PR","SD","PD"))), aes(x = `Confirmed Response_IRF`, y = OF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 26) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("CR","SD"),
                                           c("CR","PD"),
                                           c("PR","SD"),
                                           c("PR","PD")),
                        step_increase = .1)

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
boxplot_list[[2]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Response, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("R","NR")))

p0.2 <- ggplot(ABRS_metadata[sample_used,], aes(x = Response, y = ABRS_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 26) +
  ggsignif::geom_signif(comparisons = list(c("R","NR")))

p0.4 <- ggplot(ABRS_metadata[sample_used,] %>% mutate(`Confirmed Response_IRF` = factor(`Confirmed Response_IRF`, levels = c("CR","PR","SD","PD"))), aes(x = `Confirmed Response_IRF`, y = ABRS_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 26) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("CR","SD"),
                                           c("CR","PD"),
                                           c("PR","SD"),
                                           c("PR","PD")),
                        step_increase = .1)

p0.7 <- ggplot(ABRS_metadata[sample_used,] %>% mutate(`Confirmed Response_IRF` = factor(`Confirmed Response_IRF`, levels = c("CR","PR","SD","PD"))), aes(x = `Confirmed Response_IRF`, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 26) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("CR","SD"),
                                           c("CR","PD"),
                                           c("PR","SD"),
                                           c("PR","PD")),
                        step_increase = .1)

p0.8 <- ggplot(ABRS_metadata[sample_used,] %>% mutate(`Confirmed Response_IRF` = factor(`Confirmed Response_IRF`, levels = c("CR","PR","SD","PD"))), aes(x = `Confirmed Response_IRF`, y = OF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 26) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("CR","SD"),
                                           c("CR","PD"),
                                           c("PR","SD"),
                                           c("PR","PD")),
                        step_increase = .1)

# Onco-fetal score in Treg high vs Treg low----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
Treg_score_cutoff <- median(ABRS_metadata[sample_used,"Treg_score"])
ABRS_metadata[sample_used,"Treg_group"] <- "Treg_low"
ABRS_metadata[sample_used,"Treg_group"][ABRS_metadata[sample_used,"Treg_score"] > Treg_score_cutoff] <- "Treg_high"
ABRS_metadata$Treg_group <- factor(ABRS_metadata$Treg_group, levels = c("Treg_low","Treg_high"))

p <- ggplot(ABRS_metadata[sample_used,], aes(x = Treg_score, y = OF_TAM_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_GO30140_Treg_vs_OF_TAM_lm.pdf", width = 3, height = 3)

summary(lm(OF_TAM_score ~ Treg_score, ABRS_metadata[sample_used,]))
# y = 0.7892x+2.8668 adjusted.R2 = 0.3357 
cor.test(ABRS_metadata[sample_used,"OF_TAM_score"],ABRS_metadata[sample_used,"Treg_score"])
# cor 0.5860805 p-value = 3.058e-09

p <- ggplot(ABRS_metadata[sample_used,], aes(x = Treg_score, y = OF_Endo_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_GO30140_Treg_vs_OF_Endo_lm.pdf", width = 3, height = 3)

summary(lm(OF_Endo_score ~ Treg_score, ABRS_metadata[sample_used,]))
# y = 0.32009x+1.59042 adjusted.R2 = 0.1442 
cor.test(ABRS_metadata[sample_used,"OF_Endo_score"],ABRS_metadata[sample_used,"Treg_score"])
# cor 0.3927116 p-value = 0.000184

p <- ggplot(ABRS_metadata[sample_used,], aes(x = Treg_score, y = OF_CAF_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_GO30140_Treg_vs_OF_CAF_lm.pdf", width = 3, height = 3)

summary(lm(OF_CAF_score ~ Treg_score, ABRS_metadata[sample_used,]))
# y = 1.2573x+2.4436 adjusted.R2 = 0.3028
cor.test(ABRS_metadata[sample_used,"OF_CAF_score"],ABRS_metadata[sample_used,"Treg_score"])
# cor 0.557654 p-value = 2.44e-08

boxplot_list[[3]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_TAM_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
boxplot_list[[4]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_Endo_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
boxplot_list[[5]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_CAF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
boxplot_list[[6]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
boxplot_list[[7]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))

p1 <- ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))

p2.1 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
p2.2 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
p <- p2.1 + p2.2
ggsave(p, file = "./04.figures/R2_02_OFrelapsed_Treg_hi_vs_Treg_low_in_GO30140_boxplot.pdf", width = 3, height = 3)

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
Treg_score_cutoff <- median(ABRS_metadata[sample_used,"Treg_score"])
ABRS_metadata[sample_used,"Treg_group"] <- "Treg_low"
ABRS_metadata[sample_used,"Treg_group"][ABRS_metadata[sample_used,"Treg_score"] > Treg_score_cutoff] <- "Treg_high"
ABRS_metadata$Treg_group <- factor(ABRS_metadata$Treg_group, levels = c("Treg_low","Treg_high"))

p <- ggplot(ABRS_metadata[sample_used,], aes(x = Treg_score, y = OF_TAM_score)) +
  geom_point() +
  ggtitle("IMbrave_AB") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_IMbrave_Treg_vs_OF_TAM_lm.pdf", width = 3, height = 3)

summary(lm(OF_TAM_score ~ Treg_score, ABRS_metadata[sample_used,]))
# y = 0.8795x+2.8668 adjusted.R2 = 0.3776 
cor.test(ABRS_metadata[sample_used,"OF_TAM_score"],ABRS_metadata[sample_used,"Treg_score"])
# cor 0.6187989 p-value = 6.382e-14

p <- ggplot(ABRS_metadata[sample_used,], aes(x = Treg_score, y = OF_Endo_score)) +
  geom_point() +
  ggtitle("IMbrave_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_IMbrave_Treg_vs_OF_Endo_lm.pdf", width = 3, height = 3)

summary(lm(OF_Endo_score ~ Treg_score, ABRS_metadata[sample_used,]))
# y = 0.28823x+1.59042 adjusted.R2 = 0.1094 
cor.test(ABRS_metadata[sample_used,"OF_Endo_score"],ABRS_metadata[sample_used,"Treg_score"])
# cor 0.341919 p-value = 0.0001413

p <- ggplot(ABRS_metadata[sample_used,], aes(x = Treg_score, y = OF_CAF_score)) +
  geom_point() +
  ggtitle("IMbrave_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_IMbrave_Treg_vs_OF_CAF_lm.pdf", width = 3, height = 3)

summary(lm(OF_CAF_score ~ Treg_score, ABRS_metadata[sample_used,]))
# y = 1.1672x+2.2319 adjusted.R2 = 0.2768
cor.test(ABRS_metadata[sample_used,"OF_CAF_score"],ABRS_metadata[sample_used,"Treg_score"])
# cor 0.5319507 p-value = 4.799e-10

boxplot_list[[8]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_TAM_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
boxplot_list[[9]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_Endo_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
boxplot_list[[10]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_CAF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
boxplot_list[[11]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
boxplot_list[[12]] <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))

p2 <- ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))

ABRS_boxplot <- plot_grid(plotlist = boxplot_list, align = "hv", nrow = 3)
ggsave(ABRS_boxplot, file = "./04.figures/R2_02_ABRS_comparison_boxplot.pdf", width = 6, height = 8)

p <- p1 + p2
ggsave(p, file = "./04.figures/R2_02_ABRS_OF_score_in_Treg_hig_vs_Treg_low_boxplot.pdf", width = 3, height = 3)

p2.1 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
p2.2 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
p <- p2.1 + p2.2
ggsave(p, file = "./04.figures/R2_02_OFrelapsed_Treg_hi_vs_Treg_low_in_IMbrave_boxplot.pdf", width = 3, height = 3)

# Onco-fetal score in ABRS high vs ABRS low----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ABRS_score_cutoff <- median(ABRS_metadata[sample_used,"ABRS_score"])
ABRS_metadata[sample_used,"ABRS_group"] <- "ABRS_low"
ABRS_metadata[sample_used,"ABRS_group"][ABRS_metadata[sample_used,"ABRS_score"] > ABRS_score_cutoff] <- "ABRS_high"
ABRS_metadata$ABRS_group <- factor(ABRS_metadata$ABRS_group, levels = c("ABRS_low","ABRS_high"))

p <- ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = Treg_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_GO30140_ABRS_vs_Treg_lm.pdf", width = 3, height = 3)

summary(lm(ABRS_score ~ Treg_score, ABRS_metadata[sample_used,]))
# y = 1.296x+0.08113 adjusted.R2 = 0.7098
cor.test(ABRS_metadata[sample_used,"ABRS_score"],ABRS_metadata[sample_used,"Treg_score"])
# cor 0.845 p<2.2e-16

p <- ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_TAM_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_GO30140_ABRS_vs_OF_TAM_lm.pdf", width = 3, height = 3)

summary(lm(OF_TAM_score ~ ABRS_score, ABRS_metadata[sample_used,]))
# y = 0.356x+3.558 adjusted.R2 = 0.155
cor.test(ABRS_metadata[sample_used,"OF_TAM_score"],ABRS_metadata[sample_used,"ABRS_score"])
# cor 0.406 p-value = 0.0001045

p <- ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_Endo_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_GO30140_ABRS_vs_OF_Endo_lm.pdf", width = 3, height = 3)

summary(lm(OF_Endo_score ~ ABRS_score, ABRS_metadata[sample_used,]))
# y = 0.115x+1.966 adjusted.R2 = 0.03558
cor.test(ABRS_metadata[sample_used,"OF_Endo_score"],ABRS_metadata[sample_used,"ABRS_score"])
# cor 0.2166 p-value = 0.04513

p <- ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_CAF_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_GO30140_ABRS_vs_OF_CAF_lm.pdf", width = 3, height = 3)

summary(lm(OF_CAF_score ~ ABRS_score, ABRS_metadata[sample_used,]))
# y = 0.6713x+2.9644 adjusted.R2 = 0.1993
cor.test(ABRS_metadata[sample_used,"OF_CAF_score"],ABRS_metadata[sample_used,"ABRS_score"])
# cor 0.456 p-value = 9.795e-06

p1.1 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_TAM_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p1.2 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_Endo_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p1.3 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_CAF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p1.4 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p1.5 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p <- p1.1 + p1.2 + p1.3 + p1.4 + p1.5 + plot_layout(nrow = 1)
ggsave(p, file = "./04.figures/R2_02_ABRS_hi_vs_ABRS_low_in_GO30140_boxplot.pdf", width = 7, height = 3)

p2.1 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p2.2 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p <- p2.1 + p2.2
ggsave(p, file = "./04.figures/R2_02_OFrelapsed_ABRS_hi_vs_ABRS_low_in_GO30140_boxplot.pdf", width = 3, height = 3)

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ABRS_score_cutoff <- median(ABRS_metadata[sample_used,"ABRS_score"])
ABRS_metadata[sample_used,"ABRS_group"] <- "ABRS_low"
ABRS_metadata[sample_used,"ABRS_group"][ABRS_metadata[sample_used,"ABRS_score"] > ABRS_score_cutoff] <- "ABRS_high"
ABRS_metadata$ABRS_group <- factor(ABRS_metadata$ABRS_group, levels = c("ABRS_low","ABRS_high"))

p <- ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = Treg_score)) +
  geom_point() +
  ggtitle("IMbrave150_AB") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_IMbrave_ABRS_vs_Treg_lm.pdf", width = 3, height = 3)

summary(lm(ABRS_score ~ Treg_score, ABRS_metadata[sample_used,]))
# y = 1.308x+0.16589 adjusted.R2 = 0.7099
cor.test(ABRS_metadata[sample_used,"ABRS_score"],ABRS_metadata[sample_used,"Treg_score"])
# cor 0.844 p<2.2e-16

p <- ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_TAM_score)) +
  geom_point() +
  ggtitle("IMbrave_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_IMbrave_ABRS_vs_OF_TAM_lm.pdf", width = 3, height = 3)

summary(lm(OF_TAM_score ~ ABRS_score, ABRS_metadata[sample_used,]))
# y = 0.49333x+3.08179 adjusted.R2 = 0.2834 
cor.test(ABRS_metadata[sample_used,"OF_TAM_score"],ABRS_metadata[sample_used,"ABRS_score"])
# cor 0.538 p-value = 2.785e-10

p <- ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_Endo_score)) +
  geom_point() +
  ggtitle("IMbrave_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_IMbrave_ABRS_vs_OF_Endo_lm.pdf", width = 3, height = 3)

summary(lm(OF_Endo_score ~ ABRS_score, ABRS_metadata[sample_used,]))
# y = 0.13340x+1.71247 adjusted.R2 = 0.05215
cor.test(ABRS_metadata[sample_used,"OF_Endo_score"],ABRS_metadata[sample_used,"ABRS_score"])
# cor 0.2453173  p-value = 0.007165

p <- ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_CAF_score)) +
  geom_point() +
  ggtitle("IMbrave_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
ggsave(p, file = "./04.figures/R2_02_IMbrave_ABRS_vs_OF_CAF_lm.pdf", width = 3, height = 3)

summary(lm(OF_CAF_score ~ ABRS_score, ABRS_metadata[sample_used,]))
# y = 0.6351x+2.5530 adjusted.R2 = 0.1945
cor.test(ABRS_metadata[sample_used,"OF_CAF_score"],ABRS_metadata[sample_used,"ABRS_score"])
# cor 0.4486701 p-value = 3.096e-07

p1.1 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_TAM_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p1.2 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_Endo_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p1.3 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_CAF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p1.4 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p1.5 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p <- p1.1 + p1.2 + p1.3 + p1.4 + p1.5 + plot_layout(nrow = 1)
ggsave(p, file = "./04.figures/R2_02_ABRS_hi_vs_ABRS_low_in_IMbrave150_boxplot.pdf", width = 7, height = 3)

p2.1 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p2.2 <- 
  ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
p <- p2.1 + p2.2
ggsave(p, file = "./04.figures/R2_02_OFrelapsed_ABRS_hi_vs_ABRS_low_in_IMbrave_boxplot.pdf", width = 3, height = 3)

# >Survival analysis----
library(survival)
library(survminer)

sample_used <- ABRS_metadata %>% filter(Group == "GO30140_A") %>% pull(anon_sampleId)
OF_score_cutoff <- median(ABRS_metadata[sample_used,"OF_score"])
ABRS_metadata[sample_used,"OF_group"] <- "OF_low"
ABRS_metadata[sample_used,"OF_group"][ABRS_metadata[sample_used,"OF_score"] > OF_score_cutoff] <- "OF_high"
ABRS_metadata[sample_used,"OF_group"] <- factor(ABRS_metadata[sample_used,"OF_group"], levels = c("OF_low","OF_high"))
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
surv.curve <- survival::survfit(surv.obj ~ OF_group, data = ABRS_metadata[sample_used,])
cox.fit <- survival::coxph(surv.obj ~ OF_group + Gender, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)
pdf(file = "./04.figures/R2_02_ABRS_GO30140_OF_group_survival.pdf", width = 8, height = 5)
p
dev.off()

sample_used <- ABRS_metadata %>% filter(Group == "GO30140_A") %>% pull(anon_sampleId)
Treg_score_cutoff <- median(ABRS_metadata[sample_used,"Treg_score"])
ABRS_metadata[sample_used,"Treg_group"] <- "Treg_low"
ABRS_metadata[sample_used,"Treg_group"][ABRS_metadata[sample_used,"Treg_score"] > Treg_score_cutoff] <- "Treg_high"
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
surv.curve <- survival::survfit(surv.obj ~ Treg_group, data = ABRS_metadata[sample_used,])
cox.fit <- survival::coxph(surv.obj ~ Treg_group + Gender, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)
pdf(file = "./04.figures/R2_02_ABRS_GO30140_Treg_group_survival.pdf", width = 8, height = 5)
p
dev.off()

sample_used <- ABRS_metadata %>% filter(Group == "IMbrave150_AB") %>% pull(anon_sampleId)
OF_score_cutoff <- median(ABRS_metadata[sample_used,"OF_score"])
ABRS_metadata[sample_used,"OF_group"] <- "OF_low"
ABRS_metadata[sample_used,"OF_group"][ABRS_metadata[sample_used,"OF_score"] > OF_score_cutoff] <- "OF_high"
ABRS_metadata[sample_used,"OF_group"] <- factor(ABRS_metadata[sample_used,"OF_group"], levels = c("OF_low","OF_high"))
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
surv.curve <- survival::survfit(surv.obj ~ OF_group, data = ABRS_metadata[sample_used,])
cox.fit <- survival::coxph(surv.obj ~ OF_group + Gender, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)
pdf(file = "./04.figures/R2_02_ABRS_IMbrave150_OF_group_survival.pdf", width = 8, height = 5)
p
dev.off()

sample_used <- ABRS_metadata %>% filter(Group == "IMbrave150_AB") %>% pull(anon_sampleId)
Treg_score_cutoff <- median(ABRS_metadata[sample_used,"Treg_score"])
ABRS_metadata[sample_used,"Treg_group"] <- "Treg_low"
ABRS_metadata[sample_used,"Treg_group"][ABRS_metadata[sample_used,"Treg_score"] > Treg_score_cutoff] <- "Treg_high"
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
surv.curve <- survival::survfit(surv.obj ~ Treg_group, data = ABRS_metadata[sample_used,])
cox.fit <- survival::coxph(surv.obj ~ Treg_group + Gender, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)
pdf(file = "./04.figures/R2_02_ABRS_IMbrave150_Treg_group_survival.pdf", width = 8, height = 5)
p
dev.off()

sample_used <- ABRS_metadata %>% filter(Group == "IMbrave150_AB") %>% pull(anon_sampleId)
OFhiR_score_cutoff <- median(ABRS_metadata[sample_used,"OFhiR_score"])
ABRS_metadata[sample_used,"OFhiR_group"] <- "OFhiR_low"
ABRS_metadata[sample_used,"OFhiR_group"][ABRS_metadata[sample_used,"OFhiR_score"] > OFhiR_score_cutoff] <- "OFhiR_high"
ABRS_metadata[sample_used,"OFhiR_group"] <- factor(ABRS_metadata[sample_used,"OFhiR_group"], levels = c("OFhiR_low","OFhiR_high"))
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
surv.curve <- survival::survfit(surv.obj ~ OFhiR_group, data = ABRS_metadata[sample_used,])
cox.fit <- survival::coxph(surv.obj ~ OFhiR_group + Gender, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)
pdf(file = "./04.figures/R2_02_ABRS_IMbrave150_OF_group_survival.pdf", width = 8, height = 5)
p
dev.off()

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB","IMbrave150_S")) %>% pull(anon_sampleId)
OF_score_cutoff <- median(ABRS_metadata[sample_used,"OF_score"])
ABRS_metadata[sample_used,"OF_group"] <- "OF_low"
ABRS_metadata[sample_used,"OF_group"][ABRS_metadata[sample_used,"OF_score"] > OF_score_cutoff] <- "OF_high"
ABRS_metadata$OF_group <- factor(ABRS_metadata$OF_group, levels = c("OF_low","OF_high"))

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB","IMbrave150_S"), OF_group == "OF_low") %>% pull(anon_sampleId)
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
cox.fit <- survival::coxph(surv.obj ~ Treatment + Gender, data = ABRS_metadata[sample_used,])
surv.curve <- survival::survfit(surv.obj ~ Treatment, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)
pdf(file = "./04.figures/R2_02_ABRS_IMbrave150_OF_low_treatment_survival.pdf", width = 8, height = 5)
p
dev.off()

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB","IMbrave150_S"), OF_group == "OF_high") %>% pull(anon_sampleId)
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
cox.fit <- survival::coxph(surv.obj ~ Treatment + Gender, data = ABRS_metadata[sample_used,])
surv.curve <- survival::survfit(surv.obj ~ Treatment, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)
pdf(file = "./04.figures/R2_02_ABRS_IMbrave150_OF_high_treatment_survival.pdf", width = 8, height = 5)
p
dev.off()

sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_F1","GO30140_F2")) %>% pull(anon_sampleId)
OF_score_cutoff <- median(ABRS_metadata[sample_used,"OF_score"])
ABRS_metadata[sample_used,"OF_group"] <- "OF_low"
ABRS_metadata[sample_used,"OF_group"][ABRS_metadata[sample_used,"OF_score"] > OF_score_cutoff] <- "OF_high"
ABRS_metadata$OF_group <- factor(ABRS_metadata$OF_group, levels = c("OF_low","OF_high"))

sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_F1","GO30140_F2"), OF_group == "OF_low") %>% pull(anon_sampleId)
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
cox.fit <- survival::coxph(surv.obj ~ Treatment + Gender, data = ABRS_metadata[sample_used,])
surv.curve <- survival::survfit(surv.obj ~ Treatment, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)

sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_F1","GO30140_F2"), OF_group == "OF_high") %>% pull(anon_sampleId)
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
cox.fit <- survival::coxph(surv.obj ~ Treatment + Gender, data = ABRS_metadata[sample_used,])
surv.curve <- survival::survfit(surv.obj ~ Treatment, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)

sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_F1","GO30140_F2")) %>% pull(anon_sampleId)
Treg_score_cutoff <- median(ABRS_metadata[sample_used,"Treg_score"])
ABRS_metadata[sample_used,"Treg_group"] <- "Treg_low"
ABRS_metadata[sample_used,"Treg_group"][ABRS_metadata[sample_used,"Treg_score"] > Treg_score_cutoff] <- "Treg_high"
ABRS_metadata$Treg_group <- factor(ABRS_metadata$Treg_group, levels = c("Treg_low","Treg_high"))

sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_F1","GO30140_F2"), Treg_group == "Treg_low") %>% pull(anon_sampleId)
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
cox.fit <- survival::coxph(surv.obj ~ Treatment + Gender, data = ABRS_metadata[sample_used,])
surv.curve <- survival::survfit(surv.obj ~ Treatment, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)

sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_F1","GO30140_F2"), Treg_group == "Treg_high") %>% pull(anon_sampleId)
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
cox.fit <- survival::coxph(surv.obj ~ Treatment + Gender, data = ABRS_metadata[sample_used,])
surv.curve <- survival::survfit(surv.obj ~ Treatment, data = ABRS_metadata[sample_used,])
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)

sample_used <- ABRS_metadata %>% filter(`Treatment Group` %in% c("A","F1"), V2 == "ON-STUDY") %>% pull(anon_patientId)
ABRS_metadata %>% filter(anon_patientId %in% sample_used) %>% mutate(Visit = factor(Visit, levels = c("Pre-treatment","Post-treatment"))) %>%
  ggplot(aes(x = Visit, y = Treg_score, group = anon_patientId)) +
  geom_line(aes(color = anon_patientId)) + 
  labs(x = "") +
  theme_cowplot(font_size = 20)
