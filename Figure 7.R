setwd("/Volumes/ZiyiLi/PostDoc_Project/Onco-fetal Project/")
source("./onco-fetal_utils.R")
library(Seurat)
library(Polychrome)
library(dplyr)
library(dendextend)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ComplexHeatmap)

# Load data----
HCC_seu <- readRDS("./Onco-fetal CAF Identification/HCC_Ankur_release.rds")
load("./HCC immunotherapy/ABRS.rda")
ABRS_metadata <- ABRS_metadata %>% filter(!is.na(Group))
ABRS_exp <- ABRS_exp[,row.names(ABRS_metadata)]
ABRS_genes <- c("CXCR2P1","ICOS","TIMD4","CTLA4","PAX5","KLRC3","FCRL3","AIM2","GBP5","CCL4")
OF_TAM_signature <- c("FOLR2","CD163","MRC1","TREM2")
OF_Endo_signature <- c("PLVAP","ACKR1","DLL4")
OF_CAF_signature <- c("POSTN","FAP","ENG","PTGDS")
OF_high_gene <- c("CXCL5","AGR2","MUC6","CXCL1","KRT7","TEX15","CXCL6","DUOX2","KRT19","MMP7","SLC6A19","DMBT1","SYT13","MISP","MUC5B","CA9","KRT80","CLDN4","WNK2","CTNND2")
OF_low_gene <- c("PCK1","CYP2A7","FXYD1","CYP2A6","NR1I3","F9","G6PC","ZNF648","APOC4","CYP2E1","SLC25A47","CYP3A4","ADH1B","BHMT","THRSP","CYP8B1","CYP17A1","LECT2","CYP3A43","SLCO1B1")
Treg_genes <- c("FOXP3", "CTLA4", "CCR8", "BATF", "TNFRSF4", "TNFRSF18", "IKZF2", "IL2RA")
Tex_genes <- c("ENTPD1","LAYN","PDCD1","HAVCR2","TIGIT","TNFRSF9")
ABRS_metadata$ABRS_score <- colMeans(ABRS_exp[ABRS_genes,])
ABRS_metadata$OF_TAM_score <- colMeans(ABRS_exp[OF_TAM_signature,])
ABRS_metadata$OF_Endo_score <- colMeans(ABRS_exp[OF_Endo_signature,])
ABRS_metadata$OF_CAF_score <- colMeans(ABRS_exp[OF_CAF_signature,])
ABRS_metadata$OFhiR_score <- colMeans(ABRS_exp[intersect(row.names(ABRS_exp),OF_high_gene),])
ABRS_metadata$OFlowR_score <- colMeans(ABRS_exp[intersect(row.names(ABRS_exp),OF_low_gene),])
ABRS_metadata$Tex_score <- colMeans(ABRS_exp[Tex_genes,])
ABRS_metadata$Treg_score <- colMeans(ABRS_exp[Treg_genes,])

# Fig.7b----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ABRS_score_cutoff <- median(ABRS_metadata[sample_used,"ABRS_score"])
ABRS_metadata[sample_used,"ABRS_group"] <- "ABRS_low"
ABRS_metadata[sample_used,"ABRS_group"][ABRS_metadata[sample_used,"ABRS_score"] > ABRS_score_cutoff] <- "ABRS_high"
ABRS_metadata$ABRS_group <- factor(ABRS_metadata$ABRS_group, levels = c("ABRS_low","ABRS_high"))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))

# Fig.7c----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ggplot(ABRS_metadata[sample_used,], aes(x = Response, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("R","NR")))

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ggplot(ABRS_metadata[sample_used,], aes(x = Response, y = Treg_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("R","NR")))

# Fig.7d----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ABRS_score_cutoff <- median(ABRS_metadata[sample_used,"ABRS_score"])
ABRS_metadata[sample_used,"ABRS_group"] <- "ABRS_low"
ABRS_metadata[sample_used,"ABRS_group"][ABRS_metadata[sample_used,"ABRS_score"] > ABRS_score_cutoff] <- "ABRS_high"
ABRS_metadata$ABRS_group <- factor(ABRS_metadata$ABRS_group, levels = c("ABRS_low","ABRS_high"))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_TAM_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_Endo_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_CAF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ABRS_score_cutoff <- median(ABRS_metadata[sample_used,"ABRS_score"])
ABRS_metadata[sample_used,"ABRS_group"] <- "ABRS_low"
ABRS_metadata[sample_used,"ABRS_group"][ABRS_metadata[sample_used,"ABRS_score"] > ABRS_score_cutoff] <- "ABRS_high"
ABRS_metadata$ABRS_group <- factor(ABRS_metadata$ABRS_group, levels = c("ABRS_low","ABRS_high"))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_TAM_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_Endo_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_CAF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("ABRS_high","ABRS_low")))

# Extended Data Fig.10a----
HCC_seu@meta.data$ABRS_score <- colMeans(HCC_seu@assays$RNA@data[intersect(ABRS_genes,row.names(HCC_seu)),])
HCC_seu@meta.data$Treg_score <- colMeans(HCC_seu@assays$RNA@data[intersect(Treg_genes,row.names(HCC_seu)),])
HCC_seu@meta.data$temp_Cluster <- as.character(HCC_seu@meta.data$Global_Cluster)
HCC_seu@meta.data$temp_Cluster[HCC_seu@meta.data$Sub_Cluster == "CD8+ T"] <- "CD8+ T"
HCC_seu@meta.data$temp_Cluster[HCC_seu@meta.data$Sub_Cluster == "CD4+ T"] <- "CD4+ T"
HCC_seu@meta.data$temp_Cluster[HCC_seu@meta.data$Sub_Cluster == "Tregs"] <- "Tregs"
HCC_seu@meta.data$temp_Cluster <- factor(HCC_seu@meta.data$temp_Cluster, levels = c("CD8+ T","CD4+ T","Tregs","B cell","ILC","Mononuclear","Mast","Endothelium","Fibroblast","Hepatocyte","Double"))
plot.matrix <- aggregate(t(as.matrix(HCC_seu@assays$RNA@data[intersect(ABRS_genes,row.names(HCC_seu)),])),list(HCC_seu@meta.data$temp_Cluster),mean)
row.names(plot.matrix) <- plot.matrix$Group.1
plot.matrix$Group.1 <- c()
plot.matrix <- plot.matrix[1:10,]
plot.matrix.scale <- apply(plot.matrix, 2, scale)
row.names(plot.matrix.scale) <- row.names(plot.matrix)
plot.matrix_quantile <- quantile(plot.matrix.scale, c(0.05, 0.95))
plot.matrix.scale <- pmax(plot.matrix.scale, plot.matrix_quantile[1])
plot.matrix.scale <- pmin(plot.matrix.scale, plot.matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(plot.matrix.scale), max(plot.matrix.scale), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
p <- Heatmap(plot.matrix.scale,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp")

ggplot(HCC_seu@meta.data %>% filter(temp_Cluster != "Double"), aes(x = temp_Cluster, y = ABRS_score)) +
  ggbeeswarm::geom_quasirandom(col = "grey", cex = 1, alpha = 0.5) +
  geom_boxplot(outlier.size = -1, alpha = 0) +
  labs(y = "ABRS Score") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, vjust = .5),
        strip.text = element_blank(),
        axis.title.x = element_blank())

# Extended Data Fig.10b----
plot.matrix <- aggregate(t(as.matrix(HCC_seu@assays$RNA@data[intersect(Treg_genes,row.names(HCC_seu)),])),list(HCC_seu@meta.data$temp_Cluster),mean)
row.names(plot.matrix) <- plot.matrix$Group.1
plot.matrix$Group.1 <- c()
plot.matrix <- plot.matrix[1:10,]
plot.matrix.scale <- apply(plot.matrix, 2, scale)
row.names(plot.matrix.scale) <- row.names(plot.matrix)
plot.matrix_quantile <- quantile(plot.matrix.scale, c(0.05, 0.95))
plot.matrix.scale <- pmax(plot.matrix.scale, plot.matrix_quantile[1])
plot.matrix.scale <- pmin(plot.matrix.scale, plot.matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(plot.matrix.scale), max(plot.matrix.scale), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
p <- Heatmap(plot.matrix.scale,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp")

ggplot(HCC_seu@meta.data %>% filter(temp_Cluster != "Double"), aes(x = temp_Cluster, y = Treg_score)) +
  ggbeeswarm::geom_quasirandom(col = "grey", cex = 1, alpha = 0.5) +
  geom_boxplot(outlier.size = -1, alpha = 0) +
  labs(y = "Treg Score") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, vjust = .5),
        strip.text = element_blank(),
        axis.title.x = element_blank())

# Extended Data Fig.10c----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = Treg_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
summary(lm(Treg_score ~ ABRS_score, ABRS_metadata[sample_used,]))
cor.test(ABRS_metadata[sample_used,"Treg_score"],ABRS_metadata[sample_used,"ABRS_score"])

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = Treg_score)) +
  geom_point() +
  ggtitle("IMbrave150_AB") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
summary(lm(Treg_score ~ ABRS_score, ABRS_metadata[sample_used,]))
cor.test(ABRS_metadata[sample_used,"Treg_score"],ABRS_metadata[sample_used,"ABRS_score"])

# Extended Data Fig.10d----
library(survival)
library(survminer)
sample_used <- ABRS_metadata %>% filter(Group == "GO30140_A") %>% pull(anon_sampleId)
Treg_score_cutoff <- median(ABRS_metadata[sample_used,"Treg_score"])
ABRS_metadata[sample_used,"Treg_group"] <- "Treg_low"
ABRS_metadata[sample_used,"Treg_group"][ABRS_metadata[sample_used,"Treg_score"] > Treg_score_cutoff] <- "Treg_high"
surv.obj <- survival::Surv(ABRS_metadata[sample_used,"PFS \nin days\nIRF"], ABRS_metadata[sample_used,"PFS censoring\n(1=cens,0=evt)\nIRF"])
surv.curve <- survival::survfit(surv.obj ~ Treg_group, data = ABRS_metadata[sample_used,])
cox.fit <- survival::coxph(surv.obj ~ Treg_group + Gender, data = ABRS_metadata[sample_used,])
ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)

# Extended Data Fig.10e----
# see Fig.7d

# Extended Data Fig.10f----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_TAM_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
summary(lm(OF_TAM_score ~ ABRS_score, ABRS_metadata[sample_used,]))
cor.test(ABRS_metadata[sample_used,"OF_TAM_score"],ABRS_metadata[sample_used,"ABRS_score"])
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_Endo_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
summary(lm(OF_Endo_score ~ ABRS_score, ABRS_metadata[sample_used,]))
cor.test(ABRS_metadata[sample_used,"OF_Endo_score"],ABRS_metadata[sample_used,"ABRS_score"])
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_CAF_score)) +
  geom_point() +
  ggtitle("GO30140_A") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
summary(lm(OF_CAF_score ~ ABRS_score, ABRS_metadata[sample_used,]))
cor.test(ABRS_metadata[sample_used,"OF_CAF_score"],ABRS_metadata[sample_used,"ABRS_score"])

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_TAM_score)) +
  geom_point() +
  ggtitle("IMbrave150_AB") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
summary(lm(OF_TAM_score ~ ABRS_score, ABRS_metadata[sample_used,]))
cor.test(ABRS_metadata[sample_used,"OF_TAM_score"],ABRS_metadata[sample_used,"ABRS_score"])
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_Endo_score)) +
  geom_point() +
  ggtitle("IMbrave150_AB") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
summary(lm(OF_Endo_score ~ ABRS_score, ABRS_metadata[sample_used,]))
cor.test(ABRS_metadata[sample_used,"OF_Endo_score"],ABRS_metadata[sample_used,"ABRS_score"])
ggplot(ABRS_metadata[sample_used,], aes(x = ABRS_score, y = OF_CAF_score)) +
  geom_point() +
  ggtitle("IMbrave150_AB") +
  geom_smooth(method='lm', formula=y~x) +
  theme_cowplot(font_size = 7)
summary(lm(OF_CAF_score ~ ABRS_score, ABRS_metadata[sample_used,]))
cor.test(ABRS_metadata[sample_used,"OF_CAF_score"],ABRS_metadata[sample_used,"ABRS_score"])

# Extended Data Fig.10g-i----
sample_used <- ABRS_metadata %>% filter(Group %in% c("GO30140_A"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
Treg_score_cutoff <- median(ABRS_metadata[sample_used,"Treg_score"])
ABRS_metadata[sample_used,"Treg_group"] <- "Treg_low"
ABRS_metadata[sample_used,"Treg_group"][ABRS_metadata[sample_used,"Treg_score"] > Treg_score_cutoff] <- "Treg_high"
ABRS_metadata$Treg_group <- factor(ABRS_metadata$Treg_group, levels = c("Treg_low","Treg_high"))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_TAM_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_Endo_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_CAF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("GO30140_A") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))

sample_used <- ABRS_metadata %>% filter(Group %in% c("IMbrave150_AB"), Response %in% c("R","NR")) %>% pull(anon_sampleId)
Treg_score_cutoff <- median(ABRS_metadata[sample_used,"Treg_score"])
ABRS_metadata[sample_used,"Treg_group"] <- "Treg_low"
ABRS_metadata[sample_used,"Treg_group"][ABRS_metadata[sample_used,"Treg_score"] > Treg_score_cutoff] <- "Treg_high"
ABRS_metadata$Treg_group <- factor(ABRS_metadata$Treg_group, levels = c("Treg_low","Treg_high"))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_TAM_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_Endo_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_CAF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OF_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  labs(x = "") +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFhiR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
ggplot(ABRS_metadata[sample_used,], aes(x = Treg_group, y = OFlowR_score)) +
  geom_boxplot() +
  ggtitle("IMbrave150_AB") +
  theme_cowplot(font_size = 7) +
  ggsignif::geom_signif(comparisons = list(c("Treg_high","Treg_low")))
# ABRS group in Extended Fig.10i see Fig.7d

