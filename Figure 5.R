setwd("/Volumes/ZiyiLi/PostDoc_Project/Onco-fetal Project/")
source("./onco-fetal_utils.R")
library(Seurat)
library(MERINGUE)
library(Polychrome)
library(dplyr)
library(dendextend)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)

# Load data----
HCC_seu <- readRDS("./Onco-fetal CAF Identification/HCC_Ankur_release.rds")
PLANET_exp <- readRDS("./PLANET cohort/PLANET_exp.rds")
PLANET_metadata <- readRDS("./PLANET cohort/PLANET_metadata.rds")
PLANET_cibersort <- readRDS("./PLANET cohort/PLANET_cibersort.rds")
Cluster_levels <- c("B","NK","T cells","Mast","DC","pDC","SPP1+ TAM2","MT1G+ TAM3","MYH11+ CAF","Other Mononuclear","Other Endothelium","Other CAF","FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","All Hepatocyte")
PLANET_metadata[row.names(PLANET_cibersort),"OF_score"] <- rowSums(PLANET_cibersort[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")])
sample_used <- PLANET_metadata %>% filter(Bulk_sample_id %in% row.names(PLANET_cibersort), !is.na(Recurrence)) %>% pull(Bulk_sample_id)
OF_score_cutoff <- mean(PLANET_metadata[sample_used,"OF_score"]) + sd(PLANET_metadata[sample_used,"OF_score"])
PLANET_metadata[sample_used,"OF_group"] <- "OF_low"
PLANET_metadata[sample_used,"OF_group"][PLANET_metadata[sample_used,"OF_score"] > OF_score_cutoff] <- "OF_high"
PLANET_metadata$OF_group <- factor(PLANET_metadata$OF_group, levels = c("OF_low","OF_high"))

# Fig.5a----
hc <- hclust(dist(PLANET_cibersort[,Cluster_levels], method = "euclidean"), method = "ward.D2")
dend.hc <- as.dendrogram(hc)

PLANET_cibersort_bar <- PLANET_cibersort[,Cluster_levels]
PLANET_cibersort_bar$`Other Mononuclear` <- PLANET_cibersort_bar$`Other Mononuclear` + PLANET_cibersort_bar$Mast + PLANET_cibersort_bar$DC + PLANET_cibersort_bar$pDC
PLANET_cibersort_bar$`Other CAF` <- PLANET_cibersort_bar$`Other CAF` + PLANET_cibersort_bar$`MYH11+ CAF`
Cluster_levels2 <- c("B","NK","T cells","SPP1+ TAM2","MT1G+ TAM3","Other Mononuclear","Other Endothelium","Other CAF","FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","All Hepatocyte")
plot.df <- PLANET_cibersort_bar[,Cluster_levels2] %>% mutate(Sample = row.names(.)) %>% melt() %>% mutate(Sample = factor(Sample, levels = labels(dend.hc)), variable = factor(variable, levels = Cluster_levels))
ggplot(plot.df, aes(x = Sample, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", lwd = .01) +
  labs(y = "Predicted proportions", x = "") +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ha = HeatmapAnnotation(
  OF_score = anno_barplot(PLANET_metadata[labels(dend.hc),"OF_score"],
                          baseline = OF_score_cutoff, border = F,
                          gp = gpar(fill = "#808080", col = "#808080")),
  Recurrence = PLANET_metadata[labels(dend.hc),"Recurrence"],
  AFP = PLANET_metadata[labels(dend.hc),"AFPgroup"],
  RFS = PLANET_metadata[labels(dend.hc),"RFSgroup"],
  Viral.Status = PLANET_metadata[labels(dend.hc),"Viral.status"],
  TNM.Stage = PLANET_metadata[labels(dend.hc),"TNM.Stage"],
  col = list(Recurrence = c("Yes"="#F97137","No"="#4FB0AF"),
             AFP = c("AFPhi"="#6d327f","AFPlow"="#bda7cc"),
             RFS = c("0-6 month"="#d69b95","6-12 month"="#c2726d",">12 month"="#af4c4c"),
             TNM.Stage = c("Stage I"="#D1CEE8","Stage II"="#A6A4D4","Stage III"="#7575BC"),
             Viral.Status = c("Hep B carrier"="#E05D89","NBNC"="#3BAFE0")
  )
)
Heatmap(matrix = matrix(rnorm(1980), ncol = 198), 
        cluster_rows = F, cluster_columns = F,
        top_annotation = ha)

data_used <- data.frame(
  PLANET_cibersort_bar[labels(dend.hc),Cluster_levels2],
  OF_score = PLANET_metadata[labels(dend.hc),"OF_score"],
  Recurrence = PLANET_metadata[labels(dend.hc),"Recurrence"],
  AFP = PLANET_metadata[labels(dend.hc),"AFPgroup"],
  RFS = PLANET_metadata[labels(dend.hc),"RFSgroup"],
  Viral.Status = PLANET_metadata[labels(dend.hc),"Viral.status"],
  TNM.Stage = PLANET_metadata[labels(dend.hc),"TNM.Stage"])

# Fig.5b----
sample_used <- PLANET_metadata %>% filter(Bulk_sample_id %in% row.names(PLANET_cibersort), !is.na(Recurrence)) %>% pull(Bulk_sample_id)
sample_OF_table <- table(PLANET_metadata[sample_used,c("sample_ID","OF_group")])
sample_name <- row.names(sample_OF_table)
sample_OF <- sample_OF_table[,2] > 0
sample_recur_table <- table(PLANET_metadata[sample_used,c("sample_ID","Recurrence")])
sample_recur <- (sample_recur_table[,2] > 0)[sample_name]
sample_MVI_table <- table(PLANET_metadata[sample_used,c("sample_ID","MVI")])
sample_MVI <- (sample_MVI_table[,2] > 0)[sample_name]
p1 <- CrossTabPlot(table(sample_OF, sample_recur))
p2 <- CrossTabPlot(table(sample_OF, sample_MVI))
p1 + p2

# Fig.5c----
PLANET_metadata[row.names(PLANET_cibersort),c("Recurrence","OF_score")] %>% filter(!is.na(Recurrence)) %>%
  ggplot(aes(x = Recurrence, y = OF_score)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("Yes","No"))) +
  labs(x = "Recurrence", y = "Onco-fetal Score") +
  ggtitle(label = "", subtitle = "No(n=71), Yes(n=124)") +
  theme_cowplot(font_size = 8) +
  theme(strip.background = element_blank(),
        strip.switch.pad.grid = unit(1, "inch"))
cbind(Recurrence = PLANET_metadata[row.names(PLANET_cibersort),c("Recurrence")],PLANET_cibersort[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")]) %>% filter(!is.na(Recurrence)) %>% melt() %>%
  ggplot(aes(x = Recurrence, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  ggsignif::geom_signif(comparisons = list(c("Yes","No"))) +
  labs(x = "Recurrence", y = "") +
  ggtitle(label = "", subtitle = "No(n=71), Yes(n=124)") +
  theme_cowplot(font_size = 8) +
  theme(strip.background = element_blank(),
        strip.switch.pad.grid = unit(1, "inch"))

# Fig.5d----
sample_used <- PLANET_metadata %>% filter(Bulk_sample_id %in% row.names(PLANET_cibersort), !is.na(Recurrence)) %>% pull(Bulk_sample_id)
sample_OF_table <- table(PLANET_metadata[sample_used,c("sample_ID","OF_group")])
sample_name <- row.names(sample_OF_table)
sample_OF <- sample_OF_table[,2] > 0
RFS <- read.csv("./PLANET cohort/PLANET_RFS.csv", stringsAsFactors = F)
row.names(RFS) <- RFS$sample_ID
RFS["16-003799","Old.RFS"] <- 46
RFS$OFgroup <- sample_OF[row.names(RFS)]
RFS$OFgroup[RFS$OFgroup == TRUE] <- "OF high"
RFS$OFgroup[RFS$OFgroup == FALSE] <- "OF low"
RFS$OFgroup <- factor(RFS$OFgroup, levels = c("OF low","OF high"))
RFS$RFSgroup <- NA
RFS[!is.na(RFS$Old.RFS) & RFS$Old.Recurrence == "Yes" & RFS$Old.RFS <= 180,"RFSgroup"] <- "0-6 month"
RFS[!is.na(RFS$Old.RFS) & RFS$Old.Recurrence == "Yes" & RFS$Old.RFS > 180 & RFS$Old.RFS <= 360,"RFSgroup"] <- "6-12 month"
RFS[!is.na(RFS$Old.RFS) & RFS$Old.Recurrence == "Yes" & RFS$Old.RFS > 360,"RFSgroup"] <- ">12 month"
RFS$RFSgroup <- factor(RFS$RFSgroup, levels = c(">12 month","6-12 month","0-6 month"))
RFS <- RFS %>% filter(!is.na(OFgroup), Old.Recurrence == "Yes")
RFS %>% ggplot(aes(x = OFgroup, fill = RFSgroup)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("0-6 month"="#af4c4c","6-12 month"="#c2726d",">12 month"="#d69b95")) +
  theme_cowplot() +
  labs(x = "", y = "") +
  RotatedAxis()
cross_tab <- table(RFS[,c("OFgroup","RFSgroup")])
cross_tab_prop <- cross_tab/rowSums(cross_tab)

plot.data <- ROIE(table(RFS$RFSgroup,RFS$OFgroup), filter = 0)

# Fig.5e----
library(richR)
relapsed_samples <- PLANET_metadata %>% filter(Bulk_sample_id %in% row.names(PLANET_cibersort), Recurrence == "Yes") %>% pull(Bulk_sample_id)
degenes <- LIMMA(PLANET_exp[,relapsed_samples],
                 PLANET_metadata[relapsed_samples,"OF_group"])
degenes$Sig <- FALSE
degenes[degenes$adj.P.Val < 0.05 & abs(degenes$logFC) > 1,"Sig"] <- TRUE
hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
OF_high_gene_name <- degenes %>% filter(Sig == TRUE, logFC < 0) %>% pull(Symbol)
OF_high_go <- richGO(OF_high_gene_name,godata = hsago, ontology = c("BP"))
OF_low_gene_name <- degenes %>% filter(Sig == TRUE, logFC > 0) %>% pull(Symbol)
OF_low_go <- richGO(OF_low_gene_name,godata = hsago, ontology = c("BP"))
# The results will be changed based on database used
OF_high_go_show <- OF_high_go@result$Term %in% c("extracellular matrix organization","cell differentiation","cell motility","cell population proliferation","embryo development","collagen metabolic process","blood vessel development","ERK1 and ERK2 cascade","Wnt signaling pathway","negative regulation of apoptotic process","epithelial to mesenchymal transition")
OF_low_go_show <- OF_low_go@result$Term %in% c("organic acid biosynthetic process","lipid catabolic process","lipid localization","lipid transport","lipid modification","long-chain fatty acid metabolic process","cholesterol metabolic process","acute inflammatory response","fatty acid oxidation","response to cAMP")
pathway_plot.df <- rbind(
  data.frame(
    pathway = OF_high_go@result$Term[OF_high_go_show],
    pvalue = -log10(OF_high_go@result$Padj[OF_high_go_show]),
    generatio = (OF_high_go@result$Significant/OF_high_go@result$Annotated)[OF_high_go_show],
    group = "OF_high",
    database = "GO (BP)"),
  data.frame(
    pathway = OF_low_go@result$Term[OF_low_go_show],
    pvalue = log10(OF_low_go@result$Padj[OF_low_go_show]),
    generatio = (OF_low_go@result$Significant/OF_low_go@result$Annotated)[OF_low_go_show],
    group = "OF_low",
    database = "GO (BP)")
)
ggplot(pathway_plot.df, aes(x = reorder(pathway, pvalue), y = pvalue)) +
  geom_bar(stat = "identity", aes(fill = group), alpha = 0.8) +
  geom_text(data = pathway_plot.df %>% filter(pvalue < 0),
            aes(y = -sign(pvalue) * 0.5, label = pathway, hjust = 0),
            size = 2) +
  geom_text(data = pathway_plot.df %>% filter(pvalue > 0), 
            aes(y = -sign(pvalue) * 0.5, label = pathway, hjust = 1),
            size = 2) +
  scale_fill_manual(values = c("OF_high" = "#c15fa3", "OF_low" = "#f3c3b7")) +
  labs(y = "-log10 (P-Value)") +
  coord_flip() +
  theme_cowplot(font_size = 7) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "none", legend.title = element_blank())

# Fig.5f----
sample_used <- PLANET_metadata %>% filter(Bulk_sample_id %in% row.names(PLANET_cibersort), Recurrence == "Yes") %>% pull(Bulk_sample_id)
PLANET_exp_sd <- apply(PLANET_exp[,sample_used],1,sd)
genes_used <- names(sort(PLANET_exp_sd, decreasing = T))[1:1000]
pca <- princomp(PLANET_exp[genes_used,sample_used])
pca_summary <- summary(pca)
hc <- hclust(dist(pca_summary$loadings[,1:10], method = "euclidean"), method = "ward.D2")
dend.hc <- as.dendrogram(hc)
sample_group <- cutree(hc, 3)
sample_used_1_2 <- names(sample_group)[sample_group %in% c(1,2)]
degenes_1_2 <- LIMMA(PLANET_exp[,sample_used_1_2],
                     sample_group[sample_used_1_2])
sample_used_1_3 <- names(sample_group)[sample_group %in% c(1,3)]
degenes_1_3 <- LIMMA(PLANET_exp[,sample_used_1_3],
                     sample_group[sample_used_1_3])
sample_used_2_3 <- names(sample_group)[sample_group %in% c(2,3)]
degenes_2_3 <- LIMMA(PLANET_exp[,sample_used_2_3],
                     sample_group[sample_used_2_3])
degenes_OF <- LIMMA(PLANET_exp[,sample_used],
                 PLANET_metadata[sample_used,"OF_group"])
degenes_OF$Sig <- FALSE
degenes_OF[degenes_OF$adj.P.Val < 0.05 & abs(degenes_OF$logFC) > 1,"Sig"] <- TRUE
genes_used <- unique(c(
  degenes_1_2 %>% filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>% slice_max(n = 40, order_by = logFC) %>% pull(Symbol),
  degenes_1_2 %>% filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>% slice_min(n = 40, order_by = logFC) %>% pull(Symbol),
  degenes_1_3 %>% filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>% slice_max(n = 40, order_by = logFC) %>% pull(Symbol),
  degenes_1_3 %>% filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>% slice_min(n = 40, order_by = logFC) %>% pull(Symbol),
  degenes_2_3 %>% filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>% slice_max(n = 40, order_by = logFC) %>% pull(Symbol),
  degenes_2_3 %>% filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>% slice_min(n = 40, order_by = logFC) %>% pull(Symbol),
  degenes_OF %>% filter(Sig == TRUE) %>% slice_max(n = 20, order_by = logFC) %>% pull(Symbol),
  degenes_OF %>% filter(Sig == TRUE) %>% slice_min(n = 20, order_by = logFC) %>% pull(Symbol)
))
plot.matrix <- as.matrix(PLANET_exp[genes_used,sample_used])
plot.matrix <- apply(plot.matrix, 1, scale)
plot.matrix <- t(plot.matrix)
colnames(plot.matrix) <- sample_used
plot.matrix <- pmin(plot.matrix, quantile(plot.matrix,.975))
plot.matrix <- pmax(plot.matrix, quantile(plot.matrix,.025))
color_used <- circlize::colorRamp2(seq(min(plot.matrix), max(plot.matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
genes_labeled <- unique(c(
  degenes_OF %>% filter(Sig == TRUE) %>% slice_max(n = 20, order_by = logFC) %>% pull(Symbol),
  degenes_OF %>% filter(Sig == TRUE) %>% slice_min(n = 20, order_by = logFC) %>% pull(Symbol)))
ha = HeatmapAnnotation(
  OF = PLANET_metadata[sample_used,"OF_group"],
  AFP = PLANET_metadata[sample_used,"AFPgroup"],
  RFS = PLANET_metadata[sample_used,"RFSgroup"],
  col = list(OF = c("OF_low" = "#4d6ca6","OF_high" = "#eb3223"),
             AFP = c("AFPhi"="#6d327f","AFPlow"="#bda7cc"),
             RFS = c("0-6 month"="#af4c4c","6-12 month"="#c2726d",">12 month"="#d69b95")
  )
)
p <- Heatmap(plot.matrix,
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
    at = which(genes_used %in% genes_labeled),
    labels = genes_used[which(genes_used %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

# Fig.5g----
cancer_module_gene <- read.csv("./PLANET cohort/composition_of_recurrent_gene_modules.csv", stringsAsFactors = F)
cancer_module_gene_list <- list()
for(i in colnames(cancer_module_gene)){
  cancer_module_gene_list[[i]] <- cancer_module_gene[,i][cancer_module_gene[,i] %in% row.names(PLANET_exp)]
}
for(i in names(cancer_module_gene_list)){
  PLANET_metadata[colnames(PLANET_exp),i] <- colMeans(PLANET_exp[cancer_module_gene_list[[i]],])
}
PLANET_metadata$Group = paste0(PLANET_metadata$Recurrence,"_",PLANET_metadata$OF_group)
PLANET_metadata$Group[PLANET_metadata$Group %in% c("No_OF_high","No_OF_low")] <- "No"
PLANET_metadata$Group <- factor(PLANET_metadata$Group, levels = c("No","Yes_OF_low","Yes_OF_high"))
plot_list <- list()
for(i in c("Cycle","Mesenchymal","pEMT")){
  plot_list[[i]] <- 
    PLANET_metadata[!is.na(PLANET_metadata$OF_group),c("Group",i)] %>% 
    ggplot(aes_string(x = "Group", y = i)) +
    geom_boxplot() +
    ggsignif::geom_signif(comparisons = list(c("No","Yes_OF_low"),
                                             c("Yes_OF_low","Yes_OF_high")),
                          step_increase = .1, tip_length = 0) + 
    theme_cowplot(font_size = 7)
}
plot_grid(plotlist = plot_list, nrow = 1)

HRC_gene <- read.csv("../onco_fetal2/00.raw_data/coreHRC_gene_list.csv", row.names = 1, stringsAsFactors = F)
PLANET_metadata[colnames(PLANET_exp),"epiHR"] <- colMeans(PLANET_exp[intersect(HRC_gene %>% filter(epiHR == 1) %>% row.names(),row.names(PLANET_exp)),])
PLANET_metadata[colnames(PLANET_exp),"core_HRC"] <- colMeans(PLANET_exp[intersect(HRC_gene %>% filter(core_HRC == 1) %>% row.names(),row.names(PLANET_exp)),])
plot_list <- list()
for(i in c("epiHR","core_HRC")){
  plot_list[[i]] <- 
    PLANET_metadata[!is.na(PLANET_metadata$OF_group),c("Group",i)] %>% 
    ggplot(aes_string(x = "Group", y = i)) +
    geom_boxplot() +
    ggsignif::geom_signif(comparisons = list(c("No","Yes_OF_low"),
                                             c("Yes_OF_low","Yes_OF_high")),
                          step_increase = .1, tip_length = 0) + 
    theme_cowplot(font_size = 7)
}
plot_grid(plotlist = plot_list, nrow = 1)

# Fig.5h----
HCC_metabolic <- intersect(c("SLCO1B1","CYP3A43","LECT2","CYP17A1","CYP8B1","THRSP","BHMT","ADH1B","CYP3A4","SLC25A47","CYP2E1","G6PC","F9","NR1I3","CYP2A6","FXYD1","CYP2A7","PCK1"),row.names(hcc01_bin100@assays$SCT@scale.data))
HCC_EMT <- intersect(c("CTNND2","WNK2","CLDN4","KRT80","CA9","MUC5B","MISP","SYT13","DMBT1","SLC6A19","MMP7","KRT19","DUOX2","CXCL6","TEX15","KRT7","CXCL1","MUC6","AGR2","CXCL5"),row.names(hcc01_bin100@assays$SCT@scale.data))
hcc01_bin100@meta.data$HCC_metabolic_score <- colSums(hcc01_bin100@assays$SCT@scale.data[HCC_metabolic,])
hcc01_bin100@meta.data$HCC_EMT_score <- colSums(hcc01_bin100@assays$SCT@scale.data[HCC_EMT,])
hcc03_bin100@meta.data$HCC_metabolic_score <- colSums(hcc03_bin100@assays$SCT@scale.data[HCC_metabolic,])
hcc03_bin100@meta.data$HCC_EMT_score <- colSums(hcc03_bin100@assays$SCT@scale.data[HCC_EMT,])
p1 <- SpatialFeaturePlot(hcc01_bin100, features = "HCC_metabolic_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
p2 <- SpatialFeaturePlot(hcc01_bin100, features = "HCC_EMT_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
p3 <- SpatialFeaturePlot(hcc03_bin100, features = "HCC_metabolic_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')
p4 <- SpatialFeaturePlot(hcc03_bin100, features = "HCC_EMT_score", pt.size.factor = 1.5, stroke = .2, min.cutoff = "q1", max.cutoff = 'q99')

# Fig.5i----
recur_OFhigh_genes <- intersect(degenes %>% arrange(desc(logFC)) %>% head(20) %>% pull(Symbol) %>% as.character(), row.names(hcc03_bin100@assays$SCT@scale.data))
hcc03_bin100@meta.data[,"recur_OFhigh_score"] <- colSums(hcc03_bin100@assays$SCT@scale.data[recur_OFhigh_genes,])
cutoff <- quantile(hcc03_bin100@meta.data[,"recur_OFhigh_score"], .975)
high_cells <- hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data[,"recur_OFhigh_score"] > cutoff]
hcc03_bin100@meta.data[,"recur_OFhigh"] <- FALSE
hcc03_bin100@meta.data[high_cells,"recur_OFhigh"] <- TRUE
hcc03_bin100 <- CalNeighIndex(hcc03_bin100, hcc03_bin100@meta.data$cell[hcc03_bin100@meta.data$recur_OFhigh], "recur_OFhigh")
hcc03_bin100@meta.data$OF_score <- hcc03_bin100@meta.data$OF_CAF_NeighIndex + hcc03_bin100@meta.data$OF_TAM_NeighIndex + hcc03_bin100@meta.data$OF_Endo_NeighIndex
SpatialFeaturePlot(hcc03_bin100, features = "OF_score", pt.size.factor = 1.5, stroke = 0, max.cutoff = "q95") + viridis::scale_fill_viridis()
SpatialDimPlot(hcc03_bin100, group.by = "recur_OFhigh", pt.size.factor = 1.5, stroke = 0) + scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#F84141"))

# Extended Data Fig.7d----
PLANET_Music <- readRDS("./PLANET cohort/PLANET_Music.rds")
PLANET_Music <- data.frame(PLANET_Music$Est.prop.weighted, check.names = F)
PLANET_Music[,"T cells"] <- PLANET_Music[,"CD4+ T"] + PLANET_Music[,"CD8+ T"] + PLANET_Music[,"Tregs"]
PLANET_Music[,"All Hepatocyte"] <- PLANET_Music[,"Hepatocyte_P15"] + PLANET_Music[,"Hepatocyte_P7_2"] + PLANET_Music[,"Hepatocyte"]
colnames(PLANET_Music)[colnames(PLANET_Music) == "Other Endotehlium"] <- "Other Endothelium"
cluster_used <- Cluster_levels
sample_used <- row.names(PLANET_cibersort)
similarity <- c()
for(sample in row.names(PLANET_cibersort)){
  correlation <- lsa::cosine(as.vector(t(PLANET_cibersort[sample,cluster_used])),as.vector(t(PLANET_Music[sample,cluster_used])))
  similarity <- c(similarity, correlation)
}
ggplot(data.frame(cosine = similarity), aes(y = cosine)) +
  geom_boxplot(outlier.size = 1, outlier.stroke = 0) +
  labs(y = "Cosine Similarity") +
  theme_cowplot(font_size = 7) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Extended Data Fig.7e----
PLANET_metadata_temp <- PLANET_metadata %>% filter(Tissue == "T", Bulk_sample_id %in% row.names(PLANET_cibersort))
PLANET_metadata_temp <- cbind(PLANET_metadata_temp, PLANET_cibersort[row.names(PLANET_metadata_temp),])
CAF_prop_cutoff <- quantile(PLANET_metadata_temp$`PLVAP+ EC`,c(0.5,0.5))
PLANET_metadata_temp[PLANET_metadata_temp$`PLVAP+ EC` > CAF_prop_cutoff[2],"Group2"] <- "PLVAP+ EC high"
PLANET_metadata_temp[PLANET_metadata_temp$`PLVAP+ EC` <= CAF_prop_cutoff[1],"Group2"] <- "PLVAP+ EC low"
PLANET_metadata_temp$Group2 <- factor(PLANET_metadata_temp$Group2, levels = c("PLVAP+ EC low","PLVAP+ EC high",NA))
p1 <- ggplot(PLANET_metadata_temp %>% filter(!is.na(Group2)), aes(x = Group2, y = `FOLR2+ TAM1`)) +
  geom_boxplot(aes(fill = Group2), alpha = .8) +
  theme_cowplot() +
  labs(y = "Predicted proportions of\n FOLR2+ TAM1", x = "") +
  scale_fill_manual(values = c("PLVAP+ EC high" = "#F84141", "PLVAP+ EC low" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("PLVAP+ EC high","PLVAP+ EC low"))) +
  RotatedAxis() +
  theme(legend.position = "none")
p2 <- ggplot(PLANET_metadata_temp %>% filter(!is.na(Group2)), aes(x = Group2, y = `POSTN+ CAF`)) +
  geom_boxplot(aes(fill = Group2), alpha = .8) +
  theme_cowplot() +
  labs(y = "Predicted proportions of\n POSTN+ CAF", x = "") +
  scale_fill_manual(values = c("PLVAP+ EC high" = "#F84141", "PLVAP+ EC low" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("PLVAP+ EC high","PLVAP+ EC low"))) +
  RotatedAxis() +
  theme(legend.position = "none")
p1 + p2

# Extended Data Fig.7f----
OF_sample_df <- table(PLANET_metadata[,c("sample_ID","OF_group")])
OF_sample_Simpson <- 1 - apply(OF_sample_df, 1, function(x){
  (x[1]*(x[1]-1) + x[2]*(x[2]-1))/(sum(x)*(sum(x)-1))
})
OF_sample_Simpson <- table(round(OF_sample_Simpson,2)) %>% data.frame()
ggplot(OF_sample_Simpson, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, y = Freq + 1)) +
  theme_cowplot() +
  labs(x = "Simpson's Index", y = "Number of patients") +
  RotatedAxis()

# Extended Data Fig.7g----
TCGA_metadata <- read.csv("./PLANET cohort/TCGA_LIHC_phenotype.csv", sep = "\t", row.names = 1)
TCGA_cibersort <- read.csv("./PLANET cohort/TCGA_LIHC_cibersort.txt", sep = "\t", check.names = F, row.names = 1)
sample_used <- intersect(TCGA_metadata %>% filter(sample_type == "Primary Tumor", histological_type == "Hepatocellular Carcinoma", radiation_therapy == "NO", residual_tumor == "R0") %>% row.names(), row.names(TCGA_cibersort))
TCGA_metadata_temp <- TCGA_metadata[sample_used,]
TCGA_metadata_temp$Viral.status <- TCGA_metadata_temp$viral_hepatitis_serology!=""
TCGA_metadata_temp$Recurrence <- TCGA_metadata_temp$new_neoplasm_event_type %in% c("Extrahepatic Recurrence","Intrahepatic Recurrence","Locoregional Recurrence","New Primary Tumor")
TCGA_cibersort <- TCGA_cibersort[sample_used,]
TCGA_metadata_temp$OF_score <- rowSums(TCGA_cibersort[,c("FOLR2+ TAM1","POSTN+ CAF","PLVAP+ EC")])
OF_score_cutoff <- mean(TCGA_metadata_temp$OF_score) + sd(TCGA_metadata_temp$OF_score)/2
TCGA_metadata_temp$OF_group <- "OF_low"
TCGA_metadata_temp[TCGA_metadata_temp$OF_score > OF_score_cutoff,"OF_group"] <- "OF_high"
TCGA_metadata_temp$OF_group <- factor(TCGA_metadata_temp$OF_group, levels = c("OF_low","OF_high"))
OF_Recur_ct <- table(TCGA_metadata_temp[,c("OF_group","Recurrence")])
p1 <- CrossTabPlot(OF_Recur_ct)

sample_used <- PLANET_metadata %>% filter(Bulk_sample_id %in% row.names(PLANET_cibersort), !is.na(Recurrence)) %>% pull(Bulk_sample_id)
sample_OF_table <- table(PLANET_metadata[sample_used,c("sample_ID","OF_group")])
sample_name <- row.names(sample_OF_table)
sample_OF <- sample_OF_table[,2] > 0
sample_viral_table <- table(PLANET_metadata[sample_used,c("sample_ID","Viral.status")])
sample_viral <- (sample_viral_table[,1] > 0)[sample_name]
p2 <- CrossTabPlot(table(sample_OF, sample_viral))

# Extended Data Fig.7h----
OF_score <- list(
  OF_TAM_EC_CAF = rowSums(PLANET_cibersort[sample_used,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")]),
  OF_TAM_EC = rowSums(PLANET_cibersort[sample_used,c("FOLR2+ TAM1","PLVAP+ EC")]),
  OF_TAM_CAF = rowSums(PLANET_cibersort[sample_used,c("FOLR2+ TAM1","POSTN+ CAF")]),
  OF_EC_CAF = rowSums(PLANET_cibersort[sample_used,c("POSTN+ CAF","PLVAP+ EC")]),
  OF_TAM = PLANET_cibersort[sample_used,"FOLR2+ TAM1"],
  OF_EC = PLANET_cibersort[sample_used,"PLVAP+ EC"],
  OF_CAF = PLANET_cibersort[sample_used,"POSTN+ CAF"],
  OF_MVI = PLANET_metadata[sample_used,"MVI"],
  OF_Viral = PLANET_metadata[sample_used,"Viral.status"]
)
odds_ratio <- c()
for(i in names(OF_score)){
  if(i %ni% c("OF_MVI","OF_Viral","OF_AFP")){
    OF_score_cutoff <- mean(OF_score[[i]]) + sd(OF_score[[i]])
    OF_group <- OF_score[[i]] > OF_score_cutoff
    sample_OF_table <- table(PLANET_metadata[sample_used,"sample_ID"], OF_group)
    sample_name <- row.names(sample_OF_table)
    sample_OF <- sample_OF_table[,2] > 0
    sample_recur_table <- table(PLANET_metadata[sample_used,c("sample_ID","Recurrence")])
    sample_recur <- (sample_recur_table[,2] > 0)[sample_name]
    OF_Recur_ct <- table(sample_OF, sample_recur)
  }else{
    sample_recur_table <- table(PLANET_metadata[sample_used,c("sample_ID","Recurrence")])
    sample_recur <- (sample_recur_table[,2] > 0)[sample_name]
    sample_group_table <- table(PLANET_metadata[sample_used,"sample_ID"], OF_score[[i]])
    sample_group <- (sample_group_table[,2] > 0)[sample_name]
    OF_Recur_ct <- table(sample_recur, sample_group)
  }
  or <- fisher.test(OF_Recur_ct)
  odds_ratio <- c(odds_ratio, or$estimate)
}
data.frame(OF_group = names(OF_score), odds_ratio = odds_ratio) %>%
  ggplot(aes(x = 0, y = odds_ratio)) +
  geom_point() +
  lims(y = c(0, max(odds_ratio))) +
  labs(x = "", y = "Odds Ratio") +
  geom_text_repel(aes(label = OF_group),
                  nudge_x = 0.01,
                  direction = "y",
                  hjust = 0,
                  segment.size = 0.2,
                  segment.curvature = -0.1) +
  theme_cowplot(font_size = 7)

# Extended Data Fig.7i----
relapsed_samples <- PLANET_metadata %>% filter(Bulk_sample_id %in% row.names(PLANET_cibersort), Recurrence == "Yes") %>% pull(Bulk_sample_id)
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:6])
ROIE(table(PLANET_metadata[relapsed_samples,c("TNM.Stage","OF_group")])) %>% melt() %>% 
  ggplot(data = ., aes(Var2, forcats::fct_rev(Var1), fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = myPalette(100)) +
  labs(x = "", y = "") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.line = element_blank(), 
        axis.ticks = element_blank()) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 12))

# Extended Data Fig.7j----
cbind(TNM.Stage = PLANET_metadata[row.names(PLANET_cibersort),c("TNM.Stage")],PLANET_cibersort[,c("PLVAP+ EC","POSTN+ CAF","FOLR2+ TAM1")]) %>% filter(!is.na(TNM.Stage)) %>% melt() %>%
  ggplot(aes(x = TNM.Stage, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  ggsignif::geom_signif(comparisons = list(c("Stage I","Stage II"),
                                           c("Stage II","Stage III"),
                                           c("Stage I","Stage III")),
                        step_increase = .1, textsize = 2, tip_length = 0) +
  labs(x = "TNM Stage", y = "") +
  ggtitle(label = "", subtitle = "StageI(n = 87), StageII(n=63), StaegIII(n=48)") +
  theme_cowplot(font_size = 8) +
  theme(strip.background = element_blank(),
        strip.switch.pad.grid = unit(1, "inch")) +
  RotatedAxis()

# Extended Data Fig.7k----
library(survival)
library(survminer)
# data processed in Fig.5d
day_cutoff <- 720
RFS_used <- RFS %>% filter(Old.RFS <= day_cutoff)
surv.obj <- Surv(RFS_used[,"Old.RFS"], rep(1,nrow(RFS_used)))
cox.fit <- coxph(surv.obj ~ OFgroup + Gender + Age + Viral.status, data = RFS_used)
surv.curve <- survfit(surv.obj ~ OFgroup, data = RFS_used)
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)

# Extended Data Fig.7l----
relapsed_samples <- PLANET_metadata %>% filter(Bulk_sample_id %in% row.names(PLANET_cibersort), Recurrence == "Yes") %>% pull(Bulk_sample_id)
degenes <- LIMMA(PLANET_exp[,relapsed_samples],
                 PLANET_metadata[relapsed_samples,"OF_group"])
degenes$Sig <- FALSE
degenes[degenes$adj.P.Val < 0.05 & abs(degenes$logFC) > 1,"Sig"] <- TRUE
genes_labeled <- c(
  degenes %>% filter(Sig == TRUE) %>% slice_max(n = 20, order_by = logFC) %>% pull(Symbol),
  degenes %>% filter(Sig == TRUE) %>% slice_min(n = 20, order_by = logFC) %>% pull(Symbol)
)
degenes[genes_labeled,"label"] <- genes_labeled
ggplot(degenes, aes(x = -logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = Binomial_color_panel) +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  theme_cowplot(font_size = 16) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1)) +
  geom_vline(xintercept = c(-1,1), color = "grey", linetype = "dashed", lwd = 1) +
  geom_text_repel(data = degenes, aes(label = label), size = 4)

# Extended Data Fig.7m----
genes_used <- intersect(c(
  degenes %>% filter(Sig == TRUE) %>% slice_max(n = 20, order_by = logFC) %>% pull(Symbol),
  degenes %>% filter(Sig == TRUE) %>% slice_min(n = 20, order_by = logFC) %>% pull(Symbol)
),row.names(HCC_seu@assays$RNA@data))
HCC_seu@meta.data$ClusterTemp <- HCC_seu@meta.data$Sub_Cluster
HCC_seu@meta.data$ClusterTemp[HCC_seu@meta.data$Global_Cluster != "Hepatocyte"] <- "Non Malignant"
cells_used <- HCC_seu@meta.data %>% filter(Sub_Cluster != "Bi-Potent") %>% pull(CellName)
exp_matrix <- aggregate(t(as.matrix(HCC_seu@assays$RNA@data[genes_used,cells_used])),list(HCC_seu@meta.data[cells_used,"ClusterTemp"]),mean)
row.names(exp_matrix) <- exp_matrix$Group.1
exp_matrix$Group.1 <- c()
exp_matrix_scale <- apply(exp_matrix, 2, scale) %>% t()
colnames(exp_matrix_scale) <- row.names(exp_matrix)
color_used <- circlize::colorRamp2(seq(min(exp_matrix), max(exp_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
exp_matrix_quantile <- quantile(exp_matrix_scale, c(0.001, 0.999))
exp_matrix_scale <- pmax(exp_matrix_scale, exp_matrix_quantile[1])
exp_matrix_scale <- pmin(exp_matrix_scale, exp_matrix_quantile[2])
Heatmap(t(exp_matrix_scale),
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp")

# Extended Data Fig.7n----
# Analyzed in GEPIA2 (http://gepia2.cancer-pku.cn/#index)

# Extended Data Fig.7o----
# see Fig.5g

# Extended Data Fig.7q----
up_genes <- degenes %>% filter(logFC > 0, Sig == TRUE) %>% arrange(desc(logFC)) %>% pull(Symbol) %>% as.character()
up_genes <- intersect(up_genes, row.names(HCC_seu@assays$RNA@data))
up_genes <- t(HCC_seu@assays$RNA@data[up_genes,HCC_seu$Sub_ClusterNew == "Hepatocyte_P15"]) %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 1] %>% names()
sender_genes <- unique(c(
  t(HCC_seu@assays$RNA@data[,HCC_seu$Sub_ClusterNew == "POSTN+ CAF"]) %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 3] %>% names(),
  t(HCC_seu@assays$RNA@data[,HCC_seu$Sub_ClusterNew == "FOLR2+ TAM1"]) %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 3] %>% names(),
  t(HCC_seu@assays$RNA@data[,HCC_seu$Sub_ClusterNew == "PLVAP+ EC"]) %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 3] %>% names()
))
sender_genes <- sender_genes[sender_genes != "IL1B"]
nichenet_results <-
  RunNichenetr(seu = HCC_seu,
               cluster = "Hepatocyte_P15", group.by = "Sub_ClusterNew",
               receiver_genes = up_genes, sender_genes = sender_genes,
               expressed_gene_cutoff = 1,
               best_upstream_ligands_number = 50, curated = T,
               species = "homo", repo = TRUE)
genes_used <- rev(as.character(nichenet_results$p_ligand_pearson$data$y))
# genes_used <- c("FN1","VWF","COL4A1","CTGF","VCAM1","SPP1","COMP","ITGB2","ITGB1","NAMPT","INHBB","INHBA","HSPG2","COL18A1","COL1A1","HBEGF","GSTP1","JAG1","THBS2","CYR61","LAMB1","HGF","AREG","NUCB2","DLL4","THBS1","TNFSF10","NID1","ADAM15","SLIT2")
cells_used <- HCC_seu@meta.data %>% filter(Sub_ClusterNew %in% c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")) %>% row.names()
gene.mean.matrix <- 
  aggregate(t(as.matrix(HCC_seu@assays$RNA@data[genes_used,cells_used])),
            list(Cluster = HCC_seu@meta.data[cells_used,"Sub_ClusterNew"]), 
            mean)
gene.mean.df <- gene.mean.matrix %>% melt(id.vars = "Cluster", variable.name = "Gene", value.name = "Exp")
gene.mean.df$Exp <- pmin(gene.mean.df$Exp,quantile(gene.mean.df$Exp,.99))
gene.per <- 
  aggregate(t(as.matrix(HCC_seu@assays$RNA@data[genes_used,cells_used])),
            list(Group = HCC_seu@meta.data[cells_used,"Sub_ClusterNew"]),
            function(x){sum(x > 2) / length(x)}) %>% 
  melt(id.vars = "Group", variable.name = "Gene", value.name = "Per")
plot.data <- merge(gene.mean.df, gene.per)
plot.data$Gene <- as.character(plot.data$Gene)
plot.data$Group <- factor(plot.data$Group, levels = rev(c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")))
plot.data$Gene <- factor(plot.data$Gene, levels = rev(genes_used))
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplot(plot.data, aes(x = Group, y = Gene)) +
  geom_point(aes(size = Per, fill = Exp, color = Exp), shape = 21) +
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
  scale_fill_gradientn("log2(TPM + 1)",
                       colours = myColorPalette(100),
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black", barwidth = 0.8)) +
  scale_color_gradientn("log2(TPM + 1)",
                        colours = myColorPalette(100),
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks.colour = "black", barwidth = 0.8)) +
  scale_size_continuous("Percentage",
                        breaks = seq(0, 0.8, 0.2), range = c(1,3))
nichenet_results$p_ligand_pearson + theme(legend.position = "right")

# Extended Data Fig.7s----
# see Fig.5i for some pre-step calculation
hcc03_recur_dist_OF <- 
  DistFeature(seu = hcc03_bin100, 
              cells = hcc03_bin100@meta.data %>% filter(recur_OFhigh) %>% row.names(),
              feature = "OF_score",
              max.dist = 20)
hcc03_recur_dist_OF_random_list <- list()
for(i in 1:10){
  hcc03_recur_dist_OF_random_list[[i]] <- 
    DistFeature(seu = hcc03_bin100, 
                cells = sample(row.names(hcc03_bin100@meta.data), sum(hcc03_bin100@meta.data$recur_OFhigh)),
                feature = "OF_score",
                max.dist = 20)
}
hcc03_dist_OF_random %>% melt() %>% 
  ggplot(aes(x = factor(Var2), y = value)) + 
  geom_boxplot(color = "#9FA0A3", outlier.stroke = 0) + 
  geom_point(data = data.frame(value = colMeans(hcc03_recur_dist_OF$dist.feature),
                               Vars = factor(1:20)),
             aes(x = Vars, y = value), color = "#F84141") +
  labs(x = "Distance", y = "Proportions of Onco-fetal High Bins") +
  theme_cowplot(font_size = 7, line_size = .25)