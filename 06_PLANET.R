setwd("/work/lzy/project/onco_fetal/")
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

# >Generate CIBERSORTx training data----
# cells <- HCC_seu@meta.data %>% filter(Global_Cluster != "Doublet", NTF == "Tumor") %>% pull(CellName) %>% sample(30000)
# genes_used <- rowSums(HCC_seu@assays$RNA@counts[,cells]) > 100
# exp <- 2^HCC_seu@assays$RNA@counts[genes_used, cells] - 1
# exp <- apply(exp,2,function(x){signif(x,digits = 3)})
# training_exp <- exp
# colnames(training_exp) <- HCC_seu@meta.data[colnames(training_exp),"Sub_ClusterNew"]
# training_exp <- data.table(Gene = row.names(training_exp), training_exp, check.names = F)
# write_tsv(training_exp, file = "/work/lzy/project/onco_fetal/00.data/Nanostring/Input/Tumor_training2.txt", quote = F, col_names = T)
# 
# cells_used <- gsub("fetal_","",Fib_seu@meta.data %>% filter(Sub_Cluster == "ITM2C+ Fib") %>% pull(CellID))
# fetal_seu@meta.data[cells_used,"cellType"] <- "ITM2C+ Fibroblast"
# cells <- unique(c(fetal_seu@meta.data %>% row.names() %>% sample(30000),cells_used))
# cells <- cells[fetal_seu@meta.data[cells,"cellType"] != "C0"]
# genes_used <- rowSums(fetal_seu@assays$RNA@counts[,cells]) > 100
# exp <- 2^fetal_seu@assays$RNA@counts[genes_used, cells] - 1
# exp <- apply(exp,2,function(x){signif(x,digits = 3)})
# training_exp <- exp
# colnames(training_exp) <- fetal_seu@meta.data[colnames(training_exp),"cellType"]
# training_exp <- data.table(Gene = row.names(training_exp), training_exp, check.names = F)
# write_tsv(training_exp, file = "/work/lzy/project/onco_fetal/00.data/Nanostring/Input/Fetal_training.txt", quote = F, col_names = T)

# Ankur bulk----
Ankur_exp <- read.table("./00.data/HCC_bulk/Ankur_HCC_bulk.txt", sep = "\t", head = T)
Ankur_exp <- Ankur_exp[!duplicated(Ankur_exp$Gene),]
row.names(Ankur_exp) <- Ankur_exp$Gene
Ankur_exp$Gene <- c()
Ankur_sampleinfo <- read.csv("./00.data/HCC_bulk/sample_info.csv", row.names = 1, stringsAsFactors = F)
Ankur_clinical <- readr::read_tsv("./00.data/HCC_bulk/ClinicalData.tsv") %>% data.frame()
Ankur_metadata <- merge(Ankur_sampleinfo,Ankur_clinical,by.x="sample_id",by.y="Unified_ID")
row.names(Ankur_metadata) <- Ankur_metadata$Bulk_sample_id
Ankur_metadata <- Ankur_metadata[colnames(Ankur_exp),]
Ankur_metadata$Tumor <- stringr::str_split_fixed(Ankur_metadata$Bulk_sample_id,"[.]",3)[,3]
Ankur_metadata$TNM.Stage <- Ankur_metadata$TNM.Stage.V8
Ankur_metadata$TNM.Stage[Ankur_metadata$TNM.Stage.V8 %in% c("TNM Stage IA","TNM Stage IB")] <- "Stage I"
Ankur_metadata$TNM.Stage[Ankur_metadata$TNM.Stage.V8 %in% c("TNM Stage II")] <- "Stage II"
Ankur_metadata$TNM.Stage[Ankur_metadata$TNM.Stage.V8 %in% c("TNM Stage IIIA","TNM Stage IIIB")] <- "Stage III"
Ankur_metadata[,"MVI"][Ankur_metadata[,"MVI"] == "no"] <- "No"
Ankur_metadata[,"MVI"][Ankur_metadata[,"MVI"] == "yes"] <- "Yes"
Ankur_metadata$RFSgroup <- NA
Ankur_metadata[!is.na(Ankur_metadata$RFS) & Ankur_metadata$Recurrence == "Yes" & Ankur_metadata$RFS <= 180,"RFSgroup"] <- "0-6 month"
Ankur_metadata[!is.na(Ankur_metadata$RFS) & Ankur_metadata$Recurrence == "Yes" & Ankur_metadata$RFS > 180 & Ankur_metadata$RFS <= 360,"RFSgroup"] <- "6-12 month"
Ankur_metadata[!is.na(Ankur_metadata$RFS) & Ankur_metadata$Recurrence == "Yes" & Ankur_metadata$RFS > 360,"RFSgroup"] <- ">12 month"
Ankur_metadata$RFSgroup <- factor(Ankur_metadata$RFSgroup, levels = c("0-6 month","6-12 month",">12 month"))
Ankur_cibersort <- read.table("./00.data/HCC_bulk/CIBERSORTx/Output/CIBERSORTx_Results.txt", sep = "\t", head = T, check.names = F)
row.names(Ankur_cibersort) <- Ankur_cibersort$Mixture
Ankur_cibersort$Mixture <- c()
Ankur_cibersort[,"T cells"] <- Ankur_cibersort[,"CD4+ T"] + Ankur_cibersort[,"CD8+ T"] + Ankur_cibersort[,"Tregs"]
Ankur_cibersort[,"All Hepatocyte"] <- Ankur_cibersort[,"Hepatocyte_P15"] + Ankur_cibersort[,"Hepatocyte_P7_2"] + Ankur_cibersort[,"Hepatocyte"]
Ankur_cibersort$`P-value` <- Ankur_cibersort$Correlation <- Ankur_cibersort$RMSE <- c()
colnames(Ankur_cibersort)[colnames(Ankur_cibersort) == "Other Endotehlium"] <- "Other Endothelium"
Ankur_variant <- read.csv("./00.data/HCC_bulk/protein_change.csv")
row.names(Ankur_variant) <- Ankur_variant$Sample
Ankur_samples <- read.csv("./00.data/HCC_bulk/masterfile_all_dna_samples.csv")
Ankur_samples <- merge(Ankur_samples, Ankur_metadata[,c("sample_id","Bulk_sample_id","Tumor")],by.x=c("Unified_ID","Tumor"),by.y=c("sample_id","Tumor"))
row.names(Ankur_samples) <- Ankur_samples$DNA_lib
Ankur_samples <- Ankur_samples[intersect(Ankur_samples$DNA_lib ,Ankur_variant$Sample),]
Ankur_variant <- Ankur_variant[intersect(Ankur_samples$DNA_lib ,Ankur_variant$Sample),]
Ankur_variant$Sample <- c()
row.names(Ankur_variant) <- as.character(Ankur_samples[row.names(Ankur_variant),"Bulk_sample_id"])

# Music deconvoluted results----
Ankur_Music <- readRDS("./00.data/HCC_bulk/CIBERSORTx/Output/Music_deconvolution.rds")
Ankur_Music <- data.frame(Ankur_Music$Est.prop.weighted, check.names = F)
Ankur_Music[,"T cells"] <- Ankur_Music[,"CD4+ T"] + Ankur_Music[,"CD8+ T"] + Ankur_Music[,"Tregs"]
Ankur_Music[,"All Hepatocyte"] <- Ankur_Music[,"Hepatocyte_P15"] + Ankur_Music[,"Hepatocyte_P7_2"] + Ankur_Music[,"Hepatocyte"]
colnames(Ankur_Music)[colnames(Ankur_Music) == "Other Endotehlium"] <- "Other Endothelium"
cluster_used <- Cluster_levels
sample_used <- row.names(Ankur_cibersort)
similarity <- c()
for(sample in row.names(Ankur_cibersort)){
  correlation <- lsa::cosine(as.vector(t(Ankur_cibersort[sample,cluster_used])),as.vector(t(Ankur_Music[sample,cluster_used])))
  similarity <- c(similarity, correlation)
}
p <- ggplot(data.frame(cosine = similarity), aes(y = cosine)) +
  geom_boxplot(outlier.size = 1, outlier.stroke = 0) +
  labs(y = "Cosine Similarity") +
  theme_cowplot(font_size = 7) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(p, file = "./04.figures/06.HCC_bulk_predicted_proportions_cosine_similarity.pdf", width = 0.6, height = 2)

# Define Onco-fetal score----
Ankur_metadata[row.names(Ankur_cibersort),"OF_score"] <- rowSums(Ankur_cibersort[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")])
sample_used <- Ankur_metadata %>% filter(Bulk_sample_id %in% row.names(Ankur_cibersort), !is.na(Recurrence)) %>% pull(Bulk_sample_id)
OF_score_cutoff <- mean(Ankur_metadata[sample_used,"OF_score"]) + sd(Ankur_metadata[sample_used,"OF_score"])
Ankur_metadata[sample_used,"OF_group"] <- "OF_low"
Ankur_metadata[sample_used,"OF_group"][Ankur_metadata[sample_used,"OF_score"] > OF_score_cutoff] <- "OF_high"
Ankur_metadata$OF_group <- factor(Ankur_metadata$OF_group, levels = c("OF_low","OF_high"))

# Cluster tumor samples according to proportions----
hc <- hclust(dist(Ankur_cibersort[,Cluster_levels], method = "euclidean"), method = "ward.D2")
dend.hc <- as.dendrogram(hc)
Recurrence_color_panel <- c("Yes" = "#F3746C", "No" = "#69B4CE")
pdf("./04.figures/06.HCC_bulk_tumor_dendogram.pdf", width = 6, height = 3)
dend.hc %>% set("labels_col", Recurrence_color_panel[Ankur_metadata[labels(dend.hc),"Recurrence"]]) %>% plot(center = T)
dev.off()

plot.df <- Ankur_cibersort[,Cluster_levels] %>% mutate(Sample = row.names(.)) %>% melt() %>% mutate(Sample = factor(Sample, levels = labels(dend.hc)), variable = factor(variable, levels = Cluster_levels))
p <- ggplot(plot.df, aes(x = Sample, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", lwd = .01) +
  labs(y = "Predicted proportions", x = "") +
  scale_fill_manual(name = "", values = Cluster_color_panel) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(p, file = "./04.figures/06.HCC_bulk_tumor_proportions_barplot.pdf", width = 8, height = 4)

ha = HeatmapAnnotation(
  OF_score = anno_barplot(Ankur_metadata[labels(dend.hc),"OF_score"],
                          baseline = OF_score_cutoff, border = F,
                          gp = gpar(fill = "#808080", col = "#808080")),
  Recurrence = Ankur_metadata[labels(dend.hc),"Recurrence"],
  RFS = Ankur_metadata[labels(dend.hc),"RFSgroup"],
  Viral.Status = Ankur_metadata[labels(dend.hc),"Viral.status"],
  TNM.Stage = Ankur_metadata[labels(dend.hc),"TNM.Stage"],
  col = list(Recurrence = c("Yes"="#F97137","No"="#4FB0AF"),
             RFS = c("0-6 month"="#d69b95","6-12 month"="#c2726d",">12 month"="#af4c4c"),
             TNM.Stage = c("Stage I"="#D1CEE8","Stage II"="#A6A4D4","Stage III"="#7575BC"),
             Viral.Status = c("Hep B carrier"="#E05D89","NBNC"="#3BAFE0")
  )
)
pdf("./04.figures/06.HCC_bulk_sample_info_heatmap.pdf", width = 8, height = 4)
Heatmap(matrix = matrix(rnorm(1980), ncol = 198), 
        cluster_rows = F, cluster_columns = F,
        top_annotation = ha)
dev.off()

# Existence of onco-fetal niche----
Ankur_metadata_temp <- Ankur_metadata %>% filter(Tissue == "T", Bulk_sample_id %in% row.names(Ankur_cibersort))
Ankur_metadata_temp <- cbind(Ankur_metadata_temp, Ankur_cibersort[row.names(Ankur_metadata_temp),])
CAF_prop_cutoff <- quantile(Ankur_metadata_temp$`PLVAP+ EC`,c(0.5,0.5))
Ankur_metadata_temp[Ankur_metadata_temp$`PLVAP+ EC` > CAF_prop_cutoff[2],"Group2"] <- "PLVAP+ EC high"
Ankur_metadata_temp[Ankur_metadata_temp$`PLVAP+ EC` <= CAF_prop_cutoff[1],"Group2"] <- "PLVAP+ EC low"
Ankur_metadata_temp$Group2 <- factor(Ankur_metadata_temp$Group2, levels = c("PLVAP+ EC low","PLVAP+ EC high",NA))
p1 <- ggplot(Ankur_metadata_temp %>% filter(!is.na(Group2)), aes(x = Group2, y = `POSTN+ CAF`)) +
  geom_boxplot(aes(fill = Group2), alpha = .8) +
  theme_cowplot() +
  labs(y = "Predicted proportions of\n POSTN+ CAF", x = "") +
  scale_fill_manual(values = c("PLVAP+ EC high" = "#F84141", "PLVAP+ EC low" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("PLVAP+ EC high","PLVAP+ EC low"))) +
  RotatedAxis() +
  theme(legend.position = "none")
p2 <- ggplot(Ankur_metadata_temp %>% filter(!is.na(Group2)), aes(x = Group2, y = `FOLR2+ TAM1`)) +
  geom_boxplot(aes(fill = Group2), alpha = .8) +
  theme_cowplot() +
  labs(y = "Predicted proportions of\n FOLR2+ TAM1", x = "") +
  scale_fill_manual(values = c("PLVAP+ EC high" = "#F84141", "PLVAP+ EC low" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("PLVAP+ EC high","PLVAP+ EC low"))) +
  RotatedAxis() +
  theme(legend.position = "none")
p <- p1 + p2
ggsave(p, file = "./04.figures/06.HCC_bulk_onco-fetal_proportion_EChighlow_boxplot.pdf", width = 4, height = 5)

# Fisher.exact test of onco-fetal enriched patients and recurrence----
OF_Recur_ct <- table(Ankur_metadata[sample_used,c("OF_group","Recurrence")])
p1 <- CrossTabPlot(OF_Recur_ct)
OF_Viral_ct <- table(Ankur_metadata[sample_used,c("OF_group","Viral.status")])
p2 <- CrossTabPlot(OF_Viral_ct)
OF_MVI_ct <- table(Ankur_metadata[sample_used,c("OF_group","MVI")])
p3 <- CrossTabPlot(OF_MVI_ct)
p <- p1 + p2 + p3
ggsave(p, file = "./04.figures/06.HCC_bulk_Recurrence_fisher_test.pdf", width = 9, height = 3)

sample_used <- Ankur_metadata %>% filter(Bulk_sample_id %in% row.names(Ankur_cibersort), !is.na(Recurrence)) %>% pull(Bulk_sample_id)
sample_OF_table <- table(Ankur_metadata[sample_used,c("sample_ID","OF_group")])
sample_name <- row.names(sample_OF_table)
sample_OF <- sample_OF_table[,2] > 0
sample_recur_table <- table(Ankur_metadata[sample_used,c("sample_ID","Recurrence")])
sample_recur <- (sample_recur_table[,2] > 0)[sample_name]
sample_MVI_table <- table(Ankur_metadata[sample_used,c("sample_ID","MVI")])
sample_MVI <- (sample_MVI_table[,2] > 0)[sample_name]
sample_viral_table <- table(Ankur_metadata[sample_used,c("sample_ID","Viral.status")])
sample_viral <- (sample_viral_table[,1] > 0)[sample_name]
p1 <- CrossTabPlot(table(sample_OF, sample_recur))
p2 <- CrossTabPlot(table(sample_OF, sample_MVI))
p3 <- CrossTabPlot(table(sample_OF, sample_viral))
p <- p1 + p2 + p3
ggsave(p, file = "./04.figures/06.HCC_bulk_Recurrence_fisher_test_sample_level.pdf", width = 9, height = 3)

OF_score <- list(
  OF_TAM_EC_CAF = rowSums(Ankur_cibersort[sample_used,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")]),
  OF_TAM_EC = rowSums(Ankur_cibersort[sample_used,c("FOLR2+ TAM1","PLVAP+ EC")]),
  OF_TAM_CAF = rowSums(Ankur_cibersort[sample_used,c("FOLR2+ TAM1","POSTN+ CAF")]),
  OF_EC_CAF = rowSums(Ankur_cibersort[sample_used,c("POSTN+ CAF","PLVAP+ EC")]),
  OF_TAM = Ankur_cibersort[sample_used,"FOLR2+ TAM1"],
  OF_EC = Ankur_cibersort[sample_used,"PLVAP+ EC"],
  OF_CAF = Ankur_cibersort[sample_used,"POSTN+ CAF"],
  OF_MVI = Ankur_metadata[sample_used,"MVI"],
  OF_Viral = Ankur_metadata[sample_used,"Viral.status"]
)
odds_ratio <- c()
for(i in names(OF_score)){
  if(i %ni% c("OF_MVI","OF_Viral")){
    OF_score_cutoff <- mean(OF_score[[i]]) + sd(OF_score[[i]])
    OF_group <- OF_score[[i]] > OF_score_cutoff
    sample_OF_table <- table(Ankur_metadata[sample_used,"sample_ID"], OF_group)
    sample_name <- row.names(sample_OF_table)
    sample_OF <- sample_OF_table[,2] > 0
    sample_recur_table <- table(Ankur_metadata[sample_used,c("sample_ID","Recurrence")])
    sample_recur <- (sample_recur_table[,2] > 0)[sample_name]
    OF_Recur_ct <- table(sample_OF, sample_recur)
  }else{
    sample_recur_table <- table(Ankur_metadata[sample_used,c("sample_ID","Recurrence")])
    sample_recur <- (sample_recur_table[,2] > 0)[sample_name]
    sample_group_table <- table(Ankur_metadata[sample_used,"sample_ID"], OF_score[[i]])
    sample_group <- (sample_group_table[,2] > 0)[sample_name]
    OF_Recur_ct <- table(sample_recur, sample_group)
  }
  or <- fisher.test(OF_Recur_ct)
  odds_ratio <- c(odds_ratio, or$estimate)
}
# for(i in names(OF_score)){
#   OF_score_cutoff <- mean(OF_score[[i]]) + sd(OF_score[[i]])
#   OF_group <- rep("OF_low",length(sample_used))
#   OF_group[OF_score[[i]] > OF_score_cutoff] <- "OF_high"
#   OF_group <- factor(OF_group, levels = c("OF_low","OF_high"))
#   OF_Recur_ct <- table(Ankur_metadata[sample_used,"Recurrence"], OF_group)
#   or <- fisher.test(OF_Recur_ct)
#   odds_ratio <- c(odds_ratio, or$estimate)
# }
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
  theme_cowplot(font_size = 7) -> p
ggsave(p, file = "./04.figures/06.HCC_bulk_Recurrence_odds_ratio_sample_level.pdf", width = 4, height = 3)

# Consistency of onco-fetal group (shannon's index)----
OF_sample_df <- table(Ankur_metadata[,c("sample_ID","OF_group")])
OF_sample_Simpson <- 1 - apply(OF_sample_df, 1, function(x){
  (x[1]*(x[1]-1) + x[2]*(x[2]-1))/(sum(x)*(sum(x)-1))
  })
OF_sample_Simpson <- table(round(OF_sample_Simpson,2)) %>% data.frame()
p <- ggplot(OF_sample_Simpson, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, y = Freq + 1)) +
  theme_cowplot() +
  labs(x = "Simpson's Index", y = "Number of patients") +
  RotatedAxis()
ggsave(p, file = "./04.figures/06.HCC_bulk_Recurrence_patients_consistency.pdf", width = 2, height = 3)

# RFS and onco-fetal group----
library(survival)
library(survminer)
RFS <- read.csv("./00.data/HCC_bulk/67 patients with updated RFS.csv", stringsAsFactors = F)
row.names(RFS) <- RFS$sample_ID
RFS["16-003799","Old.RFS"] <- 46
RFS$OFgroup <- sample_OF[row.names(RFS)]
RFS$OFgroup[RFS$OFgroup == TRUE] <- "OF high"
RFS$OFgroup[RFS$OFgroup == FALSE] <- "OF low"
RFS$RFSgroup <- NA
RFS[!is.na(RFS$Old.RFS) & RFS$Old.Recurrence == "Yes" & RFS$Old.RFS <= 180,"RFSgroup"] <- "0-6 month"
RFS[!is.na(RFS$Old.RFS) & RFS$Old.Recurrence == "Yes" & RFS$Old.RFS > 180 & RFS$Old.RFS <= 360,"RFSgroup"] <- "6-12 month"
RFS[!is.na(RFS$Old.RFS) & RFS$Old.Recurrence == "Yes" & RFS$Old.RFS > 360,"RFSgroup"] <- ">12 month"
RFS$RFSgroup <- factor(RFS$RFSgroup, levels = c("0-6 month","6-12 month",">12 month"))
RFS$Viral.status <- plyr::revalue(
  RFS$Viral.status,
  replace = c("Hep B carrier" = "Carrier",
              "NBNC" = "NBNC",
              "Hep C carrier" = "Carrier")
)
RFS <- RFS %>% filter(!is.na(OFgroup), Old.Recurrence == "Yes")
day_cutoff <- 720
RFS_used <- RFS %>% filter(Old.RFS <= day_cutoff)
surv.obj <- Surv(RFS_used[,"Old.RFS"], rep(1,nrow(RFS_used)))
cox.fit <- coxph(surv.obj ~ OFgroup + Gender + Age + Viral.status, data = RFS_used)
surv.curve <- 
  survfit(surv.obj ~ OFgroup, data = RFS_used)
p <- ggsurvplot(
  surv.curve,  surv.median.line = "hv",
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 4),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 3)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)
ggsave(p$plot, file = "./04.figures/06.HCC_bulk_OF_high_vs_low_RFS_survival.pdf", width = 8, height = 5)

cox.fit <- coxph(surv.obj ~ Viral.status + Gender + Age , data = RFS_used)
surv.curve <- 
  survfit(surv.obj ~ Viral.status, data = RFS_used)
ggsurvplot(
  surv.curve,
  legend.title = "", font.legend = 16, legend = "right",
  pval = paste0(
    "P = ", round(summary(cox.fit)$coefficients[1, 5], 2),
    "\nHR = ", round(summary(cox.fit)$coefficients[1, 2], 2)),
  risk.table = TRUE, tables.theme = theme_cleantable()
)

RFS %>% ggplot(aes(x = RFSgroup, fill = OFgroup)) +
  geom_bar(position = "fill") +
  theme_cowplot() +
  labs(x = "", y = "") +
  RotatedAxis()

p1 <- RFS %>% ggplot(aes(x = OFgroup, fill = RFSgroup)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("0-6 month"="#d69b95","6-12 month"="#c2726d",">12 month"="#af4c4c")) +
  theme_cowplot() +
  labs(x = "", y = "") +
  RotatedAxis()
ggsave(p1, file = "./04.figures/06.HCC_bulk_OF_high_vs_low_RFS_barplot.pdf", width = 3, height = 6)

p2 <- ROIE_plot(ROIE(table(RFS$RFSgroup,RFS$OFgroup), filter = 0))
ggsave(p2, file = "./04.figures/06.HCC_bulk_OF_high_vs_low_RFS_ROIE.pdf", width = 4, height = 4)

CrossTabPlot(table(RFS$OFgroup, 
                   RFS$RFSgroup %>% 
                     plyr::revalue(replace = c("0-6 month" = "0-12 month",
                                               "6-12 month" = "0-12 month",
                                               ">12 month"= ">12 month"))
                   )) +
  labs(x = "", y = "")

# Onco-fetal niche enriched in Recurrence samples----
p1 <- cbind(TNM.Stage = Ankur_metadata[row.names(Ankur_cibersort),c("TNM.Stage")],Ankur_cibersort[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","B","T cells")]) %>% filter(!is.na(TNM.Stage)) %>% melt() %>%
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
ggsave(p1, file = "./04.figures/06.HCC_bulk_TNM_stage_boxplot.pdf", width = 5, height = 3.5)

p2 <- cbind(Recurrence = Ankur_metadata[row.names(Ankur_cibersort),c("Recurrence")],Ankur_cibersort[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","B","T cells")]) %>% filter(!is.na(Recurrence)) %>% melt() %>%
  ggplot(aes(x = Recurrence, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  ggsignif::geom_signif(comparisons = list(c("Yes","No"))) +
  labs(x = "Recurrence", y = "") +
  ggtitle(label = "", subtitle = "No(n=71), Yes(n=124)") +
  theme_cowplot(font_size = 8) +
  theme(strip.background = element_blank(),
        strip.switch.pad.grid = unit(1, "inch"))
ggsave(p2, file = "./04.figures/06.HCC_bulk_Recurrence_boxplot.pdf", width = 5, height = 3)

p3 <- Ankur_metadata[row.names(Ankur_cibersort),c("Recurrence","OF_score")] %>% filter(!is.na(Recurrence)) %>%
  ggplot(aes(x = Recurrence, y = OF_score)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("Yes","No"))) +
  labs(x = "Recurrence", y = "Onco-fetal Score") +
  ggtitle(label = "", subtitle = "No(n=71), Yes(n=124)") +
  theme_cowplot(font_size = 8) +
  theme(strip.background = element_blank(),
        strip.switch.pad.grid = unit(1, "inch"))
ggsave(p3, file = "./04.figures/06.HCC_bulk_Recurrence_onco-fetal_score_boxplot.pdf", width = 1, height = 3)

# Onco-fetal in TCGA----
TCGA_metadata <- read.csv("./00.data/TCGA_LIHC/TCGA.LIHC.phenotype.csv", sep = "\t", row.names = 1)
TCGA_cibersort <- read.csv("./00.data/TCGA_LIHC/CIBERSORTx/Output/CIBERSORTx_Results.txt", sep = "\t", check.names = F, row.names = 1)
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
p <- CrossTabPlot(OF_Recur_ct)
ggsave(p, file = "./04.figures/06.HCC_bulk_TCGA_Recurrence_fisher_test.pdf", width = 3, height = 3)

# Differences between OF-high and OF_low samples in relapsed patients----
# DEGenes
relapsed_samples <- Ankur_metadata %>% filter(Bulk_sample_id %in% row.names(Ankur_cibersort), Recurrence == "Yes") %>% pull(Bulk_sample_id)
degenes <- LIMMA(Ankur_exp[,relapsed_samples],
                 Ankur_metadata[relapsed_samples,"OF_group"])
degenes$Sig <- FALSE
degenes[degenes$adj.P.Val < 0.05 & abs(degenes$logFC) > 1,"Sig"] <- TRUE
genes_labeled <- c(
  degenes %>% filter(Sig == TRUE) %>% slice_max(n = 20, order_by = logFC) %>% pull(Symbol),
  degenes %>% filter(Sig == TRUE) %>% slice_min(n = 20, order_by = logFC) %>% pull(Symbol)
)
degenes[genes_labeled,"label"] <- genes_labeled
write.csv(degenes, file = "./03.results/DEGenes/HCC_bulk_Recurrence_OF_high_vs_OF_low.csv")
volcano_plot <- 
  ggplot(degenes, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = Binomial_color_panel) +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  theme_cowplot(font_size = 16) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))
legend <- get_legend(volcano_plot)
volcano_plot <- 
  ggAIplot(volcano_plot + theme(legend.position = "none")) +
  geom_vline(xintercept = c(-1,1), color = "grey", linetype = "dashed", lwd = 1) +
  geom_text_repel(data = degenes, aes(label = label), size = 4)
volcano_plot <- plot_grid(volcano_plot, legend, rel_widths = c(4,1))
ggsave(volcano_plot, file = "./04.figures/06.HCC_bulk_OF_high_vs_low_in_Recurrence_DEGenes_volcano.pdf", width = 7.75, height = 6)

# Expression of degenes in single cell dataset
genes_used <- intersect(c(
  degenes %>% filter(Sig == TRUE) %>% slice_max(n = 20, order_by = logFC) %>% pull(Symbol),
  degenes %>% filter(Sig == TRUE) %>% slice_min(n = 20, order_by = logFC) %>% pull(Symbol)
),row.names(HCC_seu@assays$RNA@data))
HCC_seu@meta.data$ClusterTemp <- HCC_seu@meta.data$Sub_Cluster
HCC_seu@meta.data$ClusterTemp[HCC_seu@meta.data$Global_Cluster != "Hepatocyte"] <- "Non Malignant"
cells_used <- HCC_seu@meta.data %>% filter(Sub_Cluster != "Bi-Potent") %>% pull(CellName)
exp_matrix <- aggregate(t(HCC_seu@assays$RNA@data[genes_used,cells_used]),list(HCC_seu@meta.data[cells_used,"ClusterTemp"]),mean)
row.names(exp_matrix) <- exp_matrix$Group.1
exp_matrix$Group.1 <- c()
exp_matrix <- apply(exp_matrix, 2, zscore) %>% t()
color_used <- circlize::colorRamp2(seq(min(exp_matrix), max(exp_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
exp_matrix_quantile <- quantile(exp_matrix, c(0.001, 0.999))
exp_matrix <- pmax(exp_matrix, exp_matrix_quantile[1])
exp_matrix <- pmin(exp_matrix, exp_matrix_quantile[2])
p <- Heatmap(t(exp_matrix),
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp")
pdf("./04.figures/06.HCC_bulk_OF_high_vs_low_in_Recurrence_DEGenes_heatmap.pdf", width = 4.5, height = 1.5)
draw(p)
dev.off()

# Pathway enrichment
library(richR)
hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
OF_high_gene_name <- degenes %>% filter(Sig == TRUE, logFC > 0) %>% pull(Symbol)
OF_high_go <- richGO(OF_high_gene_name,godata = hsago, ontology = c("BP"))
OF_low_gene_name <- degenes %>% filter(Sig == TRUE, logFC < 0) %>% pull(Symbol)
OF_low_go <- richGO(OF_low_gene_name,godata = hsago, ontology = c("BP"))
save(OF_high_go, OF_low_go, file = "./03.results/DEGenes/HCC_bulk_Recurrence_OF_high_vs_OF_low_pathway.rda")
OF_high_go_show <- c(3,11,37,47,143,165,216,262,346,467,616)
OF_low_go_show <- c(18,19,27,32,61,63,68,75,101,291)
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
p <- ggplot(pathway_plot.df, aes(x = reorder(pathway, pvalue), y = pvalue)) +
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
ggsave(p, file = "./04.figures/06.HCC_bulk_OF_high_vs_OF_low_in_Recurrence_pathway.pdf", width = 4, height = 3)

# Nichenet analyasis
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
save(nichenet_results, file = "./03.results/NicheNet/Bulk_OF_high_upregulated_nichenet.rda")

genes_used <- rev(as.character(nichenet_results$p_ligand_pearson$data$y))
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
p <- ggplot(plot.data, aes(x = Group, y = Gene)) +
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
ggsave(p, file = "./04.figures/06.HCC_bulk_OF_high_vs_OF_low_in_Recurrence_nichenet_bubble.pdf", width = 2.25, height = 4.5)

p <- nichenet_results$p_ligand_pearson + theme(legend.position = "right")
ggsave(p, file = "./04.figures/06.HCC_bulk_OF_high_vs_OF_low_in_Recurrence_nichenet_pearson.pdf", width = 3.5, height = 4.5)

# Depletion of onco-fetal high relapsed samples from early stage
p <- ROIE_plot(ROIE(table(Ankur_metadata[relapsed_samples,c("TNM.Stage","OF_group")])))
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:6])
p <- 
  ROIE(table(Ankur_metadata[relapsed_samples,c("TNM.Stage","OF_group")])) %>% melt() %>% 
  ggplot(data = ., aes(Var2, forcats::fct_rev(Var1), fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = myPalette(100)) +
  labs(x = "", y = "") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.line = element_blank(), 
        axis.ticks = element_blank()) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 12))
ggsave(p, file = "./04.figures/06.HCC_bulk_OF_high_vs_OF_low_in_Recurrence_ROIE.pdf", width = 4, height = 4)

# Recurrence related cell-cell communication----
HCC_seu <- readRDS("./02.processed_data/HCC_seu.rds")
HCC_seu@meta.data$Recurrence <- "Yes"
HCC_seu@meta.data$Recurrence[HCC_seu@meta.data$PatientID %in% c("P3","P6","P7","P8","P11","P12","P13")] <- "No"
HCC_seu@meta.data$Sub_ClusterNew <- factor(HCC_seu@meta.data$Sub_ClusterNew, levels = names(Cluster_color_panel_sc))

ncells_cluster <- data.frame(HCC_seu@meta.data %>% filter(NTF == "Tumor", Sub_ClusterNew != "Bi-Potent") %>% select(Recurrence, Sub_ClusterNew) %>% table())
ncells_recur <- data.frame(HCC_seu@meta.data %>% filter(NTF == "Tumor", Sub_ClusterNew != "Bi-Potent") %>% select(Recurrence) %>% table())
colnames(ncells_recur) <- c("Recurrence","Freq")
ncells_cluster <- merge(ncells_cluster, ncells_recur, by = c("Recurrence"))
colnames(ncells_cluster)[c(3,4)] <- c("ncells_cluster","ncells_recur")
ncells_cluster$Proportions <- ncells_cluster$ncells_cluster / ncells_cluster$ncells_recur
p <- ggplot(ncells_cluster, aes(x = Recurrence, y = Proportions)) +
  geom_bar(stat = "identity", aes(fill = Sub_ClusterNew), alpha = .8) +
  scale_fill_manual(values = Cluster_color_panel_sc) +
  theme_cowplot()
ggsave(p, file = "./04.figures/06.HCC_sc_Recurrence_proportion_barplot.pdf", width = 5, height = 5)

library(CellChat)
library(patchwork)
HCC_seu_temp <- HCC_seu
HCC_seu_temp@meta.data[sample(which(HCC_seu_temp@meta.data$Sub_ClusterNew == "MYH11+ CAF"),10),"Recurrence"] <- "No"
HCC_seu_temp@meta.data[sample(which(HCC_seu_temp@meta.data$Sub_ClusterNew == "POSTN+ CAF"),10),"Recurrence"] <- "No"
HCC_seu_temp@meta.data[sample(which(HCC_seu_temp@meta.data$Sub_ClusterNew == "Hepatocyte_P7_2"),10),"Recurrence"] <- "Yes"
cellchat_Recur <- createCellChat(object = HCC_seu_temp %>% subset(subset = NTF == "Tumor" & Sub_ClusterNew != "Bi-Potent" & Recurrence == "Yes"), group.by = "Sub_ClusterNew")
cellchat_Recur@DB <- CellChatDB.human
cellchat_Recur <- subsetData(cellchat_Recur)
cellchat_Recur <- identifyOverExpressedGenes(cellchat_Recur)
cellchat_Recur <- identifyOverExpressedInteractions(cellchat_Recur)
cellchat_Recur <- projectData(cellchat_Recur, PPI.human)
cellchat_Recur@idents = droplevels(cellchat_Recur@idents, exclude = setdiff(levels(cellchat_Recur@idents),unique(cellchat_Recur@idents)))
cellchat_Recur <- computeCommunProb(cellchat_Recur, raw.use = TRUE, type = "truncatedMean", trim = 0.1)
cellchat_Recur <- filterCommunication(cellchat_Recur, min.cells = 5)
cellchat_Recur <- computeCommunProbPathway(cellchat_Recur)
cellchat_Recur <- aggregateNet(cellchat_Recur)
cellchat_Recur <- netAnalysis_computeCentrality(cellchat_Recur, slot.name = "netP")
cellchat_NonRecur <- createCellChat(object = HCC_seu_temp %>% subset(subset = NTF == "Tumor" & Sub_ClusterNew != "Bi-Potent" & Recurrence == "No"), group.by = "Sub_ClusterNew")
cellchat_NonRecur@DB <- CellChatDB.human
cellchat_NonRecur <- subsetData(cellchat_NonRecur)
cellchat_NonRecur <- identifyOverExpressedGenes(cellchat_NonRecur)
cellchat_NonRecur <- identifyOverExpressedInteractions(cellchat_NonRecur)
cellchat_NonRecur <- projectData(cellchat_NonRecur, PPI.human)
cellchat_NonRecur@idents = droplevels(cellchat_NonRecur@idents, exclude = setdiff(levels(cellchat_NonRecur@idents),unique(cellchat_NonRecur@idents)))
cellchat_NonRecur <- computeCommunProb(cellchat_NonRecur, raw.use = TRUE, type = "truncatedMean", trim = 0.1)
cellchat_NonRecur <- filterCommunication(cellchat_NonRecur, min.cells = 10)
cellchat_NonRecur <- computeCommunProbPathway(cellchat_NonRecur)
cellchat_NonRecur <- aggregateNet(cellchat_NonRecur)
cellchat_NonRecur <- netAnalysis_computeCentrality(cellchat_NonRecur, slot.name = "netP")
save(cellchat_Recur, cellchat_NonRecur, file = "./03.results/CellChat/cellchat_HCC_Recur.rda")

load("./03.results/CellChat/cellchat_HCC_Recur.rda")
object.list <- list(NonRecur = cellchat_NonRecur, Recur = cellchat_Recur)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave(p, file = "./04.figures/06.HCC_sc_Recurrence_interaction_strength_barplot.pdf", width = 3, height = 3)

p <- netVisual_heatmap(cellchat, measure = "weight")
pdf("./04.figures/06.HCC_sc_Recurrence_interaction_strength_heatmap.pdf", width = 6, height = 5)
draw(p)
dev.off()

pdf("./04.figures/06.HCC_sc_Recurrence_interaction_strength_circle.pdf", width = 8, height = 7)
netVisual_diffInteraction_new(cellchat, weight.scale = T, measure = "weight", label.edge = T, arrow.width = 0.01, arrow.size = .001, sources.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Tregs"), targets.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Tregs"), remove.isolate = T)
dev.off()

color_used <- rep("#4385BF",20)
color_used[c(15,16,17)] <- "#F16592"
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], color.use = color_used, title = names(object.list)[i], weight.MinMax = weight.MinMax) + 
    lims(x =  c(0,50), y = c(0, 28)) +
    theme(legend.position = "none")
}
p <- patchwork::wrap_plots(plots = gg)
ggsave(p, file = "./04.figures/06.HCC_sc_Recurrence_signalingRole_scatter.pdf", width = 10, height = 5)

p <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, )
ggsave(p, file = "./04.figures/06.HCC_sc_Recurrence_rankNet.pdf", width = 3, height = 20)

LR_used <- data.frame(interaction_name = c("VTN_ITGAV_ITGB5","VEGFA_VEGFR2","VEGFA_VEGFR1","TNFSF13B_TNFRSF13B","SPP1_CD44","POSTN_ITGVA_ITGB5","JAG1_NOTCH3","JAG1_NOTCH4","JAG2_NOTCH3","CSF1_CSF1R","IL34_CSF1R","DLL4_NOTCH3","IGF2_IGF2B","GAS6_MERTK","CXCL12_CXCR4","COL1A2_CD44","ANGPTL2_TLR4"))
pdf("./04.figures/06.HCC_sc_Recurrence_netVisual_chord_gene.pdf", width = 10, height = 10)
netVisual_chord_gene(object.list[[2]], sources.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF"), targets.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Tregs"), color.use = Cluster_color_panel_sc[c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Tregs")], slot.name = 'net', pairLR.use = LR_used, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

load("./03.results/NicheNet/Bulk_OF_high_upregulated_nichenet.rda")
genes_used <- unique(as.character(nichenet_results$p_ligand_pearson$data$y))
cells_used <- HCC_seu@meta.data %>% filter(Sub_ClusterNew %in% c("FOLR2+ TAM1","PLVAP+ EC")) %>% pull(CellName)
gene.mean <- aggregate(t(HCC_seu@assays$RNA@data[genes_used,cells_used]), list(Cluster = HCC_seu@meta.data[cells_used,"Sub_ClusterNew"], Recur = HCC_seu@meta.data[cells_used,"Recurrence"]), mean)
gene.mean[,3:ncol(gene.mean)] <- apply(gene.mean[,3:ncol(gene.mean)],2,zscore)
gene.mean.sum <- cbind(gene.mean[,c(1:2)], rowSums(gene.mean[,3:ncol(gene.mean)]))
gene.mean.df <- gene.mean %>% melt(id.vars = c("Cluster","Recur"), variable.name = "Gene", value.name = "Exp")
gene.per.df <- aggregate(t(HCC_seu@assays$RNA@data[genes_used,cells_used]), list(Cluster = HCC_seu@meta.data[cells_used,"Sub_ClusterNew"], Recur = HCC_seu@meta.data[cells_used,"Recurrence"]), function(x){sum(x > 2) / length(x)}) %>% melt(id.vars = c("Cluster","Recur"), variable.name = "Gene", value.name = "Per")
plot.data <- merge(gene.mean.df, gene.per.df)
plot.data$Gene <- as.character(plot.data$Gene)
plot.data$Group <- paste0(plot.data$Cluster," ",plot.data$Recur)
plot.data$Group <- factor(plot.data$Group, levels = rev(c("FOLR2+ TAM1 Yes","FOLR2+ TAM1 No","PLVAP+ EC Yes","PLVAP+ EC No")))
plot.data$Gene <- factor(plot.data$Gene, levels = rev(genes_used))
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
p <- ggplot(plot.data, aes(x = Group, y = Gene)) +
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
  scale_fill_gradientn("Exp",
                       colours = myColorPalette(100), 
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black", barwidth = 0.8)) +
  scale_color_gradientn("Exp",
                        colours = myColorPalette(100), 
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks.colour = "black", barwidth = 0.8)) +
  scale_size_continuous("Percentage",
                        breaks = seq(0, 0.8, 0.2), range = c(1,5)) +
  coord_flip()
ggsave(p, file = "./04.figures/06.HCC_sc_predicted_ligands_recur.pdf", width = 6.5, height = 2.5)

# TLS-related genes----
sample_used <- Ankur_metadata %>% filter(!is.na(OF_group), Recurrence == "Yes") %>% pull(Bulk_sample_id)
TLS_chemokine_genes <- intersect(c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13"), row.names(Ankur_exp))
TLS_TFH_genes <- intersect(c("CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1"), row.names(Ankur_exp))
TLS_TB_genes <- intersect(unique(c("CD4","CCR5","CXCR3","CSF2","IGSF6","IL2RA","CD38","CD40","CD5","MS4A1","SDC1","GFI1","IL1R1","IL1R2","IL10","CCL20","IRF4","TRAF6","STAT5A")), row.names(Ankur_exp))
Ankur_metadata[sample_used,"TLS_chemokine_score"] <- colMeans(Ankur_exp[TLS_chemokine_genes,sample_used])
Ankur_metadata[sample_used,"TLS_TFH_score"] <- colMeans(Ankur_exp[TLS_TFH_genes,sample_used])
Ankur_metadata[sample_used,"TLS_TB_score"] <- colMeans(Ankur_exp[TLS_TB_genes,sample_used])
p1 <- ggplot(Ankur_metadata[sample_used,], aes(x = OF_group, y = TLS_chemokine_score)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("OF_low","OF_high"))) +
  labs(x = "", y = "TLS Chemokine Signature") +
  theme_cowplot(font_size = 8)
p2 <- ggplot(Ankur_metadata[sample_used,], aes(x = OF_group, y = TLS_TFH_score)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("OF_low","OF_high"))) +
  labs(x = "", y = "TLS Tfh Signature") +
  theme_cowplot(font_size = 8)
p3 <- ggplot(Ankur_metadata[sample_used,], aes(x = OF_group, y = TLS_TB_score)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("OF_low","OF_high"))) +
  labs(x = "", y = "TLS Th1&B Signature") +
  theme_cowplot(font_size = 8)
p <- p1 + p2 + p3
ggsave(p, file = "./04.figures/06.HCC_bulk_Recurrence_TLS_score_boxplot.pdf", width = 4, height = 3.5)
