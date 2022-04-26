setwd("/work/lzy/project/onco_fetal/")
source("../utils/utils_data_processing.R")
library(Seurat)
library(dplyr)

# >Load dataset----
CAF_HCC <- readRDS("./02.processed_data/CAF_HCC.rds")
CAF_HCC$Sub_Cluster <- as.character(CAF_HCC$Sub_Cluster)
CAF_HCC@active.ident <- factor(CAF_HCC$Sub_Cluster)

Fib_fetal <- readRDS("./02.processed_data/Fib_fetal.rds")
Fib_fetal$Sub_Cluster <- as.character(Fib_fetal$Sub_Cluster)
Fib_fetal$Sub_Cluster[Fib_fetal$Sub_Cluster == "Fib"] <- "IGFBP3+ Fib"
Fib_fetal@active.ident <- factor(Fib_fetal$Sub_Cluster)

# Table S1----
CAF_markers <- FindAllMarkers(CAF_HCC)
write.csv(CAF_markers, file = "./03.results/Tables/Table S1 A.csv")
Fib_markers <- FindAllMarkers(Fib_fetal)
write.csv(Fib_markers, file = "./03.results/Tables/Table S1 B.csv")

# Table S2----
cells_used <- CAF_HCC@meta.data %>% filter(Sub_Cluster %in% c("POSTN+ CAF", "MYH11+ CAF")) %>% pull(CellID)
degenes <- LIMMA(
  CAF_HCC@assays$RNA@data[,cells_used],
  stringr::str_split_fixed(CAF_HCC@meta.data[cells_used,"Sub_Cluster"],"\\+",2)[,1]
)
write.csv(degenes, file = "./03.results/Tables/Table S2 A.csv")

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
plot.df$Sig <- "Not significant"
plot.df[plot.df$logFC_fetal > 0 & plot.df$logFC_tumor > 0 & plot.df$p_fetal < 0.05 & plot.df$p_tumor < 0.05,"Sig"] <- "Both low up-regulated"
plot.df[plot.df$logFC_fetal > 0.5 & plot.df$logFC_tumor > 0.5 & plot.df$p_fetal < 0.05 & plot.df$p_tumor < 0.05,"Sig"] <- "Both high up-regulated"
plot.df[abs(plot.df$logFC_fetal) < 0.5 & plot.df$logFC_tumor > 0.5 & plot.df$p_fetal > 0.05 & plot.df$p_tumor < 0.05,"Sig"] <- "Tumor-specific up-regulated"
plot.df[plot.df$logFC_fetal < 0 & plot.df$p_fetal < 0.05 & plot.df$logFC_tumor < 0 & plot.df$p_tumor < 0.05,"Sig"] <- "Down-regulated"
write.csv(plot.df, file = "./03.results/Tables/Table S2 B.csv")

# Table S3----
load("./03.results/NicheNet/nichenet_results_curated.rda")
ligands_to_plot <- rbind(
  data.frame(Cluster = "FOLR2+ TAM",
             rbind(nichenet_results$`FOLR2+ TAM1_7`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% select(test_ligand, pearson),
                   nichenet_results$`FOLR2+ TAM1_1`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% select(test_ligand, pearson)) %>%
               arrange(desc(pearson)) %>% filter(!duplicated(test_ligand)) %>% top_n(50)),
  data.frame(Cluster = "POSTN+ CAF", nichenet_results$`POSTN+ CAF`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% select(test_ligand, pearson)),
  data.frame(Cluster = "PLVAP+ EC",
             rbind(nichenet_results$`PLVAP+ EC_3`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% select(test_ligand, pearson),
                   nichenet_results$`PLVAP+ EC_4`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% select(test_ligand, pearson),
                   nichenet_results$`PLVAP+ EC_9`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% select(test_ligand, pearson)) %>%
               arrange(desc(pearson)) %>% filter(!duplicated(test_ligand)) %>% top_n(50)))
colnames(ligands_to_plot)[c(2,3)] <- c("Gene","Pearson")
write.csv(ligands_to_plot, file = "./03.results/Tables/Table S3 A.csv")

library(CellChat)
cellchat <- readRDS("./03.results/CellChat/cellchat_HCC.rds")
macrophage_cluster_used <- c("FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2")
CAF_cluster_used <- c("POSTN+ CAF","MYH11+ CAF")
EC_cluster_used <- c("PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9")
Tumor_cluster_used <- c("Hepatocytes_P7_2","Hepatocytes_P9_1","Hepatocytes_P15")
T_cluster_used <- c("CD4+ T","CD8+ T","Tregs")
sources_used <- targets_used <- c(macrophage_cluster_used, EC_cluster_used, CAF_cluster_used, Tumor_cluster_used, T_cluster_used)
df.net <- subsetCommunication(cellchat, slot.name = "net", sources.use = sources_used, targets.use = targets_used)
df.net$prob[df.net$prob == 0] <- 0.00000001
df.net$prob <- -1/log(df.net$prob)
write.csv(df.net, file = "./03.results/Tables/Table S3 B.csv")

# Table S4----
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
Ankur_cibersort <- read.table("./00.data/HCC_bulk/CIBERSORTx/Output/CIBERSORTx_Results.txt", sep = "\t", head = T, check.names = F)
row.names(Ankur_cibersort) <- Ankur_cibersort$Mixture
Ankur_cibersort$Mixture <- c()
Ankur_cibersort[,"T cells"] <- Ankur_cibersort[,"CD4+ T"] + Ankur_cibersort[,"CD8+ T"] + Ankur_cibersort[,"Tregs"]
Ankur_cibersort[,"All Hepatocyte"] <- Ankur_cibersort[,"Hepatocyte_P15"] + Ankur_cibersort[,"Hepatocyte_P7_2"] + Ankur_cibersort[,"Hepatocyte"]
Ankur_cibersort$`P-value` <- Ankur_cibersort$Correlation <- Ankur_cibersort$RMSE <- c()
colnames(Ankur_cibersort)[colnames(Ankur_cibersort) == "Other Endotehlium"] <- "Other Endothelium"
bulk_info <- cbind(Ankur_metadata[row.names(Ankur_cibersort),], Ankur_cibersort)
write.csv(bulk_info, file = "./03.results/Tables/Table S4 B.csv")

# Table S5----
hcc01_bin100 <- readRDS("./02.processed_data/BGI/hcc01_bin100_seurat_processed.rds")
hcc01_bin100@active.ident <- factor(hcc01_bin100$Cluster)
hcc03_bin100 <- readRDS("./02.processed_data/BGI/hcc03_bin100_seurat_processed.rds")
hcc03_bin100@active.ident <- factor(hcc03_bin100$Cluster)

HCC01_markers <- FindAllMarkers(hcc01_bin100)
write.csv(HCC01_markers, file = "./03.results/Tables/Table S5 A.csv")
HCC03_markers <- FindAllMarkers(hcc03_bin100)
write.csv(HCC03_markers, file = "./03.results/Tables/Table S5 B.csv")

# Table S6----
Ankur_exp <- read.table("./00.data/HCC_bulk/Ankur_HCC_bulk.txt", sep = "\t", head = T)
Ankur_exp <- Ankur_exp[!duplicated(Ankur_exp$Gene),]
row.names(Ankur_exp) <- Ankur_exp$Gene
Ankur_exp$Gene <- c()
Ankur_sampleinfo <- read.csv("./00.data/HCC_bulk/sample_info.csv", row.names = 1, stringsAsFactors = F)
Ankur_clinical <- readr::read_tsv("./00.data/HCC_bulk/ClinicalData.tsv") %>% data.frame()
Ankur_metadata <- merge(Ankur_sampleinfo,Ankur_clinical,by.x="sample_id",by.y="Unified_ID")
row.names(Ankur_metadata) <- Ankur_metadata$Bulk_sample_id
Ankur_metadata <- Ankur_metadata[colnames(Ankur_exp),]
Ankur_cibersort <- read.table("./00.data/HCC_bulk/CIBERSORTx/Output/CIBERSORTx_Results.txt", sep = "\t", head = T, check.names = F)
row.names(Ankur_cibersort) <- Ankur_cibersort$Mixture
Ankur_cibersort$Mixture <- c()
Ankur_metadata[row.names(Ankur_cibersort),"OF_score"] <- rowSums(Ankur_cibersort[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")])
sample_used <- Ankur_metadata %>% filter(Bulk_sample_id %in% row.names(Ankur_cibersort), !is.na(Recurrence)) %>% pull(Bulk_sample_id)
OF_score_cutoff <- mean(Ankur_metadata[sample_used,"OF_score"]) + sd(Ankur_metadata[sample_used,"OF_score"])
Ankur_metadata[sample_used,"OF_group"] <- "OF_low"
Ankur_metadata[sample_used,"OF_group"][Ankur_metadata[sample_used,"OF_score"] > OF_score_cutoff] <- "OF_high"
Ankur_metadata$OF_group <- factor(Ankur_metadata$OF_group, levels = c("OF_low","OF_high"))
relapsed_samples <- Ankur_metadata %>% filter(Bulk_sample_id %in% row.names(Ankur_cibersort), Recurrence == "Yes") %>% pull(Bulk_sample_id)
degenes <- LIMMA(Ankur_exp[,relapsed_samples],
                 Ankur_metadata[relapsed_samples,"OF_group"])
write.csv(degenes, file = "./03.results/Tables/Table S6.csv")

# Table S7----
HCC_seu <- readRDS("./02.processed_data/HCC_seu.rds")
load("./02.processed_data/Nanostring/nanostring.rda")
Cluster_levels <- c("B","T cells","Mast","DC","pDC","SPP1+ TAM2","MT1G+ TAM3","Other Mononuclear","Other Endothelium","Other CAF","FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Hepatocyte","Hepatocyte_P7_2","Hepatocyte_P15")
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
hc <- hclust(dist(tumor_prop2, method = "euclidean"))
ns_tumor_sample$Cluster <- c()
ns_tumor_sample[names(cutree(hc,3)),"Cluster"] <- cutree(hc,3)
sample_used <- ns_tumor_sample %>% filter(Cluster %in% c(1,2)) %>% pull(NewID)
degenes_12 <- LIMMA(ns_tumor_exp_norm[,sample_used],
                    ns_tumor_sample[sample_used,"Cluster"])
write.csv(degenes_12, file = "./03.results/Tables/Table S7.csv")
