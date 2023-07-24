setwd("/Volumes/ZiyiLi/PostDoc_Project/Onco-fetal Project/")
source("./onco-fetal_utils.R")
library(Seurat)
library(CellChat)
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
HCC_seu@meta.data$Recurrence <- "Yes"
HCC_seu@meta.data$Recurrence[HCC_seu@meta.data$PatientID %in% c("P3","P6","P7","P8","P11","P12","P13")] <- "No"
load("./HCC and Fetal Liver Nanostring DSP Data/nanostring.rda")
tumor_prop <- read.table("./HCC and Fetal Liver Nanostring DSP Data/CIBERSORTx_Results.txt", sep = "\t", head = T, check.names = F, row.names = 1)
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
Cluster_levels <- c("B","T cells","Mast","DC","pDC","SPP1+ TAM2","MT1G+ TAM3","Other Mononuclear","Other Endothelium","Other CAF","FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Hepatocyte","Hepatocyte_P7_2","Hepatocyte_P15")
tumor_prop2 <- tumor_prop2[,Cluster_levels]

# Fig.6a----
ncells_cluster <- data.frame(HCC_seu@meta.data %>% filter(NTF == "Tumor", Sub_ClusterNew != "Bi-Potent") %>% dplyr::select(Recurrence, Sub_ClusterNew) %>% table())
ncells_recur <- data.frame(HCC_seu@meta.data %>% filter(NTF == "Tumor", Sub_ClusterNew != "Bi-Potent") %>% dplyr::select(Recurrence) %>% table())
colnames(ncells_recur) <- c("Recurrence","Freq")
ncells_cluster <- merge(ncells_cluster, ncells_recur, by = c("Recurrence"))
colnames(ncells_cluster)[c(3,4)] <- c("ncells_cluster","ncells_recur")
ncells_cluster$Proportions <- ncells_cluster$ncells_cluster / ncells_cluster$ncells_recur
ggplot(ncells_cluster, aes(x = Recurrence, y = Proportions)) +
  geom_bar(stat = "identity", aes(fill = Sub_ClusterNew), alpha = .8) +
  theme_cowplot()

ROIE_table <- ROIE(table(HCC_seu@meta.data %>% filter(NTF == "Tumor", Sub_ClusterNew != "Bi-Potent") %>% dplyr::select(Recurrence, Sub_ClusterNew)))

# Fig.6b----
library(CellChat)
library(patchwork)
cellchat_Recur <- createCellChat(object = HCC_seu %>% subset(subset = NTF == "Tumor" & Sub_ClusterNew != "Bi-Potent" & Recurrence == "Yes"), group.by = "Sub_ClusterNew")
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
cellchat_NonRecur <- createCellChat(object = HCC_seu %>% subset(subset = NTF == "Tumor" & Sub_ClusterNew != "Bi-Potent" & Recurrence == "No"), group.by = "Sub_ClusterNew")
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
#load("./Onco-fetal CAF Identification/HCC_cellchat.rda")
group.new = levels(cellchat_Recur@idents)
cellchat_NonRecur <- liftCellChat(cellchat_NonRecur, group.new)
object.list <- list(NonRecur = cellchat_NonRecur, Recur = cellchat_Recur)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T, arrow.width = 0.01, arrow.size = .001, sources.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Tregs"), targets.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Tregs"), remove.isolate = T)

# Fig.6c----
LR_used <- data.frame(interaction_name = c("VTN_ITGAV_ITGB1","VEGFA_VEGFR2","VEGFA_VEGFR1","TGFB3_TGFBR1_TGFBR2","GDF15_TGFBR2","PTN_NCL","SPP1_CD44","POSTN_ITGVA_ITGB5","PDGFA_PDGFRA","PDGFA_PDGFRB","JAG1_NOTCH3","JAG1_NOTCH4","JAG2_NOTCH3","MDK_SDC1","CSF1_CSF1R","IL34_CSF1R","DLL4_NOTCH3","IGF2_IGF2B","GAS6_MERTK","CXCL12_CXCR4","CXCL1_ACKR1","AREG_EGFR"))
netVisual_chord_gene(object.list[[2]], sources.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Hepatocyte_P7_2","Hepatocyte_P15"), targets.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Tregs","Hepatocyte_P7_2","Hepatocyte_P15"), slot.name = 'net', pairLR.use = LR_used, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

# Fig.6e----
plot.data <- data.frame(Sample = ns_tumor_sample[row.names(tumor_prop),"Sample"],
                        OF_score = rowSums(tumor_prop2[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")]),
                        tumor_prop2[,c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF")],
                        check.names = F)
ggplot(plot.data %>% reshape2::melt(), aes(x = Sample, y = value)) +
  geom_boxplot(aes(color = Sample), alpha = .8) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  labs(y = "Predicted proportions", x = "") +
  theme_cowplot() +
  ggsignif::geom_signif(comparisons = list(c("B015","B017")),
                        step_increase = .1, tip_length = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Fig.6f----
hc <- hclust(dist(tumor_prop2, method = "euclidean"))
dend.hc <- as.dendrogram(hc)
dend.hc %>% plot(center = T)
dend.hc %>% rect.dendrogram(k=3, border=8, lty=5, lwd=2) 

plot.df <- tumor_prop2[,colSums(tumor_prop2)!=0] %>% mutate(Sample = row.names(.)) %>% melt() %>% mutate(Sample = factor(Sample, levels = labels(dend.hc)), variable = factor(variable, levels = Cluster_levels))
ggplot(plot.df, aes(x = Sample, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity") +
  labs(y = "Predicted proportions", x = "") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

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
  Sample = ns_tumor_sample[sample_used,"Sample"]
)
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

# Fig.6g-h----
ns_tumor_sample_temp <- ns_tumor_sample %>% filter(!is.na(Cluster))
ns_tumor_sample_temp$Cluster <- paste0("C",ns_tumor_sample_temp$Cluster)
ns_tumor_sample_temp <- cbind(ns_tumor_sample_temp, tumor_prop2[row.names(ns_tumor_sample_temp), ])
ns_tumor_sample_temp$OF_score <- ns_tumor_sample_temp$`FOLR2+ TAM1` + ns_tumor_sample_temp$`PLVAP+ EC` + ns_tumor_sample_temp$`POSTN+ CAF`
p1 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = `PLVAP+ EC`)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Predicted proportions of\n PLVAP+ EC", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p2 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = `FOLR2+ TAM1`)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Predicted proportions of\n FOLR2+ TAM1", x = "") +
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
p4 <- ggplot(ns_tumor_sample_temp %>% filter(Cluster %in% c("C1","C2")), aes(x = Cluster, y = OF_score)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Onco-fetal Score", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")
p1 + p2 + p3 + p4

# Fig.6i----
prioritized_ligands <- c("LTA","ADAM17","PGF","COL2A1","COL4A1","BMP2","ITGAM","TGFB3","VEGFA","BTLA")
ns_tumor_sample_temp$prioritized_ligands_score <- colMeans(ns_tumor_exp_norm[prioritized_ligands,row.names(ns_tumor_sample_temp)])
ggplot(ns_tumor_sample_temp %>% filter(Cluster != "C3"), aes(x = Cluster, y = prioritized_ligands_score)) +
  geom_boxplot(aes(fill = Cluster), alpha = .8) +
  theme_cowplot(font_size = 8) +
  labs(y = "Shared Prioritized Ligands", x = "") +
  scale_fill_manual(values = c("C2" = "#F84141", "C1" = "#9FA0A3")) +
  ggsignif::geom_signif(comparisons = list(c("C1","C2")), textsize = 3) +
  theme(legend.position = "none")

# Fig.6j----
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
p1 + p2

# Extended Data Fig.8a----
genes_used <- c("LTA","PGF","BMP2","ITGAM","VEGFA","ADAM17","TGFB3","BTLA","TNF","TGFB1","SPP1","GAS6","VWF","ITGB2","CXCL12","VTN","TNFSF14","PLG","VCAM1","DLL4")
HCC_seu@meta.data$Sub_ClusterNew2 <- HCC_seu@meta.data$Sub_ClusterNew
HCC_seu@meta.data$Sub_ClusterNew2[HCC_seu@meta.data$Sub_ClusterNew2 %in% c("Hepatocyte","Hepatocyte_P15","Hepatocyte_P7_2")] <- "Hepatocyte"
cells_used <- HCC_seu@meta.data %>% filter(Sub_ClusterNew2 %in% c("FOLR2+ TAM1","PLVAP+ EC","SPP1+ TAM2","Hepatocyte")) %>% pull(CellName)
gene.mean <- aggregate(t(HCC_seu@assays$RNA@data[genes_used,cells_used]), list(Cluster = HCC_seu@meta.data[cells_used,"Sub_ClusterNew"], Recur = HCC_seu@meta.data[cells_used,"Recurrence"]), mean)
gene.mean[,3:ncol(gene.mean)] <- apply(gene.mean[,3:ncol(gene.mean)],2,scale)
gene.mean.sum <- cbind(gene.mean[,c(1:2)], rowSums(gene.mean[,3:ncol(gene.mean)]))
gene.mean.df <- gene.mean %>% melt(id.vars = c("Cluster","Recur"), variable.name = "Gene", value.name = "Exp")
gene.per.df <- aggregate(t(HCC_seu@assays$RNA@data[genes_used,cells_used]), list(Cluster = HCC_seu@meta.data[cells_used,"Sub_ClusterNew"], Recur = HCC_seu@meta.data[cells_used,"Recurrence"]), function(x){sum(x > 2) / length(x)}) %>% melt(id.vars = c("Cluster","Recur"), variable.name = "Gene", value.name = "Per")
plot.data <- merge(gene.mean.df, gene.per.df)
plot.data$Gene <- as.character(plot.data$Gene)
plot.data$Group <- paste0(plot.data$Cluster," ",plot.data$Recur)
plot.data$Group <- factor(plot.data$Group, levels = c("FOLR2+ TAM1 Yes","FOLR2+ TAM1 No","PLVAP+ EC Yes","PLVAP+ EC No","SPP1+ TAM2 Yes","SPP1+ TAM2 No","Hepatocyte Yes","Hepatocyte No"))
plot.data$Gene <- factor(plot.data$Gene, levels = genes_used)
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

# Extended Data Fig.8b----
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# Extended Data Fig.8c----
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
patchwork::wrap_plots(plots = gg)

# Extended Data Fig.8d----
p <- netVisual_heatmap(cellchat, measure = "weight")

# Extended Data Fig.8e----
p <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

# Extended Data Fig.8f----
p <- netVisual_bubble(cellchat, sources.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Hepatocyte_P7_2","Hepatocyte_P15"), targets.use = c("FOLR2+ TAM1","PLVAP+ EC","POSTN+ CAF","Tregs","Hepatocyte_P7_2","Hepatocyte_P15"), pairLR.use = LR_used, comparison = c(1, 2), angle.x = 45, max.quantile = .9)

# Extended Data Fig.8g----
genes_used <- c("SLIT2","ADAM15","NID1","TNFSF10","THBS1","DLL4","NUCB2","AREG","HGF","LAMB1","CYR61","THBS2","JAG1","GSTP1","HBEGF","COL1A1","COL18A1","HSPG2","INHBA","INHBB","NAMPT","ITGB1","ITGB2","COMP","SPP1","VCAM1","CTGF","COL4A1","VWF","FN1")
cells_used <- HCC_seu@meta.data %>% filter(Sub_ClusterNew %in% c("FOLR2+ TAM1","PLVAP+ EC")) %>% pull(CellName)
gene.mean <- aggregate(t(HCC_seu@assays$RNA@data[genes_used,cells_used]), list(Cluster = HCC_seu@meta.data[cells_used,"Sub_ClusterNew"], Recur = HCC_seu@meta.data[cells_used,"Recurrence"]), mean)
gene.mean[,3:ncol(gene.mean)] <- apply(gene.mean[,3:ncol(gene.mean)],2,scale)
gene.mean.sum <- cbind(gene.mean[,c(1:2)], rowSums(gene.mean[,3:ncol(gene.mean)]))
gene.mean.df <- gene.mean %>% melt(id.vars = c("Cluster","Recur"), variable.name = "Gene", value.name = "Exp")
gene.per.df <- aggregate(t(HCC_seu@assays$RNA@data[genes_used,cells_used]), list(Cluster = HCC_seu@meta.data[cells_used,"Sub_ClusterNew"], Recur = HCC_seu@meta.data[cells_used,"Recurrence"]), function(x){sum(x > 2) / length(x)}) %>% melt(id.vars = c("Cluster","Recur"), variable.name = "Gene", value.name = "Per")
plot.data <- merge(gene.mean.df, gene.per.df)
plot.data$Gene <- as.character(plot.data$Gene)
plot.data$Group <- paste0(plot.data$Cluster," ",plot.data$Recur)
plot.data$Group <- factor(plot.data$Group, levels = rev(c("FOLR2+ TAM1 Yes","FOLR2+ TAM1 No","PLVAP+ EC Yes","PLVAP+ EC No")))
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

# Extended Data Fig.9a----
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
p1 + p2

# Extended Data Fig.9b----
ns_tumor_sd <- sort(apply(ns_tumor_exp_norm,1,sd),decreasing = T)
features_used <- names(ns_tumor_sd)[1:1000]
pca_analysis <- prcomp(t(ns_tumor_exp_norm[features_used,]))
pca_analysis <- summary(pca_analysis)
pca_plot <- data.frame(
  pca_analysis$x[,1:2],
  Tissue = ns_tumor_sample$Tissue,
  Sample = ns_tumor_sample$Sample)
ggplot(pca_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Tissue, shape = Sample), size = 3, alpha = .8) +
  xlab(sprintf("PC1: %s%% variance",round(pca_analysis$importance[2,1],2) * 100)) + 
  ylab(sprintf("PC2: %s%% variance",round(pca_analysis$importance[2,2],2) * 100)) +
  theme_cowplot()

# Extended Data Fig.9c----
OF_TAM_signature <- c("FOLR2","CD163","MRC1","TREM2")
OF_Endo_signature <- c("PLVAP","ACKR1","DLL4")
OF_CAF_signature <- c("POSTN","FAP","ENG","PTGDS")
ns_tumor_sample$OF_TAM_score <- colMeans(ns_tumor_exp_norm[OF_TAM_signature,])
ns_tumor_sample$OF_Endo_score <- colMeans(ns_tumor_exp_norm[OF_Endo_signature,])
ns_tumor_sample$OF_CAF_score <- colMeans(ns_tumor_exp_norm[OF_CAF_signature,])
plot.data <- ns_tumor_sample %>% select(Tissue, Group, OF_TAM_score, OF_Endo_score, OF_CAF_score)
ggplot(plot.data %>% reshape2::melt(), aes(x = Group, y = value)) +
  geom_boxplot(aes(fill = Tissue), alpha = .8) +
  facet_wrap(~variable, scales = "free") +
  labs(y = "Expression", x = "") +
  theme_cowplot() +
  ggsignif::geom_signif(comparisons = list(c("B015_Normal","B017_Normal"),
                                           c("B015_Tumor","B017_Tumor"),
                                           c("B015_Normal","B015_Tumor"),
                                           c("B017_Normal","B017_Tumor")),
                        step_increase = .1, tip_length = 0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Extended Data Fig.9d----
degenes_12$Sig <- FALSE
degenes_12$label <- NA
degenes_12[degenes_12$adj.P.Val < 0.05 & abs(degenes_12$logFC) > 0.5,"Sig"] <- TRUE
genes_labeled <- c(
  degenes_12 %>% filter(Sig == TRUE) %>% slice_max(n = 25, order_by = logFC) %>% pull(Symbol),
  degenes_12 %>% filter(Sig == TRUE) %>% slice_min(n = 25, order_by = logFC) %>% pull(Symbol)
)
degenes_12[genes_labeled,"label"] <- genes_labeled
volcano_plot <- 
  ggplot(degenes_12, aes(x = -logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 3) +
  scale_color_manual(values = Binomial_color_panel) +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  theme_cowplot(font_size = 16) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1)) +
  geom_vline(xintercept = c(-0.5,0.5), color = "grey", linetype = "dashed", lwd = 1) +
  geom_text_repel(data = degenes_12, aes(label = label), size = 4)

# Extended Data Fig.9e----
library(AUCell)
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
p2 <- ggplot(HCC_seu@meta.data %>% filter(NTF == "Tumor") %>% arrange(C1_AUC), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = C1_AUC), size = .1) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("C1 AUC score", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = .2, barheight = 12))
p1 + p2

# Extended Data Fig.9f----
ns_tumor_sd <- sort(apply(ns_tumor_exp_norm,1,sd),decreasing = T)
ns_fetal_sd <- sort(apply(ns_fetal_exp_norm,1,sd),decreasing = T)
core_features <- unique(c(
  names(ns_tumor_sd)[1:500],
  names(ns_fetal_sd)[1:500]
))
samples_used <- ns_sample %>% pull(NewID)
hc <- hclust(dist(t(ns_exp_norm[core_features,samples_used]),method = "euclidean"))
dend.hc <- as.dendrogram(hc)
plot(dend.hc)
dend.hc %>% rect.dendrogram(k=6, border=8, lty=5, lwd=2) 

cluster <- c(rep("C1",24),rep("C2",10),rep("C3_1",9),rep("C3_2",24),rep("C4",2),rep("C5",3))
ns_sample$Cluster <- c()
ns_sample[labels(dend.hc),"Cluster"] <- cluster
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
  Sample = ns_sample[sample_used,"Sample"]
)
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

# Extended Data Fig.9g----
pca_analysis <- prcomp(t(ns_exp_norm[core_features,]))
pca_analysis <- summary(pca_analysis)
pca_plot <- data.frame(
  pca_analysis$x[,1:2],
  Tissue = ns_sample$Tissue,
  Sample = ns_sample$Sample,
  Cluster = ns_sample$Cluster,
  NewID = ns_sample$NewID)
p1 <- ggplot(pca_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Tissue, shape = Sample), size = 3, alpha = .8) +
  xlab(sprintf("PC1: %s%% variance",round(pca_analysis$importance[2,1],2) * 100)) + 
  ylab(sprintf("PC2: %s%% variance",round(pca_analysis$importance[2,2],2) * 100)) +
  theme_cowplot()
p2 <- ggplot(pca_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Cluster), size = 3, alpha = .8) +
  xlab(sprintf("PC1: %s%% variance",round(pca_analysis$importance[2,1],2) * 100)) + 
  ylab(sprintf("PC2: %s%% variance",round(pca_analysis$importance[2,2],2) * 100)) +
  theme_cowplot()
p1 + p2

# Extended Data Fig.9h----
 ggplot(ns_sample, aes(x = Cluster)) +
  geom_bar(aes(fill = Tissue), color = "white", position = "fill") +
  labs(y = "Frequency (%)") +
  theme_cowplot() + 
  RotatedAxis()

# Extended Data Fig.9i----
# see Fig.6j