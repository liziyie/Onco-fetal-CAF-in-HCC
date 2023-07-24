setwd("/Volumes/ZiyiLi/PostDoc_Project/Onco-fetal Project/")
source("./onco-fetal_utils.R")
library(Seurat)
library(Polychrome)
library(dplyr)
library(nichenetr)
library(reshape2)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)

# Load data----
HCC_seu <- readRDS("./Onco-fetal CAF Identification/HCC_Ankur_release.rds")
cluster_info <- data.frame(
  Global_Cluster = c(rep("Endothelium",8),rep("Fibroblast",9),rep("Hepatocyte",7),rep("Lymphocyte",5),rep("Myeloid",10)),
  Sub_Cluster = c("PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9","CD320+ EC","CD9+ EC","IGFBP3+ EC","PLPP3+ EC","TFF3+ EC",
                  "POSTN+ CAF","SDC2+ CAF","HSP+ CAF","MYH11+ CAF","APOA2+ CAF","CAF","Fib","ABCAB+ Fib","MT1M+ Fib",
                  "Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P8_1","Hepatocytes_P8_2","Hepatocytes_P9_1","Hepatocytes_P9_2","Hepatocytes_P15",
                  "B","NK","CD8+ T","CD4+ T","Tregs",
                  "FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2","MT1G+ TAM3","Monocyte","Mo-derived cells","DC1","DC2","pDC","Mast"
  )
)

# Fig.3a----
load("./Onco-fetal CAF Identification/nichenet_results_curated.rda")
HCC_seu$Global_Cluster <- factor(HCC_seu$Global_Cluster, levels = c("B cell","T cell","ILC","Mast","Mononuclear","Fibroblast","Endothelium","Hepatocyte","Doublet"))
cluster_info <- unique(HCC_seu@meta.data[,c("Global_Cluster","Sub_Cluster")]) %>% arrange(Global_Cluster)
row.names(cluster_info) <- cluster_info$Sub_Cluster
ligand_enrich <- data.frame(Cluster = cluster_info$Sub_Cluster)
for(cluster in names(nichenet_results)){
  genes_used <- nichenet_results[[cluster]]$ligand_activities %>% arrange(desc(pearson)) %>% pull(test_ligand)
  temp <- aggregate(colMeans(HCC_seu@assays$RNA@data[genes_used[1:20],]), list(Cluster = HCC_seu$Sub_Cluster), mean)
  row.names(temp) <- temp$Cluster
  temp <- temp[cluster_info$Sub_Cluster,]
  ligand_enrich <- cbind(ligand_enrich, temp$x)
}
colnames(ligand_enrich)[2:40] <- names(nichenet_results)
row.names(ligand_enrich) <- ligand_enrich$Cluster

cluster_used <- cluster_info %>% filter(Global_Cluster %in% c("Endothelium","Mononuclear","Hepatocyte")) %>% pull(Sub_Cluster)
plot.data1 <- ligand_enrich[,c("Cluster","POSTN+ CAF")] %>% filter(Cluster %in% cluster_used) %>% merge(.,cluster_info, by.x = "Cluster", by.y = "Sub_Cluster") %>% arrange(Global_Cluster, desc(`POSTN+ CAF`)) %>% setNames(c("Sub_Cluster","MeanExp","Global_Cluster"))
ggplot(plot.data1, aes(x = reorder(Sub_Cluster,-MeanExp), y = MeanExp, fill = Global_Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Top20 Ligands of POSTN+ CAF") +
  facet_wrap(~Global_Cluster, nrow = 1, scales = "free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p1

cluster_used <- cluster_info %>% filter(Global_Cluster %in% c("Fibroblast","Mononuclear","Hepatocyte")) %>% pull(Sub_Cluster)
plot.data2 <- ligand_enrich[,c("Cluster","PLVAP+ EC_9")] %>% filter(Cluster %in% cluster_used) %>% merge(.,cluster_info, by.x = "Cluster", by.y = "Sub_Cluster") %>% arrange(Global_Cluster, desc(`PLVAP+ EC_9`)) %>% setNames(c("Sub_Cluster","MeanExp","Global_Cluster"))
ggplot(plot.data2, aes(x = reorder(Sub_Cluster,-MeanExp), y = MeanExp, fill = Global_Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Top20 Ligands of PLVAP+ EC_9") +
  facet_wrap(~Global_Cluster, nrow = 1, scales = "free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p2

cluster_used <- cluster_info %>% filter(Global_Cluster %in% c("Fibroblast","Endothelium","Hepatocyte")) %>% pull(Sub_Cluster)
plot.data3 <- ligand_enrich[,c("Cluster","FOLR2+ TAM1_7")] %>% filter(Cluster %in% cluster_used) %>% merge(.,cluster_info, by.x = "Cluster", by.y = "Sub_Cluster") %>% arrange(Global_Cluster, desc(`FOLR2+ TAM1_7`)) %>% setNames(c("Sub_Cluster","MeanExp","Global_Cluster"))
ggplot(plot.data3, aes(x = reorder(Sub_Cluster,-MeanExp), y = MeanExp, fill = Global_Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Top20 Ligands of FOLR2+ TAM1_7") +
  facet_wrap(~Global_Cluster, nrow = 1, scales = "free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p3

cluster_used <- cluster_info %>% filter(Global_Cluster %in% c("Fibroblast","Endothelium","Hepatocyte")) %>% pull(Sub_Cluster)
plot.data4 <- ligand_enrich[,c("Cluster","FOLR2+ TAM1_1")] %>% filter(Cluster %in% cluster_used) %>% merge(.,cluster_info, by.x = "Cluster", by.y = "Sub_Cluster") %>% arrange(Global_Cluster, desc(`FOLR2+ TAM1_1`)) %>% setNames(c("Sub_Cluster","MeanExp","Global_Cluster"))
ggplot(plot.data4, aes(x = reorder(Sub_Cluster,-MeanExp), y = MeanExp, fill = Global_Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Top20 Ligands of FOLR2+ TAM1_1") +
  facet_wrap(~Global_Cluster, nrow = 1, scales = "free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p4
plot_grid(p1,p3,p2,p4)

# Fig.3b----
ligands_to_plot <- rbind(
  data.frame(Cluster = "FOLR2+ TAM",
             rbind(nichenet_results$`FOLR2+ TAM1_7`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% dplyr::select(test_ligand, pearson),
                   nichenet_results$`FOLR2+ TAM1_1`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% dplyr::select(test_ligand, pearson)) %>%
               arrange(desc(pearson)) %>% filter(!duplicated(test_ligand)) %>% top_n(50)),
  data.frame(Cluster = "POSTN+ CAF", nichenet_results$`POSTN+ CAF`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% dplyr::select(test_ligand, pearson)),
  data.frame(Cluster = "PLVAP+ EC",
             rbind(nichenet_results$`PLVAP+ EC_3`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% dplyr::select(test_ligand, pearson),
                   nichenet_results$`PLVAP+ EC_4`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% dplyr::select(test_ligand, pearson),
                   nichenet_results$`PLVAP+ EC_9`$ligand_activities %>% arrange(desc(pearson)) %>% top_n(50) %>% dplyr::select(test_ligand, pearson)) %>%
               arrange(desc(pearson)) %>% filter(!duplicated(test_ligand)) %>% top_n(50)))
colnames(ligands_to_plot)[c(2,3)] <- c("Gene","Pearson")
cells_used <- list(
  "FOLR2+ TAM" = HCC_seu@meta.data %>% filter(Sub_Cluster %in% c("POSTN+ CAF","MYH11+ CAF","DC1","DC2","PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9","Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P15")) %>% pull(CellName),
  "PLVAP+ EC" = HCC_seu@meta.data %>% filter(Sub_Cluster %in% c("POSTN+ CAF","MYH11+ CAF","DC1","DC2","FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2","Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P15")) %>% pull(CellName),
  "POSTN+ CAF" = HCC_seu@meta.data %>% filter(Sub_Cluster %in% c("DC1","DC2","FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2","PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9","Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P15")) %>% pull(CellName)
)
cluster_order <- c("FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2","DC1","DC2","POSTN+ CAF","MYH11+ CAF","PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9","Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P15")
plot.data <- list()
for(cluster in c("FOLR2+ TAM","POSTN+ CAF","PLVAP+ EC")){
  gene.per <- 
    aggregate(t(as.matrix(HCC_seu@assays$RNA@data[ligands_to_plot %>% filter(Cluster == cluster) %>% pull(Gene),cells_used[[cluster]]])),
              list(Group = HCC_seu@meta.data[cells_used[[cluster]],"Sub_Cluster"]), 
              function(x){sum(x > 0) / length(x)}) %>% 
    melt(id.vars = "Group", variable.name = "Gene", value.name = "Per")
  genes_used <- gene.per %>% filter(Per > 0.1) %>% pull(Gene) %>% as.character() %>% unique()
  gene.mean.matrix <- 
    aggregate(t(as.matrix(HCC_seu@assays$RNA@data[genes_used,cells_used[[cluster]]])),
              list(Cluster = HCC_seu@meta.data[cells_used[[cluster]],"Sub_Cluster"]), 
              mean)
  gene.mean.zscore <- apply(gene.mean.matrix[,2:ncol(gene.mean.matrix)], 2, function(x){(x - mean(x)) / sd(x)})
  row.names(gene.mean.zscore) <- gene.mean.matrix[,1]
  gene.mean.zscore.df <- 
    data.frame(Group = gene.mean.matrix[,1], gene.mean.zscore, check.names = F) %>% 
    melt(id.vars = "Group", variable.name = "Gene", value.name = "Exp")
  gene.mean.zscore <- t(gene.mean.zscore)
  plot.data[[cluster]] <- merge(merge(gene.mean.zscore.df, gene.per), ligands_to_plot %>% filter(Cluster == cluster))
  plot.data[[cluster]]$Group <- factor(plot.data[[cluster]]$Group, levels = cluster_order[cluster_order %in% plot.data[[cluster]]$Group])
}
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
p <- list()
for(cluster in c("FOLR2+ TAM","POSTN+ CAF","PLVAP+ EC")){
  p[[cluster]] <- ggplot(plot.data[[cluster]], aes(x = Group, y = reorder(Gene,Pearson))) +
    geom_point(aes(size = Per, fill = Exp, color = Exp), shape = 21) +
    scale_radius(range = c(0,4)) +
    facet_wrap(Cluster~., scales = "free", nrow = 1) +
    theme_cowplot() + 
    scale_x_discrete(limits = rev(levels(plot.data$Group))) +
    theme(
      axis.ticks = element_line(size = .5),
      axis.ticks.length = unit(0.05,"in"),
      axis.text.y = element_text(color = "black", size = 10),
      axis.text.x = element_text(color = "black", angle = 90, 
                                 hjust = 1, vjust = 0.5, size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      panel.background = element_rect(colour = "black", fill = "white")) + 
    labs(x = "", y = "") +
    scale_fill_gradientn("Exp",
                         colours = myColorPalette(100), 
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks.colour = "black", barwidth = 0.8)) +
    scale_color_gradientn("Exp",
                          colours = myColorPalette(100), 
                          guide = guide_colorbar(frame.colour = "black",
                                                 ticks.colour = "black", barwidth = 0.8))
}
p <- plot_grid(plotlist = p, nrow = 1)

# Fig.3c----
library(ggvenn)
temp <- list(
  "TAM" = ligands_to_plot %>% filter(Cluster == "FOLR2+ TAM") %>% pull(Gene) %>% unique(),
  "EC" = ligands_to_plot %>% filter(Cluster == "PLVAP+ EC") %>% pull(Gene) %>% unique(),
  "CAF" = ligands_to_plot %>% filter(Cluster == "POSTN+ CAF") %>% pull(Gene) %>% unique()
)
ggvenn(temp, fill_color = c("#E3618D","#3D98BE","#D9B72E"), fill_alpha = .7, show_percentage = TRUE)

# Fig.3d----
library(CellChat)
library(patchwork)
cellchat <- createCellChat(object = HCC_seu, group.by = "Sub_Cluster")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional", k = 10)
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural", k = 5)

color_used <- rep("#4385BF",40)
color_used[c(13,14,33,34,35,36)] <- "#F16592"
p <- netAnalysis_signalingRole_scatter(cellchat, color.use = color_used)

# Fig.3e----
df.net <- subsetCommunication(cellchat, slot.name = "net", sources.use = sources_used, targets.use = targets_used, pairLR.use = LR_used)
df.net$source.target <- paste(df.net$source, df.net$target, sep = "|")
df.net$source.target <- 
  factor(df.net$source.target, 
         levels = c(
           sort(grep("^Hepatocytes_P7",cluster_combination_used,value = T))[c(3,6,4,5,1,2,7)],
           sort(grep("^Hepatocytes_P9",cluster_combination_used,value = T))[c(3,6,4,5,1,2,7)],
           sort(grep("^Hepatocytes_P15",cluster_combination_used,value = T))[c(3,6,4,5,1,2,7)],
           sort(grep("^POSTN",cluster_combination_used,value = T))[c(5,6,3,4,7,2,1,8)],
           sort(grep("^MYH11",cluster_combination_used,value = T))[c(5,6,3,4,7,2,1,8)],
           sort(grep("^PLVAP",cluster_combination_used,value = T))[c(3,4,1,2,5,8,9,6,7,10)],
           sort(grep("^FOLR2",cluster_combination_used,value = T))[c(1,4,2,3,5,8,6,7)],
           sort(grep("^SPP1",cluster_combination_used,value = T))[c(1,4,2,3)]
         )
  )
df.net <- df.net %>% filter(!is.na(prob), !is.na(source.target))
df.net$prob[df.net$prob == 0] <- 0.00000001
df.net$prob <- -1/log(df.net$prob)
df.net$pval[df.net$pval > 0.05] = 1
df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
df.net$pval[df.net$pval > 0.001 & df.net$pval <= 0.01] = 3
df.net$pval[df.net$pval <= 0.001] = 4
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
legend.values <- c(1,2,3,4)
names(legend.values) <- c("p > 0.05","0.01 < p < 0.05","0.001 < p < 0.01","P < 0.001")
ggplot(df.net,aes(x=source.target,y=interaction_name_2)) +
  geom_point(aes(size=pval,color=prob)) +
  scale_radius(range = c(min(df.net$pval), max(df.net$pval)),
               breaks = sort(unique(df.net$pval)), labels = names(legend.values)[legend.values %in% sort(unique(df.net$pval))], name = "p-value") +
  theme_minimal() +
  theme(
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.1,"in"),
    strip.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5, size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "lightgrey",
      linetype = "dashed",
      size = 0.2
    )) +
  labs(x = "", y = "") +
  scale_color_gradientn(name = "Commun. Prob.", colours = myColorPalette(100)[c(1:29,71:99)], limits = c(min.norm, max(df.net$prob)), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barwidth = 0.8))

# Fig.3f----
LR.df <- list(
  Ligand = c("VTN","IGF1","IL34","CXCL12","GAS6","FGF7","ADGRE5","DLL4","CSF1","TNFSF10"), 
  Receptor = c("ITGB5","IGF1R","CSF1R","CXCR4","MERTK","FGFR1","CD55","NOTCH3","CSF1R","TNFRSF10B")
)
Cluster.df <- list(from = c("Hepatocytes_P7_2","Hepatocytes_P9_1","Hepatocytes_P15","POSTN+ CAF","MYH11+ CAF","PLVAP+ EC_4","PLVAP+ EC_9"), to = c("POSTN+ CAF","MYH11+ CAF","PLVAP+ EC_4","PLVAP+ EC_9","FOLR2+ TAM1_1","FOLR2+ TAM1_7"))
Ligand_mean_expression <- 
  aggregate(
    t(HCC_seu@assays$RNA@data[LR.df$Ligand,HCC_seu@meta.data %>% filter(Sub_Cluster %in% Cluster.df$from) %>% pull(CellName)]),
    list(Group = HCC_seu@meta.data[HCC_seu@meta.data %>% filter(Sub_Cluster %in% Cluster.df$from) %>% pull(CellName),"Sub_Cluster"]), 
    mean
  )
row.names(Ligand_mean_expression) <- Ligand_mean_expression$Group
Ligand_mean_expression$Group <- c()
Ligand_mean_expression_normalize <- apply(Ligand_mean_expression, 2, scale)
row.names(Ligand_mean_expression_normalize) <- row.names(Ligand_mean_expression)
Ligand.matrix <- Ligand_mean_expression_normalize
Ligand.matrix <- pmax(Ligand.matrix, -1)
Ligand.matrix <- pmin(Ligand.matrix, 2)
Ligand.matrix <- Ligand.matrix[Cluster.df$from, LR.df$Ligand]
color_used <- circlize::colorRamp2(seq(-1, 2, length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
Ligand.ht <- 
  Heatmap(t(Ligand.matrix),
          cluster_rows = FALSE, 
          cluster_columns = FALSE,
          col = color_used,
          row_names_gp = gpar(fontsize = 9),
          column_names_gp = gpar(fontsize = 9),
          heatmap_legend_param = list(
            at = c(-1, 0, 1, 2),
            labels = c("<-1", "0", "1", ">2"),
            title = "Ligand\nExp",
            legend_height = unit(3, "cm"))
  )

Receptor_mean_expression <- 
  aggregate(
    t(HCC_seu@assays$RNA@data[LR.df$Receptor,HCC_seu@meta.data %>% filter(Sub_Cluster %in% Cluster.df$to) %>% pull(CellName)]),
    list(Group = HCC_seu@meta.data[HCC_seu@meta.data %>% filter(Sub_Cluster %in% Cluster.df$to) %>% pull(CellName),"Sub_Cluster"]), 
    mean
  )
row.names(Receptor_mean_expression) <- Receptor_mean_expression$Group
Receptor_mean_expression$Group <- c()
Receptor_mean_expression_normalize <- apply(Receptor_mean_expression, 2, scale)
row.names(Receptor_mean_expression_normalize) <- row.names(Receptor_mean_expression)
Receptor.matrix <- Receptor_mean_expression_normalize
Receptor.matrix <- pmax(Receptor.matrix, -1)
Receptor.matrix <- pmin(Receptor.matrix, 1.8)
Receptor.matrix <- Receptor.matrix[Cluster.df$to, LR.df$Receptor]
color_used <- circlize::colorRamp2(seq(-1, 1.8, length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
Receptor.ht <- 
  Heatmap(t(Receptor.matrix),
          cluster_rows = FALSE, 
          cluster_columns = FALSE,
          row_names_side = "left",
          col = color_used,
          row_names_gp = gpar(fontsize = 9),
          column_names_gp = gpar(fontsize = 9),
          heatmap_legend_param = list(
            at = c(-1, 0, 0.9, 1.8),
            labels = c("<-1", "0", "0.9", ">1.8"),
            title = "Receptor\nExp",
            legend_height = unit(3, "cm"))
  )

Ligand.ht + Receptor.ht

# Extended Data Fig.4b----
cluster_info <- unique(HCC_seu@meta.data[,c("Global_Cluster","Sub_Cluster")]) %>% arrange(Global_Cluster)
row.names(cluster_info) <- cluster_info$Sub_Cluster
cluster_info$Global_Cluster <- as.character(cluster_info$Global_Cluster)
cluster_info$Global_Cluster[cluster_info$Global_Cluster %in% c("B cell","T cell","ILC")] <- "Lymphocyte"
count_network <- read.table("./Source Data/Extended Data Fig.4b.txt", sep = "\t", head = T, stringsAsFactors = F)
count_network$SOURCE[count_network$SOURCE == "MYL9+ CAF"] <- "MYH11+ CAF"
count_network$TARGET[count_network$TARGET == "MYL9+ CAF"] <- "MYH11+ CAF"
count_network$SOURCE[count_network$SOURCE == "VIM+ CAF"] <- "SDC2+ CAF"
count_network$TARGET[count_network$TARGET == "VIM+ CAF"] <- "SDC2+ CAF"
count_network <- merge(count_network, cluster_info, by.x = "TARGET", by.y = "Sub_Cluster")
count_network %>% filter(SOURCE == "POSTN+ CAF", TARGET %in% c("FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2","MT1G+ TAM3","PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9","Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P9_2","Hepatocytes_P15","CD4+ T","CD8+ T","Tregs","B","NK")) %>% 
  ggplot(.,aes(x = count, y = reorder(TARGET,count))) +
  geom_segment(aes(yend = TARGET), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(color = Global_Cluster)) +
  ggtitle("POSTN+ CAF") +
  facet_grid(Global_Cluster ~ ., scales = "free_y", space = "free_y") +
  geom_text(aes(x = count+5, label = count)) +
  labs(y = "", x = "Number of Significant Interactions") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), legend.position = "none") -> p1
count_network %>% filter(SOURCE == "PLVAP+ EC_4", TARGET %in% c("FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2","MT1G+ TAM3","POSTN+ CAF","MYH11+ CAF","Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P9_2","Hepatocytes_P15","CD4+ T","CD8+ T","Tregs","B","NK")) %>% 
  ggplot(.,aes(x = count, y = reorder(TARGET,count))) +
  geom_segment(aes(yend = TARGET), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(color = Global_Cluster)) +
  ggtitle("PLVAP+ EC_9") +
  facet_grid(Global_Cluster ~ ., scales = "free_y", space = "free_y") +
  geom_text(aes(x = count+5, label = count)) +
  labs(y = "", x = "Number of Significant Interactions") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), legend.position = "none") -> p2
count_network %>% filter(SOURCE == "FOLR2+ TAM1_1", TARGET %in% c("PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9","POSTN+ CAF","MYH11+ CAF","Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P9_2","Hepatocytes_P15","CD4+ T","CD8+ T","Tregs","B","NK")) %>% 
  ggplot(.,aes(x = count, y = reorder(TARGET,count))) +
  geom_segment(aes(yend = TARGET), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(color = Global_Cluster)) +
  ggtitle("FOLR2+ TAM1_1") +
  facet_grid(Global_Cluster ~ ., scales = "free_y", space = "free_y") +
  geom_text(aes(x = count+5, label = count)) +
  labs(y = "", x = "Number of Significant Interactions") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), legend.position = "none") -> p3
count_network %>% filter(SOURCE == "FOLR2+ TAM1_7", TARGET %in% c("PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9","POSTN+ CAF","MYH11+ CAF","Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P9_2","Hepatocytes_P15","CD4+ T","CD8+ T","Tregs","B","NK")) %>% 
  ggplot(.,aes(x = count, y = reorder(TARGET,count))) +
  geom_segment(aes(yend = TARGET), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(color = Global_Cluster)) +
  ggtitle("FOLR2+ TAM1_7") +
  facet_grid(Global_Cluster ~ ., scales = "free_y", space = "free_y") +
  geom_text(aes(x = count+5, label = count)) +
  labs(y = "", x = "Number of Significant Interactions") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), legend.position = "none") -> p4
plot_grid(p1,p2,p3,p4,nrow = 2)

# Extended Data Fig.4c----
temp <- list(
  "TAM" = ligands_to_plot %>% filter(Cluster == "FOLR2+ TAM") %>% pull(Gene) %>% unique(),
  "EC" = ligands_to_plot %>% filter(Cluster == "PLVAP+ EC") %>% pull(Gene) %>% unique(),
  "CAF" = ligands_to_plot %>% filter(Cluster == "POSTN+ CAF") %>% pull(Gene) %>% unique()
)
genes_used <- unique(c(temp$TAM, temp$EC, temp$CAF))
genes_used <- data.frame(
  Gene = genes_used,
  inTAM = genes_used %in% temp$TAM,
  inEC = genes_used %in% temp$EC,
  inCAF = genes_used %in% temp$CAF
)
genes_used$Group <- paste0(ifelse(genes_used$inTAM,"TAM_",""),
                           ifelse(genes_used$inEC,"EC_",""),
                           ifelse(genes_used$inCAF,"CAF",""))
genes_used <- genes_used %>% filter(Group %in% c("EC_CAF","TAM_CAF","TAM_EC_","TAM_EC_CAF")) %>% arrange(Group)
mean_exp <- aggregate(t(as.matrix(HCC_seu@assays$RNA@data[genes_used$Gene,HCC_seu$Global_Cluster %ni% c("Doublet")])),list(Cluster = HCC_seu$Sub_Cluster[HCC_seu$Global_Cluster %ni% c("Doublet")]),mean)
row.names(mean_exp) <- mean_exp$Cluster
mean_exp$Cluster <- c()
markers_plot_matrix <- apply(mean_exp, 2, scale) %>% t()
colnames(markers_plot_matrix) <- row.names(mean_exp)
markers_plot_matrix_quantile <- quantile(markers_plot_matrix, c(0.01, 0.98))
markers_plot_matrix <- pmax(markers_plot_matrix, markers_plot_matrix_quantile[1])
markers_plot_matrix <- pmin(markers_plot_matrix, markers_plot_matrix_quantile[2])
markers_plot_matrix <- markers_plot_matrix[,cluster_info$Sub_Cluster]
color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), max(markers_plot_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
ha <- HeatmapAnnotation(
  Group = genes_used$Group,
  col = list(Group = c("EC_CAF"="#BC3C29FF","TAM_CAF"="#0072B5FF","TAM_EC_"="#E18727FF","TAM_EC_CAF"="#20854EFF"))
  )
ha_row <-  rowAnnotation(
  df = data.frame(Global_Cluster = cluster_info$Global_Cluster),
  width = unit(0.5, "cm"))
Heatmap(t(markers_plot_matrix),
        cluster_rows = FALSE, cluster_columns = FALSE,
        top_annotation = ha,
        left_annotation = ha_row,
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 7),
        col = color_used, name = "Exp",
        column_split = genes_used$Group,
        row_split = cluster_info$Global_Cluster)

# Extended Data Fig.4d----
genes_used <- c("NOTCH3","CCL2","ITGB1","CXCL12","PTN")
gene.mean.matrix <- 
  aggregate(t(as.matrix(Fib_fetal@assays$RNA@data[genes_used,])),
            list(Cluster = Fib_fetal$Sub_Cluster), 
            mean)
gene.mean.zscore <- apply(gene.mean.matrix[,2:ncol(gene.mean.matrix)], 2, function(x){(x - mean(x)) / sd(x)})
row.names(gene.mean.zscore) <- gene.mean.matrix[,1]
gene.mean.zscore.df <- 
  data.frame(Group = gene.mean.matrix[,1], gene.mean.zscore, check.names = F) %>% 
  melt(id.vars = "Group", variable.name = "Gene", value.name = "Exp")
gene.mean.zscore <- t(gene.mean.zscore)
gene.per <- 
  aggregate(t(as.matrix(Fib_fetal@assays$RNA@data[genes_used,])),
            list(Group = Fib_fetal$Sub_Cluster),
            function(x){sum(x > 2) / length(x)}) %>% 
  melt(id.vars = "Group", variable.name = "Gene", value.name = "Per")
plot.data <- merge(gene.mean.zscore.df, gene.per)
plot.data$Gene <- factor(as.character(plot.data$Gene), levels = genes_used)
plot.data$Group <- factor(plot.data$Group, levels = rev(c("Fib","ACTA2+ Fib","COL1A1+ Fib","BGN+ Fib","APOC3+ Fib","ITM2C+ Fib","Mesenchymal","PLEK+ Fib")))
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplot(plot.data, aes(x = Group, y = Gene)) +
  geom_point(aes(size = Per, fill = Exp, color = Exp), shape = 21) +
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
                        breaks = seq(0, 0.8, 0.2), range = c(1,8))

# Extended Data Fig.4e----
# see Fig.3e

# Supplementary Fig.5a----
hepato_seu <- HCC_seu %>% subset(subset = Global_Cluster == "Hepatocyte") %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
hepato_seu <- hepato_seu %>% FindNeighbors(dims = 1:15) %>% RunUMAP(dims = 1:15)
DimPlot(hepato_seu, pt.size = .5)

# Supplementary Fig.5d----
module_genes_list <- readRDS("./Onco-fetal CAF Identification/hepatocyte_nmf_programs_genes.RDS")
for(module in names(module_genes_list)){
  cat(module,"\n")
  hepato_seu@meta.data[,paste0(module,"_Score")] <- 
    CalculateModuleScore(hepato_seu@assays$RNA@data, module_genes_list[[module]])
}
hepato_seu@meta.data[,c("Sub_Cluster","PatientID","ViralvsNonViral",paste0("M",1:8,"_Score"))] %>% melt(measure.vars = paste0("M",1:8,"_Score")) %>% 
  ggplot(aes(x = factor(Sub_Cluster), y = value)) +
  geom_boxplot(fill = "#9FA0A3", outlier.shape = 1, alpha = .8) +
  facet_wrap(~variable, scales = "free", nrow = 2) +
  theme_cowplot() +
  labs(y = "Expression", x = "") +
  theme_cowplot(font_size = 9) +
  theme(strip.text.x = element_text(size = 9),
        strip.background = element_rect(fill=NA, color=NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Supplementary Fig.5e----
umap_by_module_score <- umap::umap(hepato_seu@meta.data[,paste0("M",c(1:5,7,8),"_Score")], min_dist = .5, spread = 1)
hepato_seu@meta.data[,"UMAP_1_Module"] <- umap_by_module_score$layout[,1]
hepato_seu@meta.data[,"UMAP_2_Module"] <- umap_by_module_score$layout[,2]
ggplot(hepato_seu@meta.data, aes(x = UMAP_1_Module, y = UMAP_2_Module)) +
  geom_point(aes(color = Sub_Cluster), size = 1.2, alpha = .8) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))

# Supplementary Fig.5f----
hepato_seu.plot.data <- hepato_seu@meta.data[,c("UMAP_1_Module","UMAP_2_Module",paste0("M",1:8,"_Score"))]
hepato_seu.plot.data <- melt(hepato_seu.plot.data, id.vars = c("UMAP_1_Module","UMAP_2_Module"))
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
ggplot(hepato_seu.plot.data %>% arrange(value), aes(x = UMAP_1_Module, y = UMAP_2_Module)) +
  geom_point(aes(color = value), alpha = 0.8, size = 1.2) +
  facet_wrap(~variable, nrow = 2) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("Module Score", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))

# Supplementary Fig.5g----
require(gtable)
require(egg)
cells_used <- hepato_seu@meta.data %>% pull(CellName)
genes_used <- c("APOM","AADAC","FETUB","LARS2","FKBP2","APOA4","PGLYRP2","SPP2","FGL1","IGFBP2","AFM","CFB","SLC3A1","SDS","PROC","F10")
hepato_seu$Z2_Score <- CalculateModuleScore(hepato_seu@assays$RNA@data, genes_used)
expression.matrix <- as.matrix(HCC_seu@assays$RNA@data[genes_used,cells_used])
groupid <- factor(HCC_seu@meta.data[cells_used,"Sub_Cluster"])
# Boxplot in the left panel
boxplot.data <- data.frame(Expression = hepato_seu$Z2_Score, Group = groupid)
p_boxplot <- 
  ggplot(boxplot.data, aes(x = Group, y = Expression)) + 
  geom_boxplot(color = "#1979B5", fill = "#1979B5",
               outlier.colour = "black", outlier.shape = 1) +
  theme_bw() +
  ggsignif::geom_signif(comparisons = list(
    c("Hepatocytes_P7_2","Hepatocytes_P7_1"),
    c("Hepatocytes_P9_2","Hepatocytes_P9_1"),
    c("Hepatocytes_P7_2","Hepatocytes_P15"),
    c("Hepatocytes_P7_2","Hepatocytes_P8_1"),
    c("Hepatocytes_P7_2","Hepatocytes_P9_1"),
    c("Hepatocytes_P9_1","Hepatocytes_P8_1")
  ), step = .1) +
  theme(legend.position = "NULL",
        axis.text.y = element_text(size = 8),
        axis.title = element_blank()) +
  scale_x_discrete(limits = rev(levels(as.factor(groupid)))) +
  coord_flip()
dat <- ggplot_build(p_boxplot)$data[[1]]
dat$xsp <- 1/2 * (dat$xmax + dat$xmin) - 1/4 * (dat$xmax - dat$xmin)
dat$xep <- 1/2 * (dat$xmax + dat$xmin) + 1/4 * (dat$xmax - dat$xmin)
p_boxplot <- 
  p_boxplot + 
  geom_segment(data = dat,  aes(x = xmin, xend = xmax, y = middle, yend = middle), colour = "#2BA147", size = 2) +
  geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymin, yend = ymin), colour = "black", size = 2) +
  geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymax, yend = ymax), colour = "black", size = 2)
# Heatmap in the right panel
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
gene.median.expression <- aggregate(t(expression.matrix), list(Group = groupid), mean)
gene.median.expression[, 2:(length(genes_used) + 1)] <- apply(gene.median.expression[, 2:(length(genes_used) + 1)], 2, function(x){
  (x - mean(x))/sd(x)
})
heatmap.data <- reshape2::melt(gene.median.expression, id = "Group")
heatmap.data$Group <- as.factor(heatmap.data$Group)
heatmap.data$Group <- factor(heatmap.data$Group, levels = rev(levels(heatmap.data$Group)))
max_value <- ceiling(max(heatmap.data$value))
scale_breaks <- round(c(-max_value * 4/5, -max_value * 2/5, 0, max_value * 2/5, max_value * 4/5),2)
p_heatmap <- ggplot(heatmap.data, aes(x = factor(variable), y = Group)) +
  geom_tile(aes(fill = value), colour = "white") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        legend.title = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1, size = 8),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_gradientn(limits = c(-max_value, max_value), breaks = scale_breaks, expand = c(0,0), colours = myColorPalette(100)) +
  scale_x_discrete(breaks = unique(heatmap.data$variable), labels = unique(heatmap.data$variable)) + 
  scale_y_discrete(expand = c(0,0))
gt = ggplotGrob(p_heatmap)
leg = gtable_filter(gt, "guide-box")
leg[[1]][[1]][[1]][[1]][[1]][[2]]$height = unit(1, "npc")
pos = unit.c(unit(0.01,"npc"), unit(.25, "npc"), unit(.5, "npc"), unit(.75, "npc"), unit(.99, "npc"))
leg[[1]][[1]][[1]][[1]][[1]][[3]]$children[[1]]$y = pos
leg[[1]][[1]][[1]][[1]][[1]][[5]]$y0 = pos
leg[[1]][[1]][[1]][[1]][[1]][[5]]$y1 = pos
leg[[1]][[1]][[1]][[1]]$heights = unit.c(rep(unit(0, "mm"), 3),
                                         unit(1, "npc"),
                                         unit(0, "mm"))
leg[[1]][[1]]$heights[[3]] = sum(rep(unit(0, "mm"), 3),
                                 unit(1, "npc"),
                                 unit(0, "mm"))
p_heatmap = gtable_add_grob(gt, leg, t = 7, l = 9)
# Final plot
g1 <- ggplotGrob(p_boxplot)
g2 <- p_heatmap
fg1 <- gtable_frame(g1, width = unit(30, "null"))
fg2 <- gtable_frame(g2, width = unit(20, "null"))
p <- gtable_frame(gtable_cbind(fg1, fg2))