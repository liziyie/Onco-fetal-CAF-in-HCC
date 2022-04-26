setwd("/work/lzy/project/onco_fetal/")
source("../utils/utils_plot.R")
source("../utils/utils_data_processing.R")
source("../utils/utils_color.R")
library(Seurat)
library(Polychrome)
library(dplyr)
library(ggplot2)
library(cowplot)
library(nichenetr)
library(tidyverse)
library(ComplexHeatmap)

# >Load dataset----
HCC_seu <- readRDS("./02.processed_data/HCC_seu.rds")

Global_Cluster_color_panel <- c("Lymphocyte" = "#69B4CE","B cell" = "#69B4CE", "T cell" = "#ECAFCF", "ILC" = "#A0D7C9", "Mast" = "#BB4A94", "Mononuclear" = "#F3746C", "Myeloid" = "#F3746C", "Hepatocyte" = "#CAA57D", "Fibroblast" = "#C35338", "Double" = "#C35338", "Endothelium" = "#EAA944")

Hepatocyte_color_panel <- c("Hepatocytes_P7_1"="#BB0021FF","Hepatocytes_P7_2"="#EE0000FF","Hepatocytes_P8_1"="#008280FF","Hepatocytes_P8_2"="#008B45FF","Hepatocytes_P9_1"="#5F559BFF","Hepatocytes_P9_2"="#631879FF","Hepatocytes_P15"="#3B4992FF")

TAM_color_panel <- c(
  "FOLR2+ TAM1_1" = "#6A3D9A", "FOLR2+ TAM1_7" = "#CAB2D6",
  "SPP1+ TAM2" = "#EA6846"
)

cluster_info <- data.frame(
  Global_Cluster = c(rep("Endothelium",8),rep("Fibroblast",9),rep("Hepatocyte",7),rep("Lymphocyte",5),rep("Myeloid",10)),
  Sub_Cluster = c("PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9","CD320+ EC","CD9+ EC","IGFBP3+ EC","PLPP3+ EC","TFF3+ EC",
                  "POSTN+ CAF","SDC2+ CAF","HSP+ CAF","MYH11+ CAF","APOA2+ CAF","CAF","Fib","ABCAB+ Fib","MT1M+ Fib",
                  "Hepatocytes_P7_1","Hepatocytes_P7_2","Hepatocytes_P8_1","Hepatocytes_P8_2","Hepatocytes_P9_1","Hepatocytes_P9_2","Hepatocytes_P15",
                  "B","NK","CD8+ T","CD4+ T","Tregs",
                  "FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2","MT1G+ TAM3","Monocyte","Mo-derived cells","DC1","DC2","pDC","Mast"
  )
)

# >Tumor cells----
# InferCNV----
library(infercnv)
infercnv_obj <- readRDS("../Analysis/Results/InferCNV/EC_as_ref_sampled_20201127/21_denoiseHMMi6.NF_NA.SD_1.5.NL_FALSE.infercnv_obj")
plot_cnv(infercnv_obj, out_dir = "./Results/InferCNV/tmp", title = "InferCNV", output_format = "pdf")

# Expression of Zone2 hepatocytes markers----
require(gtable)
require(egg)
cells_used <- HCC_seu@meta.data %>% filter(Global_Cluster == "Hepatocyte") %>% pull(CellName)
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
pdf("./04.figures/04.Zone2_Hepatocytes_score_heatmap.pdf", width = 10, height = 3)
grid.draw(p)
dev.off()

# UMAP plot----
hepato_seu <- HCC_seu %>% subset(subset = Global_Cluster == "Hepatocyte") %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
hepato_seu <- hepato_seu %>% FindNeighbors(dims = 1:15) %>% RunUMAP(dims = 1:15)
p1 <- DimPlot(hepato_seu, pt.size = .5) + scale_color_manual(values = Hepatocyte_color_panel)
p1 <- ggAIplot(p1)
p2 <- DimPlot(hepato_seu, group.by = "PatientID", pt.size = .5) + scale_color_manual(values = c54[c(1:3,5:16)])
p2 <- ggAIplot(p2)
p <- p1 + p2
ggsave(p, file = "./04.figures/04.Hepatocytes_umap.pdf", width = 11, height = 4)

# NMF analysis----
source("../utils/utils_nmf_programs.R")
expr_tumor <- list()
for(cluster in HCC_seu@meta.data %>% filter(Global_Cluster == "Hepatocyte") %>% pull(Sub_Cluster) %>% unique){
  cells_used <- HCC_seu@meta.data %>% filter(Sub_Cluster == cluster) %>% pull(CellName)
  expr_tumor[[cluster]] <- as.matrix(HCC_seu@assays$RNA@data[,cells_used])
}
w_basis_tumor <- h_coef_tumor <- list()
for(i in names(expr_tumor)) {
  w <- NULL
  h <- NULL
  for(j in 6:9) {
    n <- nmf_programs(expr_tumor[[i]], is.log=T, rank=j)
    colnames(n$w_basis) <- paste0(i, "_", j, ".", 1:j)
    colnames(n$h_coef) <- paste0(i, "_", j, ".", 1:j)
    w <- cbind(w, n$w_basis)
    h <- cbind(h, n$h_coef)
  }
  w_basis_tumor[[i]] <- w
  h_coef_tumor[[i]] <- h
}
saveRDS(w_basis_tumor, "./02.processed_data/hepatocyte_nmf_w_basis.RDS")
saveRDS(h_coef_tumor, "./02.processed_data/hepatocyte_nmf_h_coef.RDS")   

nmf_programs_genes_tumor <- readRDS("./02.processed_data/hepatocyte_nmf_w_basis.RDS")
nmf_programs_cells_tumor <- readRDS("./02.processed_data/hepatocyte_nmf_h_coef.RDS")

# get gene programs (top 50 genes by NMF score)
nmf_programs_sig_tumor <- lapply(nmf_programs_genes_tumor, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
# select robust NMF programs
nmf_filter_tumor <- robust_nmf_programs(nmf_programs_sig_tumor, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)
nmf_programs_sig_tumor <- lapply(nmf_programs_sig_tumor, function(x) x[, is.element(colnames(x), nmf_filter_tumor),drop=F])
nmf_programs_sig_tumor <- do.call(cbind, nmf_programs_sig_tumor)
# calculate similarity between programs
nmf_intersect_tumor <- apply(nmf_programs_sig_tumor , 2, function(x) apply(nmf_programs_sig_tumor , 2, function(y) length(intersect(x,y)))) 
# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_tumor <- hclust(as.dist(50-nmf_intersect_tumor), method="average")
nmf_intersect_hc_tumor <- reorder(as.dendrogram(nmf_intersect_hc_tumor), colMeans(nmf_intersect_tumor))
nmf_intersect_tumor <- nmf_intersect_tumor[order.dendrogram(nmf_intersect_hc_tumor), order.dendrogram(nmf_intersect_hc_tumor)]
# plot similarity matrix heatmap     
nmf_intersect_meltI_tumor <- reshape2::melt(nmf_intersect_tumor) 
p1 <- ggplot(data = nmf_intersect_meltI_tumor, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  scale_x_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltI_tumor$Var1)[seq(1, length(unique(nmf_intersect_meltI_tumor$Var1)), by=5)], labels= seq(1, length(unique(nmf_intersect_meltI_tumor$Var1)), by=5)) + 
  scale_y_discrete(name="\nPrograms", breaks=unique(nmf_intersect_meltI_tumor$Var2)[seq(1, length(unique(nmf_intersect_meltI_tumor$Var2)), by=5)], labels= seq(1, length(unique(nmf_intersect_meltI_tumor$Var2)), by=5))  +
  geom_vline(xintercept=c(424,504), linetype="longdash", size=0.6)+
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
# for each program, plot correlation between nmf cell scores and cell complexity
complexity_tumor <- lapply(expr_tumor, function(x) apply(x, 2, function(y) length(which(y != 0))))
nmf_corr_complexity_tumor <- data.frame(matrix(ncol = 1, nrow = ncol(nmf_intersect_tumor)), row.names=colnames(nmf_intersect_tumor))
colnames(nmf_corr_complexity_tumor) <- "corr"                                     
for(i in rownames(nmf_corr_complexity_tumor)) {
  a <- gsub(".{4}$", "", i)
  b <- nmf_programs_cells_tumor[[a]][,i]
  c <- complexity_tumor[[a]]
  nmf_corr_complexity_tumor[i,1] <- cor(b,c)
}
p2 <- ggplot(nmf_corr_complexity_tumor, aes(y=corr, x = 1:nrow(nmf_corr_complexity_tumor))) +
  geom_smooth(span=0.1, se = FALSE, color="gray36", size= 0.8, method = "loess") + 
  geom_point(alpha=0.3, size = 1.5) + 
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "gray97"), panel.grid = element_blank(), axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text.y = element_text(size = 10), plot.margin=unit(c(1,1,-0.6,1), "cm") ) + 
  labs(y="Cell complexity\ncorrelation", x="") + 
  scale_y_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-0.8, 0.8, 0.4)) + 
  scale_x_continuous(expand = c(0,0))
pdf("./04.figures/04.Hepatocytes_nmf_metaprograms_complexity.pdf", height = 6, width = 8, onefile=FALSE) 
egg::ggarrange(p2,p1, nrow=2, ncol = 1, heights = c(1,5))
dev.off()
# manually remove programs associated with cell complexity and save output
nmf_intersect_tumor <- nmf_intersect_tumor[-c(1:3),-c(1:3)]
saveRDS(nmf_intersect_tumor, file = "./02.processed_data/hepatocyte_nmf_intersect_tumor.RDS")
saveRDS(nmf_programs_sig_tumor[,colnames(nmf_intersect_tumor)], file = "./02.processed_data/hepatocyte_nmf_programs_sig.RDS")
# hierarchical clustering of the filtered similarity matrix 
nmf_intersect_hc_tumor <- hclust(as.dist(50-nmf_intersect_tumor), method="average")
nmf_intersect_hc_tumor <- reorder(as.dendrogram(nmf_intersect_hc_tumor), colMeans(nmf_intersect_tumor))
nmf_intersect_tumor <- nmf_intersect_tumor[order.dendrogram(nmf_intersect_hc_tumor), order.dendrogram(nmf_intersect_hc_tumor)]
nmf_intersect_melt_tumor <- melt(nmf_intersect_tumor)
p1 <- ggplot(data = nmf_intersect_melt_tumor, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")+  
  theme(axis.ticks = element_blank(), panel.background = element_blank(), 
        panel.border=element_rect(colour = "gray40", size = 0.4, fill=F),  
        axis.line = element_blank(), axis.title.x=element_blank(), 
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10), 
        legend.title = element_text(size=10), legend.text = element_text(size = 9),
        legend.text.align = 0.5, legend.direction = "horizontal", 
        legend.justification="bottom", plot.margin = unit(c(0.5,3,0.5,0.5), "cm")) + 
  scale_x_discrete(name="\nPrograms", breaks=unique(nmf_intersect_melt_tumor$Var2)[seq(1, length(unique(nmf_intersect_melt_tumor$Var2)), by=5)], labels= seq(1, length(unique(nmf_intersect_melt_tumor$Var2)), by=5)) +
  scale_y_discrete(name="\nPrograms", breaks=unique(nmf_intersect_melt_tumor$Var2)[seq(1, length(unique(nmf_intersect_melt_tumor$Var2)), by=5)], labels= seq(1, length(unique(nmf_intersect_melt_tumor$Var2)), by=5)) +
  guides(color = guide_colourbar(barheight = 0.6, barwidth = 4.4, title.position = "top", title.hjust = 0.5)) 
# plot similarity matrix cancer type annotation
nmf_annot_type_tumor <- data.frame("cell_lines" =gsub(".{4}$", "",rownames(nmf_intersect_tumor)), stringsAsFactors = F)
nmf_annot_type_tumor$Sample <- factor(stringr::str_split_fixed(nmf_annot_type_tumor$cell_lines, "_", 3)[,2])
sample_color <- data.frame(
  "Sample" = c("P7","P8","P9","P15"), 
  "Color"= c("#CAA57D","#9E6BAB","#AFB2B7","#588198"),
  stringsAsFactors = F
)
row.names(sample_color) <- sample_color$Sample
p2 <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=Sample, color=Sample)) +
  geom_tile() +
  scale_color_manual(values = sample_color[match(levels(nmf_annot_type_tumor$Sample), sample_color$Sample), "Color"] , name="") +
  scale_fill_manual(values = sample_color[match(levels(nmf_annot_type_tumor$Sample), sample_color$Sample), "Color"] , name="") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank(), legend.position=c(1.3,-4), plot.margin = unit(c(0.5,3,-0.6,0.5), "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = guide_legend(title = "", override.aes = list(size = 5), ncol = 1), color=FALSE) 
pdf("./04.figures/04.Hepatocytes_nmf_metaprograms_cluster.pdf", height = 6, width = 8.5, onefile=FALSE)
egg::ggarrange(p2,p1, nrow = 2, heights = c(1,50))
dev.off()
# GO pathway enrich and annotation for each module 
module_genes <- list(
  M1 = sort(table(nmf_programs_sig_tumor[,as.character(unique(nmf_intersect_melt_tumor$Var1))[1:4]])),
  M2 = sort(table(nmf_programs_sig_tumor[,as.character(unique(nmf_intersect_melt_tumor$Var1))[5:10]])),
  M3 = sort(table(nmf_programs_sig_tumor[,as.character(unique(nmf_intersect_melt_tumor$Var1))[11:15]])),
  M4 = sort(table(nmf_programs_sig_tumor[,as.character(unique(nmf_intersect_melt_tumor$Var1))[16:19]])),
  M5 = sort(table(nmf_programs_sig_tumor[,as.character(unique(nmf_intersect_melt_tumor$Var1))[20:23]])),
  M6 = sort(table(nmf_programs_sig_tumor[,as.character(unique(nmf_intersect_melt_tumor$Var1))[24:30]])),
  M7 = sort(table(nmf_programs_sig_tumor[,as.character(unique(nmf_intersect_melt_tumor$Var1))[31:34]])),
  M8 = sort(table(nmf_programs_sig_tumor[,as.character(unique(nmf_intersect_melt_tumor$Var1))[35:39]]))
)
genes_used_list <- list()
for(i in names(module_genes)){
  genes_used_list[[i]] <- names(module_genes[[i]])[module_genes[[i]] >= 2]
}
save(genes_used_list, file = "./02.processed_data/hepatocyte_nmf_programs_genes.RDS")
library(gprofiler2)
gostres <- gost(query = genes_used_list,
  organism = "hsapiens", ordered_query = FALSE,
  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
  measure_underrepresentation = FALSE, evcodes = TRUE,
  user_threshold = 0.05, correction_method = "g_SCS",
  domain_scope = "annotated", custom_bg = NULL,
  numeric_ns = "", sources = NULL, as_short_link = FALSE)
plot_list <- list()
for(query_used in unique(gostres$result$query)){
  plot_list[[query_used]] <- 
    gostres$result %>% filter(query == query_used) %>% arrange(p_value) %>% head(15) %>% 
    ggplot(aes(x = -log10(p_value), y = reorder(term_name, -p_value))) +
    geom_bar(stat = "identity", aes(fill = source)) +
    labs(y = "") + ggtitle(query_used) +
    scale_y_discrete(guide = ggplot2::guide_axis(n.dodge = 2), 
                     labels = function(x) stringr::str_wrap(x, width = 20)) +
    scale_fill_nejm() +
    theme_cowplot()
}
p <- plot_grid(plotlist = plot_list)
ggsave(p, file = "./04.figures/04.Hepatocytes_program_enrich.pdf", width = 25, height = 15)
# Calculate the score of each module in tumor cells
tumor_cells <- HCC_seu@meta.data %>% filter(Global_Cluster == "Hepatocyte") %>% pull(CellName)
HCC_tumor <- HCC_seu %>% subset(CellName %in% tumor_cells)
for(i in names(genes_used_list)){
  HCC_tumor@meta.data[,paste0(i,"_Sig")] <- colMeans(HCC_tumor@assays$RNA@data[genes_used_list[[i]],])
}
p <- 
  HCC_tumor@meta.data[,c("Sub_Cluster","PatientID","ViralvsNonViral",paste0("M",1:8,"_Sig"))] %>% melt(measure.vars = paste0("M",1:8,"_Sig")) %>% 
  ggplot(aes(x = factor(Sub_Cluster), y = value)) +
  geom_boxplot(fill = "#9FA0A3", outlier.shape = 1, alpha = .8) +
  facet_wrap(~variable, scales = "free", nrow = 2) +
  theme_cowplot() +
  labs(y = "Expression", x = "") +
  theme_cowplot(font_size = 9) +
  theme(strip.text.x = element_text(size = 9),
        strip.background = element_rect(fill=NA, color=NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(p, file = "./04.figures/04.Hepatocytes_nmf_signature_score.pdf", width = 8, height = 6)

# Cell cycle score----
hepato_seu@meta.data$M1_sig <- colMeans(hepato_seu@assays$RNA@data[genes_used_list$M1,])
hepato_seu@meta.data$hepato_UMAP_1 <- hepato_seu@reductions$umap@cell.embeddings[,1]
hepato_seu@meta.data$hepato_UMAP_2 <- hepato_seu@reductions$umap@cell.embeddings[,2]
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p <- ggplot(hepato_seu@meta.data %>% arrange(M1_sig), aes(x = hepato_UMAP_1, y = hepato_UMAP_2)) +
  geom_point(aes(color = M1_sig), alpha = 0.8, size = 1) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_colour_gradientn("M1 Signature", colors = myColorPalette(100)) +
  theme_cowplot(font_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 12))
p <- ggAIplot(p)
ggsave(p, file = "./04.figures/04.Hepatocytes_M1_signature_umap.pdf", width = 8.5, height = 6.5)

# Calculate Module score and cluster----
module_genes_list <- readRDS("./02.processed_data/hepatocyte_nmf_programs_genes.RDS")
for(module in names(module_genes_list)){
  cat(module,"\n")
  hepato_seu@meta.data[,paste0(module,"_Score")] <- 
    CalculateModuleScore(hepato_seu@assays$RNA@data, module_genes_list[[module]])
}

p <- 
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
ggsave(p, file = "./04.figures/04.Hepatocyte_module_score.pdf", width = 8, height = 6)

umap_by_module_score <- umap::umap(hepato_seu@meta.data[,paste0("M",c(1:5,7,8),"_Score")], min_dist = .5, spread = 1)
hepato_seu@meta.data[,"UMAP_1_Module"] <- umap_by_module_score$layout[,1]
hepato_seu@meta.data[,"UMAP_2_Module"] <- umap_by_module_score$layout[,2]
p <- 
  ggplot(hepato_seu@meta.data, aes(x = UMAP_1_Module, y = UMAP_2_Module)) +
  geom_point(aes(color = Sub_Cluster), size = 1.2, alpha = .8) +
  scale_color_manual(name = "", values = Hepatocyte_color_panel) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))
p.pdf <- ggAIplot(p)
ggsave(p.pdf, file = "./04.figures/04.Hepatocyte_module_based_UMAP.pdf", width = 8, height = 6.5)

hepato_seu.plot.data <- hepato_seu@meta.data[,c("UMAP_1_Module","UMAP_2_Module",paste0("M",1:8,"_Score"))]
hepato_seu.plot.data <- melt(hepato_seu.plot.data, id.vars = c("UMAP_1_Module","UMAP_2_Module"))
myColorPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
p <- 
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
p.pdf <- ggAIplot.grid(p, "variable")
ggsave(p.pdf, file = "./04.figures/04.Hepatocyte_module_score_UMAP.pdf", width = 22, height = 10.5)

saveRDS(hepato_seu, file = "./02.processed_data/hepatocyte_seu.rds")

# # >Calcualte cluster proportions----
# sub_cluster_num <- data.frame(table(HCC_seu$Sub_Cluster, HCC_seu$PatientID, HCC_seu$NTF))
# colnames(sub_cluster_num) <- c("Sub_Cluster", "PatientID", "Source", "nSub")
# global_cluster_num <- data.frame(table(HCC_seu$Global_Cluster, HCC_seu$PatientID, HCC_seu$NTF))
# colnames(global_cluster_num) <- c("Global_Cluster", "PatientID", "Source", "nGlobal")
# CD45_num <- data.frame(table(HCC_seu$CD45, HCC_seu$PatientID, HCC_seu$NTF))
# colnames(CD45_num) <- c("CD45", "PatientID", "Source", "nCD45")
# patient_num <- data.frame(table(HCC_seu$PatientID, HCC_seu$NTF))
# colnames(patient_num) <- c("PatientID", "Source", "nPatient")
# cluster_info <- unique(HCC_seu@meta.data[,c("Sub_Cluster","Global_Cluster")])
# cluster_info$CD45 <- plyr::revalue(
#   cluster_info$Global_Cluster,
#   replace = c("B cell" = "CD45pos", "Doublet" = "CD45pos", "Endothelium" = "CD45neg",
#               "Fibroblast" = "CD45neg", "Hepatocyte" = "CD45neg", "ILC" = "CD45pos",
#               "Mast" = "CD45pos", "Mononuclear" = "CD45pos", "T cell" = "CD45pos")
# )
# cluster_prop_df <- merge(sub_cluster_num, cluster_info)
# cluster_prop_df <- merge(cluster_prop_df, global_cluster_num)
# cluster_prop_df <- merge(cluster_prop_df, CD45_num)
# cluster_prop_df <- merge(cluster_prop_df, patient_num)
# cluster_prop_df$pGlobal <- cluster_prop_df$nSub/cluster_prop_df$nGlobal
# cluster_prop_df$pCD45 <- cluster_prop_df$nSub/cluster_prop_df$nCD45
# cluster_prop_df$pPatient <- cluster_prop_df$nSub/cluster_prop_df$nPatient
# cluster_info$id <- paste0(cluster_info$Sub_Cluster,cluster_info$Global_Cluster)
# cluster_prop_df$id <- paste0(cluster_prop_df$Sub_Cluster, cluster_prop_df$Global_Cluster)
# cluster_prop_df <- cluster_prop_df[cluster_prop_df$id %in% cluster_info$id,]
# cluster_prop_df$id <- c()
# save(cluster_prop_df, file = "./Data/cluster_prop_df.rda")
# 
# # >Correlation of cluster proportions----
# cluster.prop.df <- cluster_prop_df %>% filter(Source == "Tumor", Global_Cluster %in% c("Endothelium","Fibroblast","Mononuclear"), PatientID != "HN", Sub_Cluster != "ABCAB+ Fib", Sub_Cluster != "MT1M+ Fib") %>% select(PatientID, Sub_Cluster, pPatient) %>% tidyr::spread(PatientID, pPatient)
# row.names(cluster.prop.df) <- cluster.prop.df$Sub_Cluster
# cluster.prop.df$Sub_Cluster <- c()
# cluster.prop.df$P10 <- c()
# cluster_used <- row.names(cluster.prop.df)[apply(cluster.prop.df,1,sd) != 0]
# cluster.prop.df <- cluster.prop.df[cluster_used,]
# cluster.prop.cor <- cor(t(cluster.prop.df))
# color_used <- circlize::colorRamp2(seq(min(cluster.prop.cor), max(cluster.prop.cor), length = 9), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:10]))
# pdf("./Figures/04.Cluster_proportion_corr.pdf", width = 7.5, height = 7)
# Heatmap(cluster.prop.cor,col = color_used,name = "corr")
# dev.off()

# >Nichenetr----
# Analysis----
# 1. Calculate DEGenes
HCC_seu$Global_Cluster_new <- HCC_seu$Global_Cluster
HCC_seu$Global_Cluster_new[HCC_seu$Global_Cluster_new %in% c("B cell","ILC","Mast")] <- "Other"
degenes <- list()
for(global_cluster_new in c("Endothelium","Fibroblast","Mononuclear","T cell","Other","Hepatocyte")){
  cells_used <- HCC_seu@meta.data %>% filter(Global_Cluster_new == global_cluster_new) %>% pull(CellName)
  HCC_seu_subset <- HCC_seu[,cells_used]
  FindDEGenes(
    expression_matrix = HCC_seu_subset@assays$RNA@data,
    groupid = HCC_seu_subset$Sub_Cluster,
    cutoff = 1,
    logFC = 0.5,
    out.prefix = paste0("./03.results/DEGenes/All/",global_cluster_new))
  for(sub_cluster in unique(HCC_seu_subset$Sub_Cluster)){
    cat(sub_cluster,"\n")
    degenes[[sub_cluster]] <- read.csv(paste0("./Results/DEGenes/All/",global_cluster_new,"_cutoff1_",sub_cluster,"_de_genes.csv"))
    # degenes[[sub_cluster]] <- FindMarkers(HCC_seu_subset, ident.1 = sub_cluster, group.by = "Sub_Cluster")
  }
}
saveRDS(degenes, file = "./03.results/DEGenes/all_cluster_degenes2.rds")

# 2. Run Nichenetr (Server)
degenes <- readRDS("./03.results/DEGenes/all_cluster_degenes2.rds")
nichenet_results <- list()
sub_clusters <- unique(HCC_seu$Sub_Cluster)
sub_clusters <- sub_clusters[sub_clusters != "Bi-Potent"]

for(cluster in sub_clusters){
  cat(cluster,"\n")
  genes <- degenes[[cluster]] %>% filter(Exp.Per.In > 0.3, Exp.Per.Out < 0.8, AUC > 0.5) %>% arrange(desc(AUC)) %>% pull(Symbol)
  if(length(genes) > 350) genes <- genes[1:350]
  nichenet_results[[cluster]] <-
    RunNichenetr(seu = HCC_seu,
                 cluster = cluster, group.by = "Sub_Cluster",
                 receiver_genes = genes, sender_genes = NULL,
                 best_upstream_ligands_number = 30, curated = T,
                 species = "homo", repo = TRUE)
  ggsave(nichenet_results[[cluster]]$combined_plot, file = paste0("./03.results/NicheNet/",cluster,"_curated_summary.pdf"), width = 15, height = 12)
}
save(nichenet_results, file = "./03.results/NicheNet/nichenet_results_curated.rda")

# Build network----
load("./03.results/NicheNet/nichenet_results_curated.rda")
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
ligand_enrich[,c("Cluster","POSTN+ CAF")] %>% filter(Cluster %in% cluster_used) %>% merge(.,cluster_info, by.x = "Cluster", by.y = "Sub_Cluster") %>% arrange(Global_Cluster, desc(`POSTN+ CAF`)) %>% setNames(c("Sub_Cluster","MeanExp","Global_Cluster")) %>% 
  ggplot(., aes(x = reorder(Sub_Cluster,-MeanExp), y = MeanExp, fill = Global_Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Top20 Ligands of POSTN+ CAF") +
  scale_fill_manual(values = Global_Cluster_color_panel) +
  facet_wrap(~Global_Cluster, nrow = 1, scales = "free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p1

cluster_used <- cluster_info %>% filter(Global_Cluster %in% c("Fibroblast","Mononuclear","Hepatocyte")) %>% pull(Sub_Cluster)
ligand_enrich[,c("Cluster","PLVAP+ EC_9")] %>% filter(Cluster %in% cluster_used) %>% merge(.,cluster_info, by.x = "Cluster", by.y = "Sub_Cluster") %>% arrange(Global_Cluster, desc(`PLVAP+ EC_9`)) %>% setNames(c("Sub_Cluster","MeanExp","Global_Cluster")) %>% 
  ggplot(., aes(x = reorder(Sub_Cluster,-MeanExp), y = MeanExp, fill = Global_Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Top20 Ligands of PLVAP+ EC_9") +
  scale_fill_manual(values = Global_Cluster_color_panel) +
  facet_wrap(~Global_Cluster, nrow = 1, scales = "free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p2

cluster_used <- cluster_info %>% filter(Global_Cluster %in% c("Fibroblast","Endothelium","Hepatocyte")) %>% pull(Sub_Cluster)
ligand_enrich[,c("Cluster","FOLR2+ TAM1_7")] %>% filter(Cluster %in% cluster_used) %>% merge(.,cluster_info, by.x = "Cluster", by.y = "Sub_Cluster") %>% arrange(Global_Cluster, desc(`FOLR2+ TAM1_7`)) %>% setNames(c("Sub_Cluster","MeanExp","Global_Cluster")) %>% 
  ggplot(., aes(x = reorder(Sub_Cluster,-MeanExp), y = MeanExp, fill = Global_Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Top20 Ligands of FOLR2+ TAM1_7") +
  scale_fill_manual(values = Global_Cluster_color_panel) +
  facet_wrap(~Global_Cluster, nrow = 1, scales = "free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p3

cluster_used <- cluster_info %>% filter(Global_Cluster %in% c("Fibroblast","Endothelium","Hepatocyte")) %>% pull(Sub_Cluster)
ligand_enrich[,c("Cluster","FOLR2+ TAM1_1")] %>% filter(Cluster %in% cluster_used) %>% merge(.,cluster_info, by.x = "Cluster", by.y = "Sub_Cluster") %>% arrange(Global_Cluster, desc(`FOLR2+ TAM1_1`)) %>% setNames(c("Sub_Cluster","MeanExp","Global_Cluster")) %>% 
  ggplot(., aes(x = reorder(Sub_Cluster,-MeanExp), y = MeanExp, fill = Global_Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Top20 Ligands of FOLR2+ TAM1_1") +
  scale_fill_manual(values = Global_Cluster_color_panel) +
  facet_wrap(~Global_Cluster, nrow = 1, scales = "free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> p4

p <- plot_grid(p1,p3,p2,p4)
ggsave(p, file = "./04.figures/04.Nichenet_Top20_ligands_enrichment.pdf", width = 17, height = 8)

# Expression of prioritized ligands----
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
ggsave(p, file = "./04.figures/04.Prioritized_ligands_bubble_heatmap.pdf", width = 12, height = 6)

# Shared ligands----
library(ggvenn)
temp <- list(
  "TAM" = ligands_to_plot %>% filter(Cluster == "FOLR2+ TAM") %>% pull(Gene) %>% unique(),
  "EC" = ligands_to_plot %>% filter(Cluster == "PLVAP+ EC") %>% pull(Gene) %>% unique(),
  "CAF" = ligands_to_plot %>% filter(Cluster == "POSTN+ CAF") %>% pull(Gene) %>% unique()
)
p <- ggvenn(temp, fill_color = c("#E3618D","#3D98BE","#D9B72E"), fill_alpha = .7, show_percentage = TRUE)
ggsave(p, file = "./04.figures/04.Common_ligands_venn.pdf", width = 10, height = 10)

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
genes_used %>% filter(Group == "TAM_EC_CAF")
genes_used_list <- list(
  common = genes_used %>% filter(Group %in% c("EC_CAF","TAM_CAF","TAM_EC_","TAM_EC_CAF")) %>% arrange(Group),
  specific = genes_used %>% filter(Group %in% c("EC_","TAM_","CAF")) %>% arrange(Group)
)
figure_list <- list()
for(genes_type in names(genes_used_list)){
  genes_used <- genes_used_list[[genes_type]]
  mean_exp <- aggregate(t(HCC_seu@assays$RNA@data[genes_used$Gene,HCC_seu$Global_Cluster %ni% c("Doublet")]),list(Cluster = HCC_seu$Sub_Cluster[HCC_seu$Global_Cluster %ni% c("Doublet")]),mean)
  row.names(mean_exp) <- mean_exp$Cluster
  mean_exp$Cluster <- c()
  markers_plot_matrix <- apply(mean_exp, 2, zscore) %>% t()
  markers_plot_matrix_quantile <- quantile(markers_plot_matrix, c(0.01, 0.98))
  markers_plot_matrix <- pmax(markers_plot_matrix, markers_plot_matrix_quantile[1])
  markers_plot_matrix <- pmin(markers_plot_matrix, markers_plot_matrix_quantile[2])
  
  markers_plot_matrix <- markers_plot_matrix[,cluster_info$Sub_Cluster]
  color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), max(markers_plot_matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
  ha <- HeatmapAnnotation(
    Group = genes_used$Group,
    col = list(Group = c("CAF"="#7876B1FF","EC_"="#6F99ADFF","TAM_"="#FFDC91FF","EC_CAF"="#BC3C29FF","TAM_CAF"="#0072B5FF","TAM_EC_"="#E18727FF","TAM_EC_CAF"="#20854EFF"))
  )
  ha_row <-  rowAnnotation(
    df = data.frame(Global_Cluster = cluster_info$Global_Cluster), 
    col = list(Global_Cluster = Global_Cluster_color_panel), 
    width = unit(0.5, "cm"))
  figure_list[[genes_type]] <- 
    Heatmap(t(markers_plot_matrix),
            cluster_rows = FALSE, cluster_columns = FALSE,
            top_annotation = ha,
            left_annotation = ha_row,
            column_names_gp = gpar(fontsize = 7),
            row_names_gp = gpar(fontsize = 7),
            col = color_used, name = "Exp",
            column_split = genes_used$Group,
            row_split = cluster_info$Global_Cluster)
}
pdf("./04.figures/04.Nichenet_Top_ligands_expression.pdf", width = 20, height = 6)
draw(figure_list$specific + figure_list$common, ht_gap = unit(1,"cm"))
dev.off()

# Shared signatures----
TAM_cells_used <- HCC_seu@meta.data %>% filter(Sub_Cluster %in% c("FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2","MT1G+ TAM3")) %>% pull(CellName)
TAM_DEGenes <- LIMMA(
  expression_matrix = HCC_seu@assays$RNA@data[,TAM_cells_used],
  groupid = HCC_seu@meta.data[TAM_cells_used,"Sub_Cluster"] %in% c("FOLR2+ TAM1_1","FOLR2+ TAM1_7")
  )
CAF_cells_used <- HCC_seu@meta.data %>% filter(Sub_Cluster %in% c("POSTN+ CAF","MYH11+ CAF")) %>% pull(CellName)
CAF_DEGenes <- LIMMA(
  expression_matrix = HCC_seu@assays$RNA@data[,CAF_cells_used],
  groupid = HCC_seu@meta.data[CAF_cells_used,"Sub_Cluster"] == "POSTN+ CAF"
)
EC_cells_used <- HCC_seu@meta.data %>% filter(Sub_Cluster %in% c("PLVAP+ EC_3","PLPP3+ EC","PLVAP+ EC_9","PLVAP+ EC_4","IGFBP3+ EC")) %>% pull(CellName)
EC_DEGenes <- LIMMA(
  expression_matrix = HCC_seu@assays$RNA@data[,EC_cells_used],
  groupid = HCC_seu@meta.data[EC_cells_used,"Sub_Cluster"] %in% c("PLVAP+ EC_3","PLVAP+ EC_4","PLVAP+ EC_9")
)
intersect(CAF_degenes, EC_degenes)
intersect(CAF_degenes, TAM_degenes)
intersect(TAM_degenes, EC_degenes)
intersect(TAM_degenes, intersect(CAF_degenes, EC_degenes))
temp <- list(
  "TAM" = TAM_DEGenes %>% filter(Grp == "C_TRUE", adj.P.Val < 0.05, logFC < -0.5) %>% arrange(logFC) %>% pull(Symbol),
  "EC" = EC_DEGenes %>% filter(Grp == "C_TRUE", adj.P.Val < 0.05, logFC < -0.5) %>% arrange(logFC) %>% pull(Symbol),
  "CAF" = CAF_DEGenes %>% filter(Grp == "C_TRUE", adj.P.Val < 0.05, logFC < -0.5) %>% arrange(logFC) %>% pull(Symbol)
)
p <- ggvenn(temp, fill_color = c("#E3618D","#3D98BE","#D9B72E"), fill_alpha = .7, show_percentage = TRUE)
ggsave(p, file = "./04.figures/04.Common_upregulation_signatures_venn.pdf", width = 10, height = 10)

# >CellChat----
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
saveRDS(cellchat, file = "./03.results/CellChat/cellchat_HCC.rds")
cellchat <- readRDS("./03.results/CellChat/cellchat_HCC.rds")

# Interaction Strength-----
color_used <- rep("#4385BF",40)
color_used[c(13,14,33,34,35,36)] <- "#F16592"
p <- netAnalysis_signalingRole_scatter(cellchat, color.use = color_used)
ggsave(p, file = "./04.figures/04.CellChat_Interaction_strength.pdf", width = 7, height = 6)

# Singaling clustering----
p <- netVisual_embedding(cellchat, type = "structural", label.size = 3.5, color.use = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF"))
ggsave(p, file = "./04.figures/04.CellChat_Signaling_structural_cluster.pdf", width = 7, height = 6)

# Bubble heatmap----
macrophage_cluster_used <- c("FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2")
CAF_cluster_used <- c("POSTN+ CAF","MYH11+ CAF")
EC_cluster_used <- c("PLVAP+ EC_4","PLVAP+ EC_9")
Tumor_cluster_used <- c("Hepatocytes_P7_2","Hepatocytes_P9_1","Hepatocytes_P15")
T_cluster_used <- c("CD4+ T","CD8+ T","Tregs")
cluster_combination_used <- c()
for(i in Tumor_cluster_used){
  for(j in c(CAF_cluster_used,EC_cluster_used,macrophage_cluster_used)){
    cluster_combination1 <- paste0(i,"|",j)
    cluster_combination_used <- c(cluster_combination_used, cluster_combination1)
  }
}
for(i in CAF_cluster_used){
  for(j in c(EC_cluster_used,macrophage_cluster_used,T_cluster_used)){
    cluster_combination1 <- paste0(i,"|",j)
    cluster_combination2 <- paste0(j,"|",i)
    cluster_combination_used <- c(cluster_combination_used, cluster_combination1,cluster_combination2)
  }
}
for(i in EC_cluster_used){
  for(j in c(CAF_cluster_used,macrophage_cluster_used)){
    cluster_combination1 <- paste0(i,"|",j)
    cluster_combination2 <- paste0(j,"|",i)
    cluster_combination_used <- c(cluster_combination_used, cluster_combination1,cluster_combination2)
  }
}
for(i in macrophage_cluster_used){
  for(j in c(CAF_cluster_used,EC_cluster_used)){
    cluster_combination1 <- paste0(i,"|",j)
    cluster_combination2 <- paste0(j,"|",i)
    cluster_combination_used <- c(cluster_combination_used, cluster_combination1,cluster_combination2)
  }
}
cluster_combination_used <- unique(cluster_combination_used)
sources_used <- targets_used <- c(macrophage_cluster_used, EC_cluster_used, CAF_cluster_used, Tumor_cluster_used, T_cluster_used)
LR_used <- data.frame(interaction_name = c("ADGRE5_CD55","ANGPTL4_SDC2","AREG_EGFR","CCL14_CCR1","CCL19_CCR7","CCL21_CCR7","CCL8_CCR1","CD99_CD99L2","CSF1_CSF1R","CXCL12_CXCR4","CXCL16_CXCR6","CXCL2_ACKR1","CXCL3_ACKR1","CXCL8_ACKR1","DLL4_NOTCH3","EREG_EGFR","FGF7_FGFR1","GAS6_MERTK","HBEGF_EGFR","IGF1_IGF1R","IGF2_IGF2R","IL34_CSF1R","ITGB2_ICAM1","JAG1_NOTCH3","JAG1_NOTCH4","JAG2_NOTCH3","LAMB2_CD44","LGALS9_HAVCR2","MDK_SDC2","MIF_CD74_CD44","PGF_VEGFR1","PTN_SDC2","SEMA3C_NRP2_PLXNA2","SPP1_CD44","SPP1_CD44","THBS1_CD47","THBS1_SDC4","TNFSF10_TNFRSF10B","TNFSF12_TNFRSF12A","TNXB_SDC4","VEGFA_VEGFR2","VTN_ITGAV_ITGB5"))
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
min.norm <- 0
for(LRpair in as.character(unique(df.net$interaction_name))){
  prob.sig <- df.net[df.net$interaction_name == LRpair,"prob"]
  prob.sig.normalized <- scale(c(prob.sig,rep(0,length(cluster_combination_used)-length(prob.sig))))[1:length(prob.sig)]
  if(min.norm > min(prob.sig.normalized)){min.norm <- min(prob.sig.normalized)}
  df.net[df.net$interaction_name == LRpair,"prob"] <- prob.sig.normalized
}
df.net$pval[df.net$pval > 0.05] = 1
df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
df.net$pval[df.net$pval > 0.001 & df.net$pval <= 0.01] = 3
df.net$pval[df.net$pval <= 0.001] = 4
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
legend.values <- c(1,2,3,4)
names(legend.values) <- c("p > 0.05","0.01 < p < 0.05","0.001 < p < 0.01","P < 0.001")
p <- 
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
  scale_color_gradientn(name = "Commun. Prob.", colours = myColorPalette(100)[c(11:31,71:91)], limits = c(min.norm, max(df.net$prob)), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barwidth = 0.8))
ggsave(p, file = "./04.figures/04.CellChat_LR_dotplot.pdf", width = 15, height = 11)

# Ligand-receptor heatmap----
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
p <- Ligand.ht + Receptor.ht
pdf("./04.figures/04.LR_heatmap.pdf", width = 6, height = 3)
draw(p, ht_gap = unit(0.5, "cm"), auto_adjust = FALSE, merge_legend = TRUE)
dev.off()

# Chrod plot----
sub_cluster_used.list <- list(
  "VTN_ITGAV_ITGB5" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Hepatocyte","Fibroblast")) %>% select(Global_Cluster, Sub_Cluster) %>% arrange(Global_Cluster) %>% pull(Sub_Cluster) %>% unique(),
  "IL34_CSF1R" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Endothelium","Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% arrange(Global_Cluster) %>% pull(Sub_Cluster) %>% unique(),
  "CSF1_CSF1R" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Endothelium","Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% arrange(Global_Cluster) %>% pull(Sub_Cluster) %>% unique(),
  "CXCL12_CXCR4" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% arrange(Global_Cluster) %>% pull(Sub_Cluster) %>% unique(),
  "GAS6_MERTK" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% arrange(Global_Cluster) %>% pull(Sub_Cluster) %>% unique(),
  "DLL4_NOTCH3" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Endothelium","Fibroblast")) %>% select(Global_Cluster, Sub_Cluster) %>% arrange(Global_Cluster) %>% pull(Sub_Cluster) %>% unique(),
  "TNFSF10_TNFRSF10B" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Hepatocyte","Endothelium","Fibroblast")) %>% select(Global_Cluster, Sub_Cluster) %>% arrange(Global_Cluster) %>% pull(Sub_Cluster) %>% unique()
)
pdf("./04.figures/04.LR_chord_plot.pdf", width = 6.5, height = 6.5)
plot.list <- list()
for(LRpair in names(sub_cluster_used.list)){
  pathway <- cellchat@DB$interaction[LRpair,"pathway_name"]
  plot.list[[LRpair]] <- netVisual_individual_new(cellchat, signaling = pathway, pairLR.use = LRpair, clusters.use = sub_cluster_used.list[[LRpair]], layout = "chord")
}
dev.off()

# >CellPhoneDB----
# Significant LR pairs counts----
cluster_info <- unique(HCC_seu@meta.data[,c("Global_Cluster","Sub_Cluster")]) %>% arrange(Global_Cluster)
row.names(cluster_info) <- cluster_info$Sub_Cluster
cluster_info$Global_Cluster <- as.character(cluster_info$Global_Cluster)
cluster_info$Global_Cluster[cluster_info$Global_Cluster %in% c("B cell","T cell","ILC")] <- "Lymphocyte"
count_network <- read.table("./03.results/CellPhoneDB/Results/HCC_tumor_statistical_analysis_201127/out/count_network.txt", sep = "\t", head = T, stringsAsFactors = F)
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
  scale_color_manual(values = Global_Cluster_color_panel) +
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
  scale_color_manual(values = Global_Cluster_color_panel) +
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
  scale_color_manual(values = Global_Cluster_color_panel) +
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
  scale_color_manual(values = Global_Cluster_color_panel) +
  labs(y = "", x = "Number of Significant Interactions") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), legend.position = "none") -> p4
p <- plot_grid(p1,p2,p3,p4,nrow = 2)
ggsave(p, file = "./04.figures/04.CellPhoneDB_Significant_LR.pdf", width = 8, height = 8)

# Bubble heatmap----
macrophage_cluster_used <- c("FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2")
CAF_cluster_used <- c("POSTN+ CAF","MYL9+ CAF")
EC_cluster_used <- c("PLVAP+ EC_4","PLVAP+ EC_9")
Tumor_cluster_used <- c("Hepatocytes_P7_2","Hepatocytes_P9_1")
cluster_combination_used <- c()
for(i in Tumor_cluster_used){
  for(j in c(CAF_cluster_used,EC_cluster_used)){
    cluster_combination1 <- paste0(i,"|",j)
    cluster_combination_used <- c(cluster_combination_used, cluster_combination1)
  }
}
for(i in CAF_cluster_used){
  for(j in c(EC_cluster_used,macrophage_cluster_used)){
    cluster_combination1 <- paste0(i,"|",j)
    cluster_combination2 <- paste0(j,"|",i)
    cluster_combination_used <- c(cluster_combination_used, cluster_combination1,cluster_combination2)
  }
}
for(i in EC_cluster_used){
  for(j in c(CAF_cluster_used,macrophage_cluster_used)){
    cluster_combination1 <- paste0(i,"|",j)
    cluster_combination2 <- paste0(j,"|",i)
    cluster_combination_used <- c(cluster_combination_used, cluster_combination1,cluster_combination2)
  }
}
for(i in macrophage_cluster_used){
  for(j in c(CAF_cluster_used,EC_cluster_used)){
    cluster_combination1 <- paste0(i,"|",j)
    cluster_combination2 <- paste0(j,"|",i)
    cluster_combination_used <- c(cluster_combination_used, cluster_combination1,cluster_combination2)
  }
}
cluster_combination_used <- unique(cluster_combination_used)
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
mean <- read.table("./Results/CellPhoneDB/Results/HCC_tumor_statistical_analysis_201127/means.txt", sep = "\t", head = T, stringsAsFactors = F, check.names = F)
LR_info <- mean[,1:11]
cluster_combination_used <- intersect(cluster_combination_used, colnames(mean))
mean <- mean[,cluster_combination_used]
pvalue <- read.table("./Results/CellPhoneDB/Results/HCC_tumor_statistical_analysis_201127/pvalues.txt", sep = "\t", head = T, stringsAsFactors = F, check.names = F)
pvalue <- pvalue[,cluster_combination_used]
# LR_used <- rowSums(pvalue < 0.01) > 0
# LR_used <- LR_used & (apply(mean,1,max) > 0.5)
LR_used <- c("VTN_aVb1 complex","VEGFA_FLT1","TTR_NGFR","TNFSF10_TNFRSF10D","TNFSF10_RIPK1","THY1_aXb2 complex","TGFB3_TGFBR3","SPP1_CD44","PROS1_AXL","PLXNB2_PTN","PECAM1_CD38","NRP2_VEGFA","NOTCH3_JAG2","NOTCH2_DLL4","MERTK_GAS6","LGALS9_SLC1A5","LGALS9_HAVCR2","LAMC1_a6b1 complex","JAG1_NOTCH4","HLA-F_LILRB2","HLA-DRB1_OGN","HGF_CD44","FLT1_VEGFB","FBN1_a5b1 complex","EGFR_GRN","DPP4_CXCL12","DLL4_NOTCH3","CXCL12_CXCR4","CSF1R_IL34","CSF1R_CSF1","CD99_PILRA","CD74_MIF","CD46_JAG1","CD44_SELE","CD40_TNFSF13B","CCRL2_CCL19","CCR1_CCL14","CCL2_ACKR1","AXL_IL15RA","ANXA1_FPR1","ANGPT2_TEK","ADORA3_ENTPD1","ACVR_1B2A receptor_INHBB")
row.names(mean) <- row.names(pvalue) <- LR_info$interacting_pair
mean_norm <- apply(mean[LR_used,],1,scale)
row.names(mean_norm) <- colnames(mean)
mean_vec <- unlist(as.data.frame(t(mean_norm)))
pvalue_vec <- unlist(as.data.frame(pvalue[LR_used,]))
pvalue_vec[pvalue_vec==0] = 0.0009
# plot.data <- data.frame(
#   pair = rep(LR_info$interacting_pair[LR_used], length(cluster_combination_used)),
#   clusters = rep(cluster_combination_used, each = sum(LR_used)),
#   pvalue = pvalue_vec, mean = mean_vec
# )
plot.data <- data.frame(
  pair = rep(LR_used, length(cluster_combination_used)),
  clusters = rep(cluster_combination_used, each = length(LR_used)),
  pvalue = pvalue_vec, mean = mean_vec
)
plot.data$clusters <- 
  factor(plot.data$clusters, 
         levels = c(
           sort(grep("^Hepa",cluster_combination_used,value = T)),
           sort(grep("^POSTN",cluster_combination_used,value = T)),
           sort(grep("^MYL9",cluster_combination_used,value = T)),
           sort(grep("^PLVAP",cluster_combination_used,value = T)),
           sort(grep("^FOLR2",cluster_combination_used,value = T)),
           sort(grep("^SPP1",cluster_combination_used,value = T))
           )[c(1,4,2,3,5,8,6,7,9,10,13,11,12,14,15,18,16,17,19,20,23,21,22,24,25,28,26,27,29,32,30,31,33,36,34,35,37,40,38,39)]
         )
p <- 
  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_radius(range = c(0.5,4)) +
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
  scale_fill_gradientn(colours = c(myColorPalette(100)[1:30], myColorPalette(100)[71:100]), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barwidth = 0.8)) +
  scale_color_gradientn(colours = c(myColorPalette(100)[1:30], myColorPalette(100)[71:100]), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barwidth = 0.8))
ggsave(p, file = "./Figures/04.CellPhoneDB_LR_dotplot.pdf", width = 11, height = 10)

# Ligand-receptor circos plot----
LR.df <- data.frame(Ligand = c("TTR","NOTCH3","NOTCH2","IL34","CSF1","SLC1A5","CCL19","CXCL12"), Receptor = c("NGFR","JAG2","DLL4","CSF1R","CSF1R","LGALS9","CCRL2","CXCR4"))
row.names(LR.df) <- paste0(LR.df$Ligand,"-",LR.df$Receptor)
subcluster_num <- table(HCC_seu@meta.data$Sub_Cluster) %>% data.frame() %>% setNames(c("Sub_Cluster","Freq"))
row.names(subcluster_num) <- subcluster_num$Sub_Cluster
sub_cluster_used.list <- list(
  "TTR-NGFR" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Hepatocyte","Fibroblast")) %>% select(Global_Cluster, Sub_Cluster) %>% unique() %>% arrange(Global_Cluster) %>% pull(Sub_Cluster),
  "NOTCH3-JAG2" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Endothelium","Fibroblast")) %>% select(Global_Cluster, Sub_Cluster) %>% unique() %>% arrange(Global_Cluster) %>% pull(Sub_Cluster),
  "NOTCH2-DLL4" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Endothelium","Fibroblast")) %>% select(Global_Cluster, Sub_Cluster) %>% unique() %>% arrange(Global_Cluster) %>% pull(Sub_Cluster),
  "IL34-CSF1R" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Endothelium","Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% unique() %>% arrange(Global_Cluster) %>% pull(Sub_Cluster),
  "CSF1-CSF1R" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Endothelium","Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% unique() %>% arrange(Global_Cluster) %>% pull(Sub_Cluster),
  "SLC1A5-LGALS9" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Endothelium","Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% unique() %>% arrange(Global_Cluster) %>% pull(Sub_Cluster),
  "CCL19-CCRL2" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% unique() %>% arrange(Global_Cluster) %>% pull(Sub_Cluster),
  "CXCL12-CXCR4" = HCC_seu@meta.data %>% filter(Global_Cluster %in% c("Fibroblast", "Mononuclear")) %>% select(Global_Cluster, Sub_Cluster) %>% unique() %>% arrange(Global_Cluster) %>% pull(Sub_Cluster)
)
# Calculate interaction strength
interaction_strength.list <- list()
for(i in 1:nrow(LR.df)){
  ligand_exp <- aggregate(HCC_seu@assays$RNA@data[LR.df$Ligand[i],], list(Sub_Cluster = HCC_seu$Sub_Cluster), mean)
  ligand_pro <- aggregate(HCC_seu@assays$RNA@data[LR.df$Ligand[i],], list(Sub_Cluster = HCC_seu$Sub_Cluster), function(x){sum(x > 0)/length(x)})
  ligand_exp$x[ligand_pro$x < 0.1] <- 0
  receptor_exp <- aggregate(HCC_seu@assays$RNA@data[LR.df$Receptor[i],], list(Sub_Cluster = HCC_seu$Sub_Cluster), mean)
  receptor_pro <- aggregate(HCC_seu@assays$RNA@data[LR.df$Receptor[i],], list(Sub_Cluster = HCC_seu$Sub_Cluster), function(x){sum(x > 0)/length(x)})
  receptor_exp$x[receptor_pro$x < 0.1] <- 0
  interaction_strength <- ligand_exp$x %*% t(receptor_exp$x)
  row.names(interaction_strength) <- colnames(interaction_strength) <- ligand_exp$Sub_Cluster
  interaction_strength.list[[paste0(LR.df$Ligand[i],"-",LR.df$Receptor[i])]] <- interaction_strength
}
# data for plot
cluster_highlight.list <- list(
  "NOTCH3-JAG2" = c("POSTN+ CAF","MYH11+ CAF","PLVAP+ EC_9"),
  "NOTCH2-DLL4" = c("POSTN+ CAF","PLVAP+ EC_4","PLVAP+ EC_9"),
  "IL34-CSF1R" = c("POSTN+ CAF","FOLR2+ TAM1_1","FOLR2+ TAM1_7"),
  "CSF1-CSF1R" = c("PLVAP+ EC_9","FOLR2+ TAM1_1","FOLR2+ TAM1_7"),
  "SLC1A5-LGALS9" = c("PLVAP+ EC_9","FOLR2+ TAM1_1","FOLR2+ TAM1_7"),
  "CCL19-CCRL2" = c("POSTN+ CAF","FOLR2+ TAM1_1","FOLR2+ TAM1_7"),
  "CXCL12-CXCR4" = c("POSTN+ CAF","FOLR2+ TAM1_1","FOLR2+ TAM1_7")
)
from_cluster_highlight.list <- list(
  "NOTCH3-JAG2" = c("POSTN+ CAF","MYH11+ CAF"),
  "NOTCH2-DLL4" = c("POSTN+ CAF"),
  "IL34-CSF1R" = c("POSTN+ CAF"),
  "CSF1-CSF1R" = c("PLVAP+ EC_9"),
  "SLC1A5-LGALS9" = c("PLVAP+ EC_9"),
  "CCL19-CCRL2" = c("POSTN+ CAF"),
  "CXCL12-CXCR4" = c("POSTN+ CAF")
)
to_cluster_highlight.list <- list(
  "NOTCH3-JAG2" = c("PLVAP+ EC_9"),
  "NOTCH2-DLL4" = c("PLVAP+ EC_4","PLVAP+ EC_9"),
  "IL34-CSF1R" = c("FOLR2+ TAM1_1","FOLR2+ TAM1_7"),
  "CSF1-CSF1R" = c("FOLR2+ TAM1_1","FOLR2+ TAM1_7"),
  "SLC1A5-LGALS9" = c("FOLR2+ TAM1_1","FOLR2+ TAM1_7"),
  "CCL19-CCRL2" = c("FOLR2+ TAM1_1","FOLR2+ TAM1_7"),
  "CXCL12-CXCR4" = c("FOLR2+ TAM1_1","FOLR2+ TAM1_7")
)
links.df <- initialize.df <- list()
for(LRpair in names(cluster_highlight.list)){
  cluster_used <- sub_cluster_used.list[[LRpair]]
  interaction_used <- interaction_strength.list[[LRpair]][cluster_used,cluster_used]
  normalize_factor <- min(subcluster_num[cluster_used,"Freq"] / rowSums(interaction_used),
                          t(subcluster_num[cluster_used,"Freq"] / t(interaction_used)))
  interaction_used_normalized <- interaction_used * normalize_factor
  interaction_used_normalized <- interaction_used_normalized[rowSums(interaction_used_normalized) != 0, colSums(interaction_used_normalized) != 0]
  interaction_used_normalized <- t(apply(interaction_used_normalized, 1, cumsum))
  interaction_used_normalized <- cbind(Start = 0, interaction_used_normalized)
  
  links.df[[LRpair]] <- data.frame(from = NA, from.start = NA, from.end = NA, to = NA, to.start = NA, to.end = NA, stringsAsFactors = F)
  for(i in 1:nrow(interaction_used_normalized)){
    from = row.names(interaction_used_normalized)[i]
    for(j in 2:ncol(interaction_used_normalized)){
      to = colnames(interaction_used_normalized)[j]
      links.df[[LRpair]] <- rbind(links.df[[LRpair]], c(from, interaction_used_normalized[i,(j-1)], interaction_used_normalized[i,j], to, 0, interaction_used_normalized[i,j] - interaction_used_normalized[i,(j-1)]))
    }
  }
  links.df[[LRpair]]$from.factors <- links.df[[LRpair]]$from
  links.df[[LRpair]]$to.factors <- links.df[[LRpair]]$to
  links.df[[LRpair]] <- links.df[[LRpair]][-1,]
  links.df[[LRpair]][,c(2,3,5,6)] <- apply(links.df[[LRpair]][,c(2,3,5,6)], 2, as.numeric)
  links.df[[LRpair]]$color <- "lightgrey"
  links.df[[LRpair]][links.df[[LRpair]]$to.factors %in% to_cluster_highlight.list[[LRpair]] & links.df[[LRpair]]$from.factors %in% from_cluster_highlight.list[[LRpair]],"color"] <- "#C04146"
  links.df[[LRpair]] <- links.df[[LRpair]] %>% arrange(desc(color))
  
  initialize.df[[LRpair]] <- data.frame(Sub_Cluster = cluster_used, xmin = 0, xmax = subcluster_num[cluster_used,"Freq"], color = "lightgrey", stringsAsFactors = F) %>% mutate(factors = Sub_Cluster)
  initialize.df[[LRpair]][initialize.df[[LRpair]]$factors %in% cluster_highlight.list[[LRpair]],"color"] <- "#C04146"
  ###
  for(i in initialize.df[[LRpair]]$Sub_Cluster){
    if(sum(links.df[[LRpair]]$from == i) > 0){
      links.df[[LRpair]][links.df[[LRpair]]$from == i,"from.start"] <- links.df[[LRpair]][links.df[[LRpair]]$from == i,"from.start"]/max(links.df[[LRpair]][links.df[[LRpair]]$from == i,"from.start"])*initialize.df[[LRpair]][initialize.df[[LRpair]]$Sub_Cluster == i,"xmax"]
      links.df[[LRpair]][links.df[[LRpair]]$from == i,"from.end"] <- links.df[[LRpair]][links.df[[LRpair]]$from == i,"from.end"]/max(links.df[[LRpair]][links.df[[LRpair]]$from == i,"from.end"])*initialize.df[[LRpair]][initialize.df[[LRpair]]$Sub_Cluster == i,"xmax"]
    }
  }
}
save(links.df, initialize.df, file = "./Data/circos_plot_data.rda")

# Circos plot
library(circlize)
pdf("./Figures/04.LR_circos_plot.pdf", width = 10, height = 10)
for(LRpair in names(links.df)){
  circos.par(cell.padding = c(0.02, 0, 0.02, 0))
  circos.initialize(factors = initialize.df[[LRpair]]$Sub_Cluster, xlim = initialize.df[[LRpair]][,c(2,3)])
  for(i in 1:3){
    circos.track(factors = initialize.df[[LRpair]]$Sub_Cluster, ylim = c(0,1), bg.border = NA)
  }
  circos.track(factors = initialize.df[[LRpair]]$Sub_Cluster, ylim = c(0,1), 
               bg.border = NA, bg.col = initialize.df[[LRpair]]$color, track.height = 0.05,
               panel.fun = function(x, y) {
                 circos.text(x = CELL_META$xcenter, 
                             y = CELL_META$cell.ylim[2] + uy(5, "mm"),
                             labels = CELL_META$sector.index, facing = "clockwise", cex = 1,
                             niceFacing = T, adj = c(0,0.5))
               }
  )
  for(i in 1:nrow(links.df[[LRpair]])){
    circos.link(links.df[[LRpair]]$from[i],
                c(links.df[[LRpair]]$from.start[i], links.df[[LRpair]]$from.end[i]),
                links.df[[LRpair]]$to[i],
                c(links.df[[LRpair]]$to.start[i], links.df[[LRpair]]$to.end[i]),
                col = links.df[[LRpair]]$color[i],
                directional = 1, arr.length = 0.1
    )
  }
  title(LRpair, cex.main = 1.5)
  circos.clear()
}
dev.off()

# >IL34 vs CSF1 induced siganutres----
bulk_degenes <- read.table("../Published_data/IL34_bulk/FinalDataMatrix.txt", sep = "\t", head = T)
bulk_degenes <- bulk_degenes[-1,]
bulk_degenes <- bulk_degenes[!duplicated(bulk_degenes$Scan.REF),]
row.names(bulk_degenes) <- bulk_degenes$Scan.REF
bulk_degenes$Scan.REF <- c()
colnames(bulk_degenes) <- c("IL34_untreated_d3","mCSF_untreated_d1","IL34_untreated_d7","mCSF_untreated_d3","mCSF_untreated_d5","mCSF_untreated_d7","mCSF_treated_d7","IL34_untreated_d1","IL34_untreated_d5","IL34_treated_d7")
gene_symbol <- read.table("../Published_data/IL34_bulk/A-AGIL-28.adf.txt", sep = "\t")
gene_symbol <- gene_symbol[!duplicated(gene_symbol$V1) & !duplicated(gene_symbol$V2),]
row.names(gene_symbol) <- gene_symbol$V1
bulk_degenes <- bulk_degenes[row.names(gene_symbol),]
row.names(bulk_degenes) <- gene_symbol$V2
bulk_degenes[,] <- apply(bulk_degenes[,],2,function(x){as.numeric(x)})
bulk_degenes <- bulk_degenes[,c("IL34_untreated_d1","IL34_untreated_d3","IL34_untreated_d5","IL34_untreated_d7","IL34_treated_d7","mCSF_untreated_d1","mCSF_untreated_d3","mCSF_untreated_d5","mCSF_untreated_d7","mCSF_treated_d7")]
degenes <- data.frame(
  IL34_minus_mCSF = bulk_degenes$IL34_untreated_d7 - bulk_degenes$mCSF_untreated_d7,
  IL34_treat = bulk_degenes$IL34_untreated_d7 - bulk_degenes$IL34_treated_d7,
  mCSF_treat = bulk_degenes$mCSF_untreated_d7 - bulk_degenes$mCSF_treated_d7,
  Gene = row.names(bulk_degenes),
  row.names = row.names(bulk_degenes)
)
IL34_induced_genes <- degenes %>% filter(IL34_minus_mCSF > 2, abs(IL34_treat) > 1) %>% arrange(desc(IL34_minus_mCSF)) %>% pull(Gene)
mCSF_induced_genes <- degenes %>% filter(IL34_minus_mCSF < -2, abs(mCSF_treat) > 1) %>% arrange(IL34_minus_mCSF) %>% pull(Gene)

bulk_degenes2 <- read.csv("../Published_data/IL34_bulk/E0079_lfc0.5.csv")
bulk_degenes2 <- bulk_degenes2[!duplicated(bulk_degenes2$ILMN_GENE),2:4]
row.names(bulk_degenes2) <- bulk_degenes2$ILMN_GENE
IL34_induced_genes2 <- bulk_degenes2 %>% filter(CSF.IL34.lfc < 0, CSF.IL34.pv < 0.05) %>% arrange(CSF.IL34.lfc) %>% pull(ILMN_GENE)
mCSF_induced_genes2 <- bulk_degenes2 %>% filter(CSF.IL34.lfc > 0, CSF.IL34.pv < 0.05) %>% arrange(desc(CSF.IL34.lfc)) %>% pull(ILMN_GENE)

IL34_induced_markers <- intersect(intersect(IL34_induced_genes, IL34_induced_genes2), row.names(HCC_seu))
mCSF_induced_markers <- intersect(intersect(mCSF_induced_genes, mCSF_induced_genes2), row.names(HCC_seu))
macrophage_meta <- HCC_seu@meta.data %>% filter(Global_Cluster == "Mononuclear")
macrophage_meta$IL34_induced_score <- colMeans(HCC_seu@assays$RNA@data[IL34_induced_markers,macrophage_meta$CellName])
macrophage_meta$mCSF_induced_score <- colMeans(HCC_seu@assays$RNA@data[mCSF_induced_markers,macrophage_meta$CellName])
macrophage_meta[,c("Sub_Cluster","IL34_induced_score","mCSF_induced_score")] %>% filter(Sub_Cluster %in% c("FOLR2+ TAM1_1", "FOLR2+ TAM1_7", "SPP1+ TAM2")) %>% melt(.,id.vars = "Sub_Cluster") %>% 
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot(aes(color = Sub_Cluster, fill = Sub_Cluster), alpha = .4, outlier.size = -1) +
  facet_wrap(~Sub_Cluster, nrow = 1) +
  scale_color_manual(values = TAM_color_panel) +
  scale_fill_manual(values = TAM_color_panel) +
  ggsignif::geom_signif(comparisons = list(
    c("IL34_induced_score","mCSF_induced_score"))
  ) + 
  theme_cowplot(font_size = 12) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") -> p
ggsave(p, file = "./Figures/04.IL34_mCSF_induced_score_boxplot.pdf", width = 4, height = 6)
macrophage_meta %>% filter(Sub_Cluster %in% c("FOLR2+ TAM1_1","FOLR2+ TAM1_7","SPP1+ TAM2")) %>% ggplot(., aes(x = IL34_induced_score, y = mCSF_induced_score)) +
  geom_point(aes(color = Sub_Cluster), size = 3, alpha = .6) +
  scale_color_manual(name = "", values = TAM_color_panel) +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1)) -> p
p.pdf <- ggAIplot(p) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0.2, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -0.2, linetype = "dashed")
ggsave(p.pdf, file = "./Figures/04.IL34_mCSF_induced_score_scatter.pdf", width = 6, height = 4.5)
