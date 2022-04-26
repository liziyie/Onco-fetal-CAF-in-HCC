# Cirrhosis and healthy liver_Ramachandran_Nature----
setwd("/work/lzy/project/onco_fetal/00.data/cirrhotic_liver/")
library(Seurat)
library(dplyr)
# Healthy----
CAF_healthy1A_seu <- CreateSeuratObject(counts = Read10X("./healthy/healthy1-A"), project = "healthy1-A")
CAF_healthy1B_seu <- CreateSeuratObject(counts = Read10X("./healthy/healthy1-B"), project = "healthy1-B")
CAF_healthy2_seu <- CreateSeuratObject(counts = Read10X("./healthy/healthy2"), project = "healthy2")
CAF_healthy3A_seu <- CreateSeuratObject(counts = Read10X("./healthy/healthy3-A"), project = "healthy3-A")
CAF_healthy3B_seu <- CreateSeuratObject(counts = Read10X("./healthy/healthy3-B"), project = "healthy3-B")
CAF_healthy4_seu <- CreateSeuratObject(counts = Read10X("./healthy/healthy4"), project = "healthy4")
CAF_healthy_NatureRam <- merge(CAF_healthy1A_seu, c(CAF_healthy1B_seu, CAF_healthy2_seu, CAF_healthy3A_seu, CAF_healthy3B_seu, CAF_healthy4_seu), add.cells.id = c("H1A","H1B","H2","H3A","H3B","H4"), project = "NatureRam")
rm(CAF_healthy1A_seu, CAF_healthy1B_seu, CAF_healthy2_seu, CAF_healthy3A_seu, CAF_healthy3B_seu, CAF_healthy4_seu);gc()
CAF_healthy_NatureRam <- CAF_healthy_NatureRam %>% NormalizeData() %>% ScaleData() %>% FindVariableFeatures(., nfeatures = 2000) %>% RunPCA(., npcs = 30) %>% FindNeighbors(., dims = 1:15) %>% FindClusters(., resolution = 0.5)
CAF_healthy_NatureRam <- RunUMAP(CAF_healthy_NatureRam, dims = 1:15)
DimPlot(CAF_healthy_NatureRam, reduction = "umap", label = T)
DotPlot(CAF_healthy_NatureRam, features = c("PTPRC","EPCAM","CDH5","FLT1","MGP","LIFR","CD3D","IL7R","CD4","CD8A","NKG7","CD79A","CD14","FCGR3A","CD86","TPSAB1","BATF3","FCER1A","LILRA4","CCR7","ACTA2","IGHG1","S100A8","HBB","CD163","MAFB","CYP3A4","CYP3A7","MKI67","FCRL5","JSRP1","KRT7","KRT19","APOC2","HSD17B6","ALB","ORM1","COL1A1","FAP","SPARC","MYL9","PDGFRB","VIM","DCN","THY1","PDGFRA")) + RotatedAxis()
CAF_healthy_NatureRam$MajorCluster <- plyr::revalue(
  CAF_healthy_NatureRam$RNA_snn_res.0.5,
  replace = c("0" = "CD4 T", "1" = "Mono/Macro", "2" = "Hepatocyte", "3" = "CD8 T", "4" = "Endothelial", 
              "5" = "NK", "6" = "Endothelial", "7" = "Fibroblast", "8" = "NK", "9" = "Endothelial", 
              "10" = "Fibroblast", "11" = "Mono/Macro", "12" = "Hepatocyte", "13" = "DC", "14" = "Plasma", 
              "15" = "Mono/Macro", "16" = "Endothelial", "17" = "Cycling", "18" = "B cell", "19" = "Hepatocyte"))
DimPlot(CAF_healthy_NatureRam, reduction = "umap", group.by = "MajorCluster", label = T) + 
  theme(legend.position = "none")
CAF_healthy_Ramachandran_Nature_2019 <- subset(CAF_healthy_NatureRam, subset = RNA_snn_res.0.5 %in% c(7,10))
CAF_healthy_Ramachandran_Nature_2019 <- CAF_healthy_Ramachandran_Nature_2019 %>% FindVariableFeatures(nfeatures = 1000) %>% RunPCA(npcs = 30) %>% FindNeighbors(dims = 1:15) %>% FindClusters(resolution = 0.3)
CAF_healthy_Ramachandran_Nature_2019 <- RunUMAP(CAF_healthy_Ramachandran_Nature_2019, dims = 1:10)
DimPlot(CAF_healthy_Ramachandran_Nature_2019, label = T)
DotPlot(CAF_healthy_Ramachandran_Nature_2019, features = c("PTPRC","EPCAM","FLT1","ESAM","COL1A1","FAP","SPARC","MYL9","ACTA2","PDGFRB","VIM","DCN","ECM1","TAGLN","RGS5","TIMP1","PDGFRA")) + RotatedAxis()
CAF_healthy_Ramachandran_Nature_2019$Sub_Cluster <- plyr::revalue(
  CAF_healthy_Ramachandran_Nature_2019$seurat_clusters,
  replace = c("0" = "VSMC",
              "1" = "HSC",
              "2" = "HSC",
              "3" = "SAMes")
)
save(CAF_healthy_Ramachandran_Nature_2019, file = "./healthy/CAF_healthy_Ramachandran_Nature_2019.rda")

# Cirrhosis----
CAF_cirrhosis1A_seu <- CreateSeuratObject(counts = Read10X("./cirrhotic/cirrhotic1-A/"), project = "cirrhotic1-A")
CAF_cirrhosis1B_seu <- CreateSeuratObject(counts = Read10X("./cirrhotic/cirrhotic1-B/"), project = "cirrhotic1-B")
CAF_cirrhosis2_seu <- CreateSeuratObject(counts = Read10X("./cirrhotic/cirrhotic2/"), project = "cirrhotic2")
CAF_cirrhosis3_seu <- CreateSeuratObject(counts = Read10X("./cirrhotic/cirrhotic3/"), project = "cirrhotic3")
CAF_cirrhosis_NatureRam <- merge(CAF_cirrhosis1A_seu, c(CAF_cirrhosis1B_seu, CAF_cirrhosis2_seu, CAF_cirrhosis3_seu), add.cells.id = c("C1A","C1B","C2","C3"), project = "NatureRam")
rm(CAF_cirrhosis1A_seu, CAF_cirrhosis1B_seu, CAF_cirrhosis2_seu, CAF_cirrhosis3_seu);gc()
CAF_cirrhosis_NatureRam <- CAF_cirrhosis_NatureRam %>% NormalizeData() %>% ScaleData() %>% FindVariableFeatures(., nfeatures = 2000) %>% RunPCA(., npcs = 30) %>% FindNeighbors(., dims = 1:15) %>% FindClusters(., resolution = 0.5)
CAF_cirrhosis_NatureRam <- RunUMAP(CAF_cirrhosis_NatureRam, dims = 1:15)
DimPlot(CAF_cirrhosis_NatureRam, reduction = "umap", label = T)
DotPlot(CAF_cirrhosis_NatureRam, features = c("PTPRC","EPCAM","CDH5","FLT1","MGP","LIFR","CD3D","IL7R","CD4","CD8A","NKG7","CD79A","CD14","FCGR3A","CD86","TPSAB1","BATF3","FCER1A","LILRA4","CCR7","ACTA2","IGHG1","S100A8","HBB","CD163","MAFB","CYP3A4","CYP3A7","MKI67","FCRL5","JSRP1","KRT7","KRT19","APOC2","HSD17B6","ALB","ORM1","COL1A1","FAP","SPARC","MYL9","PDGFRB","VIM","DCN","THY1","PDGFRA")) + RotatedAxis()
CAF_cirrhosis_NatureRam$MajorCluster <- plyr::revalue(
  CAF_cirrhosis_NatureRam$RNA_snn_res.0.5,
  replace = c("0" = "Endothelial", "1" = "Endothelial", "2" = "NK", "3" = "T", "4" = "Epithelial", 
              "5" = "Mono/Macro", "6" = "Epithelial", "7" = "Mono/Macro", "8" = "Mono/Macro", "9" = "Endothelial", 
              "10" = "Epithelial", "11" = "Plasma", "12" = "Fibroblast", "13" = "Cycling", "14" = "Endothelial", 
              "15" = "Plasma", "16" = "Fibroblast", "17" = "Endothelial"))
DimPlot(CAF_cirrhosis_NatureRam, reduction = "umap", group.by = "MajorCluster", label = T) + 
  theme(legend.position = "none")
CAF_cirrhosis_Ramachandran_Nature_2019 <- subset(CAF_cirrhosis_NatureRam, subset = RNA_snn_res.0.5 %in% c(12,16))
CAF_cirrhosis_Ramachandran_Nature_2019 <- CAF_cirrhosis_Ramachandran_Nature_2019 %>% FindVariableFeatures(nfeatures = 1000) %>% RunPCA(npcs = 30) %>% FindNeighbors(dims = 1:15) %>% FindClusters(resolution = 0.6)
CAF_cirrhosis_Ramachandran_Nature_2019 <- RunUMAP(CAF_cirrhosis_Ramachandran_Nature_2019, dims = 1:10)
DimPlot(CAF_cirrhosis_Ramachandran_Nature_2019, label = T)
DotPlot(CAF_cirrhosis_Ramachandran_Nature_2019, features = c("PTPRC","EPCAM","FLT1","ESAM","COL1A1","FAP","SPARC","MYL9","ACTA2","PDGFRB","VIM","DCN","ECM1","TAGLN","RGS5","TIMP1","PDGFRA")) + RotatedAxis()
CAF_cirrhosis_Ramachandran_Nature_2019$Sub_Cluster <- plyr::revalue(
  CAF_cirrhosis_Ramachandran_Nature_2019$seurat_clusters,
  replace = c("0" = "HSC",
              "1" = "VSMC",
              "2" = "SAMes",
              "3" = "Mesothelia")
)
save(CAF_cirrhosis_Ramachandran_Nature_2019, file = "./cirrhotic/CAF_cirrhosis_Ramachandran_Nature_2019.rda")

# Integrate CAF in HCC, healthy liver and cirrhosis----
setwd("/work/lzy/project/onco_fetal/")
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
CAF_HCC <- readRDS("./02.processed_data/CAF_HCC.rds")
CAF_HCC$Sub_Cluster <- as.character(CAF_HCC$Sub_Cluster)
CAF_HCC$Sub_Cluster <- factor(CAF_HCC$Sub_Cluster, levels = c("CAF","HSP+ CAF","MYH11+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF","Fib","MT1M+ Fib","ABCAB+ Fib"))
CAF_HCC$Source <- CAF_HCC$NTF
load("./00.data/cirrhotic_liver/healthy/CAF_healthy_Ramachandran_Nature_2019.rda")
CAF_healthy_Ramachandran_Nature_2019$Source <- "Healthy"
load("./00.data/cirrhotic_liver/cirrhotic/CAF_cirrhosis_Ramachandran_Nature_2019.rda")
CAF_cirrhosis_Ramachandran_Nature_2019$Source <- "Cirrhosis"
CAF.list <- list(CAF_HCC, CAF_healthy_Ramachandran_Nature_2019, CAF_cirrhosis_Ramachandran_Nature_2019)
for (i in 1:length(CAF.list)) {
  CAF.list[[i]] <- NormalizeData(CAF.list[[i]], verbose = FALSE)
  CAF.list[[i]] <- FindVariableFeatures(CAF.list[[i]], nfeatures = 1500, verbose = FALSE)
}
CAF.anchors <- FindIntegrationAnchors(object.list = CAF.list, dims = 1:30, k.filter = 50)
CAF.integrated <- IntegrateData(anchorset = CAF.anchors, dims = 1:30)
DefaultAssay(CAF.integrated) <- "integrated"
CAF.integrated <- ScaleData(CAF.integrated, verbose = FALSE)
CAF.integrated <- RunPCA(CAF.integrated, npcs = 30, verbose = FALSE)
CAF.integrated <- FindNeighbors(CAF.integrated, reduction = "pca", dims = 1:15)
CAF.integrated <- FindClusters(CAF.integrated, resolution = .6)
CAF.integrated <- RunUMAP(CAF.integrated, reduction = "pca", dims = 1:15, verbose = FALSE, min.dist = .1, spread = 1)
CAF.integrated$Integrated_cluster <- CAF.integrated$integrated_snn_res.0.6

CAF.integrated$RNA_snn_res.0.3 <- CAF.integrated$RNA_snn_res.0.5 <- CAF.integrated$RNA_snn_res.0.6 <- CAF.integrated$integrated_snn_res.0.6 <- CAF.integrated$seurat_clusters <- c()
CAF.integrated$integrated_UMAP_1 <- CAF.integrated@reductions$umap@cell.embeddings[,1]
CAF.integrated$integrated_UMAP_2 <- CAF.integrated@reductions$umap@cell.embeddings[,2]
saveRDS(CAF.integrated, file = "../02.processed_data/CAF_integration_CAF.rds")

# Integrate CAF in HCC, fetal and healthy liver----
setwd("/work/lzy/project/onco_fetal/")
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
CAF_HCC <- readRDS("./02.processed_data/CAF_HCC.rds")
CAF_HCC$Sub_Cluster <- as.character(CAF_HCC$Sub_Cluster)
CAF_HCC$Sub_Cluster[CAF_HCC$Sub_Cluster == "Fib"] <- "Fib(Tumor)"
CAF_HCC$Sub_Cluster <- factor(CAF_HCC$Sub_Cluster, levels = c( "CAF","HSP+ CAF","MYH11+ CAF","SDC2+ CAF","POSTN+ CAF","APOA2+ CAF","Fib(Tumor)","MT1M+ Fib","ABCAB+ Fib"))
CAF_HCC$Source <- CAF_HCC$NTF
Fib_fetal <- readRDS("./02.processed_data/Fib_fetal.rds")
Fib_fetal$Sub_Cluster <- as.character(Fib_fetal$Sub_Cluster)
Fib_fetal$Sub_Cluster[Fib_fetal$Sub_Cluster == "Fib"] <- "Fib(Fetal)"
Fib_fetal$Sub_Cluster <- factor(Fib_fetal$Sub_Cluster, levels = c("Fib(Fetal)","ACTA2+ Fib","COL1A1+ Fib","BGN+ Fib","APOC3+ Fib","ITM2C+ Fib","Mesenchymal","PLEK+ Fib"))
Fib_fetal$Source <- "Fetal"
load("./00.data/cirrhotic_liver/healthy/CAF_healthy_Ramachandran_Nature_2019.rda")
CAF_healthy_Ramachandran_Nature_2019$Source <- "Healthy"
CAF.list <- list(CAF_HCC, Fib_fetal %>% subset(Sub_Cluster != "Mesenchymal"), CAF_healthy_Ramachandran_Nature_2019)
for (i in 1:length(CAF.list)) {
  CAF.list[[i]] <- NormalizeData(CAF.list[[i]], verbose = FALSE)
  CAF.list[[i]] <- FindVariableFeatures(CAF.list[[i]], nfeatures = 1500, verbose = FALSE)
}
CAF.anchors <- FindIntegrationAnchors(object.list = CAF.list, dims = 1:30, k.filter = 50)
CAF.integrated <- IntegrateData(anchorset = CAF.anchors, dims = 1:30)
DefaultAssay(CAF.integrated) <- "integrated"
CAF.integrated <- ScaleData(CAF.integrated, verbose = FALSE)
CAF.integrated <- RunPCA(CAF.integrated, npcs = 30, verbose = FALSE)
CAF.integrated <- FindNeighbors(CAF.integrated, reduction = "pca", dims = 1:15)
CAF.integrated <- FindClusters(CAF.integrated, resolution = .6)
CAF.integrated <- RunUMAP(CAF.integrated, reduction = "pca", dims = 1:15, verbose = FALSE, min.dist = .1, spread = 1)
CAF.integrated$Integrated_cluster <- CAF.integrated$integrated_snn_res.0.6
CAF.integrated$RNA_snn_res.0.3 <- CAF.integrated$RNA_snn_res.0.5 <- CAF.integrated$integrated_snn_res.0.6 <- CAF.integrated$seurat_clusters <- c()
CAF.integrated$Integrated_cluster[CAF.integrated$Integrated_cluster == 12] <- 4
CAF.integrated$integrated_UMAP_1 <- CAF.integrated@reductions$umap@cell.embeddings[,1]
CAF.integrated$integrated_UMAP_2 <- CAF.integrated@reductions$umap@cell.embeddings[,2]
saveRDS(CAF.integrated, file = "./02.processed_data/CAF_integration_onco-fetal.rds")

# Music deconvolution----
library(MuSiC)
library(Biobase)
setwd("/work/lzy/project/onco_fetal/")
Ankur_exp <- read.table("./00.data/HCC_bulk/Ankur_HCC_bulk.txt", sep = "\t", head = T)
Ankur_exp <- Ankur_exp[!duplicated(Ankur_exp$Gene),]
row.names(Ankur_exp) <- Ankur_exp$Gene
Ankur_exp$Gene <- c()
Ankur_sampleinfo <- read.csv("./00.data/HCC_bulk/sample_info.csv", row.names = 1, stringsAsFactors = F)
Ankur_clinical <- readr::read_tsv("./00.data/HCC_bulk/ClinicalData.tsv") %>% data.frame()
Ankur_metadata <- merge(Ankur_sampleinfo,Ankur_clinical,by.x="sample_id",by.y="Unified_ID")
row.names(Ankur_metadata) <- Ankur_metadata$Bulk_sample_id
Ankur_metadata <- Ankur_metadata[colnames(Ankur_exp),]
Ankur_exp_temp <- Ankur_exp[,Ankur_metadata %>% filter(Tissue == "T") %>% pull(Bulk_sample_id)]
bulk.eset <- ExpressionSet(as.matrix(Ankur_exp), phenoData = AnnotatedDataFrame(Ankur_metadata))

sc_exp <- read.table("./00.data/HCC_bulk/CIBERSORTx/Input/Tumor_training2.txt", sep = "\t", head = T, stringsAsFactors = F, row.names = 1, check.names = F)
sc_metadata <- data.frame(Cluster = colnames(sc_exp), row.names = paste0("Cell",1:ncol(sc_exp)))
sc_metadata$sampleID <- row.names(sc_metadata)
colnames(sc_exp) <- paste0("Cell",1:ncol(sc_exp))
sc.eset <- ExpressionSet(as.matrix(sc_exp), phenoData = AnnotatedDataFrame(sc_metadata))
Ankur_Music_prop <- music_prop(
  bulk.eset = bulk.eset, sc.eset = sc.eset, 
  clusters = 'Cluster', samples = 'sampleID', 
  verbose = F)
saveRDS(Ankur_Music_prop, file = "./00.data/HCC_bulk/CIBERSORTx/Output/Music_deconvolution.rds")
