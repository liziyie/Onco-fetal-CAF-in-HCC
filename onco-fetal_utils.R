# Color panel----
Binomial_color_panel <- c("TRUE" = "#E64B35", "FALSE" = "lightgrey")

# >>Not include in----
`%ni%` <- Negate(`%in%`)

# >>LIMMA----
LIMMA <- function(expression_matrix, groupid, doAUC = FALSE) {
  ## Differential expressed genes in two groups by Limma.
  ##
  ## Args:
  #' @expression_matrix: Gene*cell matrix.
  #' @groupid: The groupid of each cell, there should be only two groups.
  ##
  ## Returns:
  ## A dataframe with the output of limma and expression percentage.
  library(limma)
  library(ROCR)
  expression_matrix <- as.matrix(expression_matrix)
  groupid <- as.character(groupid)
  groupid_raw <- unique(groupid)
  groupid[groupid == groupid_raw[1]] <- "GroupA"
  groupid[groupid == groupid_raw[2]] <- "GroupB"
  contrast <<- paste0(levels(factor(groupid)), collapse = "-")
  design <- model.matrix( ~ 0 + factor(groupid))
  colnames(design) <- levels(factor(groupid))
  rownames(design) <-
    colnames(expression_matrix)  # design data used in limma
  contrast.matrix <- makeContrasts(contrast, levels = design)
  fit <- lmFit(expression_matrix, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput <-  topTable(fit2, coef = 1, n = Inf)
  nrDEG <-  na.omit(tempOutput)
  nrDEG$Symbol <- row.names(nrDEG)
  positive_group <-
    row.names(fit2$contrasts)[fit2$contrasts == 1]  # high expression when logFC > 0
  negative_group <-
    row.names(fit2$contrasts)[fit2$contrasts == -1]  # low expression when logFC < 0
  nrDEG$Grp <-
    c(negative_group, positive_group)[as.numeric(nrDEG$logFC > 0) + 1]
  nrDEG$Grp[nrDEG$Grp == "GroupA"] <- groupid_raw[1]
  nrDEG$Grp[nrDEG$Grp == "GroupB"] <- groupid_raw[2]
  cell.Grp1 <- which(groupid == levels(as.factor(groupid))[1])
  cell.Grp2 <- which(groupid == levels(as.factor(groupid))[2])
  Exp.Mean.Grp1 <-
    rowMeans(expression_matrix[nrDEG$Symbol, cell.Grp1])  # mean expression in the group
  Exp.Mean.Grp2 <-
    rowMeans(expression_matrix[nrDEG$Symbol, cell.Grp2])  # mean expression out the group
  Exp.Per.Grp1 <-
    apply(expression_matrix[nrDEG$Symbol, cell.Grp1], 1, Expression.Per)  # expression percentage in the group
  Exp.Per.Grp2 <-
    apply(expression_matrix[nrDEG$Symbol, cell.Grp2], 1, Expression.Per)  # expression percentage out the group
  de.genes.all <- cbind(nrDEG, Exp.Mean.Grp1, Exp.Per.Grp1, Exp.Mean.Grp2, Exp.Per.Grp2)
  colnames(de.genes.all)[9:12] <-
    c(paste0(c("Exp.Mean.", "Exp.Per."), groupid_raw[1]),
      paste0(c("Exp.Mean.", "Exp.Per."), groupid_raw[2]))
  if (!is.na(de.genes.all[1, 1])) {
    if(doAUC){
      for (k in 1:nrow(de.genes.all)) {
        category <- as.numeric(groupid_new == de.genes.all$Grp[k])
        pred <-
          prediction(expression_matrix[de.genes.all$Symbol[k], ], category)
        pauc <- performance(pred, measure = "auc")
        de.genes.all[k, "AUC"] <- pauc@y.values[[1]]
      }  # calculate the AUC for each gene
    }
    de.genes.all <- de.genes.all %>% arrange(desc(logFC))
    de.genes.all
    return(de.genes.all)
  } else{
    print("No significant genes!")
  }
  
}

Expression.Per <- function(x, cutoff = 0.1) {
  # percent of gene-expressed cell, cutoff = 0.1
  return(sum(x > cutoff) / length(x))
}

# >>BatchEntropy----
BatchEntropy <- function(input.data, group.id, k.used = 30, dimension.used = "tSNE") {
  ## Calculate the cell entropy according to the cell group.
  ##
  ## Args:
  #' @input.data: A matrix with each cell in a row. The column could be
  #' tSNE coordinate, expression matrix or even a dist format.
  #' @group.id: A vector with the same length as the row of input.data, or
  #' a list of vectors.
  #' @k.used: The k used to build a kNN-graph.
  #' @dimension.used: The method to reduce the dimension, tSNE by default,
  #' could also be PCA or raw.
  ##
  ## Returns:
  ## A vector with each cell's entropy.
  library(dbscan)
  library(entropy)
  if (dimension.used == "raw") {
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = input.data, k = k.used, sort = FALSE, search = "dist")
  }
  if (dimension.used == "tSNE") {
    cat(paste("Reduce the dimension using", dimension.used, "at", Sys.time(), "\n"))
    tSNE.coor <- Rtsne::Rtsne(input.data)
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = tSNE.coor$Y, k = k.used, sort = FALSE, search = "dist")
  }
  if (dimension.used == "PCA") {
    cat(paste("Reduce the dimension using", dimension.used, "at", Sys.time(), "\n"))
    PCA.coor <- prcomp(input.data, center = FALSE)
    PCA.cumsd <- cumsum(PCA.coor$sdev) / sum(PCA.coor$sdev)
    nPCs.used <- which(PCA.cumsd > 0.9)[1]
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = PCA.coor$x[, 1:nPCs.used], k = k.used, sort = FALSE, search = "dist")
  }
  if (!is.list(group.id)) {
    group.id <- list(default = group.id)
  }
  cell_entropy <- list()
  for (i in names(group.id)) {
    knn_group <- matrix(group.id[[i]][knn_graph$id],
                        nrow = nrow(input.data),
                        byrow = FALSE)
    row.names(knn_group) <- row.names(input.data)
    colnames(knn_group) <- 1:k.used
    cat(paste("Calculate the cell entropy of", i, "at", Sys.time(), "\n"))
    cell_entropy[[i]] <- apply(knn_group, 1, function(x) {
      entropy(table(x))
    })
  }
  return(cell_entropy)
}

# >>Confusion heatmap----
Confusion_heatmap_new <- function (ori, prd, color = NULL) 
{
  cross.validation.filt <- tibble(ori = ori, prd = prd) %>% 
    dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
  cross.validation.filt[is.na(cross.validation.filt)] = 0
  cross.validation.filt[, -1] <- round(cross.validation.filt[, 
                                                             -1]/rowSums(cross.validation.filt[, -1]), 2)
  cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", 
                                                                   value = "Prob", -ori)
  p <- cross.validation.filt %>% ggplot(aes(ori, prd, fill = Prob)) + 
    geom_tile() + theme(axis.title = element_text(size = 0)) + 
    theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
    theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                         panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                         axis.ticks = element_blank(), axis.title = element_blank()) + 
    theme(axis.text.y = element_text(color = "black"), 
          axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1))
  if(is.null(color)){
    p <- p + scale_fill_viridis()
  }else{
    p <- p + scale_fill_gradientn(colors = color)
  }
  return(p)
}

# Chi-square test plot---
CrossTabPlot <- function(crosstab){
  pvalue <- fisher.test(crosstab)
  ct <- crosstab %>% data.frame()
  ct$Prop <- round(ct$Freq/sum(ct$Freq) * 100,2)
  item1 <- unique(ct[,1])[1]
  item2 <- unique(ct[,2])[1]
  ct[(ct[,1] == item1) == (ct[,2] == item2),"Group"] <- "Group1"
  ct[(ct[,1] == item1) != (ct[,2] == item2),"Group"] <- "Group2"
  p <- ggplot(ct, aes_string(x = colnames(ct)[1], y = colnames(ct)[2])) +
    geom_tile(aes(fill = Group)) +
    ggtitle(paste0("p-value = ",format(pvalue$p.value,digits=3))) +
    geom_text(aes(label = paste0(Freq,"\n(",Prop,"%)"))) +
    scale_fill_manual(values = c("Group1" = "#9bc0cc", "Group2" = "#e1e2e3")) +
    theme_minimal() +
    theme(legend.position = "none")
  return(p)
}

# >>RO/E----
ROIE <- function(crosstab, filter = NULL){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  if(is.null(filter)){filter = 10}
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  rowsum.matrix[rowsum.matrix <= filter] <- 0
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  colsum.matrix[colsum.matrix <= filter] <- 0
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  roie <- roie[rowSums(roie)!=0,]
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        if(m2[i,j] == 0){
          div.result[i,j] <- 0
        }else{
          div.result[i,j] <- m1[i,j] / m2[i,j]
        }
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

# >>ROIE_plot----
ROIE_plot <- function(ROIE_table, max = 2.5){
  ROIE_table <- melt(ROIE_table)
  ROIE_table$text <- round(ROIE_table$value,2)
  ROIE_table$value[ROIE_table$value > max] <- max
  ROIE_table$value[ROIE_table$value == 0] <- NA
  ROIE_table$text[ROIE_table$text == 0] <- "-"
  p <- 
    ggplot(data = ROIE_table, aes(Var2, Var1, fill = value)) +
    geom_tile() + 
    scale_fill_gradientn(name = "Ro/e", colours = colorRampPalette(brewer.pal(4, "YlOrRd"))(100), na.value = "lightgrey") +
    geom_text(aes(label = text)) +
    labs(x = "", y = "") +
    theme_cowplot(font_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          axis.line = element_blank(), 
          axis.ticks = element_blank()) +
    guides(fill = guide_colourbar(barwidth = 1, barheight = 10))
  return(p)
}

# >>Run Nichenetr----
RunNichenetr <- function(seu, cluster, group.by = "Sub_Cluster", receiver_genes, sender_genes = NULL, expressed_genes_receiver = NULL, background_expressed_genes = NULL, expressed_gene_cutoff = 2, best_upstream_ligands_number = 30, species = "homo", curated = FALSE, repo = TRUE, only_genes_mode = FALSE){
  library(nichenetr)
  if(species == "homo"){
    if(exists("ligand_target_matrix_human") & exists("weighted_networks_human") & exists("lr_network_human")){
      if(repo) cat("Required database already exists.\n")
    }else{
      if(repo) cat("Load required database.\n")
      load("/Volumes/ZiyiLi/Utils/nichenetr_human.rda")
    }
    nichenetr_data <- list(
      ligand_target_matrix = ligand_target_matrix_human,
      lr_network = lr_network_human,
      weighted_networks = weighted_networks_human
    )
  }else if(species == "mouse"){
    if(exists("ligand_target_matrix_mouse") & exists("weighted_networks_mouse") & exists("lr_network_mouse")){
      if(repo) cat("Required database already exists.\n")
    }else{
      if(repo) cat("Load required database.\n")
      load("/Volumes/ZiyiLi/Utils/nichenetr_mouse.rda")
    }
    nichenetr_data <- list(
      ligand_target_matrix = ligand_target_matrix_mouse,
      lr_network = lr_network_mouse,
      weighted_networks = weighted_networks_mouse
    )
  }
  if(curated){
    nichenetr_data$lr_network <- nichenetr_data$lr_network %>% filter(database %ni% c("ppi_prediction_go","ppi_prediction"))
  }
  if(repo) cat("Extract expressed genes.\n")
  geneset_oi = receiver_genes
  if(!only_genes_mode){
    cells_used <- row.names(seu@meta.data)[seu@meta.data[,group.by] == cluster]
    expression = t(seu@assays$RNA@data[,cells_used])
    if(is.null(expressed_genes_receiver)){
      expressed_genes_receiver = expression %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= expressed_gene_cutoff] %>% names()
    }
    background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(nichenetr_data[["ligand_target_matrix"]])]
    if(is.null(sender_genes)){
      # expressed_genes_sender = colnames(expression)
      expressed_genes_sender = expression %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= expressed_gene_cutoff] %>% names()
    }else{
      expressed_genes_sender = sender_genes
    }
  }else{
    if(is.null(expressed_genes_receiver)){
      expressed_genes_receiver = background_expressed_genes 
    }
    background_expressed_genes = background_expressed_genes %>% .[. %in% rownames(nichenetr_data[["ligand_target_matrix"]])]
    expressed_genes_sender = sender_genes
  }
  if(repo) cat("Calculate ligand activities.\n")
  ligands = nichenetr_data[["lr_network"]] %>% pull(from) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  receptors = nichenetr_data[["lr_network"]] %>% pull(to) %>% unique()
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  lr_network_expressed = nichenetr_data[["lr_network"]] %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
  potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
  ligand_activities = predict_ligand_activities(
    geneset = geneset_oi, 
    background_expressed_genes = background_expressed_genes, 
    ligand_target_matrix = nichenetr_data[["ligand_target_matrix"]],
    potential_ligands = potential_ligands)
  best_upstream_ligands = ligand_activities %>% top_n(best_upstream_ligands_number, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
  
  if(repo) cat("Infer top-predicted target genes.\n")
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = nichenetr_data[["ligand_target_matrix"]], n = 250) %>% bind_rows()
  active_ligand_target_links_df <- active_ligand_target_links_df[!is.na(active_ligand_target_links_df$weight),]
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = nichenetr_data[["ligand_target_matrix"]], cutoff = 0.25)
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets = active_ligand_target_links_df$target %>% unique()
  order_targets <- order_targets[order_targets %in% row.names(active_ligand_target_links)]
  
  if(repo) cat("Infer top-predicted receptors.\n")
  lr_network_top = nichenetr_data[["lr_network"]] %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  lr_network_top_df = nichenetr_data[["weighted_networks"]]$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  lr_network_top_df = lr_network_top_df %>% tidyr::spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  if(repo) cat("Plot figures.\n")
  ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = 
    vis_ligand_pearson %>% 
    make_heatmap_ggplot(
      "Prioritized ligands","", 
      color = "darkorange", x_axis_position = "top", 
      legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + 
    theme(legend.text = element_text(size = 9),
          axis.ticks = element_blank(),
          axis.title.x = element_text(),
          axis.text.y = element_text(face = "italic"))
  
  if(only_genes_mode) {cluster = NULL}
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = 
    vis_ligand_target %>% 
    make_heatmap_ggplot(
      "Potential ligands",paste0(cluster," upregulated genes"), 
      color = "purple",legend_position = "top", 
      x_axis_position = "top",
      legend_title = "Regulatory potential") + 
    scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + 
    theme(axis.text.x = element_text(face = "italic"),
          axis.ticks = element_blank(),
          axis.title.y = element_blank())
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands]
  p_ligand_receptor_network = 
    vis_ligand_receptor_network %>% t() %>% 
    make_heatmap_ggplot(
      "Potential ligands",paste0("Receptors expressed by ",cluster), 
      color = "mediumvioletred", 
      x_axis_position = "top",legend_title = "Prior interaction potential") + 
    theme(axis.ticks = element_blank(), 
          axis.title.y = element_blank(),
          axis.text = element_text(face = "italic"))
  
  figures_without_legend = cowplot::plot_grid(
    p_ligand_pearson + theme(legend.position = "none"),
    NULL,
    p_ligand_target_network + theme(legend.position = "none", axis.text.y = element_blank()),
    NULL, NULL, 
    p_ligand_receptor_network + theme(legend.position = "none"),
    align = "hv", nrow = 2,
    rel_widths = c(ncol(vis_ligand_pearson) + 6, -5, ncol(vis_ligand_target) - 6)
  )
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
    NULL,
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    NULL,
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor_network)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, -.3, 1, -.2, 1))
  combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,1), nrow = 2, align = "hv")
  
  results <- list(
    ligand_activities = ligand_activities,
    best_upstream_ligands = best_upstream_ligands,
    active_ligand_target_links_df = active_ligand_target_links_df,
    active_ligand_target_links = active_ligand_target_links,
    lr_network_top_matrix = lr_network_top_matrix,
    p_ligand_pearson = p_ligand_pearson,
    p_ligand_target_network = p_ligand_target_network,
    p_ligand_receptor_network = p_ligand_receptor_network,
    combined_plot = combined_plot
  )
  return(results)
}

# >>Define module score----
CalculateModuleScore <- function(exp, features, iter = 1000, nbin = 20, random_add = 1e+30){
  exp <- t(t(exp) - colMeans(exp))
  features_mean_exp <- colMeans(exp[features,])
  gene.avg <- rowMeans(exp)
  gene.avg <- gene.avg[order(gene.avg)]
  gene.cut <- cut_number(x = gene.avg + rnorm(length(gene.avg))/random_add, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(gene.cut) <- names(gene.avg)
  features.cut <- data.frame(table(gene.cut[features]))
  ctrl_mean_exp <- pbapply::pblapply(1:iter, function(k){
    ctrl.use <- c()
    for(i in 1:nrow(features.cut)){
      ctrl.use <- c(ctrl.use, sample(names(gene.cut)[which(gene.cut == features.cut$Var1[i])],
                                     features.cut$Freq[i], replace = F))
    }
    return(colMeans(exp[ctrl.use,]))
  })
  ctrl_mean_exp <- do.call(rbind, ctrl_mean_exp)
  p_value <- colSums(t(t(ctrl_mean_exp) > features_mean_exp)) / iter
  p_value[p_value == 0] <- 0.1/iter
  module_score <- -log10(p_value)
  module_score <- (module_score - min(module_score)) / max(module_score)
  return(module_score)
}

# >>Find neighbor spots----
FindNeighborSpots <- function(seu, cells, datatype = "Stereo-seq"){
  results <- c()
  for(cell in cells){
    results <- c(results, FindNeighborSpot(seu, cell, dist = 1, self = T, datatype = datatype))
  }
  results <- unique(results)
  return(results)
}

FindNeighborSpot <- function(seu, cell, dist = 1, self = T, datatype = "Stereo-seq"){
  if(datatype == "Stereo-seq"){
    imagerow0 <- seu@images$slice1@coordinates[cell,"imagerow"]
    imagecol0 <- seu@images$slice1@coordinates[cell,"imagecol"]
    row_range <- (imagerow0-dist):(imagerow0+dist)
    col_range <- (imagecol0-dist):(imagecol0+dist)
    cells <- seu@images$slice1@coordinates %>% filter(imagerow %in% row_range, imagecol %in% col_range) %>% row.names()
  } else if(datatype == "Visium"){
    row0 <- seu@images$slice1@coordinates[cell,"row"]
    col0 <- seu@images$slice1@coordinates[cell,"col"]
    row_range <- (row0-dist):(row0+dist)
    col_range <- (col0-dist):(col0+dist)
    cells <- seu@images$slice1@coordinates %>% filter(row %in% row_range, col %in% col_range) %>% row.names()
    cells <- c(cells, seu@images$slice1@coordinates %>% filter(row == row0, col %in% c(col0+2, col0-2)) %>% row.names())
  }
  if(self == FALSE){
    cells <- setdiff(cells, cell)
  }
  return(cells)
}

# >>Calculate Neighborhood Index----
CalNeighIndex <- function(seu, cells, name = ""){
  seu@meta.data[,paste0(name,"_NeighIndex")] <- 0
  for(cell in cells){
    neighbor_spots <- FindNeighborSpot(seu, cell, dist = 1, self = T)
    seu@meta.data[neighbor_spots,paste0(name,"_NeighIndex")] <- seu@meta.data[neighbor_spots,paste0(name,"_NeighIndex")] + 1/9
  }
  return(seu)
}

# >>Calculate Distance Expression----
DistFeature <- function(seu, cells, feature, max.dist = 20){
  seu@images$slice1@coordinates <- seu@images$slice1@coordinates %>% mutate(coord = paste0(imagerow,"_",imagecol))
  ncells <- length(cells)
  dist.df <- data.frame(x_coord = 0, y_coord = 0, dist = 0)
  for(i in -max.dist:max.dist){
    for(j in -max.dist:max.dist){
      spot.dist <- sqrt(i^2 + j^2)
      dist.df <- rbind(dist.df, c(i,j,spot.dist))
    }
  }
  dist.df <- dist.df[-1,]
  dist.list <- list()
  for(i in 1:max.dist){
    dist.list[[i]] <- dist.df %>% filter(dist > i-1, dist <= i)
  }
  dist.feature.matrix <- matrix(numeric(ncells*max.dist), nrow = ncells, ncol = max.dist)
  row.names(dist.feature.matrix) <- cells
  colnames(dist.feature.matrix) <- 1:max.dist
  for(i in 1:ncells){
    cell <- cells[i]
    imagerow <- seu@images$slice1@coordinates[cell,"imagerow"]
    imagecol <- seu@images$slice1@coordinates[cell,"imagecol"]
    dist.feature <- c()
    for(j in 1:max.dist){
      spots_dist_coord <- t(c(imagerow,imagecol) + t(dist.list[[j]][,c(1,2)]))
      spots_dist_coord <- paste0(spots_dist_coord[,1],"_",spots_dist_coord[,2])
      cells.distj <- seu@images$slice1@coordinates %>% filter(coord %in% spots_dist_coord) %>% row.names()
      mean_feature <- sum(seu@meta.data[cells.distj,feature] > 0.2) / length(cells.distj)
      dist.feature <- c(dist.feature, mean_feature)
    }
    dist.feature.matrix[i,] <- dist.feature
  }
  results <- list(
    cells_used = cells,
    dist.list = dist.list,
    dist.feature = dist.feature.matrix
  )
  return(results)
}
