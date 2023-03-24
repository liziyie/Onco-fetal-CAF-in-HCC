# Load Packages----
library(Matrix)

# >>Not include in----
`%ni%` <- Negate(`%in%`)

# >>RemoveBatch----
RemoveBatch <- function(raw.counts, batch, min.cells = 20) {
  ## Remove batch effect in counts level.
  ##
  ## Args:
  #' @raw_counts: gene*cell matrix, Matrix(sparse) format.
  #' @batch: Cells batch vector, should be in the same length as the columns
  #' of raw.counts.
  #' @min.cells: The minimum cells in a batch.
  ##
  ## Returns:
  ## A sacled counts matrix.
  batch <- as.character(batch)
  elements.in.batch <- sort(unique(batch))
  for (element in elements.in.batch) {
    n.element <- sum(batch == element)
    if (n.element <= min.cells) {
      batch[batch == element] = sample(batch[!(batch == element)], n.element)
    }
  }  # Replace the batch having cells less than min.cells with other batch
  batch <- factor(batch)
  cell.depth <- Matrix::colSums(raw.counts)
  batch.depth <- tapply(cell.depth, batch, sum, simplify = T)
  batch.depth <-
    setNames(as.integer(batch.depth), names(batch.depth))  # Depth for each cell and each batch
  gene.average <-
    (Matrix::rowSums(raw.counts) + length(levels(batch))) / (sum(cell.depth) + length(levels(batch)))  # Average counts for each gene
  gene_counts.batch <-
    Matrix(RowSumsdgcMatrix(raw.counts, batch), sparse = T)
  gene_percentage.batch <-
    t(log(gene_counts.batch + 1) - log(batch.depth + 1))  # Total counts percentage for each gene in each batch
  batch_factor <-
    exp(gene_percentage.batch - log(gene.average))  # Batch_factor for each gene in each batch
  cells.batch = split(x = colnames(raw.counts), f = batch)  # Split all cells by batch
  scaled_counts <- lapply(seq_along(cells.batch), function(idx) {
    batch.name = names(cells.batch)[idx]
    cells_in_batch = cells.batch[[idx]]
    return(raw.counts[, cells_in_batch, drop = F] / batch_factor[, batch.name])
  })
  scaled_counts = do.call("cbind", scaled_counts)
  scaled_counts = scaled_counts[, colnames(raw.counts)]  # Calculated the scaled counts
  return(scaled_counts)
}

RowSumsdgcMatrix <- function(x, group, reorder = TRUE) {
  ## Calculate row sums for dgcMatrix(sparse) by group.
  ##
  ## Args:
  #' @x: Sparse matrix.
  #' @group: The group to calculate row sums, should be in the same length as
  #' x's columns.
  ##
  ## Returns:
  ## A matrix has the same columns as x, the rows of the matrix is same as the
  ## levels of group.
  result <- c()
  if (reorder == TRUE) {
    order <- sort(unique(group))
  } else{
    order <- unique(group)
  }
  for (i in order) {
    rowsumi <- Matrix::rowSums(x[, group == i])
    result <- rbind(result, rowsumi)
  }
  row.names(result) <- order
  return(result)
}

# >> Counts2TPM----
Counts2TPM <- function(counts, gene.length) {
  ## Transform the counts matrix into tpm matrix. Be very slow when the
  ## counts matrix is big.
  ##
  ## Args:
  #' @counts: gene*cell matrix.
  #' @gene.length: the gene length for each gene in each cell, should be in the same
  #' dimension as counts matrix, calculated from kallisto.
  ##
  ## Returns:
  ## the tpm matrix.
  if(is.vector(gene.length) & length(counts) == length(gene.length)){
    scaled.counts <- counts / gene.length
    tpm <- 10 ^ 6 * scaled.counts / sum(scaled.counts)
    return(tpm)
  }
  if (sum(dim(counts) == dim(gene.length)) == 2) {
    scaled.counts <- counts / gene.length
    tpm <- 10 ^ 6 * t(t(scaled.counts) / colSums(scaled.counts))
    return(tpm)
  }
  else{
    warning("The dimensions of counts and genelength matrix are different")
  }
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

# >>glm.predict----
glm.predict <-
  function(train.data, train.group, downsample = FALSE, sample.cells = 0, genes.used = NA, test.data, test.group, alpha = 0.99, nfolds = 10, seed = 903417) {
    ## Calculate the similarities of the train data and test data.
    ##
    ## Args:
    #' @train.data: A train data matrix with each cell in a column and each gene
    #' in a row.
    #' @train.group: A vector with the same length as the column of train.data.
    #' @downsample: Whether to sample cells in each cluster to the minimum cluster size.
    #' @sample.cells: Sample cells in each group of cells in train data, if 0 do not
    #' sample cells.
    #' @genes.used: Use a subset of genes in both the train and test data.
    #' @test.data: A test data matrix with each cell in a column and each gene
    #' in a row.
    #' @test.group: A vector with the same length as the column of train.data.
    #' @alpha: The elasticnet mixing parameter, with 0≤α≤1, passed to cv.glmnet.
    #' @nfolds: Number of folds, passed to cv.glmnet.
    ##
    ## Returns:
    ## The probability of each cell in the test.data to be predicted as each group.
    library(glmnet)
    library(ComplexHeatmap)
    set.seed(seed)
    glm.fits <- list()
    glm.predict <- list()
    if (length(genes.used) > 1) {
      train.data <- train.data[genes.used,]
      test.data <- test.data[genes.used,]
      if (length(genes.used) <= 50) {
        cat("There were less than 50 features used in the training data!\n")
      }
    }
    if (sample.cells == 0 & downsample) {
      sample.cells <- max(50, min(table(train.group)))
    }
    if (sample.cells > 0) {
      ngroup <- length(unique(train.group))
      if (ncol(train.data) >= sample.cells * ngroup) {
        cells_used <- c()
        for (groupi in sort(unique(train.group))) {
          if (length(which(train.group == groupi)) > sample.cells) {
            cells_used <-
              c(cells_used, sample(which(train.group == groupi), sample.cells))
          } else{
            cells_used <- c(cells_used, which(train.group == groupi))
          }
        }
        train.data <- train.data[, cells_used]
        train.group <- train.group[cells_used]
      }
    }
    for (groupi in sort(unique(train.group))) {
      fac <-  factor(train.group == groupi)
      glm.fits[[groupi]] <-
        cv.glmnet(x = t(train.data), fac, offset = getPopulationOffset(fac), 
                  family = 'binomial', intercept = FALSE, 
                  alpha = alpha, nfolds = nfolds, type.measure = 'class'
        )
      glm.predict[[groupi]] <-
        predict(
          object = glm.fits[[groupi]],
          newx = t(test.data),
          newoffset = rep(0, ncol(test.data)),
          s = 'lambda.min'
        )
    }
    glm.predict.df <- data.frame(do.call(cbind, glm.predict))
    colnames(glm.predict.df) <- sort(unique(train.group))
    glm.predict.df.prob <- (1 + exp(-glm.predict.df)) ** -1
    glm.cluster <-
      colnames(glm.predict.df.prob)[apply(glm.predict.df.prob, 1, which.max)]
    glm.predict.mean <-
      apply(glm.predict.df, 2, function(e)
        sapply(split(e, test.group), mean))
    glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
    heatmap <- Heatmap(
      t(glm.predict.mean.prob),
      name = 'Predicted\nSimilarity',
      column_title = 'test data',
      row_title = 'train data',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )
    return(
      list(
        test.group = test.group,
        logits = glm.predict.df,
        probability = glm.predict.df.prob,
        cluster = glm.cluster,
        heatmap = heatmap
      )
    )
  }

getPopulationOffset = function(y) {
  ## Calculate the offset value used in glm.predict.
  if (!is.factor(y))
    y = factor(y)
  if (length(levels(y)) != 2)
    stop("y must be a two-level factor")
  off = sum(y == levels(y)[2]) / length(y)
  off = log(off / (1 - off))
  return(rep(off, length(y)))
}

# >>Correlative network----
Corrnetwork <- function(bulk.data, scRNA.data, cluster, signatures.list = NULL, bulk.abundance.matrix = NULL, sc.avg = NULL, sc.proportion = NULL, out.prefix = NULL, exp.cutoff = 1, pro.cutoff = 0.2, ngenes.enrich = 13, enrich.cutoff = 1.96, force = FALSE){
  #' @bulk.data: A gene*sample bulk data expression matrix.
  #' @scRNA.data: A gene*cell scRNA-seq data expression matrix.
  #' @cluster: A vector of the clusters of scRNA-seq data.
  #' @signatures.list: A list of the signature genes of each cluster in scRNA-seq data, the names of
  #' signature list should be consistent with the categories of cluster.
  #' @bulk.abundance.matrix: A sample*proportion matrix.
  #' @out.prefix: The path to save the output.
  #' expression of marker genes, could also be "Cibersortx".
  #' @exp.cutoff: The cutoff of mean expression used in removing the cluster self-expression genes.
  #' @pro.cutoff: The proportion of expressed cells used in removing the cluster self-expression genes.
  #' @ngenes.enrich: The number of top correlation genes used to enrich the clusters.
  
  set.seed(903417)
  dir.create(out.prefix, showWarnings = F, recursive = T)
  
  # Filter genes
  cat("Get shared genes between the scRNA-seq data and the bulk data.\n")
  genes.used <- intersect(row.names(scRNA.data), row.names(bulk.data))
  scRNA.data <- scRNA.data[genes.used,]
  bulk.data <- bulk.data[genes.used,]
  
  # Check the clusters' signatures provided or not
  if(!is.null(bulk.abundance.matrix)){
    bulk.samples.used <- intersect(row.names(bulk.abundance.matrix),colnames(bulk.data))
    bulk.data <- bulk.data[,bulk.samples.used]
    bulk.abundance.matrix <- bulk.abundance.matrix[bulk.samples.used,]
    if(sum(!(cluster %in% colnames(bulk.abundance.matrix)))){
      warning("Not all clusters' signatures were provided, and only the clusters with provided signatures would be used!")
      cells.used <- cluster %in% colnames(bulk.abundance.matrix)
      scRNA.data <- scRNA.data[,cells.used]
      cluster <- cluster[cells.used]
    }
    bulk.abundance.matrix <- bulk.abundance.matrix[,unique(cluster)]
    bulk.abundance <- as.list(as.data.frame(bulk.abundance.matrix))
  }else{
    cat("Filter the signature list with provided clusters and genes in scRNA-seq.\n")
    if(sum(!(cluster %in% names(signatures.list)))){
      warning("Not all clusters' signatures were provided, and only the clusters with provided signatures would be used!")
      cells.used <- cluster %in% names(signatures.list)
      scRNA.data <- scRNA.data[,cells.used]
      cluster <- cluster[cells.used]
    }
    signatures.list <- signatures.list[unique(cluster)]
    signature_used <- lapply(signatures.list, function(x){intersect(x, row.names(bulk.data))})
    if(!is.null(signature_used)) write.csv(t(plyr::ldply(signature_used, rbind)), file = sprintf("%s/signature_used.csv", out.prefix))
    
    # Calculate the cluster abundance in bulk data
    cat("Calculate the abundance of each cluster in bulk data.\n")
    bulk.abundance <- list()
    for(group in names(signatures.list)){
      bulk.abundance[[group]] <- rowMeans(apply(bulk.data[signature_used[[group]],], 1, zscore))
    }
    bulk.abundance.matrix <- do.call(cbind.data.frame, bulk.abundance)
    row.names(bulk.abundance.matrix) <- names(bulk.abundance[[1]])
  }
  write.csv(bulk.abundance.matrix, file = sprintf("%s/bulk_abundance.csv", out.prefix))
  
  # Calculate the average gene expression across clusters in scRNA-seq
  cat("Calculate each gene's mean expression across clusters in scRNA-seq.\n")
  if(!force & file.exists(sprintf("%s/scRNA-seq_mean_expression.csv", out.prefix))){
    cat(sprintf("Read existed file in %s/scRNA-seq_mean_expression.csv\n", out.prefix))
    scRNA_cluster_avg <- read.csv(sprintf("%s/scRNA-seq_mean_expression.csv", out.prefix), check.names = F, row.names = 1)
  }else if(!is.null(sc.avg)){
    cat("Use provided scRNA_cluster_avg.\n")
    scRNA_cluster_avg <- sc.avg
  }else{
    scRNA_cluster_avg <- aggregate(t(as.matrix(scRNA.data)), list(cluster), mean)
    row.names(scRNA_cluster_avg) <- scRNA_cluster_avg$Group.1
    scRNA_cluster_avg$Group.1 <- c()
    # Tranform the mean expression to zscore
    scRNA_cluster_avg <- t(scRNA_cluster_avg)
    write.csv(scRNA_cluster_avg, file = sprintf("%s/scRNA-seq_mean_expression.csv", out.prefix))
  }
  scRNA_cluster_zscore <- t(apply(scRNA_cluster_avg, 1, zscore))
  
  # Calculate the expression proportion of each gene in each cluster
  cat("Calculate each gene's expressed proportion across clusters in scRNA-seq.\n")
  if(!force & file.exists(sprintf("%s/scRNA-seq_proportion.csv", out.prefix))){
    cat(sprintf("Read existed file in %s/scRNA-seq_proportion.csv\n", out.prefix))
    scRNA_cluster_pro <- read.csv(sprintf("%s/scRNA-seq_proportion.csv", out.prefix), check.names = F, row.names = 1)
  }else if(!is.null(sc.proportion)){
    cat("Use provided scRNA_cluster_pro.\n")
    scRNA_cluster_pro <- sc.proportion
  }else{
    scRNA_cluster_pro <- aggregate(t(as.matrix(scRNA.data)), list(cluster), function(x){sum(x > 0) / length(x)})
    row.names(scRNA_cluster_pro) <- scRNA_cluster_pro$Group.1
    scRNA_cluster_pro$Group.1 <- c()
    scRNA_cluster_pro <- t(scRNA_cluster_pro)
    write.csv(scRNA_cluster_pro, file = sprintf("%s/scRNA-seq_proportion.csv", out.prefix))
  }

  # Calculate the correlation matrix and adjust it to remove the background
  cat("Calculate the gene-cluster correlation matrix and remove the self-exprssed genes in each cluster.\n")
  gene2group.corr.raw <- gene2group.corr <- list()
  if(!force & file.exists(sprintf("%s/correlation_matrix_exp_%s_pro_%s.csv",out.prefix,exp.cutoff,pro.cutoff))){
    cat(sprintf("Read existed file in %s/correlation_matrix_exp_%s_pro_%s.csv\n",out.prefix,exp.cutoff,pro.cutoff))
    gene2group.corr.matrix <- read.csv(sprintf("%s/correlation_matrix_exp_%s_pro_%s.csv",out.prefix,exp.cutoff,pro.cutoff), check.names = F, row.names = 1)
  }else{
    for(group in names(bulk.abundance)){
      gene2group.corr.raw[[group]] <- gene2group.corr[[group]] <- apply(bulk.data, 1, function(x){cor(x, bulk.abundance[[group]], method = "pearson")})
      gene2group.corr[[group]][scRNA_cluster_avg[,group] >= exp.cutoff & scRNA_cluster_pro[,group] >= pro.cutoff] = 0
    }
    gene2group.corr.raw.matrix <- do.call(cbind.data.frame, gene2group.corr.raw)
    row.names(gene2group.corr.raw.matrix) <- names(gene2group.corr.raw[[1]])
    gene2group.corr.matrix <- do.call(cbind.data.frame, gene2group.corr)
    row.names(gene2group.corr.matrix) <- names(gene2group.corr[[1]])
    write.csv(gene2group.corr.raw.matrix, file = sprintf("%s/correlation_raw_matrix_exp_%s_pro_%s.csv",out.prefix,exp.cutoff,pro.cutoff))
    write.csv(gene2group.corr.matrix, file = sprintf("%s/correlation_matrix_exp_%s_pro_%s.csv",out.prefix,exp.cutoff,pro.cutoff))
    
  }

  # Calculate the enrich score of each cluster
  cat("Calculate the enrichment score of each cluster.\n")
  if(!force & file.exists(sprintf("%s/enrichment_zscore_exp_%s_pro_%s_ngenes_%s.csv",out.prefix,exp.cutoff,pro.cutoff,ngenes.enrich))){
    cat(sprintf("Read existed file in %s/enrichment_zscore_exp_%s_pro_%s_ngenes_%s.csv\n",out.prefix,exp.cutoff,pro.cutoff,ngenes.enrich))
    enrich.zscore.matrix <- read.csv(sprintf("%s/enrichment_zscore_exp_%s_pro_%s_ngenes_%s.csv",out.prefix,exp.cutoff,pro.cutoff,ngenes.enrich), check.names = F, row.names = 1)
  }else{
    enrich.score <- list()
    for(i in 1:ncol(gene2group.corr.matrix)){
      top.corr.gene <- row.names(gene2group.corr.matrix)[order(gene2group.corr.matrix[,i], decreasing = TRUE)][1:ngenes.enrich]
      top.corr.gene.exp.matrix <- scRNA_cluster_zscore[top.corr.gene,]
      enrich.score[[colnames(gene2group.corr.matrix)[i]]] <- colMeans(top.corr.gene.exp.matrix)
    }
    enrich.score.matrix <- do.call(cbind.data.frame, enrich.score)
    enrich.zscore.matrix <- apply(enrich.score.matrix, 2, zscore)
    write.csv(enrich.zscore.matrix, file = sprintf("%s/enrichment_zscore_exp_%s_pro_%s_ngenes_%s.csv",out.prefix,exp.cutoff,pro.cutoff,ngenes.enrich))
  }
  
  # Plot the network
  library(igraph)
  Nodes <- apply(enrich.zscore.matrix, 2, function(x){
    row.names(enrich.zscore.matrix)[which(x > enrich.cutoff)] %>% data.frame(stringsAsFactors = F)
  })
  network.m <- matrix(data = 0, nrow = length(unique(cluster)), ncol = length(unique(cluster)), dimnames = list(unique(cluster), unique(cluster)))
  for(i in names(Nodes)){
    network.m[Nodes[[i]]$.,i] <- 1
  }
  diag(network.m) <- 0
  nodes <- data.frame(Sub_Cluster = row.names(network.m),
                      size = rowSums(network.m) + colSums(network.m)) %>% filter(size != 0)
  links <- melt(network.m) %>% filter(value != 0)
  net <- graph_from_data_frame(d = links, vertices = nodes, directed = T)
  V(net)$size <- nodes$size * 2
  save(nodes, links, net, file = sprintf("%s/network_cutoff_exp_%s_pro_%s_ngenes_%s_cutoff_%s.rda",out.prefix,exp.cutoff,pro.cutoff,ngenes.enrich,enrich.cutoff))
  pdf(sprintf("%s/network_exp_%s_pro_%s_ngenes_%s_cutoff_%s.pdf",out.prefix,exp.cutoff,pro.cutoff,ngenes.enrich,enrich.cutoff), width = 10, height = 10)
  plot(net, edge.arrow.size=.5, vertex.color="gold", 
       vertex.frame.color="lightgrey", vertex.label.color="black", 
       vertex.label.cex=0.8, vertex.label.dist=0.5)
  dev.off()
}


# >>Ligand-receptor pairs----
# Sub Function
C2Cmatrix <- function(x){
  ## Transform a long LR-pairs matrix to a wide matrix
  L.group <- qstrsplit(names(x), pattern = " > ", select.region = 1)
  R.group <- qstrsplit(names(x), pattern = " > ", select.region = 2)
  x.df <- data.frame(x, L.group, R.group, stringsAsFactors = F)
  x.group <- sort(unique(x.df$L.group))
  ngroup <- length(x.group)
  c2c.matrix <- matrix(rep(0, ngroup * ngroup), nrow = ngroup)
  row.names(c2c.matrix) <- colnames(c2c.matrix) <- x.group
  for(i in 1:nrow(c2c.matrix)){
    c2c.matrix[i,] <- as.vector(x.df %>% dplyr::filter(L.group == x.group[i]) %>% dplyr::arrange(R.group) %>% dplyr::select(x) %>% t())
  }
  return(c2c.matrix)
}

# Main Function
CCommunication <- function(Expression.matrix, groupid, interaction.used, use.permutation = F, permutation.mean, permutation.sd, exp.cutoff, per.cutoff, show.progress = T) {
  ## LR-based Cell-cell communication analysis
  ##
  ## Args:
  #' @Expression.matrix: A Gene * Cell matrix.
  #' @groupid: The groupid of each cell.
  #' @interaction.used: The ligand and receptor used to calculate.
  #' @use.permutation: Whether to calculate the p-value for LR pairs exp.
  #' @permutation.mean: Permutation parameter used to calcualte the p-value, 
  #' only used when use.permutation is TRUE.
  #' @permutation.sd: Permutation parameter used to calcualte the p-value, 
  #' only used when use.permutation is TRUE.
  #' @exp.cutoff: The cutoff used for whether a gene is expressed in a cell.
  #' @per.cutoff: The cutoff used for the percentage of cells with a gene expressed
  #' in a cluster.
  #' @show.progress: Whether to show the progress.
  #'
  #' Returns:
  #' A list includes 7 matrix: LR.per, LR.exp, LR.pairs.exp, LR.pairs.exp.pvalue,
  #' LR.pairs.counts, LR.pairs.attraction, LR.pairs.affinity. The last 4 when 
  #' use.permutation is TRUE.
  
  library(data.table)
  # Genes and ligand receptor pairs used
  LRName.used <- interaction.used$LRName
  nInteraction <- length(LRName.used)
  ngroup <- length(unique(groupid))
  
  # The percentage of cells with LR genes expressed in each cluster
  if (show.progress) {
    cat("Calculate the percentage of cells with LR genes expressed in each cluster\n")
  }
  
  Expression.matrix.DT <- 
    setDT(data.frame(t(Expression.matrix), groupid = groupid, check.names = F))
  
  LR.per <- data.frame(Expression.matrix.DT[, lapply(.SD, function(x)
    sum(x > exp.cutoff) / length(x)), keyby = .(groupid)], check.names = F)
  rownames(LR.per) <- LR.per$groupid
  LR.per$groupid <- c()
  
  # The mean expression of LR genes in each cluster
  if (show.progress) {
    cat("Calculate the mean expression of LR genes in each cluster\n")
  }
  
  LR.exp <- data.frame(Expression.matrix.DT[, lapply(.SD, mean), keyby = .(groupid)], check.names = F)
  rownames(LR.exp) <- LR.exp$groupid
  LR.exp$groupid <- c()
  LR.exp.norm <- LR.exp
  LR.exp.norm[LR.per < per.cutoff] <- 0
  
  # The expression of LR pairs in each cluster pair
  if (show.progress) {
    cat("Calculate the expression of LR pairs in each cluster pair\n")
  }
  
  LR.pairs.exp <-
    data.frame(matrix(rep(0, nInteraction * ngroup * ngroup), nrow = nInteraction))
  row.names(LR.pairs.exp) <- LRName.used
  
  for (i in 1:ngroup) {
    column.matched <- (ngroup * (i - 1) + 1):(ngroup * i)
    colnames(LR.pairs.exp)[column.matched] <-
      paste0(row.names(LR.exp)[i], " > ", row.names(LR.exp))
    LR.pairs.exp[LRName.used, column.matched] <-
      t(LR.exp.norm[i,interaction.used$Ligand])[,1] * t(LR.exp.norm[, interaction.used$Receptor])
  }
  LR.pairs.exp <- Matrix(as.matrix(LR.pairs.exp), sparse = T)
  
  # return result
  result <- list(
    LR.per = LR.per,
    LR.exp = LR.exp,
    LR.pairs.exp = LR.pairs.exp
  )
  if(use.permutation){
    permutation.mean.vector <- as.vector(permutation.mean)
    permutation.sd.vector <- as.vector(permutation.sd)
    LR.pairs.exp.vector <- as.vector(LR.pairs.exp)
    LR.pairs.exp.pvalue.vector <- mapply(pnorm, LR.pairs.exp.vector, permutation.mean.vector, permutation.sd.vector, lower.tail = FALSE)
    LR.pairs.exp.pvalue.vector[LR.pairs.exp.vector == 0] = 1
    LR.pairs.exp.pvalue <- matrix(LR.pairs.exp.pvalue.vector, nrow = nrow(LR.pairs.exp))
    row.names(LR.pairs.exp.pvalue) <- row.names(LR.pairs.exp)
    colnames(LR.pairs.exp.pvalue) <- colnames(LR.pairs.exp)
    LR.pairs.counts.vector <- colSums(LR.pairs.exp.pvalue < 0.05)
    LR.pairs.counts <- C2Cmatrix(LR.pairs.counts.vector)
    LR.pairs.exp.temp <- LR.pairs.exp
    LR.pairs.exp.temp[LR.pairs.exp.pvalue >= 0.05] <-  0
    LR.pairs.attraction.vector <- colSums(LR.pairs.exp.temp)
    LR.pairs.attraction <- C2Cmatrix(LR.pairs.attraction.vector)
    LR.pairs.affinity <- LR.pairs.attraction + t(LR.pairs.attraction)
    result <- c(result, list(LR.pairs.exp.pvalue = LR.pairs.exp.pvalue,
                             LR.pairs.counts = LR.pairs.counts,
                             LR.pairs.attraction = LR.pairs.attraction,
                             LR.pairs.affinity = LR.pairs.affinity))
  }
  return(result)
}

DoPermutation <- function(Expression.matrix, groupid, permutation.mode = "LR", use.permutation = F, n.permutation = 1000, exp.cutoff = 3, per.cutoff = 0.2, show.progress = T, out.prefix = "", do.return = F){
  ## Permutation step in LR-based cell-cell communication analysis, "LR" mode should be run first
  ##
  ## Args:
  #' @Expression.matrix: A Gene * Cell matrix.
  #' @groupid: The groupid of each cell.
  #' @permutation.mode: "LR" for ligand-receptor pairs expression permutation 
  #' and "C2C" for cluster-cluster counts/attraction/affinity permutation. 
  #' @use.permutation: Whether to calculate the p-value for LR pairs exp, when use 
  #' the permutation in "C2C" mode, should be given as the path of the permutation 
  #' paramter rda file created in "LR" mode.
  #' @n.permutation: Premutation times.
  #' @exp.cutoff: The cutoff used for whether a gene is expressed in a cell.
  #' @per.cutoff: The cutoff used for the percentage of cells with a gene expressed
  #' in a cluster.
  #' @show.progress: Whether to show the progress.
  #' @out.prefix: The path to save the output.
  #' @do.return: Whether to return a list as result.
  #'
  #' Returns:
  #' Two rda file storing the permutation lists and permutation parmaters in both two
  #' modes. A list with the p-value of LR.pairs counts/attraction/affinity in "C2C" mode. 
  
  library(pbapply)
  # Load ligand-receptor pairs
  load("./raw_data/interaction_190619.rda")
  
  # Create the output directory
  Expression.matrix <- as.matrix(Expression.matrix)
  out.dir <- sprintf("%s_e%s_p%s", out.prefix, exp.cutoff, per.cutoff)
  dir.name <- dirname(out.dir)
  dir.create(dir.name, showWarnings = F, recursive = T)
  
  # Cells used
  group.used <- names(table(groupid))[(table(groupid) >= 5 / per.cutoff)]
  cells.used <- groupid %in% group.used
  Expression.matrix <- Expression.matrix[, cells.used]
  groupid <- groupid[cells.used]
  
  # Genes and ligand receptor pairs used
  genes.in.expression <- row.names(Expression.matrix)
  interaction.used <- interaction_all %>% filter(Ligand %in% genes.in.expression & Receptor %in% genes.in.expression)
  interaction.used[] <- lapply(interaction.used, as.character)
  genes.used <- unique(c(interaction.used$Ligand, interaction.used$Receptor))
  Expression.matrix <- Expression.matrix[genes.used, ]
  
  if(permutation.mode == "LR"){
    
    permutation.groupid <- sample(groupid)
    time_used <- system.time(
      CCommunication(Expression.matrix, permutation.groupid, interaction.used, use.permutation = F, exp.cutoff = exp.cutoff, per.cutoff = per.cutoff, show.progress = F)
    )
    cat(paste0("The time need to do the permutation is: ", lubridate::seconds_to_period(time_used[3] * n.permutation), "\n"))
    
    if(show.progress){
      cat("Do the permutation\n")
    }
    permutation.list <- pblapply(1:n.permutation, function(i){
      permutation.groupid <- sample(groupid)
      return(CCommunication(Expression.matrix, permutation.groupid, interaction.used, use.permutation = F, exp.cutoff = exp.cutoff, per.cutoff = per.cutoff, show.progress = F)$LR.pairs.exp)
    })
    
    if(show.progress){
      cat("Save the permutation result\n")
    }
    save(permutation.list, file = sprintf("%s_n%s_%s_mode_permutation_list.rda", out.dir, n.permutation, permutation.mode))
    
    if(show.progress){
      cat("Calculate the mean and sd of each LR pairs in each cluster pairs\n")
    }
    permutation.summary <- pblapply(1:nrow(permutation.list[[1]]), function(i){
      temp <- do.call(rbind, lapply(permutation.list, function(x) x[i,]))
      permutation.mean <- apply(temp, 2, mean)
      permutation.sd <- apply(temp, 2, sd)
      return(list(mean = permutation.mean, sd = permutation.sd))
    })
    permutation.mean <- do.call(rbind, lapply(permutation.summary, function(x) x$mean))
    permutation.sd <- do.call(rbind, lapply(permutation.summary, function(x) x$sd))
    save(permutation.mean, permutation.sd, file = sprintf("%s_n%s_%s_mode_permutation_parameter.rda", out.dir, n.permutation, permutation.mode))
    
    if(do.return){
      return(list(permutation.mean = permutation.mean, permutation.sd = permutation.sd))
    }
    
  }else if(permutation.mode == "C2C"){
    load(use.permutation)
    time_used <- system.time(
      result <- CCommunication(Expression.matrix, groupid, interaction.used, use.permutation = T, permutation.mean = permutation.mean, permutation.sd = permutation.sd, exp.cutoff = exp.cutoff, per.cutoff = per.cutoff, show.progress = T)
    )
    write.csv(result$LR.per, file = sprintf("%s_LR_percentage.csv", out.dir))
    write.csv(result$LR.exp, file = sprintf("%s_LR_expression.csv", out.dir))
    write.csv(as.matrix(result$LR.pairs.exp), file = sprintf("%s_LR_pairs_expression.csv", out.dir))
    write.csv(result$LR.pairs.exp.pvalue, file = sprintf("%s_LR_pairs_expression_pvalue.csv", out.dir))
    write.csv(result$LR.pairs.counts, file = sprintf("%s_LR_pairs_counts.csv", out.dir))
    write.csv(result$LR.pairs.attraction, file = sprintf("%s_LR_pairs_attraction.csv", out.dir))
    write.csv(result$LR.pairs.affinity, file = sprintf("%s_LR_pairs_affinity.csv", out.dir))
    assign(paste0(basename(out.dir), "_interaction"), result)
    save(list = paste0(basename(out.dir), "_interaction"), file = sprintf("%s_interaction.rda", out.dir))
    if(do.return){return(result)}else{return()}
  }
}

# >>DEGenes----
DEGenes <- function(expression_matrix, groupid, out.prefix) {
  ## Differential expressed genes in multiple groups calculated by ANOVA.
  ##
  ## Args:
  #' @expression_matrix: Gene*cell matrix.
  #' @groupid: The groupid of each cell.
  #' @out.prefix: The output filepath with the file name prefix.
  ##
  ## Returns:
  ## A csv file names as outname_all_de_genes.csv in OUTDIR.
  library(pbapply)
  library(ROCR)
  library(stringr)
  library(genefilter)
  out.dir <- dirname(out.prefix)
  dir.create(out.dir, showWarnings = F, recursive = T)
  Expression.matrix <- as.matrix(expression_matrix)  # matrix format
  ngroup <- length(unique(groupid))  # number of groups
  group_names <-
    GroupPair(exprs = Expression.matrix[1, ], group = factor(groupid))
  pboptions(type = "timer",
            style = 1,
            char = "=")
  de.genes_all <- pblapply(1:nrow(Expression.matrix), function(i) {
    temp <-
      data.frame(exprs = Expression.matrix[i, ], group = factor(groupid))
    fml <- aov(exprs ~ group, data = temp)
    if (!is.na(summary(fml)[[1]][1, 'Pr(>F)'])) {
      # make sure that there is fml output
      if (summary(fml)[[1]][1, 'Pr(>F)'] < 0.05) {
        # p-value < 0.05
        diff_gene_info <- c(row.names(Expression.matrix)[i],
                            summary(fml)[[1]][1, 'F value'],
                            summary(fml)[[1]][1, 'Pr(>F)'])  # geneid F-value p-value
        out.tukeyHSD <- TukeyHSD(fml)$group  # ANOVA result
        exp.percent <-
          aggregate(temp$exprs, by = list(temp$group), Expression.Per) # expression percentage
        exp.avg <-
          aggregate(temp$exprs, by = list(temp$group), mean)
        exp.sd <- aggregate(temp$exprs, by = list(temp$group), sd)
        diff_gene_full_info <-
          c(
            diff_gene_info,
            out.tukeyHSD[, "diff"],
            out.tukeyHSD[, "p adj"],
            exp.percent$x,
            exp.avg$x,
            exp.sd$x
          )
        return(diff_gene_full_info)
      }
    }
  })
  de.genes_all <-
    data.frame(do.call(rbind, de.genes_all), check.names = F)
  colnames(de.genes_all) <-
    c(
      "Symbol",
      "Fvalue",
      "Pvalue",
      paste0("HSD.diff.", group_names$Group.pair),
      paste0("HSD.padj.", group_names$Group.pair),
      paste0("Exp.Per.", group_names$Group),
      paste0("Exp.avg.", group_names$Group),
      paste0("Exp.sd.", group_names$Group)
    )  # set the colnames of the final output
  de.genes_all <-
    dplyr::arrange(de.genes_all, desc(Fvalue))  # sort by F-value
  row.names(de.genes_all) <-
    1:nrow(de.genes_all)  # set the rownames of the final output
  write.csv(de.genes_all, sprintf("%s_all_de_genes.csv", out.prefix))  # save in the csv format
  
}

# >>GrpDEGenes----
GrpDEGenes <- function(out.prefix, cutoff = "ALL", expression_matrix, groupid, logFC = 1, adj.p = 0.01, doAUC = FALSE) {
  ## Filter each group marker genes form ANOVA result.
  ##
  ## Args:
  #' @out.prefix: The output filepath with the file name prefix.
  #' @cutoff: The gene should be highly expressed than how many other groups. Could
  #' be a number, a vector or "ALL".By default is "ALL" which will output all the
  #' possible cutoff.
  #' @expression_matrix: Gene*cell matrix.
  #' @groupid: The groupid of each cell.
  #' @logFC: The logFC cutoff.
  #' @adj.p: The adjusted p-value cutoff.
  ##
  ## Returns:
  ## A list of csv file names as outname_cutoff_de_genes.csv and
  ## outname_cutoff_group_de_genes.csv in OUTDIR.
  library(ROCR)
  degenes <- read.csv(sprintf("%s_all_de_genes.csv", out.prefix), row.names = 1, stringsAsFactors = F, check.names = F)  # read the output of DEGenes
  ngroup <- length(unique(groupid))  # number of groups
  expression_matrix <- as.matrix(expression_matrix)
  group_names <- GroupPair(exprs = expression_matrix[1, ], group = factor(groupid))
  if (cutoff[1] == "ALL") {
    # calculate the result of all possible cutoff
    for (cutoffi in 1:(ngroup - 1)) {
      sprintf("Cutoff = %s", cutoffi)
      GrpDEGenes(out.prefix, cutoffi, expression_matrix, groupid)
    }
  } else if (length(cutoff) > 1) {
    for (cutoffi in cutoff) {
      sprintf("Cutoff = %s", cutoffi)
      GrpDEGenes(out.prefix, cutoffi, expression_matrix, groupid)
    }
  } else{
    tukeyHSD_trans <- list()
    out.tukeyHSD.name <- list()
    de.genes <- c()
    de.genes.grp <- c()
    for (i in 1:ngroup) {
      tukeyHSD_trans[[group_names$Group[i]]] <- c(rep(1, i - 1), rep(-1, ngroup - i))  # a 1/-1 vector to unifying the logFC in tukeyHSD output to group i / other groups
      out.tukeyHSD.name[[group_names$Group[i]]] <- GrepGroupPair(group_names$Group[i], group_names$Group.pair)
    }
    for (j in group_names$Group) {
      cell.In <- which(groupid == j)
      cell.Out <- which(groupid != j)
      out.tukeyHSD.group <- degenes[, c(paste0("HSD.diff.", out.tukeyHSD.name[[j]]),
                                        paste0("HSD.padj.", out.tukeyHSD.name[[j]]))]
      Diffgroupn <- rowSums(as.matrix(out.tukeyHSD.group[, 1:(ngroup - 1)]) %*% diag(tukeyHSD_trans[[j]]) > logFC & out.tukeyHSD.group[, ngroup:(2 * ngroup - 2)] < adj.p)  # a gene is considered highly expressed in one group than another when logFC > 1 and adj.p-value < 0.01
      de.genes <-
          degenes[Diffgroupn >= cutoff, c(
            "Symbol", "Fvalue", "Pvalue",
            paste0("HSD.diff.", out.tukeyHSD.name[[j]]),
            paste0("HSD.padj.", out.tukeyHSD.name[[j]]),
            paste0("Exp.Per.", group_names$Group),
            paste0("Exp.avg.", group_names$Group),
            paste0("Exp.sd.", group_names$Group))]
      if (!is.na(de.genes[1, 3])) {
        de.genes[, paste0("HSD.diff.", out.tukeyHSD.name[[j]])] <-
          as.matrix(de.genes[, paste0("HSD.diff.", out.tukeyHSD.name[[j]])]) %*% diag(tukeyHSD_trans[[j]])  # %*% is to unify the logFC
        if (nrow(de.genes) > 1) {
          de.genes$Exp.Mean.In <-
            rowMeans(expression_matrix[de.genes$Symbol, cell.In])  # mean expression in the group
          de.genes$Exp.Mean.Out <-
            rowMeans(expression_matrix[de.genes$Symbol, cell.Out])  # mean expression out the group
          de.genes$Exp.Sd.In <-
            rowSds(expression_matrix[de.genes$Symbol, cell.In])  # standard deviation in the group
          de.genes$Exp.Sd.Out <-
            rowSds(expression_matrix[de.genes$Symbol, cell.Out])  # standard deviation in the group
          de.genes$Exp.Per.In <-
            apply(expression_matrix[de.genes$Symbol, cell.In], 1, Expression.Per)  # expression percentage in the group
          de.genes$Exp.Per.Out <-
            apply(expression_matrix[de.genes$Symbol, cell.Out], 1, Expression.Per)  # expression percentage out the group
        } else{
          de.genes$Exp.Mean.In <-
            mean(expression_matrix[de.genes$Symbol, cell.In])  # mean expression in the group
          de.genes$Exp.Mean.Out <-
            mean(expression_matrix[de.genes$Symbol, cell.Out]) # mean expression out the group
          de.genes$Exp.Sd.In <-
            sd(expression_matrix[de.genes$Symbol, cell.In])  # standard deviation in the group
          de.genes$Exp.Sd.Out <-
            sd(expression_matrix[de.genes$Symbol, cell.Out])  # standard deviation in the group
          de.genes$Exp.Per.In <-
            Expression.Per(expression_matrix[de.genes$Symbol, cell.In])  # expression percentage in the group
          de.genes$Exp.Per.Out <-
            Expression.Per(expression_matrix[de.genes$Symbol, cell.Out])  # expression percentage out the group
          }
        category <- as.numeric(groupid == j)
        AUC <- pblapply(1:nrow(de.genes), function(k) {
          pred <- prediction(expression_matrix[de.genes$Symbol[k], ], category)
          pauc <- performance(pred, measure = "auc")
          return(pauc@y.values[[1]])
        })
        de.genes <- cbind(de.genes, AUC = do.call(rbind, AUC))
        de.genes <- data.frame(de.genes, stringsAsFactors = F, check.names = F) %>% dplyr::arrange(desc(AUC))
        write.csv(de.genes, sprintf("%s_cutoff%g_%s_de_genes.csv", out.prefix, cutoff, j))  # save each group result in csv format
        column_used <-
          !(colnames(de.genes) %in% c(
            paste0("HSD.diff.", out.tukeyHSD.name[[j]]),
            paste0("HSD.padj.", out.tukeyHSD.name[[j]])))
        de.genes.grp <- rbind(de.genes.grp, cbind(de.genes[, column_used], Group = j))
      }
    }
    de.genes.grp <- arrange.vars(de.genes.grp, vars = c("Group" = 2, "AUC" = 3))
    write.csv(de.genes.grp, sprintf("%s_cutoff%g_all_de_genes.csv", out.prefix, cutoff))  # save all groups result in csv format
  }
}

# >> Find Differential Genes----
FindDEGenes <- function(expression_matrix, groupid, out.prefix, cutoff = "ALL", logFC = 1, adj.p = 0.01, ncell = 1, doAUC = FALSE){
  groupid_n <- table(groupid)
  groupname_kept <- names(groupid_n)[groupid_n > ncell]
  cells_kept <- groupid %in% groupname_kept
  expression_matrix <- expression_matrix[,cells_kept]
  groupid <- groupid[cells_kept]
  cat("Call DEGenes function\n")
  DEGenes(expression_matrix = expression_matrix, 
          groupid = groupid,
          out.prefix = out.prefix)
  cat("Call GrpDEGenes function\n")
  GrpDEGenes(out.prefix = out.prefix,
             cutoff = cutoff,
             logFC = logFC,
             adj.p = adj.p,
             expression_matrix = expression_matrix,
             groupid = groupid,
             doAUC = doAUC)
}

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

GroupPair <- function(exprs, group) {
  # the group order and the tukeyHSD output pairs order
  group <- as.factor(group)
  pre.temp <- data.frame(exprs = exprs, group = group)
  pre.fml <- aov(exprs ~ group, data = pre.temp)
  pre.out.tukeyHSD <- TukeyHSD(pre.fml)$group
  return(list(
    Group = levels(group),
    Group.pair = row.names(pre.out.tukeyHSD)
  ))
}

GrepGroupPair <- function(group, group.pair) {
  if (!grepl("-", group)) {
    group.pair.df <- str_split_fixed(group.pair, "-", 4)
    group.in.group.pair.df <- group.pair.df == group
    group_name <-
      group.pair[group.in.group.pair.df[, 1] |
                   group.in.group.pair.df[, 2] |
                   group.in.group.pair.df[, 3] | group.in.group.pair.df[, 4]]
  } else {
    group_name <- grep(group, group.pair, value = T)
  }
  return(group_name)
}

# >> Arrange variables in a data frame----
arrange.vars <- function(data, vars) {
  ## Arrange df vars by position
  ##
  ## Args:
  #' @data: a data frame
  #' @vars: must be a named vector, e.g. c("var.name"=1)
  ##
  ## Returns:
  ## new data frame
  
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  
  ##sanity checks
  stopifnot(!any(duplicated(var.nms)),!any(duplicated(var.pos)))
  stopifnot(is.character(var.nms),
            is.numeric(var.pos))
  stopifnot(all(var.nms %in% data.nms))
  stopifnot(all(var.pos > 0),
            all(var.pos <= var.nr))
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[!(data.nms %in% var.nms)]
  stopifnot(length(out.vec) == var.nr)
  
  ##re-arrange vars by position
  data <- data[, out.vec]
  return(data)
}

# >> Summary data by group----
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                      conf.interval = .95, .drop = TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2 (xx[[col]], na.rm=na.rm),
                     mean = mean    (xx[[col]], na.rm=na.rm),
                     median = median(xx[[col]], na.rm=na.rm),
                     sd   = sd      (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
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

# >>Infer CNV----
#' Calculate smooth gene expression
#' @param gexp.mat assay matrix, should be scaled column is cell , row is gene
#' @param gene.range genomic range for the assays ,names is symbol
#' @param chrs Chromosomes to be used (default: paste0('chr', c(1:22, 'X')))
#' @param windows.size Window size for sliding window mean. Must be odd number. (default: 101)
#' @param return.list Boolean for whether to return splitted matrix list or a combined matrix

calc_smoothed_exprs  <- function(gexp.mat, gene.range, chrs = paste0("chr", c(1:22, "X")), 
                                 window.size = 101, verbose = F,  return.list = T) {
  # split 
  tl <- split_assay_by_chromosome(gexp.mat, gene.range)
  tl <- tl[chrs]
  # smoothing
  tlsmooth <- lapply(names(tl), function(chrom) {
    if(verbose) cat(chrom, "\n")
    ws = window.size
    d <- tl[[chrom]]
    
    if ( nrow(d) < window.size) {
      ws = nrow(d)
      if ( (ws %% 2) == 0) ws = ws - 1
    }
    d.smooth <- caTools::runmean(d, k = ws, endrule = "mean", align = "center")
    .getSmoothRange <- function(genes, ws, verbose = F) {
      out <- rep(NA, length(genes))
      k2 <- ws %/% 2
      
      for(i in 1:k2) {
        out[i] = paste0(strsplit(genes[1], "_|:")[[1]][3], "-", strsplit(genes[i+k2], "_|:|-")[[1]][4])
      }
      
      for(i in (k2+1):(length(genes)-k2 )) {
        
        out[i] = paste0(strsplit(genes[i-k2], "_|:")[[1]][3], "-", strsplit(genes[(i-k2)+ws-1], "_|:|-")[[1]][4])
      }
      for(i in (length(genes)-k2+1):length(genes)) {
        out[i] = paste0(strsplit(genes[i], "_|:")[[1]][3], "-", strsplit(genes[length(genes)], "_|:|-")[[1]][4])
      }
      
      out
    }
    smoothRanges = .getSmoothRange(rownames(d), ws = ws )
    rownames(d.smooth) = paste0(chrom,":",smoothRanges)
    colnames(d.smooth) = colnames(d)
    d.smooth
    
  })
  names(tlsmooth) <- names(tl)
  
  if( return.list) {
    return(tlsmooth)
  } else {
    return(do.call("rbind", tlsmooth))
  }
}

#' split assay matrix by chromosome 
#' @param gexp.mat assay matrix, should be scaled column is cell , row is gene
#' @param gene.range genomic range for the assays ,names is symbol
split_assay_by_chromosome <- function(gexp.mat, gene.range) {
  ov.genes = intersect(rownames(gexp.mat), names(gene.range))
  if(length(ov.genes) < 10) {
    info =sprintf("There are less than %s genes overlap, check your gene name", length(ov.genes))
    cat(info, "\n")
  }
  
  gexp.mat = gexp.mat[ov.genes,]
  gene.range = gene.range[ov.genes]
  #--- Turn to data.frame
  gene.range = GenomicRanges::sort(gene.range, ignore.strand = TRUE) # Sort first
  gos = GenomicRanges::as.data.frame(gene.range)
  gexp.mat = gexp.mat[rownames(gos),]
  #--- add position info to genes
  rownames(gos) <- paste0(names(gene.range),"_", paste0(gos$seqnames, ":", gos$start, ":", gos$end))
  rownames(gexp.mat) <- paste0(rownames(gexp.mat), "_", paste0(gos$seqnames, ":", gos$start, ":", gos$end))
  tl <- lapply(split(gos, gos$seqnames) , function(genes.at.chr) {
    gexp.mat[rownames(genes.at.chr),]  
  })
  
  return(tl)
}

# >>Run Nichenetr----
RunNichenetr <- function(seu, cluster, group.by = "Sub_Cluster", receiver_genes, sender_genes = NULL, expressed_genes_receiver = NULL, background_expressed_genes = NULL, expressed_gene_cutoff = 2, best_upstream_ligands_number = 30, species = "homo", curated = FALSE, repo = TRUE, only_genes_mode = FALSE){
  library(nichenetr)
  if(species == "homo"){
    if(exists("ligand_target_matrix_human") & exists("weighted_networks_human") & exists("lr_network_human")){
      if(repo) cat("Required database already exists.\n")
    }else{
      if(repo) cat("Load required database.\n")
      load("/work/lzy/project/utils/nichenetr_human.rda")
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
      load("/work/lzy/project/utils/nichenetr_mouse.rda")
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

RunNichenetr_exp <- function(exp, receiver_oi_genes, receiver_expressed_genes = NULL, sender_genes, receiver_name = "Receiver", sender_name = "Sender", background_expressed_genes = NULL, best_upstream_ligands_number = 30, species = "homo", curated = FALSE, repo = TRUE){
  library(nichenetr)
  if(species == "homo"){
    if(exists("ligand_target_matrix_human") & exists("weighted_networks_human") & exists("lr_network_human")){
      if(repo) cat("Required database already exists.\n")
    }else{
      if(repo) cat("Load required database.\n")
      load("/work/lzy/project/utils/nichenetr_human.rda")
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
      load("/work/lzy/project/utils/nichenetr_mouse.rda")
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
  geneset_oi = receiver_oi_genes
  if(is.null(background_expressed_genes)){
    background_expressed_genes <- row.names(exp)[rowMeans(exp) > 1]
  }
  background_expressed_genes = background_expressed_genes %>% .[. %in% rownames(nichenetr_data[["ligand_target_matrix"]])]
  if(is.null(receiver_expressed_genes)) {receiver_expressed_genes = receiver_oi_genes}
  expressed_genes_receiver = receiver_expressed_genes
  expressed_genes_sender = sender_genes
  
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
      paste0("Prioritized ligands from ", sender_name), "", 
      color = "darkorange", x_axis_position = "top", 
      legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + 
    theme(legend.text = element_text(size = 9),
          axis.ticks = element_blank(),
          axis.title.x = element_text(),
          axis.text.y = element_text(face = "italic"))
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = 
    vis_ligand_target %>% 
    make_heatmap_ggplot(
      "Potential ligands",paste0(receiver_name," upregulated genes"), 
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
      "Potential ligands",paste0("Receptors expressed by ",receiver_name), 
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

# >>SCENIC----
auc_thresh_kmeans_df <-  function (regulonAUC) 
{
  regulons <- rownames(regulonAUC)
  kmeans_thresholds <- list()
  print("Processing regulon distributions...")
  for (regulon_no in 1:length(regulons)) {
    svMisc::progress(regulon_no)
    regulon <- regulons[regulon_no]
    df <- data.frame(auc = regulonAUC[regulon, 
    ], cells = names(regulonAUC[regulon, 
    ]), regulon = regulon)
    df <- df %>% subset(auc > 0)
    km <- kmeans(df$auc, centers = 2)
    df$cluster <- as.factor(km$cluster)
    cluster1_max <- max(subset(df, cluster == 1)$auc)
    cluster2_max <- max(subset(df, cluster == 2)$auc)
    if (cluster1_max > cluster2_max) {
      df <- df %>% mutate(cluster = gsub(2, 3, cluster)) %>% 
        mutate(cluster = gsub(1, 2, cluster)) %>% mutate(cluster = gsub(3, 
                                                                        1, cluster))
    }
    df <- df %>% arrange(desc(auc))
    df_sub <- df %>% subset(cluster == 1)
    auc_thresholds <- df_sub[1, ]$auc
    kmeans_thresholds[[regulon]] <- auc_thresholds
  }
  print("Done evaluating thresholds...")
  return(kmeans_thresholds)
}

binarize_regulons_df <- function(regulonAUC, thresholds) 
{
  binary_regulon_list <- list()
  for (regulon_no in 1:length(names(thresholds))) {
    svMisc::progress(regulon_no)
    regulon <- names(thresholds)[regulon_no]
    auc_df <- data.frame(auc = regulonAUC[regulon, 
    ], cells = names(regulonAUC[regulon, 
    ]))
    auc_df <- auc_df %>% mutate(regulon = if_else(auc >= thresholds[regulon], 1, 0)) %>% select(-auc)
    colnames(auc_df) <- c("cells", regulon)
    binary_regulon_list[[regulon]] <- auc_df
  }
  return(binary_regulon_list)
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

# >>Calculate Proximity Index----
CalProxiIndex <- function(seu, cells, name = ""){
  seu@meta.data[,paste0(name,"_ProxiIndex")] <- 0
  for(cell in cells){
    neighbor_spots <- FindNeighborSpot(seu, cell, dist = 10, self = F)
    for(neighbor_spot in neighbor_spots){
      distance <- dist(rbind(seu@images$slice1@coordinates[neighbor_spot,c("imagerow","imagecol")], seu@images$slice1@coordinates[cell,c("imagerow","imagecol")]))
      seu@meta.data[neighbor_spot,paste0(name,"_ProxiIndex")] <- max(seu@meta.data[neighbor_spots,paste0(name,"_ProxiIndex")], 1/(1+distance))
    }
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

# >>Generate meta clusters----
MetaCluster <- function(exp, cluster, n_cells = NA, n_metacells = NA, max.n_metacells = NA){
  MetaExp <- data.frame(exp = 0)
  for(clusteri in unique(cluster)){
    MetaExp_cluster <- data.frame(exp = 0)
    cells_used <- cluster == clusteri
    exp_used <- exp[,cells_used]
    n_total_cells <- ncol(exp_used)
    if(is.na(n_cells)){
      n_cells <- ceiling(n_total_cells/n_metacells)
    }else if(is.na(n_metacells)){
      n_metacells <- ceiling(n_total_cells/n_cells)
      if(!is.na(max.n_metacells) & n_metacells > max.n_metacells){
        n_metacells = max.n_metacells
      }
    }
    n_total_cells_for_meta <- n_cells * n_metacells
    if(n_total_cells_for_meta <= n_total_cells){
      cells_used_for_meta <- sample(1:n_total_cells, n_total_cells_for_meta)
    }else if(n_total_cells_for_meta > n_total_cells){
      cells_used_for_meta <- c(1:n_total_cells, sample(1:n_total_cells, n_total_cells_for_meta - n_total_cells, replace = T))
    }
    cell_group_for_meta <- split(cells_used_for_meta, sample(rep(1:n_metacells, n_cells)))
    for(i in 1:n_metacells){
      MetaExp_tmp <- rowMeans(exp_used[,cell_group_for_meta[[i]]])
      MetaExp_cluster <- cbind(MetaExp_cluster, MetaExp_tmp)
    }
    MetaExp_cluster <- MetaExp_cluster[,-1]
    colnames(MetaExp_cluster) <- paste0(clusteri,"_",1:n_metacells)
    MetaExp <- cbind(MetaExp, MetaExp_cluster)
  }
  MetaExp <- MetaExp[,-1]
}

# >>Define module score----
CalculateModuleScore <- function(exp, features, iter = 1000, nbin = 20){
  exp <- t(t(exp) - colMeans(exp))
  features_mean_exp <- colMeans(exp[features,])
  gene.avg <- rowMeans(exp)
  gene.avg <- gene.avg[order(gene.avg)]
  gene.cut <- cut_number(x = gene.avg + rnorm(length(gene.avg))/1e+30, 
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

# >>Compare clusters in different datasets----
Compare2Clusters <- function(exp1, cluster1, exp2, cluster2, ngenes = 1000, upper_gene = F, homo_gene = F, name1 = "D1", name2 = "D2"){
  if(homo_gene == T){
    gene_homo <- read.csv("/work/lzy/project/utils/homomuris.csv", row.names = 3, stringsAsFactors = F)
    exp1_gene <- row.names(exp1)
    exp2_gene <- row.names(exp2)
    exp1_gene[exp1_gene %in% row.names(gene_homo)] <- gene_homo[exp1_gene[exp1_gene %in% row.names(gene_homo)] ,"symbol"]
    exp2_gene[exp2_gene %in% row.names(gene_homo)] <- gene_homo[exp2_gene[exp2_gene %in% row.names(gene_homo)] ,"symbol"]
    exp1 <- SetRowNames(exp1,exp1_gene)
    exp2 <- SetRowNames(exp2,exp2_gene)
  }
  if(homo_gene == T | upper_gene == T){
    exp1 <- SetRowNames(exp1,toupper(row.names(exp1)))
    exp2 <- SetRowNames(exp2,toupper(row.names(exp2)))
  }
  genes_used <- intersect(row.names(exp1), row.names(exp2))
  exp1_avg <- aggregate(t(exp1[genes_used,]), list(Cluster = cluster1), mean)
  row.names(exp1_avg) <- exp1_avg$Cluster
  exp1_avg$Cluster <- c()
  exp1_avg <- t(exp1_avg)
  colnames(exp1_avg) <- paste0(name1,"_",colnames(exp1_avg))
  exp2_avg <- aggregate(t(exp2[genes_used,]), list(Cluster = cluster2), mean)
  row.names(exp2_avg) <- exp2_avg$Cluster
  exp2_avg$Cluster <- c()
  exp2_avg <- t(exp2_avg)
  colnames(exp2_avg) <- paste0(name2,"_",colnames(exp2_avg))
  merged_exp <- cbind(exp1_avg, exp2_avg)
  merged_exp <- merged_exp[rowSums(merged_exp) > 0,]
  merged_exp_norm <- preprocessCore::normalize.quantiles(as.matrix(merged_exp))
  row.names(merged_exp_norm) <- row.names(merged_exp)
  colnames(merged_exp_norm) <- c(colnames(exp1_avg),colnames(exp2_avg))
  merged_exp_norm <- sva::ComBat(as.matrix(merged_exp_norm), c(rep(1,ncol(exp1_avg)), rep(2,ncol(exp2_avg))), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  if(ngenes > nrow(merged_exp_norm)){
    genes_used <- row.names(merged_exp_norm)
  }else{
    genes_used <- names(sort(apply(merged_exp_norm,1,sd),decreasing = T))[1:ngenes] 
  }
  pca_analysis <- prcomp(t(merged_exp_norm[genes_used,]))
  pca_analysis <- summary(pca_analysis)
  fit <- pvclust::pvclust(t(pca_analysis$x), method.hclust = "ward.D", method.dist = "euclidean")
  return(list(
    merged_exp = merged_exp_norm, 
    genes_used = genes_used, 
    pca_analysis = pca_analysis, 
    hc_fit = fit)
    )
}

# >>Run force-directed graph----
# a function to compute force-directed graph; return coordinates as a data frame
# arguments:
# pca.df: pca data as a data frame
# snn   : a nearest-neighbor graph as a sparse data matrix
runFDG = function(pca.df, snn, iterations = 600, working.dir, python.addr){
  current.wd = getwd()
  setwd(working.dir)
  # generate unique name for pca data file
  pca.data.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  pca.data.fname = paste(pca.data.fname, ".csv", sep = "")
  # generate unique name for snn file
  snn.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  snn.fname = paste(snn.fname, ".smm", sep = "")
  # generate unique name for fdg coordinates
  fdg.coordinates.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  fdg.coordinates.fname = paste(fdg.coordinates.fname, ".csv", sep = "")
  write.csv(pca.df, pca.data.fname)
  writeMM(obj=snn, file=snn.fname)
  command = gsub(pattern="ITER", replacement=as.character(iterations), paste(python.addr, "/work/lzy/project/utils/make_fdg.py ITER", sep = " "))
  command = paste(command, paste(c(pca.data.fname, snn.fname, fdg.coordinates.fname), collapse = " "), sep = " ")
  system(command, wait = T)
  fdg_coordinates = read.csv(fdg.coordinates.fname, header = FALSE)
  colnames(fdg_coordinates) = c("X", "Y")
  rownames(fdg_coordinates) = rownames(pca.df)
  file.remove(c(pca.data.fname, snn.fname, fdg.coordinates.fname))
  setwd(current.wd)
  return(fdg_coordinates)
}

# >>Convert human and mouse genes----
convertHumanGeneList <- function(x){
  require("biomaRt")
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
  mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 105)
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = x, 
                   mart = human, 
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

convertMouseGeneList <- function(x){
  require("biomaRt")
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
  mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 105)
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = x, 
                   mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

# >>String Split----
qstrsplit <- function(x, pattern = "_", n = 2, select.region = 1){
  library(stringr)
  ## Split the string and return the column you need
  x <- as.character(x)
  splited_entry <- str_split_fixed(x, pattern, n)
  return(splited_entry[,select.region])
}

# >>zscore----
zscore <- function(x){
  return( (x - mean(x))/sd(x) )
}

# >>one-tail test----
t.test.onetail <- function(x, y){
  t1 <- t.test(x, y, alternative = "less")
  t2 <- t.test(x, y, alternative = "greater")
  if(t1$p.value < t2$p.value){
    result <- t1
  }else{
    result <- t2
  }
  return(result)
}
wilcox.test.onetail <- function(x, y){
  w1 <- wilcox.test(x, y, alternative = "less")
  w2 <- wilcox.test(x, y, alternative = "greater")
  if(w1$p.value < w2$p.value){
    result <- w1
  }else{
    result <- w2
  }
  return(result)
}

# >>Use `%>%`----
`%>%` <- dplyr::`%>%`

# >>Set row names and remove duplicates----
SetRowNames <- function(exp, row_names){
  if(sum(duplicated(row_names)) > 0){
    is_duplicates <- duplicated(row_names)
    exp <- exp[!is_duplicates,]
    row_names <- row_names[!is_duplicates]
  }
  row.names(exp) <- row_names
  return(exp)
}

# >>Generate cluster combination list----
ClustComb <- function(cluster_vector){
  cluster_combination_list <- list()
  for(i in 1:(length(cluster_vector)-1)){
    for(j in (i+1):length(cluster_vector)){
      temp <- c(cluster_vector[i],cluster_vector[j])
      cluster_combination_list <- append(cluster_combination_list, list(temp))
    }
  }
  return(cluster_combination_list)
}

# >>Merge Seurat object----
MergeSeurat <- function(seu_list){
  nfile = length(seu_list)
  if(nfile == 2){
    seu = merge(seu_list[[1]],seu_list[[2]])
  }else if(nfile >2){
    seu = merge(seu_list[[1]],seu_list[2:nfile])
  }
}

# >>ReadMtx----
ReadMtx <- function(mtx, cells, features, cell.column = 1, feature.column = 2, cell.sep = "\t", feature.sep = "\t", skip.cell = 0, skip.feature = 0, mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE) {
  library(httr)
  all.files <- list(
    "expression matrix" = mtx,
    "barcode list" = cells,
    "feature list" = features
  )
  for (i in seq_along(along.with = all.files)) {
    uri <- tryCatch(
      expr = {
        con <- url(description = all.files[[i]])
        close(con = con)
        all.files[[i]]
      },
      error = function(...) {
        return(normalizePath(path = all.files[[i]], winslash = '/'))
      }
    )
    err <- paste("Cannot find", names(x = all.files)[i], "at", uri)
    uri <- build_url(url = parse_url(url = uri))
    if (grepl(pattern = '^[A-Z]?:///', x = uri)) {
      uri <- gsub(pattern = '^://', replacement = '', x = uri)
      if (!file.exists(uri)) {
        stop(err, call. = FALSE)
      }
    } else {
      if (!Online(url = uri, seconds = 2L)) {
        stop(err, call. = FALSE)
      }
      if (file_ext(uri) == 'gz') {
        con <- url(description = uri)
        uri <- gzcon(con = con, text = TRUE)
      }
    }
    all.files[[i]] <- uri
  }
  cell.barcodes <- read.table(
    file = all.files[['barcode list']],
    header = FALSE,
    sep = cell.sep,
    row.names = NULL,
    skip = skip.cell,
    stringsAsFactors = F
  )
  feature.names <- read.table(
    file = all.files[['feature list']],
    header = FALSE,
    sep = feature.sep,
    row.names = NULL,
    skip = skip.feature,
    stringsAsFactors = F
  )
  # read barcodes
  bcols <- ncol(x = cell.barcodes)
  if (bcols < cell.column) {
    stop(
      "cell.column was set to ",
      cell.column,
      " but ",
      cells,
      " only has ",
      bcols,
      " columns.",
      " Try setting the cell.column argument to a value <= to ",
      bcols,
      "."
    )
  }
  cell.names <- cell.barcodes[, cell.column]
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  # read features
  fcols <- ncol(x = feature.names)
  if (fcols < feature.column) {
    stop(
      "feature.column was set to ",
      feature.column,
      " but ",
      features,
      " only has ",
      fcols, " column(s).",
      " Try setting the feature.column argument to a value <= to ",
      fcols,
      "."
    )
  }
  if (any(is.na(x = feature.names[, feature.column]))) {
    na.features <- which(x = is.na(x = feature.names[, feature.column]))
    replacement.column <- ifelse(test = feature.column == 2, yes = 1, no = 2)
    if (replacement.column > fcols) {
      stop(
        "Some features names are NA in column ",
        feature.column,
        ". Try specifiying a different column.",
        call. = FALSE
      )
    } else {
      warning(
        "Some features names are NA in column ",
        feature.column,
        ". Replacing NA names with ID from column ",
        replacement.column,
        ".",
        call. = FALSE
      )
    }
    feature.names[na.features, feature.column] <- feature.names[na.features, replacement.column]
  }
  feature.names <- feature.names[, feature.column]
  if (unique.features) {
    feature.names <- make.unique(names = as.character(feature.names))
  }
  data <- readMM(file = all.files[['expression matrix']])
  if (mtx.transpose) {
    data <- t(x = data)
  }
  if (length(x = cell.names) != ncol(x = data)) {
    stop(
      "Matrix has ",
      ncol(data),
      " columns but found ", length(cell.names),
      " barcodes. ",
      ifelse(
        test = length(x = cell.names) > ncol(x = data),
        yes = "Try increasing `skip.cell`. ",
        no = ""
      ),
      call. = FALSE
    )
  }
  if (length(x = feature.names) != nrow(x = data)) {
    stop(
      "Matrix has ",
      nrow(data),
      " rows but found ", length(feature.names),
      " features. ",
      ifelse(
        test = length(x = feature.names) > nrow(x = data),
        yes = "Try increasing `skip.feature`. ",
        no = ""
      ),
      call. = FALSE
    )
  }
  
  colnames(x = data) <- cell.names
  rownames(x = data) <- feature.names
  data <- as.sparse(x = data)
  return(data)
}
