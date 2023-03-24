# Load Packages----
require(ggplot2)
require(ggforce)
library(ComplexHeatmap)
require(ggrepel)
require(gtable)
require(grid)
require(egg)
require(patchwork)
require(cowplot)
require(png)
require(dplyr)
require(reshape2)
require(RColorBrewer)

# Functions----
NoAxesTheme <- function(){
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "mm")
  )
}
theme_mine <- function(base_size = 18, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = base_size),
      strip.text.y = element_text(size = base_size),
      axis.text.x = element_text(size = base_size - 4),
      axis.text.y = element_text(size = base_size - 4, hjust = 1),
      axis.ticks =  element_line(colour = "black"), 
      axis.title.x= element_text(size = base_size - 2),
      axis.title.y= element_text(size = base_size - 2, angle = 90),
      panel.background = element_blank(), 
      panel.border =element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.spacing = unit(1.0, "lines"), 
      plot.background = element_blank(), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1)
    )
}
RotatedAxis <- function(angle = 45, ...){
  rotated.theme <- theme(axis.text.x = element_text(angle = angle, hjust = 1), validate = TRUE, ...)
  return(rotated.theme)
}

calcCenters <- function(cluster, reduction) {
  df <- data.frame(Cluster = as.factor(cluster), Dim1 = reduction[,1], Dim2 = reduction[,2])
  centers <- df %>%
    group_by(Cluster) %>%
    summarise(mean_x = median(Dim1),
              mean_y = median(Dim2))
  return(centers)
  
}

addLabels <- function(centers, label_size = 3, label_short = FALSE) {
  if (label_short) centers <- suppressWarnings(
    tidyr::separate(centers, Cluster, into = c("Cluster", "Cluster_long"), extra = "drop"))
  ggrepel::geom_text_repel(data = centers,
                              aes(x = mean_x, y = mean_y),
                              label = centers$Cluster,
                              size = label_size,
                              alpha = 0.8,
                              segment.alpha = 0.8,
                              force = 2,
                              segment.size = 0.5,
                              arrow = arrow(length = unit(0.01, 'npc')))
}

# AI friendly plot----
require(ggplot2)
require(png)
require(ggedit)

deepcopy <- function(p) {
  unserialize(serialize(p, NULL))
}

gg0point <- function(plot){
  newplot <- deepcopy(plot)
  newplot$layers[[1]]$aes_params$size <- -1
  return(newplot)
}

ggAIplot <- function(plot, calibration = c(0,0,0,0), width = 10, height = 10, dpi = 300){
  plot.png <- plot + theme_nothing()
  ggsave(plot.png, filename = "./temp.png", width = width, height = height, dpi = dpi)
  img <- readPNG("./temp.png")
  file.remove("./temp.png")
  blank.plot <- gg0point(plot)
  range.values <- c(
    ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$y.range
  )
  plot.pdf <- blank.plot +
    annotation_raster(img, 
                      xmin = range.values[1] + calibration[1], 
                      xmax = range.values[2] + calibration[2],
                      ymin = range.values[3] + calibration[3], 
                      ymax = range.values[4] + calibration[4])
  return(plot.pdf)
}

annotation_raster.grid <- 
  function (raster, xmin, xmax, ymin, ymax, interpolate = FALSE, data) {
    raster <- grDevices::as.raster(raster)
    layer(data = data, mapping = NULL, stat = StatIdentity, 
          position = PositionIdentity, geom = GeomRasterAnn, inherit.aes = FALSE, 
          params = list(raster = raster, xmin = xmin, xmax = xmax, 
                        ymin = ymin, ymax = ymax, interpolate = interpolate))
  }

ggAIplot.grid <- function(plot, facet.by, calibration = c(0,0,0,0), width = 10, height = 10, dpi = 300){
  plot.data <- plot$data
  plot.pdf <- gg0point(plot)
  range.values <- c(
    ggplot_build(plot = plot.pdf)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = plot.pdf)$layout$panel_params[[1]]$y.range
  )
  for(variable in unique(plot.data[,facet.by])){
    temp.plot.png <- plot + theme_nothing()
    temp.plot.png$data <- temp.plot.png$data[temp.plot.png$data[,facet.by] == variable,]
    ggsave(temp.plot.png, filename = "./temp.png", width = width, height = height, dpi = dpi)
    img <- readPNG("./temp.png")
    file.remove("./temp.png")
    ggimg <- annotation_raster.grid(img, 
                                    xmin = range.values[1] + calibration[1], 
                                    xmax = range.values[2] + calibration[2],
                                    ymin = range.values[3] + calibration[3], 
                                    ymax = range.values[4] + calibration[4],
                                    data = temp.plot.png$data[1,])
    plot.pdf <- plot.pdf + ggimg
  }
  return(plot.pdf)
}

ggAISpatial <- function(plot, width = 10, height = 10, dpi = 300, add.legend = TRUE){
  plot <- plot + theme(legend.position = "right")
  legend <- get_legend(plot)
  plot.png <- plot + theme_nothing()
  ggsave(plot.png, filename = "./temp.png", width = width, height = height, dpi = dpi)
  img <- readPNG("./temp.png")
  file.remove("./temp.png")
  range.values <- c(
    ggplot_build(plot = plot)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = plot)$layout$panel_params[[1]]$y.range
  )
  p <- ggplot() +
    lims(x = c(range.values[1], range.values[2]), 
        y = c(range.values[3], range.values[4])) + 
    theme_nothing() +
    annotation_raster(img,
                      xmin = range.values[1], xmax = range.values[2],
                      ymin = range.values[3], ymax = range.values[4])
  if(add.legend){
    p <- p + legend + plot_layout(widths = c(5,1))
  }
  return(p)
}

# Calculate linear regression model and label the function in ggplot2----
# devtools::source_gist("524eade46135f6348140")
# load("/data2/lzy/CRC_CD45_project/PAPER/final_release/scripts/smooth_func.rda")

# >>Summary the meatdata plot----
SummaryMatadataPlot <- function(metadata, plot.type = "Violin", plot.item, color.by, label = TRUE, pt.size = 1.5, facet.nrow = NULL, facet.ncol = NULL, legend.nrow = NULL, legend.ncol = NULL, reduction.used = "tSNE", font.size = 6, show.cutoff = 10, vector.friendly = TRUE, png.arguments = c(10,10,300), add.ggplot = NULL){
  ## Plot the summary of metadata.
  ##
  ## Args:
  metadata <- as.data.frame(metadata)
  if(plot.type == "Violin"){
    p <- ggplot(metadata, aes_string(x = color.by, y = plot.item)) +
      geom_violin(aes_string(fill = color.by, color = color.by), alpha = 0.6) + 
      geom_jitter(alpha = 0.1, color = "lightgrey") +
      geom_boxplot(width = 0.1, outlier.colour = NA, aes_string(color = color.by)) +
      theme_cowplot(font_size = 12) +
      theme(legend.position = "NULL",
            axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
      guides(fill = guide_legend(ncol = legend.ncol, nrow = legend.nrow)) +
      add.ggplot
  }else if(plot.type == "Pie"){
    require(scales)
    metadata[,plot.item] <- factor(metadata[,plot.item])
    metadata[,color.by] <- factor(metadata[,color.by])
    pie_df <- list()
    for(i in sort(unique(metadata[,plot.item]))){
      pie_df[[i]] <- metadata[metadata[,plot.item] == i,] %>% select(plot.item, color.by) %>% table() %>% melt()
      pie_df[[i]]$per <- round(pie_df[[i]]$value/sum(pie_df[[i]]$value) * 100, 2)
      pie_df[[i]]$Label <- paste0(pie_df[[i]]$per, "%")
      pie_df[[i]]$Label[pie_df[[i]]$per <= show.cutoff] <- NA
      pie_df[[i]] <- pie_df[[i]] %>% 
        mutate(end = 2 * pi * cumsum(value)/sum(value),
               start = lag(end, default = 0),
               middle = 0.5 * (start + end),
               hjust = ifelse(middle > pi, 1, 0),
               vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
    }
    pie_df <- do.call(rbind, pie_df)
    p <- ggplot(pie_df, aes_string(fill = color.by)) + 
      geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                       start = start, end = end), color = "white") +
      geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = Label,
                    hjust = hjust, vjust = vjust), 
                size = font.size, na.rm = TRUE) +
      coord_fixed() +
      facet_wrap(paste0(".~", plot.item), nrow = facet.nrow, ncol = facet.ncol) +
      scale_x_continuous(limits = c(-2, 2),
                         name = "", breaks = NULL, labels = NULL) +
      scale_y_continuous(limits = c(-2, 2),
                         name = "", breaks = NULL, labels = NULL) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.justification = "right",
        plot.margin = unit(c(-0.15,0,-0.2,0),"in"),
        strip.text = element_text(size = font.size * 3),
        legend.title = element_blank(),
        legend.text = element_text(size = font.size * 2.5)
      ) +
      guides(fill = guide_legend(ncol = legend.ncol, nrow = legend.nrow)) +
      add.ggplot
  }else if(plot.type == "Bar"){
    metadata[,plot.item] <- factor(metadata[,plot.item])
    metadata[,color.by] <- factor(metadata[,color.by])
    p <- ggplot(metadata, aes_string(x = plot.item, fill = color.by)) +
      geom_bar(position = "fill", alpha = 0.7, color = "white") +
      labs(x = "", y = "Percentage") + theme_cowplot(font_size = font.size * 3) + 
      scale_x_discrete(limits = rev(sort(unique(metadata[,plot.item])))) +
      theme(legend.title = element_blank(), legend.position = "top") + coord_flip() +
      guides(fill = guide_legend(override.aes = list(size = font.size), ncol = legend.ncol, nrow = legend.nrow)) +
      add.ggplot
  }else if(plot.type == "Scatter"){
    x.char <- paste0(plot.item, "_", reduction.used, "_1")
    y.char <- paste0(plot.item, "_", reduction.used, "_2")
    p <- ggplot(metadata, aes_string(x = x.char, y = y.char)) + 
      geom_point(aes_string(color = color.by), size = pt.size, alpha = 0.6) +
      xlab(paste0(reduction.used, "1")) + 
      ylab(paste0(reduction.used, "2")) + 
      theme_cowplot(font_size = font.size * 3) +
      theme(legend.title = element_blank()) +
      guides(colour = guide_legend(override.aes = list(size = font.size))) +
      add.ggplot
    if(vector.friendly){
      png.plot <- p + NoAxesTheme() + theme(legend.position = "NULL")
      ggsave(png.plot, filename = "./temp.png", width = png.arguments[1], height = png.arguments[2], dpi = png.arguments[3])
      img <- readPNG("./temp.png")
      file.remove("./temp.png")
      blank.plot <-  ggplot(metadata, aes_string(x = x.char, y = y.char)) + 
        geom_point(aes_string(color = color.by), size = -1, alpha = 0.6) +
        xlab(paste0(reduction.used, "1")) + 
        ylab(paste0(reduction.used, "2")) + 
        theme_cowplot(font_size = font.size * 3) +
        theme(legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = font.size))) +
        add.ggplot
      range.values <- c(
        ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$x.range,
        ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$y.range
      )
      p <- blank.plot +
        annotation_raster(img, xmin = range.values[1], xmax = range.values[2],
                          ymin = range.values[3], ymax = range.values[4])
    }
    if(label){
      cluster_center_tsne <- aggregate(metadata[,c(x.char, y.char)], list(metadata[,color.by]), median)
      p <- p + geom_label_repel(data = cluster_center_tsne, 
                         aes_string(x = x.char, y = y.char, label = "Group.1"), 
                         size = font.size/1.5)
    }
  }
  return(p)
}

# >>Gene Plot----
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x)
  ymax <- m + sd(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

gene_filter <- function(genes.input, genes.all){
  genes.used <- c()
  for(gene in genes.input){
    if(!(gene %in% genes.all)){
      cat(sprintf("Gene --%s-- not in the expression matrix!", gene), "\n")
    }else{
      genes.used <- c(genes.used, gene)
    }
  }
  return(genes.used)
}

GenePlot <- function(expression.matrix, metadata, genelist, order.genes = F, plot.type = "Violin", color.palette = "YlOrRd", brewer.pal.num = 8, color.direction = 1, font.size = 16, no.axes = FALSE, no.legend = FALSE, ncol = NULL, group.by, reduction.tag = "Global", reduction.used = "tSNE", pt.size = 1.5, vector.friendly = TRUE, png.arguments = c(10,10,100), color.by.group = F, size.range = 10, percent.exp.cutoff = 1, zscore.transformed = T, exp.filter.per = 0, add.ggplot = NULL){
  ## Plot the gene violin plot or scatter plot.
  ##
  ## Args:
  #' @expression.matrix: Gene*cell matrix.
  #' @metadata: The columns used are group plot by and tSNE.tag.
  #' @genelist: A vector including all the genes or a list of genes to plot.
  #' @plot.type: Should be "Violin" or "Scatter" or "TwoPanel" or "Bubble", default by "Violin".
  #' @color.palette: The color palette in Rcolorbrewer.
  #' @color.direction: The order of the color palette, should be 1 or -1.
  #' @font.size: The size of the text in the plot.
  #' @no.axes: Whether to plot the axises.
  #' @no.legend: Whether to plot the legend.
  #' @ncol: Used in "Scatter" and "Violin". The column of the plot when multiple genes in genelist.
  #' @group.by: Used in "Violin", "TwoPanel" and "Bubble" plot type. Which group should the plot grouped by.
  #' @reduction.tag: Which tSNE or UMAP coordinate should be used in the metadata, the counterpart column is "reduction.tag_tSNE_1".
  #' @reduction.used: Which type of reduction should be used in the metadata, "tSNE" or "UMAP".
  #' @pt.size: Used in "Scatter". The point size.
  #' @vector.friendly: Used in "Scatter". If TRUE, the plot will be saved as an AI-friendly form.
  #' @png.arguments: Used in "Scatter". If vector.friendly is TRUE, this parameter sets the size of the temporary png file saved, recommended in agreement with the final pdf/png size.
  #' @color.by.group: Used in "Violin" plot type. Should each group colored by cluster type or mean expression.
  #' @size.range: Used in "Bubble". The max size of bubble.
  #' @percent.exp.cutoff: Used in "Bubble". Only cells with gene expression higher than this cutoff will be considered as a gene-expressed cell.
  #' @zscore.transformed: Used in "Bubble". Whether to use z-score transformed expression in bubble plot.
  #' @exp.filter.per: Used in "Bubble". The gene expression matrix quantile [exp.filter.per, 1 - exp.filter.per] will be calculated first. And all the gene expression will be limitted to this expression quantile.
  #' @add.ggplot: Used in "Scatter", "Violin" and "Bubble". A ggplot command added behind the plot function with "+". 
  if(!is.list(genelist)){
    genes.used <- gene_filter(genelist, row.names(expression.matrix))
  }else{
    genes.used <- lapply(genelist, function(x){gene_filter(x, row.names(expression.matrix))})
    genes.used <- stack(genes.used)
    colnames(genes.used) <- c("Gene", "Gene.group")
  }
  
  expression.matrix <- as.matrix(expression.matrix)
  
  if(color.direction == -1){
    myColorPalette <- colorRampPalette(rev(brewer.pal(8, color.palette)))
  }else{
    myColorPalette <- colorRampPalette(brewer.pal(8, color.palette))
  }
  
  if(plot.type == "Scatter"){
    plotlist <- c()
    for(gene in genes.used){
      gene_plot.df <- data.frame(Expression = expression.matrix[gene,],
                                 tSNE_1 = metadata[,paste0(reduction.tag,"_",reduction.used,"_1")],
                                 tSNE_2 = metadata[,paste0(reduction.tag,"_",reduction.used,"_2")])
      plotlist[[gene]] <- 
        ggplot(gene_plot.df %>% arrange(Expression), aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(color = Expression), alpha = 0.7, size = pt.size) +
        scale_colour_gradientn(colors = myColorPalette(100)) + 
        labs(x = "", y = "") + ggtitle(gene) +
        theme_cowplot(font_size = font.size) +
        add.ggplot
      if(vector.friendly){
        png.plot <- plotlist[[gene]] + NoAxesTheme() + theme(legend.position = "NULL", title = element_blank())
        ggsave(png.plot, filename = "./temp.png", width = png.arguments[1], height = png.arguments[2], dpi = png.arguments[3])
        img <- readPNG("./temp.png")
        file.remove("./temp.png")
        blank.plot <-  ggplot(gene_plot.df, aes(x = tSNE_1, y = tSNE_2)) +
          geom_point(aes(color = Expression), alpha = 0.7, size = -1) +
          scale_colour_gradientn(colors = myColorPalette(100)) + 
          labs(x = "", y = "") + ggtitle(gene) +
          theme_cowplot(font_size = font.size) +
          add.ggplot
        range.values <- c(
          ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$x.range,
          ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$y.range
        )
        plotlist[[gene]] <- blank.plot +
          annotation_raster(img, xmin = range.values[1], xmax = range.values[2],
                            ymin = range.values[3], ymax = range.values[4])
      }
      if(no.axes){
        plotlist[[gene]] <- plotlist[[gene]] + NoAxesTheme()
      }
      if(no.legend){
        plotlist[[gene]] <- plotlist[[gene]] + theme(legend.position = "NULL")
      }
    }
    p <- plot_grid(plotlist = plotlist, ncol = ncol, align = "hv")
  }
  
  if(plot.type == "Violin"){
    plotlist <- c()
    groupid <- factor(metadata[,group.by])
    first.gene <- genes.used[1]
    if(!is.null(ncol)){
      last.gene <- tail(genes.used, ncol)
    }else{
      last.gene <- tail(genes.used, 1)
    }
    for(gene in genes.used){
      gene_ave_group <- aggregate(expression.matrix[gene,], list(groupid), mean)
      gene_plot.df <- data.frame(Expression = expression.matrix[gene,],
                                 Group = groupid,
                                 Ave_group = gene_ave_group$x[groupid])
      if(!color.by.group){
        plotlist[[gene]] <- 
          ggplot(gene_plot.df, aes(x = Group, y = Expression)) + 
          geom_violin(aes(color = Ave_group, fill = Ave_group), scale = "width") +
          scale_fill_gradientn(colors = myColorPalette(100)) + 
          scale_color_gradientn(colors = myColorPalette(100)) +
          labs(x = "", y = gene) + add.ggplot
      }else{
        plotlist[[gene]] <- 
          ggplot(gene_plot.df, aes(x = Group, y = Expression)) + 
          geom_violin(aes(fill = Group), color = "black", scale = "width", alpha = 0.8) +
          labs(x = "", y = gene) + add.ggplot
      }
      if(!(gene %in% last.gene)){
        plotlist[[gene]] <-
          plotlist[[gene]] + 
          theme_cowplot(font_size = font.size) +
          NoAxesTheme() +
          theme(axis.title.y = element_text(angle = 0, size = font.size))
      }else{
        if(length(genes.used) == 1){
          if(!color.by.group){
            plotlist[[gene]] <- 
              ggplot(gene_plot.df, aes(x = Group, y = Expression)) + 
              geom_violin(aes(fill = Ave_group), color = "black", scale = "width") +
              scale_fill_gradientn(colors = myColorPalette(100)) +
              add.ggplot
          }else{
            plotlist[[gene]] <- 
              ggplot(gene_plot.df, aes(x = Group, y = Expression)) + 
              geom_violin(aes(fill = Group), color = "black", scale = "width") +
              add.ggplot
          }
          plotlist[[gene]] <- 
            plotlist[[gene]] +
            stat_summary(fun.y = "mean", geom = "point", size = 2, color = "black") +
            labs(x = "", y = "Exp") + ggtitle(gene) +
            theme_cowplot(font_size = font.size) +
            theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
        }else{
          plotlist[[gene]] <- 
            plotlist[[gene]] +
            theme_cowplot(font_size = font.size) + NoAxesTheme() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = font.size), axis.title.y = element_text(angle = 0, size = font.size))
        }
      }
      if(no.axes){
        plotlist[[gene]] <- plotlist[[gene]] + NoAxesTheme()
      }
      if(no.legend){
        plotlist[[gene]] <- plotlist[[gene]] + theme(legend.position = "NULL")
      }
      if(gene == first.gene){
        p <- plotlist[[gene]]
      }else{
        p <- p + plotlist[[gene]]
      }
    }
    if(length(genes.used) > 1){
      p <- p + plot_layout(ncol = ncol)
    }
  }
  
  if(plot.type == "TwoPanel"){
    groupid <- factor(metadata[,group.by])
    # Boxplot in the left panel----
    signature.mean.expression <- colMeans(expression.matrix[genes.used,])
    boxplot.data <- data.frame(Expression = signature.mean.expression, Group = groupid)
    p_boxplot <- 
      ggplot(boxplot.data, aes(x = Group, y = Expression)) + 
      geom_boxplot(color = "#1979B5", fill = "#1979B5", 
                   outlier.colour = "black", outlier.shape = 1) + 
      theme_bw() +
      theme(legend.position = "NULL",
            axis.text.y = element_text(size = font.size),
            axis.title = element_blank()) +
      scale_x_discrete(limits = rev(levels(as.factor(groupid)))) +
      coord_flip() +
      add.ggplot
    dat <- ggplot_build(p_boxplot)$data[[1]]
    dat$xsp <- 1/2 * (dat$xmax + dat$xmin) - 1/4 * (dat$xmax - dat$xmin)
    dat$xep <- 1/2 * (dat$xmax + dat$xmin) + 1/4 * (dat$xmax - dat$xmin)
    p_boxplot <- 
      p_boxplot + 
      geom_segment(data = dat,  aes(x = xmin, xend = xmax, y = middle, yend = middle), colour = "#2BA147", size = 2) +
      geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymin, yend = ymin), colour = "black", size = 2) +
      geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymax, yend = ymax), colour = "black", size = 2)
    # Heatmap in the right panel----
    gene.median.expression <- aggregate(t(expression.matrix[genes.used,]), list(Group = groupid), mean)
    gene.median.expression[, 2:(length(genes.used) + 1)] <- apply(gene.median.expression[, 2:(length(genes.used) + 1)], 2, function(x){
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
                                       hjust = 1, size = font.size),
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
    # Final plot----
    g1 <- ggplotGrob(p_boxplot)
    g2 <- p_heatmap
    fg1 <- gtable_frame(g1, height = unit(10, "null"))
    fg2 <- gtable_frame(g2, height = unit(40, "null"))
    p <- gtable_frame(gtable_cbind(fg1, fg2))
    grid.draw(p)
  }
  if(plot.type == "Bubble"){
    if(!is.list(genelist)){
      genes.used <- data.frame(Gene = genes.used, Gene.group = "A")
    }
    groupid <- factor(metadata[,group.by])
    if(zscore.transformed){
      gene.mean.matrix <- aggregate(t(expression.matrix[genes.used$Gene,]), list(Group = groupid), function(x){mean(x[x!=0])})
      gene.mean.matrix[is.na(gene.mean.matrix)] <- 0
      gene.mean.zscore <- apply(gene.mean.matrix[,2:ncol(gene.mean.matrix)], 2, function(x){(x - mean(x)) / sd(x)})
      gene.mean.zscore.df <- data.frame(Group = gene.mean.matrix[,1], gene.mean.zscore)
      gene.mean <- melt(gene.mean.zscore.df, id.vars = "Group", variable.name = "Gene", value.name = "Exp")
    }else{
      gene.mean <- melt(aggregate(t(expression.matrix[genes.used$Gene,]), list(Group = groupid), function(x){mean(x[x!=0])}), id.vars = "Group", variable.name = "Gene", value.name = "Exp")
    }
    gene.per <- melt(aggregate(t(expression.matrix[genes.used$Gene,]), list(Group = groupid), function(x){sum(x > percent.exp.cutoff) / length(x)}), id.vars = "Group", variable.name = "Gene", value.name = "Per")
    plot.data <- merge(merge(gene.mean, gene.per), genes.used)
    if(order.genes){
      plot.data.order <- plot.data %>% filter(Group %in% Gene.group)
      plot.data.order <- plot.data.order[order(plot.data.order$Gene.group, -plot.data.order$Exp),]
      plot.data$Gene <- factor(plot.data$Gene, levels = unique(plot.data.order$Gene))
    }else{
      plot.data$Gene <- factor(plot.data$Gene, levels = unique(plot.data$Gene))
    }
    gene.mean.quantile <- quantile(plot.data$Exp, c(exp.filter.per, 1 - exp.filter.per))
    plot.data$Exp <- pmax(plot.data$Exp, gene.mean.quantile[1])
    plot.data$Exp <- pmin(plot.data$Exp, gene.mean.quantile[2])
    p <- ggplot(plot.data, aes(x = Gene, y = Group)) +
      geom_point(aes(size = Per, fill = Exp), shape = 21) +
      facet_grid(~Gene.group, scales = "free", space = "free") +
      theme_minimal() + 
      scale_y_discrete(limits = rev(levels(plot.data$Group))) +
      theme(
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.1,"in"),
        strip.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = font.size),
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5, size = font.size),
        legend.text = element_text(size = font.size),
        legend.title = element_text(size = font.size),
        panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "lightgrey",
          linetype = "dashed",
          size = 0.2
        )) + 
      labs(x = "", y = "") +
      scale_fill_gradientn(colours = myColorPalette(100), 
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barwidth = 0.8)) +
      scale_size_continuous(breaks = seq(0, 0.8, 0.2), range = c(1,size.range)) +
      add.ggplot
  }
  
  return(p)
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

# >>CellChat----
netVisual_aggregate_new <- function (object, signaling, signaling.name = NULL, color.use = NULL, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, clusters.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, vertex.weight.max = NULL, vertex.size.max = 15, weight.scale = TRUE, edge.weight.max = NULL, edge.width.max = 8, layout = c("hierarchy", "circle", "chord"), thresh = 0.05, from = NULL, to = NULL, bidirection = NULL, vertex.size = NULL, pt.title = 12, title.space = 6, vertex.label.cex = 0.8, group = NULL, cell.order = NULL, small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, ...) {
  layout <- match.arg(layout)
  if (is.null(clusters.use)){
    clusters.use <- unique(object@idents)
  }
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents)[clusters.use])
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 3, sum) != 0]
  }
  else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 0]
  }
  if (length(pairLR.name.use) == 0) {
    stop(paste0("There is no significant communication of ", signaling.name))
  }
  else {
    pairLR <- pairLR[pairLR.name.use, ]
  }
  nRow <- length(pairLR.name.use)
  prob <- prob[, , pairLR.name.use]
  pval <- pval[, , pairLR.name.use]
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
    pval <- replicate(1, pval, simplify = "array")
  }
  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1, 2), sum)
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow = c(1, 2), ps = pt.title)
    netVisual_hierarchy1(prob.sum[clusters.use, clusters.use], sources.use = sources.use, targets.use = targets.use, vertex.receiver = vertex.receiver, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex, ...)
    netVisual_hierarchy2(prob.sum[clusters.use, clusters.use], sources.use = sources.use, targets.use = targets.use, vertex.receiver = setdiff(1:nrow(prob.sum[clusters.use, clusters.use]), vertex.receiver), remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex, ...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    gg <- recordPlot()
  }
  else if (layout == "circle") {
    prob.sum <- apply(prob, c(1, 2), sum)
    gg <- netVisual_circle(prob.sum[clusters.use, clusters.use], remove.isolate = remove.isolate, top = top, sources.use = sources.use, targets.use = targets.use, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex, ...)
  }
  else if (layout == "chord") {
    prob.sum <- apply(prob, c(1, 2), sum)
    gg <- netVisual_chord_cell_internal(prob.sum[clusters.use, clusters.use], sources.use = sources.use, targets.use = targets.use, color.use = color.use, remove.isolate = remove.isolate, group = group, cell.order = cell.order, lab.cex = vertex.label.cex, small.gap = small.gap, big.gap = big.gap, scale = scale, reduce = reduce, title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
  }
  return(gg)
}

netVisual_individual_new <- function (object, signaling, signaling.name = NULL, pairLR.use = NULL, color.use = NULL, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, clusters.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex = 0.8, weight.scale = TRUE, edge.weight.max = NULL, edge.width.max = 8, layout = c("hierarchy", "circle", "chord"), height = 5, thresh = 0.05, group = NULL, cell.order = NULL, small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, nCol = NULL, ...) {
  layout <- match.arg(layout)
  if (is.null(clusters.use)){
    clusters.use <- unique(object@idents)
  }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents)[clusters.use])
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = F)
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  if (!is.null(pairLR.use)) {
    if (is.data.frame(pairLR.use)) {
      pairLR.name <- intersect(pairLR.name, as.character(pairLR.use$interaction_name))
    } else {
      pairLR.name <- intersect(pairLR.name, as.character(pairLR.use))
    }
    if (length(pairLR.name) == 0) {
      stop("There is no significant communication for the input L-R pairs!")
    }
  }
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 3, sum) != 0]
  }else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 0]
  }
  if (length(pairLR.name.use) == 0) {
    stop(paste0("There is no significant communication of ", signaling.name))
  }else {
    pairLR <- pairLR[pairLR.name.use, ]
  }
  nRow <- length(pairLR.name.use)
  prob <- prob[, , pairLR.name.use]
  pval <- pval[, , pairLR.name.use]
  if (is.null(nCol)) {
    nCol <- min(length(pairLR.name.use), 2)
  }
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
    pval <- replicate(1, pval, simplify = "array")
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max = max(prob)
  }
  if (layout == "hierarchy") {
    par(mfrow = c(nRow, 2), mar = c(5, 4, 4, 2) + 0.1)
    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[, , i][clusters.use,clusters.use]
      netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, title.name = signalName_i, ...)
      netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i), vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, title.name = signalName_i, ...)
    }
    gg <- recordPlot()
  }
  else if (layout == "circle") {
    par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd = TRUE)
    gg <- vector("list", length(pairLR.name.use))
    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[, , i][clusters.use, clusters.use]
      gg[[i]] <- netVisual_circle(prob.i, remove.isolate = remove.isolate, sources.use = sources.use, targets.use = targets.use, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, title.name = signalName_i)
    }
  }
  else if (layout == "chord") {
    par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd = TRUE)
    gg <- vector("list", length(pairLR.name.use))
    for (i in 1:length(pairLR.name.use)) {
      title.name <- pairLR$interaction_name_2[i]
      net <- prob[, , i][clusters.use, clusters.use]
      gg[[i]] <- netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, group = group, cell.order = cell.order, lab.cex = vertex.label.cex, small.gap = small.gap, big.gap = big.gap, scale = scale, reduce = reduce, title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
    }
  }
  return(gg)
}

netVisual_diffInteraction_new <- function (object, comparison = c(1, 2), measure = c("count", "weight", "count.merged", "weight.merged"), color.use = NULL, 
          title.name = NULL, sources.use = NULL, targets.use = NULL, 
          remove.isolate = FALSE, top = 1, weight.scale = FALSE, vertex.weight = 20, 
          vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex = 1, 
          vertex.label.color = "black", edge.weight.max = NULL, edge.width.max = 8, 
          alpha.edge = 0.6, label.edge = FALSE, edge.label.color = "black", 
          edge.label.cex = 0.8, edge.curved = 0.2, shape = "circle", 
          layout = in_circle(), margin = 0.2, arrow.width = 1, arrow.size = 0.2) 
{
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, "#b2182b", 
                               "#2166ac")
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                                1])] <- loop.angle[edge.start[which(edge.start[, 
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}
# >>Monocle3_heatmap----
heatmap_sm = function (obj, assay.name = "exprs", out.prefix = NULL, ncell.downsample = NULL, 
                       ave.by = NULL, columns = NULL, columns.order = NULL, gene.desc = NULL, 
                       colSet = list(), pdf.width = 16, pdf.height = 15, do.scale = TRUE, 
                       z.lo = NULL, z.hi = NULL, z.step = 1, exp.title = "Exp", n_smooth = 50,
                       do.clustering.row = T, do.clustering.col = T, dend.col = FALSE, 
                       dend.row = FALSE, clustering.distance = "spearman", clustering.method = "complete", 
                       k.row = 1, k.col = 1, palette.name = NULL, annotation_legend_param = list(), 
                       ann.bar.height = 1.5, mytitle = "", ...) 
{
  library(sscClust)
  library(gridBase)
  library(circlize)
  library(zoo)
  if (!is.null(gene.desc) && ("Group" %in% colnames(gene.desc)) && 
      ("geneID" %in% colnames(gene.desc))) {
    obj <- obj[gene.desc$geneID, ]
  }
  if (!is.null(ncell.downsample) && ncell.downsample < ncol(obj)) {
    obj <- obj[, sample(seq_len(ncol(obj)), ncell.downsample)]
  }
  n <- nrow(obj)
  m <- ncol(obj)
  if (n < 3) {
    loginfo(sprintf("Too few genes: n=%s", n))
    return(NULL)
  }
  if (m < 3) {
    loginfo(sprintf("Too few samples: m=%s", m))
    return(NULL)
  }
  if (is.null(ave.by)) {
    obj <- ssc.assay.hclust(
      obj, assay.name = assay.name,
      order.col = 
        if(is.logical(dend.col) && FALSE == dend.col) do.clustering.col else FALSE, 
      order.row = 
        if(is.logical(dend.row) && FALSE == dend.row) do.clustering.row else FALSE, 
      clustering.distance = "spearman", clustering.method = "complete",
      k.row = 1, k.col = 1)
  }
  else {
    obj <- ssc.average.cell(obj, assay.name = assay.name, 
                            column = ave.by, ret.type = "sce")
    columns <- intersect(ave.by, columns)
    columns.order <- intersect(ave.by, columns.order)
  }
  ha.col <- NULL
  annDF <- data.frame()
  if (!is.null(columns)) {
    if (!is.null(columns.order)) {
      obj <- ssc.order(obj, columns.order = columns.order)
    }
    annDF <- as.data.frame(colData(obj)[columns])
    if (length(colSet) == 0) {
      for (i in seq_along(columns)) {
        x <- columns[i]
        if (class(colData(obj)[, x]) == "numeric") {
          if (all(colData(obj)[, x] <= 1) && all(colData(obj)[, x] >= 0)) {
            Y.level <- c(0, 1)
          }
          else {
            Y.level <- pretty(colData(obj)[, x], n = 8)
          }
          colSet[[x]] <- colorRamp2(seq(Y.level[1], 
                                        Y.level[length(Y.level)], length = 7), 
                                    rev(brewer.pal(n = 7, "RdYlBu")), space = "LAB")
          annotation_legend_param[[x]] <- list(color_bar = "continuous", 
                                               legend_direction = "horizontal", 
                                               legend_width = unit(4, "cm"), 
                                               legend_height = unit(2, "cm"))
        }
        else {
          group.value <- sort(unique(colData(obj)[, x]))
          colSet[[x]] <- structure(auto.colSet(length(group.value), name = "Accent"), names = group.value)
        }
      }
    }
    g.show.legend <- T
    ha.col <- ComplexHeatmap::HeatmapAnnotation(
      df = annDF, col = colSet, 
      show_legend = g.show.legend, 
      simple_anno_size = unit(ann.bar.height, "cm"), annotation_legend_param = annotation_legend_param)
  }
  obj <- ssc.order(obj, columns.order = NULL, gene.desc = gene.desc)
  dat.plot <- as.matrix(assay(obj, assay.name))
  rownames(dat.plot) <- unname(rowData(obj)$display.name)
  if (do.scale) {
    rowM <- rowMeans(dat.plot, na.rm = T)
    rowSD <- apply(dat.plot, 1, sd, na.rm = T)
    dat.plot <- sweep(dat.plot, 1, rowM)
    dat.plot <- sweep(dat.plot, 1, rowSD, "/")
    if (!is.null(z.lo)) {
      dat.plot[dat.plot < z.lo] <- z.lo
    }
    if (!is.null(z.hi)) {
      dat.plot[dat.plot > z.hi] <- z.hi
    }
  }
  else {
    tmp.var <- pretty((dat.plot), n = 8)
    if (is.null(z.lo)) {
      z.lo <- tmp.var[1]
    }
    if (is.null(z.hi)) {
      z.hi <- tmp.var[length(tmp.var)]
    }
    if (is.null(z.step)) {
      z.step <- tmp.var[2] - tmp.var[1]
    }
  }
  if (!is.null(out.prefix)) {
    pdf(sprintf("%s.pdf", out.prefix), width = pdf.width, 
        height = pdf.height)
  }
  par(mar = c(4, 12, 4, 4))
  plot.new()
  title(main = mytitle, cex.main = 2)
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  if (is.null(palette.name)) {
    exp.palette <- rev(brewer.pal(n = 7, name = ifelse(do.scale, 
                                                       "RdBu", "RdYlBu")))
  }
  else {
    if (palette.name=="blue2green2red"){
      bks <- seq(-3.1, 3.1, by=0.1)
      exp.palette <- colorRamps::blue2green2red(length(bks) - 1)
    }
    else{
      exp.palette <- rev(brewer.pal(n = 7, name = palette.name))
    }
  }
  # smooth value
  dat.plot = t(apply(dat.plot, 1, function(x){rollmean(x, n_smooth, fill="extend")}))
  if(is.null(z.lo)){z.lo <- min(dat.plot)}
  if(is.null(z.hi)){z.hi <- max(dat.plot)}
  set.seed(1)
  ht <- ComplexHeatmap::Heatmap(
    dat.plot, name = exp.title, 
    col = colorRamp2(seq(z.lo, z.hi, length = 100), colorRampPalette(exp.palette)(100), space = "LAB"),
    column_dend_height = unit(6, "cm"), row_dend_width = unit(1, "cm"), 
    column_names_gp = grid::gpar(fontsize = 12 * 28/max(m, 32)), 
    row_names_gp = grid::gpar(fontsize = 18 *28/max(n, 32)), 
    show_heatmap_legend = T, 
    row_names_max_width = unit(10, "cm"), 
    cluster_columns = dend.col, cluster_rows = dend.row,
    row_dend_reorder = FALSE, column_dend_reorder = FALSE,
    heatmap_legend_param = list(grid_width = unit(0.8, "cm"),
                                grid_height = unit(0.8, "cm"), 
                                at = seq(z.lo, z.hi, z.step), 
                                title_gp = grid::gpar(fontsize = 14, fontface = "bold"), 
                                label_gp = grid::gpar(fontsize = 12),
                                color_bar = "continuous"), top_annotation = ha.col, 
                                ...)
  ComplexHeatmap::draw(ht, newpage = FALSE)
  if (!is.null(ha.col)) {
    for (i in seq_along(names(ha.col@anno_list))) {
      ComplexHeatmap::decorate_annotation(
        names(ha.col@anno_list)[i],
        {
          grid.text(names(ha.col@anno_list)[i], unit(-4, "mm"), 
                    gp = grid::gpar(fontsize = 14), just = "right")
        })
    }
  }
  if (!is.null(out.prefix)) {
    dev.off()
  }
  return(ht)
}

# >>Regression function in ggplot2----
stat_smooth_func <- function(mapping = NULL, data = NULL, geom = "smooth", position = "identity", ..., method = "auto", formula = y ~ x, se = TRUE, n = 80, span = 0.75, fullrange = FALSE, level = 0.95, method.args = list(), na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, xpos = NULL, ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto(
  "StatSmooth",
  Stat,
  
  setup_params = function(data, params) {
    # Figure out what type of smoothing to do: loess for small datasets,
    # gam with a cubic regression basis for large data
    # This is based on the size of the _largest_ group.
    if (identical(params$method, "auto")) {
      max_group <- max(table(data$group))
      
      if (max_group < 1000) {
        params$method <- "loess"
      } else {
        params$method <- "gam"
        params$formula <-
          y ~ s(x, bs = "cs")
      }
    }
    if (identical(params$method, "gam")) {
      params$method <- mgcv::gam
    }
    
    params
  },
  
  compute_group = function(data, scales, method = "auto", formula = y ~ x, se = TRUE, n = 80, span = 0.75, fullrange = FALSE, xseq = NULL, level = 0.95, method.args = list(), na.rm = FALSE, xpos = NULL, ypos = NULL) {
    if (length(unique(data$x)) < 2) {
      # Not enough data to perform fit
      return(data.frame())
    }
    
    if (is.null(data$weight))
      data$weight <- 1
    
    if (is.null(xseq)) {
      if (is.integer(data$x)) {
        if (fullrange) {
          xseq <- scales$x$dimension()
        } else {
          xseq <- sort(unique(data$x))
        }
      } else {
        if (fullrange) {
          range <- scales$x$dimension()
        } else {
          range <- range(data$x, na.rm = TRUE)
        }
        xseq <-
          seq(range[1], range[2], length.out = n)
      }
    }
    # Special case span because it's the most commonly used model argument
    if (identical(method, "loess")) {
      method.args$span <- span
    }
    
    if (is.character(method))
      method <- match.fun(method)
    
    base.args <-
      list(quote(formula),
           data = quote(data),
           weights = quote(weight))
    model <-
      do.call(method, c(base.args, method.args))
    
    m = model
    eq <-
      substitute(
        italic(y) == a + b %.% italic(x) * "," ~  ~ italic(R) ~ "=" ~ r,
        list(
          a = format(coef(m)[[1]], digits = 3),
          b = format(coef(m)[[2]], digits = 3),
          r = format(sqrt(summary(m)$r.squared), digits = 3)
        )
      )
    func_string = as.character(as.expression(eq))
    
    if (is.null(xpos))
      xpos = min(data$x) * 0.9
    if (is.null(ypos))
      ypos = max(data$y) * 0.9
    data.frame(x = xpos, y = ypos, label = func_string)
    
  },
  
  required_aes = c("x", "y")
)

# >>Confusion heatmap new----
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

# >>Boxplot + heatmap----
BoxHeatmap <- function(expression.matrix, groupid){
  # Boxplot in the left panel
  signature.mean.expression <- colMeans(expression.matrix)
  boxplot.data <- data.frame(Expression = signature.mean.expression, Group = groupid)
  p_boxplot <- 
    ggplot(boxplot.data, aes(x = Group, y = Expression)) + 
    geom_boxplot(color = "#1979B5", fill = "#1979B5", 
                 outlier.colour = "black", outlier.shape = 1,
                 outlier.stroke = .25, outlier.alpha = .5) +
    theme_bw() +
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
  gene.median.expression[, 2:(nrow(expression.matrix) + 1)] <- apply(gene.median.expression[, 2:(nrow(expression.matrix) + 1)], 2, function(x){
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
  return(p)
}
