#' Define the number of composition clusters
#'
#' This function iterates through defined number of clusters and plot the total withinss of kmeans clustering
#'
#' @param Celltype_deconvolved_probs probability matrix of deconvoluted spatial spots. Rows corresponds to spatial spots, columns corresponds to cell types. can be calculated with or without single-cell reference
#' @param max_k Maximum number of expected clusters to iterate
#' @return an elbow plot of the k values to decide for the best k
#' @export

Composition.cluster.k <- function(Celltype_deconvolved_probs, max_k){
  sample.wss <- NULL
  for (i in 1:max_k){
    sample.fit=kmeans(Celltype_deconvolved_probs, centers = i)
    sample.wss <- c(sample.wss, sample.fit$tot.withinss)
  }
  plot(1:max_k, sample.wss, type = "o")
}


#' Clustering spatial spots based on their deconvoluted celltype composition
#'
#' This function clusters spatial spots based on their deconvoluted celltype composition
#' @param spatial.object processed seurat object of the spatial sample with all annotation
#' @param Celltype_deconvolved_probs probability matrix of deconvoluted spatial spots. Rows corresponds to spatial spots, columns corresponds to cell types. can be calculated with or without single-cell reference
#' @param k number of expected clusters
#' @return a seurat spatial object with an added column to metadata called CompositionCluster_CC
#' @export

Composition.cluster <- function(spatial.object, Celltype_deconvolved_probs, k){
  sample.spot.clusters <- kmeans(Celltype_deconvolved_probs, k)
  sample.spot.clusters$cluster <- paste("CC", sample.spot.clusters$cluster, sep = "")
  spatial.object <- AddMetaData(spatial.object, sample.spot.clusters$cluster, "CompositionCluster_CC")
  return(spatial.object)
}



#' Umap of spatial spots based on their cell type composition profile
#'
#' This function provides umap coordinates and plots of spatial spots based on celltype composition profiles
#' @param spatial.object processed seurat object of the spatial sample with all annotation
#' @param Celltype_deconvolved_probs probability matrix of deconvoluted spatial spots. Rows corresponds to spatial spots, columns corresponds to cell types. can be calculated with or without single-cell reference
#' @return a list containing UMAP coordinates of spots, ggplot of umap color coded by comosition clusters, ggplot of umap with deconvolution pie charts
#' @export

Composition_cluster_umap <- function(spatial.object, Celltype_deconvolved_probs){
  sample.spotDecomp.umap <- umap(Celltype_deconvolved_probs,labels=as.factor(spatial.object@meta.data[rownames(Celltype_deconvolved_probs), "CompositionCluster_CC"]),
                                 dotsize = 1)
  sample.spotDecomp.umap.mat <- sample.spotDecomp.umap$layout; colnames(sample.spotDecomp.umap.mat) <- c("x", "y")
  sample.spotDecomp.umap.mat <- cbind(as.data.frame(sample.spotDecomp.umap.mat), CompositionCluster_CC=as.factor(spatial.object@meta.data[rownames(Celltype_deconvolved_probs), "CompositionCluster_CC"]))
  sample.spotDecomp.umap.mat <- cbind(sample.spotDecomp.umap.mat, Slide=as.factor(spatial.object@meta.data[rownames(Celltype_deconvolved_probs), "orig.ident"]))
  sample.spotDecomp.umap.mat <- as.data.frame(sample.spotDecomp.umap.mat)
  umap.cluster.gg <- ggplot(sample.spotDecomp.umap.mat, aes(x=x, y=y, color=CompositionCluster_CC)) + geom_point() +   theme_bw() +theme (
      axis.text.x = element_text(angle = 30, vjust = 0.7, size = 15, face = "bold", colour = "black"),
      axis.text.y = element_text( hjust = 1, size = 20, face = "bold", colour = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.text=element_text(size=15))
  sample.pos.umap <- sample.spotDecomp.umap.mat[, c("x", "y")]
  sample.pos.umap <- sample.pos.umap*100
  umap.deconv.gg <- spot.piecharts(Celltype_deconvolved_probs, sample.pos.umap,r=10, lwd = 0.05)
  output.list <- list(umap.table=sample.spotDecomp.umap.mat, umap.cluster.gg=umap.cluster.gg, umap.deconv.gg=umap.deconv.gg)
}



#' Defining Condition-specific and Condition-common composition clusters
#'
#' This function provides Defining Condition-specific and Condition-common composition clusters
#' @param spatial.object processed seurat object of the spatial sample with all annotation
#' @param Cluster_Column Column from Metadata that contains Composition cluster labels
#' @param Condition Column from Metadata that contains Condition labels
#' @return a list containing Condition-specific composition clusters, Condition-common composition clusters and ggplot of the percentages
#' @export

Composition_cluster_group_pct <- function(Spatial.object, Cluster_Column, Condition){
  sample.cond.compclust.table <- data.frame(unclass(table(Spatial.object@meta.data[,Cluster_Column], Spatial.object@meta.data[,Condition])))
  sample.cond.compclust.table.pct <- round(prop.table(as.matrix(sample.cond.compclust.table),1),2)
  #colnames(sample.cond.compclust.table.pct) <- paste(colnames(sample.cond.compclust.table.pct), "pct", sep = "_")
  rownames(sample.cond.compclust.table.pct) <- paste("CompCluster", rownames(sample.cond.compclust.table.pct),sep = "_")
  #sample.CondSpecific.compClusts <- c()
  #sample.CondComm.compClusts <- c()
  #for (sample.clust in rownames(sample.cond.compclust.table)){
  #  sample.ratio <- sample.cond.compclust.table[sample.clust, "pct.g1"]
  #  if (sample.ratio >= 75 | sample.ratio <= 25) {sample.CondSpecific.compClusts <- c(sample.CondSpecific.compClusts, sample.clust) }
  #  else { sample.CondComm.compClusts <- c(sample.CondComm.compClusts, sample.clust)}
  #}

  sample.cond.compclust.table.melt <- reshape2::melt(sample.cond.compclust.table.pct)
  sample.cond.compclust.table.melt$Var1 <- factor(sample.cond.compclust.table.melt$Var1)
  plot <- ggplot(sample.cond.compclust.table.melt, aes(x=Var2, y=Var1, fill= value)) + geom_tile()+
    scale_fill_viridis(discrete=FALSE) + theme_classic(base_size = 17) + theme (
      axis.text.x = element_text(angle = 30, vjust = 0.7, size = 17, face = "bold", colour = "black"),
      axis.text.y = element_text( hjust = 1, size = 17, face = "bold", colour = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank())
  #output.list <- list(CondSpecific= sample.CondSpecific.compClusts, CondComm = sample.CondComm.compClusts, plot = plot)
  return(plot)
}


#' Defining driver celltypes of each composition cluster
#'
#' This function Defining driver celltypes of each composition cluster
#' @param spatial.object processed seurat object of the spatial sample with all annotation
#' @param COI Composition Cluster of Interest
#' @param Celltype_deconvolved_probs probability matrix of deconvoluted spatial spots. Rows corresponds to spatial spots, columns corresponds to cell types. can be calculated with or without single-cell reference
#' @return a boxplot of cell type probabilities for cluster of interest
#' @export


Composition_cluster_enrichedCelltypes <- function(Spatial.object, COI, Celltype_deconvolved_probs){
  COI.spots <- rownames(Spatial.object@meta.data[which(Spatial.object@meta.data$CompositionCluster_CC == COI),])
  COI.topic.probs <- Celltype_deconvolved_probs[COI.spots,]
  COI.topic.probs.medians <- sort(colMedians(COI.topic.probs),decreasing = T)
  COI.leading.topics <- names(COI.topic.probs.medians[1:2])
  COI.topic.probs.melt <- melt(COI.topic.probs); colnames(COI.topic.probs.melt) <- c("spot_id", "Celltype", "prob")
  Topic.prob.all.df <- COI.topic.probs.melt
  Topic.prob.all.df$Topic <- factor(Topic.prob.all.df$Celltype)
  p <- ggplot(Topic.prob.all.df, aes(x=Topic, y=prob, fill=Topic)) + geom_boxplot(show.legend = FALSE) + theme_bw() + theme (
    axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1,size = 10, face = "bold", colour = "black"),
    axis.text.y = element_text( hjust = 1, size = 17, face = "bold", colour = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) + ggtitle(paste("Compositionn Cluster", COI, sep = " "))
  return(p)
}




spot.piecharts <- function (theta, pos, topicOrder = seq(ncol(theta)), topicCols = rainbow(ncol(theta)),
          groups = NA, group_cols = NA, r = max(0.4, max(pos)/nrow(pos) *
                                                  4), lwd = 0.5, showLegend = TRUE, plotTitle = NA, overlay = NA)
{
  if (!is.matrix(theta) & !is.data.frame(theta)) {
    stop("`theta` must be a matrix or data.frame.")
  }
  if (!is.matrix(pos) & !is.data.frame(pos)) {
    stop("`pos` must be a matrix or data.frame with exactly 2 columns named `x` and `y`.")
  }
  if ((any(!colnames(pos) %in% c("x", "y")) == TRUE) | (dim(pos)[2] !=
                                                        2)) {
    stop("`pos` must have exactly 2 columns named `x` and `y`.")
  }
  theta_ordered <- theta[, topicOrder]
  theta_ordered <- as.data.frame(theta_ordered)
  pixels <- intersect(rownames(theta_ordered), rownames(pos))
  pixels <- rownames(theta_ordered)[which(rownames(theta_ordered) %in%
                                            pixels)]
  theta_ordered_pos <- merge(data.frame(theta_ordered), data.frame(pos),
                             by = 0)
  rownames(theta_ordered_pos) <- theta_ordered_pos[, "Row.names"]
  theta_ordered_pos <- theta_ordered_pos[pixels, ]
  topicColumns <- colnames(theta_ordered_pos)[2:(dim(theta_ordered_pos)[2] -
                                                   2)]
  if (is.na(groups[1])) {
    groups <- rep("0", dim(theta_ordered_pos)[1])
    theta_ordered_pos$Pixel.Groups <- groups
  }
  else {
    theta_ordered_pos$Pixel.Groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
    group_cols <- c(`0` = "gray")
  }
  message("Plotting scatterpies for ", dim(theta_ordered_pos)[1],
          " pixels with ", length(topicColumns), " cell-types...this could take a while if the dataset is large.",
          "\n")
  if (!is.na(overlay[1])) {
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = 0:dim(overlay)[2],
                                                y = 0:dim(overlay)[1])) + ggplot2::coord_equal(xlim = c(0,
                                                                                                        dim(overlay)[2]), ylim = c(0, dim(overlay)[1]), expand = FALSE) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.line = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 12,
                                                                                                     colour = "black"), legend.title = ggplot2::element_text(size = 12,
                                                                                                                                                             colour = "black")) + ggplot2::annotation_raster(overlay,
                                                                                                                                                                                                             xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      scatterpie::geom_scatterpie(ggplot2::aes(x = x, y = y,
                                               group = Row.names, r = r, color = Pixel.Groups),
                                  lwd = lwd, data = theta_ordered_pos, cols = topicColumns,
                                  legend_name = "Topics") + ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  }
  else {
    p <- ggplot2::ggplot() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                            axis.line = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                                            axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                            axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                                            panel.background = ggplot2::element_blank(), plot.background = ggplot2::element_blank(),
                                            legend.text = ggplot2::element_text(size = 12, colour = "black"),
                                            legend.title = ggplot2::element_text(size = 12, colour = "black")) +
      scatterpie::geom_scatterpie(ggplot2::aes(x = x, y = y,
                                               group = Row.names, r = r, color = Pixel.Groups),
                                  lwd = lwd, data = theta_ordered_pos, cols = topicColumns,
                                  legend_name = "CellTypes") + ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  }
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  p <- p + ggplot2::coord_equal()
  return(p)
}
