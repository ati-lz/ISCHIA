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
  umap.deconv.gg <- vizAllTopics(Celltype_deconvolved_probs, sample.pos.umap,r=10, lwd = 0.05,
                                 topicCols = c("red", "yellowgreen", "salmon","coral4","blue", "green", "yellow", "purple", "pink", "black", "cyan", "darkcyan"))
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

