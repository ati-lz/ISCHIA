#' Differentially activated cell-type network between two conditions
#'
#' This function calculates Differentially activated cell type network between two conditions
#'
#' @param spatial.object processed seurat object of the spatial sample with all annotation
#' @param Celltype_deconvolved_probs probability matrix of deconvoluted spatial spots. Rows corresponds to spatial spots, columns corresponds to cell types. can be calculated with or without single-cell reference
#' @param COI1 First Composition cluster of interest for the differential celltype interaction network analysis
#' @param COI2 Second Composition cluster of interest for the differential celltype interaction network analysis
#' @return a list containing the estimated difference of condition specific precision matrices. celltype pairs diff activity scores, The estimated differential network over only the connected nodes
#' @export


DiffNet.Celltype <- function(spatial.object, Celltype_deconvolved_probs, COI1, COI2){
  Celltype_deconvolved_probs <- as.matrix(Celltype_deconvolved_probs)
  COI1.spots <- rownames(spatial.object@meta.data[which(spatial.object@meta.data$TopicComposition_cluster == COI1),])
  COI2.spots <- rownames(spatial.object@meta.data[which(spatial.object@meta.data$TopicComposition_cluster == COI2),])

  Celltype_deconvolved_probs.COI1 <- Celltype_deconvolved_probs[COI1.spots,]
  Celltype_deconvolved_probs.COI2 <- Celltype_deconvolved_probs[COI2.spots,]

  input.mat.list.celltypes <- list(Celltype_deconvolved_probs.COI1, Celltype_deconvolved_probs.COI2)
  print("Dtrace running...")
  dtrace.results.celltypes= Dtrace(input.mat.list.celltypes, 0.45, covType = "spearman")
  print("Dtrace Finished")
  net.dtrace.celltypes = dtrace.results.celltypes$Delta.graph.connected
  tkid.celltypes <- tkplot(net.dtrace.celltypes, vertex.size= degree(net.dtrace.celltypes)*1.5, layout =layout_with_fr,
                           vertex.color="red", vertex.label.cex=0.8, edge.width =1.5, edge.color="orange")

  dtrace.results.celltypes.Delta <- dtrace.results.celltypes$Delta
  Delta.graph = (abs(dtrace.results.celltypes$Delta) > 1e-05) * 1
  diag(dtrace.results.celltypes.Delta) <- 0
  dtrace.results.celltypes.nonzero.ind <- which(dtrace.results.celltypes.Delta!=0,arr.ind = T)
  dtrace.results.nonzeroz <- data.frame(pair=NA, value=0)
  for (i in 1:nrow(dtrace.results.celltypes.Delta)){
    for (j in 1:ncol(dtrace.results.celltypes.Delta)){
      if(dtrace.results.celltypes.Delta[i,j] != 0){
        val <- dtrace.results.celltypes.Delta[i,j]
        rowname <- rownames(dtrace.results.celltypes.Delta)[i]
        colname <- colnames(dtrace.results.celltypes.Delta)[j]
        entry.line <- data.frame(pair= paste(rowname, colname, sep = "_"), value=val)
        dtrace.results.nonzeroz <- rbind(dtrace.results.nonzeroz, entry.line)}}}
  celltypes.diff.scores <- dtrace.results.nonzeroz[-1,]
  celltypes.diff.scores <- celltypes.diff.scores[order(celltypes.diff.scores$value, decreasing = T),]

  final.output <- list(dtrace.results.celltypes.Delta, celltypes.diff.scores, net.dtrace.celltypes)
  return(final.output)
}



#' Differentially activated ligand-receptor network between two conditions
#'
#' This function calculates Differentially activated ligand-receptor network between two conditions
#'
#' @param spatial.object processed seurat object of the spatial sample with all annotation
#' @param LR.list list of interested ligands and receptors to be explored. Common genes with spatial expression matrix
#' @param COI1 First Composition cluster of interest for the differential activated l-R network analysis
#' @param Cond1 The consition that COI1 spots would be subsetted from.
#' @param COI2 Second Composition cluster of interest for the differential activated l-R network analysis
#' @param Cond2 The consition that COI2 spots would be subsetted from.
#' @return a list containing the estimated difference of condition specific precision matrices. L-R pairs diff activity scores, The estimated differential network over only the connected nodes
#' @export

DiffNet.LR <- function(spatial.object, LR.list, COI1, Cond1, COI2, Cond2){
  spatial.object.exp <- as.matrix(spatial.object@assays$Spatial@counts)
  spatial.object.exp.LR.subset.raw <- spatial.object.exp[LR.list,]
  spatial.object.exp.LR.subset.raw[which(spatial.object.exp.LR.subset.raw[,] != 0)]
  spatial.object.exp.LR.subset.raw <- spatial.object.exp.LR.subset.raw[which(rowSums(spatial.object.exp.LR.subset.raw) > 0),]
  spatial.object.exp.LR.subset.raw <- spatial.object.exp.LR.subset.raw[,which(colSums(spatial.object.exp.LR.subset.raw) > 0)]
  LR.spatial.exp.mat <- spatial.object.exp.LR.subset.raw

  LR.spatial.exp.mat.t <- as.matrix(t(LR.spatial.exp.mat))
  COI1.spots <- rownames(spatial.object@meta.data[which(spatial.object@meta.data$TopicComposition_cluster %in% COI1 & spatial.object@meta.data$orig.ident %in% Cond1),])
  COI2.spots <- rownames(spatial.object@meta.data[which(spatial.object@meta.data$TopicComposition_cluster %in% COI2 & spatial.object@meta.data$orig.ident %in% Cond2),])

  LR.spatial.exp.mat.t.COI1 <- LR.spatial.exp.mat.t[COI1.spots,]
  LR.spatial.exp.mat.t.COI2 <- LR.spatial.exp.mat.t[COI2.spots,]

  input.mat.list <- list(LR.spatial.exp.mat.t.COI1, LR.spatial.exp.mat.t.COI2)
  print("Dtrace running...")
  dtrace.results.LR= Dtrace(input.mat.list, 0.45, covType = "spearman")
  print("Dtrace Finished")
  net.dtrace.LR = dtrace.results.LR$Delta.graph.connected
  tkid.LR <- tkplot(net.dtrace.LR, vertex.size= degree(net.dtrace.LR)*1.5, layout =layout_with_fr,
                    vertex.color="red", vertex.label.cex=0.8, edge.width =1.5, edge.color="orange")

  dtrace.results.LR.Delta <- dtrace.results.LR$Delta
  Delta.graph = (abs(dtrace.results.LR$Delta) > 1e-05) * 1
  diag(dtrace.results.LR.Delta) <- 0
  dtrace.results.LR.nonzero.ind <- which(dtrace.results.LR.Delta!=0,arr.ind = T)
  dtrace.results.nonzeroz <- data.frame(pair=NA, value=0)
  for (i in 1:nrow(dtrace.results.LR.Delta)){
    for (j in 1:ncol(dtrace.results.LR.Delta)){
      if(dtrace.results.LR.Delta[i,j] != 0){
        val <- dtrace.results.LR.Delta[i,j]
        rowname <- rownames(dtrace.results.LR.Delta)[i]
        colname <- colnames(dtrace.results.LR.Delta)[j]
        entry.line <- data.frame(pair= paste(rowname, colname, sep = "_"), value=val)
        dtrace.results.nonzeroz <- rbind(dtrace.results.nonzeroz, entry.line)}}}
  LR.diff.scores <- dtrace.results.nonzeroz[-1,]
  LR.diff.scores <- LR.diff.scores[order(LR.diff.scores$value, decreasing = T),]

  final.output <- list(dtrace.results.LR.Delta, LR.diff.scores, net.dtrace.LR)
  return(final.output)
}

