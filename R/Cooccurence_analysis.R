#' Calculate sigificantly cooccuring celltypes in composition clusters of interest
#'
#' This function calculates the significant cooccuring cell types based on deconvoluted celltype probability matrix
#'
#' @param spatial.object processed seurat object of the spatial sample with all annotation
#' @param deconv.prob.mat probability matrix of deconvoluted spatial spots. Rows corresponds to spatial spots, columns corresponds to cell types. can be calculated with or without single-cell reference
#' @param COI Composition cluster of interest, a vector of Cluster IDs, of spatial spots which we would like to calculate significant cell type cooccurences in
#' @param Condition Condition of interest,a vector of Conditions, of spatial spots which we would like to calculate significant cell type cooccurences in
#' @param prob.th Probability threshold to convert the deconvolution probability matrix to binary presence/absence matrix
#' @return a list containing the results table of cell type coocuurences
#' @export
spatial.celltype.cooccurence <- function(spatial.object, deconv.prob.mat, COI, Condition, prob.th){
  COI.spots <- rownames(spatial.object@meta.data[which(spatial.object@meta.data$TopicComposition_cluster %in% COI & spatial.object@meta.data$orig.ident %in% Condition),])
  coocur.COI.exp <- t(deconv.prob.mat[COI.spots,])
  coocur.COI.exp <- biclust::binarize(coocur.COI.exp, threshold=prob.th)
  cooccur.COI.res <- cooccur(mat = coocur.COI.exp, type = "spp_site", spp_names = TRUE, prob = "comb", thresh = FALSE)
  summary(cooccur.COI.res)
  prob.table(cooccur.COI.res)
  cooccur.COI.res$results <- cbind(cooccur.COI.res$results, pair = paste(cooccur.COI.res$results$sp1_name, cooccur.COI.res$results$sp2_name,sep = "_"))
  rownames(cooccur.COI.res$results) <- cooccur.COI.res$results$pair
  p1 <- plot(cooccur.COI.res)
  return(cooccur.COI.res)
}



#' Calculate sigificantly cooccuring Ligand-receptor pairs in composition clusters of interest
#'
#' This function calculates the significant cooccuring Ligand-Receprot pairs based on Ligand-receptor expression matrix
#'
#' @param spatial.object processed seurat object of the spatial sample with all annotation
#' @param COI vector including Composition clusters of interest of spatial spots which we would like to calculate significant Ligand-receptor cooccurences in
#' @param Condition vector including conditions from orig.ident to select COI from
#' @param LR.list list of interested ligands and receptors to be explored. Common genes with spatial expression matrix
#' @param LR.pairs list of interested ligands and receptors pairs, separated by Underline
#' @param exp.th expression threshold for raw counts of ligand-receptors to convert gene expression matrix to binary presence/absence matrix
#' @param corr.th correlation threshold between ligand and receptor pair expression
#' @return a list of significantly cooccuring and correlating ligand annd receptors in composition cluster of interest and result table of cooccurence analysis
#' @export

Enriched.LRs <- function(spatial.object, COI, Condition, LR.list, LR.pairs, exp.th, corr.th){
  print("preparing L-R presence/absence matrix")
  #Idents(spatial.object) <- "TopicComposition_cluster"
  #COI.DEGs <- FindMarkers(spatial.object, ident.1 = COI, min.pct = 0.1, only.pos = T)
  #LR.list.DE <- intersect(LR.list, rownames(COI.DEGs))

  spatial.object.exp <- spatial.object@assays$Spatial@counts
  spatial.object.exp.norm <- spatial.object@assays$Spatial@data
  spatial.obj.exp.LR.subset.raw <- spatial.object.exp[LR.list,]
  #spatial.obj.exp.LR.subset.raw <- spatial.object.exp[LR.list.DE,]
  spatial.obj.exp.LR.subset.raw.binary <- as.matrix(binarize(spatial.obj.exp.LR.subset.raw, threshold = exp.th))
  spatial.obj.exp.LR.subset.raw.binary <- spatial.obj.exp.LR.subset.raw.binary[which(rowSums(spatial.obj.exp.LR.subset.raw.binary) > 0),]
  spatial.obj.exp.LR.subset.raw.binary <- spatial.obj.exp.LR.subset.raw.binary[,which(colSums(spatial.obj.exp.LR.subset.raw.binary) > 0)]
  LR.presence.absence.mat <- spatial.obj.exp.LR.subset.raw.binary

  COI.spots <- names(spatial.object$TopicComposition_cluster[which(spatial.object$TopicComposition_cluster %in% COI & spatial.object$orig.ident %in% Condition)])
  rest.of.spots <- setdiff(rownames(spatial.object@meta.data), COI.spots)
  #Correlation
  print("Calculating L-R pairs correlation")
  COI.cors <- cor(t(as.matrix(spatial.object.exp[LR.list,COI.spots])))
  #COI.cors <- cor(t(as.matrix(spatial.object.exp[LR.list.DE,COI.spots])))
  COI.cors[which(is.na(COI.cors))] <- 0

  #Coocurance
  print("preparing for cooccurence")
  common.spots <- intersect(colnames(LR.presence.absence.mat), COI.spots)
  coocur.COI.exp <- as.data.frame(as.matrix(LR.presence.absence.mat[,common.spots]))
  print("cooccurence calculation starts...")
  cooccur.COI.res <- cooccur(mat = coocur.COI.exp, type = "spp_site", thresh = TRUE, spp_names = TRUE)
  print("cooccurence calculation Ended")
  #save(cooccur.COI.res, file = "/Users/ati/Documents/Projects/Visium/Visium_IBD/HumanIBD_P4_cooccur_COI5_res_2.RData")
  summary(cooccur.COI.res)
  prob.table(cooccur.COI.res)
  cooccur.COI.res$results <- cbind(cooccur.COI.res$results, pair = paste(cooccur.COI.res$results$sp1_name, cooccur.COI.res$results$sp2_name,sep = "_"))
  rownames(cooccur.COI.res$results) <- cooccur.COI.res$results$pair
  plot(cooccur.COI.res)
  common.pairs <- intersect(LR.pairs, rownames(cooccur.COI.res$results))

  #enriched LRs
  COI.enrcihed.LRs <- data.frame(from=NA, to=NA, correlation=NA, ligand_FC= NA, Receptor_FC = NA)#, Freq.pct.1 =NA, Freq.pct.2= NA)
  pair.count <- 0
  for (pair in common.pairs){
    pair.count <- pair.count+1
    print(paste(pair.count, length(common.pairs), sep = "_"))
    LR.pair.words <- unlist(strsplit(pair, split = "_"))
    LR.pair.ligand <- LR.pair.words[1]
    LR.pair.Receptor <- LR.pair.words[2]
    ligand.exp.COI.mean <- mean(spatial.object.exp.norm[LR.pair.ligand, COI.spots])
    ligand.exp.otherspots.mean <- mean(spatial.object.exp.norm[LR.pair.ligand, rest.of.spots])
    ligand.FC <- round(ligand.exp.COI.mean/ligand.exp.otherspots.mean,4)
    Receptor.exp.COI.mean <- mean(spatial.object.exp.norm[LR.pair.Receptor, COI.spots])
    Receptor.exp.otherspots.mean <- mean(spatial.object.exp.norm[LR.pair.Receptor, rest.of.spots])
    Receptor.FC <- round(Receptor.exp.COI.mean/Receptor.exp.otherspots.mean,4)
    #if (LR.pair.ligand %in% rownames(P4.merge.LRs) & LR.pair.Receptor%in% rownames(P4.merge.LRs)){
    #if (LR.logFCs[LR.pair.ligand, "avg_log2FC"] > 0.5 & LR.logFCs[LR.pair.Receptor, "avg_log2FC"] > 0.5 & COI.cors[LR.pair.ligand, LR.pair.Receptor] > 0.5){
    #if ( COI.LR.pair.group.freq.pct[pair] > 50 & Villus.LR.pair.group.freq.pct[pair] < 50){
    if ( cooccur.COI.res$results[pair, "p_gt"] < 0.05 & COI.cors[LR.pair.ligand, LR.pair.Receptor] > corr.th){
      added.row <- data.frame(from=LR.pair.ligand, to=LR.pair.Receptor, correlation=COI.cors[LR.pair.ligand, LR.pair.Receptor], ligand_FC= ligand.FC, Receptor_FC = Receptor.FC)#, Freq.pct.1 =COI.LR.pair.group.freq.pct[pair], Freq.pct.2= Villus.LR.pair.group.freq.pct[pair])
      COI.enrcihed.LRs <- rbind(COI.enrcihed.LRs, added.row)
    }
  }
  COI.enrcihed.LRs <- COI.enrcihed.LRs[-1,]
  COI.enrcihed.LRs <- COI.enrcihed.LRs[order(COI.enrcihed.LRs$correlation, decreasing = T),]
  COI.enrcihed.LRs <- cbind(COI.enrcihed.LRs, pair= paste(COI.enrcihed.LRs$from, COI.enrcihed.LRs$to, sep = "_"))
  Output.list <- list(COI.enrcihed.LRs=COI.enrcihed.LRs, cooccurence.table=cooccur.COI.res)
  return(Output.list)
}


#' Find LRs that are significantly cooccuring in one group and not in the other group
#'
#' This function Find LRs that are significantly cooccuring in one group and not in the other group
#'
#' @param Enriched.LRs.G1 Output of Enriched.LRs function for Group1
#' @param Enriched.LRs.G2 Output of Enriched.LRs function for Group2
#' @param G1.max.pth Max pvalue threshold for significancy levels of cooccuring LRs in Group1
#' @param G2.min.pth Min pvalue threshold for non-significancy levels of cooccuring LRs in Group2
#' @return List of LRs enriched in Group1 and not in Group2
#' @export

Diff.cooc.LRs <- function(Enriched.LRs.G1, Enriched.LRs.G2, G1.max.pth, G2.min.pth){
  G1.cooc.res <- Enriched.LRs.G1$cooccurence.table$results
  rownames(G1.cooc.res) <- G1.cooc.res$pair
  G2.cooc.res <- Enriched.LRs.G2$cooccurence.table$results
  rownames(G2.cooc.res) <- G2.cooc.res$pair

  G1.enriched.LRs <- data.frame(pair = NA, G1_pval= NA, G2_pval = NA, pval_dif =NA, obs_cooc = NA)
  for (lr_pair in rownames(G1.cooc.res)){
    if(lr_pair %in% rownames(G2.cooc.res)){
      G1.cooc.res.pgt <- as.numeric(G1.cooc.res[lr_pair,"p_gt"])
      G2.cooc.res.pgt <- as.numeric(G2.cooc.res[lr_pair,"p_gt"])
      G1.cooc.obs <- as.numeric(G1.cooc.res[lr_pair,"obs_cooccur"])
      G2.cooc.obs <- as.numeric(G2.cooc.res[lr_pair,"obs_cooccur"])
      G1.cooc.exp <- as.numeric(G1.cooc.res[lr_pair,"exp_cooccur"])
      G2.cooc.exp <- as.numeric(G2.cooc.res[lr_pair,"exp_cooccur"])
      pval.diff <- G2.cooc.res.pgt-G1.cooc.res.pgt
      if (G1.cooc.res.pgt < G1.max.pth & G2.cooc.res.pgt > G2.min.pth & G1.cooc.obs > 10 & G1.cooc.obs!= G1.cooc.exp & G2.cooc.obs!=G2.cooc.exp & G2.cooc.obs < G1.cooc.obs){
        print(lr_pair)
        out <- c(pair = lr_pair, G1_pval= G1.cooc.res.pgt, G2_pval = G2.cooc.res.pgt, pval_dif= pval.diff, obs_cooc=G1.cooc.obs)
        G1.enriched.LRs <- rbind(G1.enriched.LRs, out)
      }
    }
  }
  G1.enriched.LRs <- G1.enriched.LRs[-1,]
  rownames(G1.enriched.LRs) <- G1.enriched.LRs$pair
  G1.enriched.LRs.sorted <- G1.enriched.LRs[order(G1.enriched.LRs[,"obs_cooc"], decreasing = T),]
  G1.enriched.LRs.sorted <- G1.enriched.LRs.sorted[,-1]
  #G1.enriched.LRs.sorted <- G1.enriched.LRs[order(G1.enriched.LRs[,"obs_cooc"], decreasing = T),]

  return(G1.enriched.LRs.sorted)
}


#' Plot Enriched LRs in
#'
#' This function plots Enriched LR interactions in the data frame provided by you
#'
#' @param selected.LRs Selected list of LR interactions
#' @return Chord Diagram of LR interactions
#' @export

ChordPlot.Enriched.LRs <- function(selected.LRs){
  adjacencyData <- with(selected.LRs, table(from, to))
  p <- chordDiagram(adjacencyData, transparency = 0.5,annotationTrack = c("name","grid"))
  return(p)
}


#' Plot group specific LRs
#'
#' This function Plots group specific LRs
#' @param Group1.enrcihed.LRs data.frame of selected LRs with "pair" column
#' @param Group2.enrcihed.LRs data.frame of selected LRs with "pair" column
#' @return Chord Diagram of LR interactions
#' @export

SankeyPlot.Diff.LRs <- function(Group1.enrcihed.LRs,Group2.enrcihed.LRs){
  connection.df.all.pre <- rbind(Group1.enrcihed.LRs, Group2.enrcihed.LRs)
  connection.df.all.pre.dups <- duplicated(connection.df.all.pre$pair)
  connection.df.all <- data.frame(from=NA, to=NA, correlation=NA,ligand_FC=NA, Receptor_FC=NA, pair=NA,group=NA)
  for(i in 1:nrow(connection.df.all.pre)){
    print(i)
    if(connection.df.all.pre[i, "pair"] %in% Group1.enrcihed.LRs$pair & connection.df.all.pre[i, "pair"] %in% Group2.enrcihed.LRs$pair){
      add.row <- cbind(connection.df.all.pre[i,], group="Common")
      connection.df.all[i,] <- add.row}
    else if(connection.df.all.pre[i, "pair"] %in% Group1.enrcihed.LRs$pair){
      add.row <- cbind(connection.df.all.pre[i,], group="Condition1")
      connection.df.all[i,] <- add.row}
    else if(connection.df.all.pre[i, "pair"] %in% Group2.enrcihed.LRs$pair){
      add.row <- cbind(connection.df.all.pre[i,], group="Condition2")
      connection.df.all[i,] <- add.row}
  }
  connection.df.all <- na.omit(connection.df.all)
  connection.df.all$group <- factor(connection.df.all$group)

  #Visualizing enriched LRs
  links <- connection.df.all
  #links$group2<- c(rep("Con1", 20), rep("Con2",10), rep("Common", 3))
  nodes <- data.frame(name=c(as.character(links$from), as.character(links$to)) %>% unique())
  links$IDsource <- match(links$from, nodes$name)-1
  links$IDtarget <- match(links$to, nodes$name)-1
  nodes$group <- as.factor(c("my_unique_group"))
  my_color <- 'd3.scaleOrdinal() .domain(["Condition1", "Condition2", "Common", "my_unique_group"]) .range(["grey", "red","#69b3a2", "grey"])'
  p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
                Value = "correlation", NodeID = "name",
                colourScale=my_color, LinkGroup="group",fontSize= 15)
  return(p)
}


#' Find differentially expressed genes between Spatial spots Double positive and Double negative for an LR pair
#'
#' This function Finds differentially expressed genes between Spatial spots Double positive and Double negative for an LR pair
#' @param SeuratObj Spatial Seurat Object
#' @param Enriched.LRs.df data.frame of selected LR, output of Diff.cooc.LRs function
#' @param COI Composition Cluster of Interest
#' @param Condition vector specifying conditions of interest from orig.ident metadata column
#' @return list of Differentially expressed genes between Spatial spots Double positive and Double negative for LRs in Enriched.LRs.df
#' @export

Binary.LR.Diffexp <- function(SeuratObj, Enriched.LRs.df, COI, Condition){
  output.list <- list()
  for(LR.pair in rownames(Enriched.LRs.df)){
    print(LR.pair)
    ligand <- unlist(strsplit(LR.pair, split = "_"))[1]
    receptor <- unlist(strsplit(LR.pair, split = "_"))[2]
    count.mat <- SeuratObj@assays$Spatial@data
    COI.cells <- rownames(SeuratObj@meta.data[which(SeuratObj@meta.data$TopicComposition_cluster %in% COI & SeuratObj@meta.data$orig.ident %in% Condition),])
    count.mat.subset <- t(as.matrix(count.mat[c(ligand,receptor),COI.cells]))
    LR.pair.exp.vec <- rep("None", nrow(SeuratObj@meta.data))
    names(LR.pair.exp.vec) <- rownames(SeuratObj@meta.data)
    for (cell in rownames(count.mat.subset)){
      if(sum(count.mat.subset[cell,] !=0) == 2){
        LR.pair.exp.vec[cell] <- paste(COI,"DoublePositive", sep = "_")
      }
      else{
        LR.pair.exp.vec[cell] <- paste(COI,"DoubleNegative", sep = "_")
      }
    }
    SeuratObj <- AddMetaData(SeuratObj, LR.pair.exp.vec, "LR_pair_exp_vec")
    Idents(SeuratObj) <- "LR_pair_exp_vec"
    COI.LR.pos.neg.diffexp <- FindMarkers(SeuratObj, ident.1 = c(paste(COI,"DoublePositive", sep = "_")), ident.2 = c(paste(COI,"DoubleNegative", sep = "_")), only.pos = T)
    output.list[[LR.pair]] <- COI.LR.pos.neg.diffexp
  }
  return(output.list)
}



#' Function for producing a heatmap co-occurrence visualization.
#'
#' Heatmap visualization of the pairwise species associations revealed by a cooccur analysis
#' @param Cooccure.resulst Output of spatial.celltype.cooccurence funcion
#' @return heatmap plot of cooccurence results
#' @export

plot.celltype.cooccurence <- function(Cooccure.resulst){

    x <- Cooccure.resulst
    dim <- x$species
    comat_pos <- comat_neg <- matrix(nrow=dim,ncol=dim)
    co_tab <- x$result
    for (i in 1:nrow(co_tab)){
      comat_pos[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_gt"]
      comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_gt"]

      row.names(comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]])

    }
    for (i in 1:nrow(co_tab)){
      comat_neg[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_lt"]
      comat_neg[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_lt"]
    }
    comat <- ifelse(comat_pos>=0.05,0,1) + ifelse(comat_neg>=0.05,0,-1)
    colnames(comat) <- 1:dim
    row.names(comat) <- 1:dim

    if ("spp_key" %in% names(x)){

      sp1_name <- merge(x=data.frame(order=1:length(colnames(comat)),sp1=colnames(comat)),y=x$spp_key,by.x="sp1",by.y="num",all.x=T)
      sp2_name <- merge(x=data.frame(order=1:length(row.names(comat)),sp2=row.names(comat)),y=x$spp_key,by.x="sp2",by.y="num",all.x=T)

      colnames(comat) <- sp1_name[with(sp1_name,order(order)),"spp"]
      row.names(comat) <- sp2_name[with(sp2_name,order(order)),"spp"]

    }

    comat[is.na(comat)] <- 0

    origN <- nrow(comat)

    postN <- nrow(comat)

    ind <- apply(comat, 1, function(x) all(x==0))
    comat <- comat[names(sort(ind)),]
    ind <- apply(comat, 2, function(x) all(x==0))
    comat <- comat[,names(sort(ind))]

    #comat
    data.m = melt(comat)
    colnames(data.m) <- c("X1","X2","value")
    data.m$X1 <- as.character(data.m$X1)
    data.m$X2 <- as.character(data.m$X2)

    meas <- as.character(unique(data.m$X2))

    dfids <- subset(data.m, X1 == X2)

    X1 <- data.m$X1
    X2 <- data.m$X2

    df.lower = subset(data.m[lower.tri(comat),],X1 != X2)

    X1 <- df.lower$X1
    X2 <- df.lower$X2
    value <- df.lower$value

    p <- ggplot(df.lower, aes(X1, X2)) + geom_tile(aes(fill = factor(value,levels=c(-1,0,1))), colour ="white")
    p <- p + scale_fill_manual(values = c("light blue","dark gray","red"), name = "", labels = c("negative","random","positive"),drop=FALSE) +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),plot.title = element_text(vjust=-4,size=20, face="bold"),panel.background = element_rect(fill='white', colour='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.9, 0.5),legend.text=element_text(size=18)) +
      ggtitle("Celltype Co-occurrence Matrix") +
      xlab("") + ylab("") +
      scale_x_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) +
      scale_y_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE)
    p <- p + geom_text(data=dfids,aes(label=X1),hjust=1,vjust=0,angle = -22.5, size=5, fontface="bold")#, color="dark gray")


    return(p)
}



