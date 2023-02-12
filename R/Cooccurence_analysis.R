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
  COI.spots <- rownames(spatial.object@meta.data[which(spatial.object@meta.data$CompositionCluster_CC %in% COI & spatial.object@meta.data$orig.ident %in% Condition),])
  coocur.COI.exp <- t(deconv.prob.mat[COI.spots,])
  coocur.COI.exp <- biclust::binarize(coocur.COI.exp, threshold=prob.th)
  cooccur.COI.res <- ISCHIA.cooccur(mat = coocur.COI.exp, type = "spp_site", spp_names = TRUE, prob = "comb", thresh = FALSE)
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
  #Idents(spatial.object) <- "CompositionCluster_CC"
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

  COI.spots <- names(spatial.object$CompositionCluster_CC[which(spatial.object$CompositionCluster_CC %in% COI & spatial.object$orig.ident %in% Condition)])
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
  cooccur.COI.res <- ISCHIA.cooccur(mat = coocur.COI.exp, type = "spp_site", thresh = TRUE, spp_names = TRUE)
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
    COI.cells <- rownames(SeuratObj@meta.data[which(SeuratObj@meta.data$CompositionCluster_CC %in% COI & SeuratObj@meta.data$orig.ident %in% Condition),])
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
    data.m = reshape2::melt(comat)
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




ISCHIA.cooccur <- function (mat, type = "spp_site", thresh = TRUE, spp_names = FALSE,
          true_rand_classifier = 0.1, prob = "hyper", site_mask = NULL,
          only_effects = FALSE, eff_standard = TRUE, eff_matrix = FALSE)
{
  if (type == "spp_site") {
    spp_site_mat <- mat
  }
  if (type == "site_spp") {
    spp_site_mat <- t(mat)
  }
  if (spp_names == TRUE) {
    spp_key <- data.frame(num = 1:nrow(spp_site_mat), spp = row.names(spp_site_mat))
  }
  if (!is.null(site_mask)) {
    if (nrow(site_mask) == nrow(spp_site_mat) & ncol(site_mask) ==
        ncol(spp_site_mat)) {
      N_matrix <- create.N.matrix(site_mask)
    }
    else {
      stop("Incorrect dimensions for site_mask, aborting.")
    }
  }
  else {
    site_mask <- matrix(data = 1, nrow = nrow(spp_site_mat),
                        ncol = ncol(spp_site_mat))
    N_matrix <- matrix(data = ncol(spp_site_mat), nrow = nrow(spp_site_mat),
                       ncol = nrow(spp_site_mat))
  }
  spp_site_mat[spp_site_mat > 0] <- 1
  tsites <- ncol(spp_site_mat)
  nspp <- nrow(spp_site_mat)
  spp_pairs <- choose(nspp, 2)
  incidence <- prob_occur <- obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow = spp_pairs,
                                                                                  ncol = 3)
  incidence <- prob_occur <- matrix(nrow = nrow(N_matrix),
                                    ncol = ncol(N_matrix))
  for (spp in 1:nspp) {
    if (spp < nspp) {
      for (spp_next in (spp + 1):nspp) {
        incidence[spp, spp_next] <- sum(site_mask[spp,
        ] * site_mask[spp_next, ] * mat[spp, ])
        incidence[spp_next, spp] <- sum(site_mask[spp,
        ] * site_mask[spp_next, ] * mat[spp_next, ])
      }
    }
  }
  prob_occur <- incidence/N_matrix
  pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)),
                       style = 3)
  row <- 0
  for (spp in 1:nspp) {
    if (spp < nspp) {
      for (spp_next in (spp + 1):nspp) {
        pairs <- sum(as.numeric(mat[spp, site_mask[spp,
        ] * site_mask[spp_next, ] == 1] == 1 & mat[spp_next,
                                                   site_mask[spp, ] * site_mask[spp_next, ] ==
                                                     1] == 1))
        row <- row + 1
        obs_cooccur[row, 1] <- spp
        obs_cooccur[row, 2] <- spp_next
        obs_cooccur[row, 3] <- pairs
        prob_cooccur[row, 1] <- spp
        prob_cooccur[row, 2] <- spp_next
        prob_cooccur[row, 3] <- prob_occur[spp, spp_next] *
          prob_occur[spp_next, spp]
        exp_cooccur[row, 1] <- spp
        exp_cooccur[row, 2] <- spp_next
        exp_cooccur[row, 3] <- prob_cooccur[row, 3] *
          N_matrix[spp, spp_next]
      }
    }
    setTxtProgressBar(pb, spp)
  }
  if (thresh == TRUE) {
    n_pairs <- nrow(prob_cooccur)
    prob_cooccur <- prob_cooccur[exp_cooccur[, 3] >= 1, ]
    obs_cooccur <- obs_cooccur[exp_cooccur[, 3] >= 1, ]
    exp_cooccur <- exp_cooccur[exp_cooccur[, 3] >= 1, ]
    n_omitted <- n_pairs - nrow(prob_cooccur)
    pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)),
                         style = 3)
  }
  output <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(output) <- c("sp1", "sp2", "sp1_inc", "sp2_inc",
                        "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt",
                        "p_gt")
  for (row in 1:nrow(obs_cooccur)) {
    sp1 <- obs_cooccur[row, 1]
    sp2 <- obs_cooccur[row, 2]
    sp1_inc <- incidence[sp1, sp2]
    sp2_inc <- incidence[sp2, sp1]
    max_inc <- max(sp1_inc, sp2_inc)
    min_inc <- min(sp1_inc, sp2_inc)
    nsite <- N_matrix[sp1, sp2]
    psite <- as.numeric(nsite + 1)
    prob_share_site <- rep(x = 0, times = psite)
    if (prob == "hyper") {
      if (only_effects == FALSE) {
        all.probs <- phyper(0:min_inc, min_inc, nsite -
                              min_inc, max_inc)
        prob_share_site[1] <- all.probs[1]
        for (j in 2:length(all.probs)) {
          prob_share_site[j] <- all.probs[j] - all.probs[j -
                                                           1]
        }
      }
      else {
        for (j in 0:nsite) {
          if ((sp1_inc + sp2_inc) <= (nsite + j)) {
            if (j <= min_inc) {
              prob_share_site[(j + 1)] <- 1
            }
          }
        }
      }
    }
    if (prob == "comb") {
      if (only_effects == FALSE) {
        for (j in 0:nsite) {
          if ((sp1_inc + sp2_inc) <= (nsite + j)) {
            if (j <= min_inc) {
              prob_share_site[(j + 1)] <- coprob(max_inc = max_inc,
                                                 j = j, min_inc = min_inc, nsite = nsite)
            }
          }
        }
      }
      else {
        for (j in 0:nsite) {
          if ((sp1_inc + sp2_inc) <= (nsite + j)) {
            if (j <= min_inc) {
              prob_share_site[(j + 1)] <- 1
            }
          }
        }
      }
    }
    p_lt <- 0
    p_gt <- 0
    for (j in 0:nsite) {
      if (j <= obs_cooccur[row, 3]) {
        p_lt <- prob_share_site[(j + 1)] + p_lt
      }
      if (j >= obs_cooccur[row, 3]) {
        p_gt <- prob_share_site[(j + 1)] + p_gt
      }
      if (j == obs_cooccur[row, 3]) {
        p_exactly_obs <- prob_share_site[(j + 1)]
      }
    }
    p_lt <- round(p_lt, 5)
    p_gt <- round(p_gt, 5)
    p_exactly_obs <- round(p_exactly_obs, 5)
    prob_cooccur[row, 3] <- round(prob_cooccur[row, 3], 3)
    exp_cooccur[row, 3] <- round(exp_cooccur[row, 3], 1)
    output[row, ] <- c(sp1, sp2, sp1_inc, sp2_inc, obs_cooccur[row,
                                                               3], prob_cooccur[row, 3], exp_cooccur[row, 3], p_lt,
                       p_gt)
    setTxtProgressBar(pb, nspp + row)
  }
  close(pb)
  if (spp_names == TRUE) {
    sp1_name <- merge(x = data.frame(order = 1:length(output$sp1),
                                     sp1 = output$sp1), y = spp_key, by.x = "sp1", by.y = "num",
                      all.x = T, sort = FALSE)
    sp2_name <- merge(x = data.frame(order = 1:length(output$sp2),
                                     sp2 = output$sp2), y = spp_key, by.x = "sp2", by.y = "num",
                      all.x = T, sort = FALSE)
    output$sp1_name <- sp1_name[with(sp1_name, order(order)),
                                "spp"]
    output$sp2_name <- sp2_name[with(sp2_name, order(order)),
                                "spp"]
  }
  true_rand <- (nrow(output[(output$p_gt >= 0.05 & output$p_lt >=
                               0.05) & (abs(output$obs_cooccur - output$exp_cooccur) <=
                                          (tsites * true_rand_classifier)), ]))
  output_list <- list(call = match.call(), results = output,
                      positive = nrow(output[output$p_gt < 0.05, ]), negative = nrow(output[output$p_lt <
                                                                                              0.05, ]), co_occurrences = (nrow(output[output$p_gt <
                                                                                                                                        0.05 | output$p_lt < 0.05, ])), pairs = nrow(output),
                      random = true_rand, unclassifiable = nrow(output) - (true_rand +
                                                                             nrow(output[output$p_gt < 0.05, ]) + nrow(output[output$p_lt <
                                                                                                                                0.05, ])), sites = N_matrix, species = nspp, percent_sig = (((nrow(output[output$p_gt <
                                                                                                                                                                                                            0.05 | output$p_lt < 0.05, ])))/(nrow(output))) *
                        100, true_rand_classifier = true_rand_classifier)
  if (spp_names == TRUE) {
    output_list$spp_key <- spp_key
    output_list$spp.names = row.names(spp_site_mat)
  }
  else {
    output_list$spp.names = c(1:nrow(spp_site_mat))
  }
  if (thresh == TRUE) {
    output_list$omitted <- n_omitted
    output_list$pot_pairs <- n_pairs
  }
  class(output_list) <- "cooccur"
  if (only_effects == F) {
    output_list
  }
  else {
    effect.sizes(mod = output_list, standardized = eff_standard,
                 matrix = eff_matrix)
  }
}


coprob <-
  function(max_inc,j,min_inc,nsite){
    as.numeric(round(chooseZ(max_inc,j) * chooseZ(nsite - max_inc, min_inc - j),0) / round(chooseZ(nsite,min_inc),0))
  }

