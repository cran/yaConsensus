plot.yaConsensus <-
function (x, G = 2, 
                               annotation = NULL,
                               annotation.colorCode = NULL, 
                               matching_clustering = NULL, 
                               consensus_colors = NULL,
                               reduce_to = -1,
                               main = NULL,
                               ...) {
  if(!is.null(annotation) & is.null(annotation.colorCode))
    stop("annotation.colorCode required.")
  
  x <- get_consensus_dissimilarity(x, G = G) 
  hhc <- hclust(x$consensus.diss, method = "complete")
  

  is_reduced <- FALSE
  #    if(length(attr(x$consensus.diss, "Labels")) > K) {
  if(reduce_to > 0) {
    K <- max(reduce_to, 50)
    is_reduced <- TRUE
    
    rclust <- as.factor(cutree(hhc, k = K))
    hhc <- cut_k_leaves(hhc, length(hhc$labels) - K)
    #    rclust <- blast_clustering(hhc)
    levels(rclust) <- hhc$labels <- paste0("r", 1:K)
    annotation.tmp <- annotation
    annotation.tmp$consensus.col <- annotation.tmp$consensus <- as.character(rclust)
    
    x$consensus.diss <- reduce_dist(x$consensus.diss, given = rclust)
    #hhc$height <- (hhc$height - min(hhc$height))/(max(hhc$height) - min(hhc$height))
    #x$consensus.diss <- cophenetic(hhc)
    
    if(!is.null(annotation)) {
      tmp.df <- data.frame(row.names = levels(rclust), tmp = levels(rclust))
      for(j in 1:ncol(annotation)) {
        if(ncol(annotation) > 1) tmp.df$tmp <- NA
        colnames(tmp.df)[j] <- colnames(annotation)[j]
        for(i in 1:nrow(tmp.df)) {
          tt <- table(as.character(annotation[which(rclust == levels(rclust)[i]), j]))
          tmp.df[i, j] <- ifelse(length(tt) == 1, names(tt), NA)
        }
      }
      annotation <- tmp.df
    }
  }
  
  clust <- factor(cutree(hhc, k = G))
  levels(clust) <- paste0("cc", levels(clust))
  
  annotation_col <- annotation_row <- data.frame(consensus = clust, row.names = names(clust))
  
  if(!is.null(consensus_colors)) clust.col <- consensus_colors else
    if(is.null(matching_clustering)) clust.col <- ccolors(G) else {
      if(matching_clustering %in% colnames(annotation)) {
        wwhich <- which(colnames(annotation) == matching_clustering)
        if(is_reduced) {
          annotation.tmp$consensus <- clust[annotation.tmp$consensus]
          clust.col <- match_colors(annotation.tmp[, wwhich], annotation.tmp$consensus)
        } else clust.col <- match_colors(annotation[, wwhich], clust)
        clust.col <- annotation.colorCode[clust.col]
      } else stop("The matching_color provided does not match any of the column in the annotation." )
    }
  names(clust.col) <- levels(clust)
  
  ann_colors <- list(consensus = clust.col)
  
  if(!is.null(annotation)) {
    
    for(k in 1:ncol(annotation)) {
      annotation_row$tmp <- factor(annotation[,k])
      ann_colors$tmp <- annotation.colorCode[levels(annotation_row$tmp)]
      colnames(annotation_row)[ncol(annotation_row)] <- colnames(annotation)[k]
      names(ann_colors)[length(ann_colors)] <- colnames(annotation)[k]
    }
  }
  
  this_color <- ifelse(x$bySample, "red2", "royalblue")
  
  tmp <- 1- as.matrix(x$consensus.diss)[hhc$order,]
  
  if(is.null(main) & !is_reduced) 
    main <- paste0("yaConsensus - ", length(x$labels), " items/", 
                    round(100*x$epsilon), "% sampling rate/",x$runs, " samplings")
  
  if(is.null(main) & is_reduced)  
    main <- paste0("yaConsensus - ", K, " out of ", length(x$labels), " items/",
                    round(100*x$epsilon), "% sampling rate/",x$runs, " samplings")
  
  pheatmap(tmp, 
           cluster_rows = FALSE, show_rownames = FALSE,
           cluster_cols = hhc, cutree_cols = G, show_colnames = FALSE,
           annotation_col=annotation_col, 
           annotation_row=annotation_row, 
           annotation_colors = ann_colors,
           color = colorRampPalette(c("white", this_color))(50),
           main = main) 
  
  annotation_row$consensus.col <- clust
  levels(annotation_row$consensus.col) <- clust.col[annotation_row$consensus.col]
  
  k <- 1
  if(is_reduced) {
    annotation.tmp$consensus.col <- annotation.tmp$consensus <- factor(annotation.tmp$consensus)
    levels(annotation.tmp$consensus.col) <- clust.col[levels(annotation.tmp$consensus.col)]
    annotation_row <- annotation.tmp
  }
  
  ans <- list(annotation = annotation_row, ann_colors = ann_colors, 
              hclust = hhc, statistics = print.yaConsensus(x, verbose = FALSE))
  class(ans) <- "yaConsensus_plot"
  invisible(ans)
}
