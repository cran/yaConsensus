plot.yaConsensus <-
function (x, G = 2, 
                              annotation = NULL,
                              annotation.colorCode = NULL, 
                              matching_clustering = NULL, 
                              consensus_colors = NULL,
                              ...) {
  if(!is.null(annotation) & is.null(annotation.colorCode))
    stop("annotation.colorCode required.")
  
  
  x <- get_consensus_dissimilarity(x, G = G) 
  
#  suppressPackageStartupMessages(require(pheatmap))
  
  hhc <- hclust(x$consensus.diss, method = "complete")
  clust <- factor(cutree(hhc, k = G))
  levels(clust) <- paste0("cc", levels(clust))
  
  annotation_col <- annotation_row <- data.frame(consensus = clust, row.names = names(clust))
  
  if(!is.null(consensus_colors)) clust.col <- consensus_colors else
    if(is.null(matching_clustering)) clust.col <- ccolors(G) else {
      if(matching_clustering %in% colnames(annotation)) {
        wwhich <- which(colnames(annotation) == matching_clustering)
        clust.col <- match_colors(annotation[, wwhich], clust)
        clust.col <- annotation.colorCode[clust.col]
      } else stop("The matching_color provided does not match any of the column in the annotation." )
    }
  names(clust.col) <- levels(clust)
  
  ann_colors <- list(consensus = clust.col)
  
  if(!is.null(annotation)) 
    for(k in 1:ncol(annotation)) {
      annotation_row$tmp <- factor(annotation[,k])
      ann_colors$tmp <- annotation.colorCode[levels(annotation_row$tmp)]
      colnames(annotation_row)[ncol(annotation_row)] <- colnames(annotation)[k]
      names(ann_colors)[length(ann_colors)] <- colnames(annotation)[k]
    }
  
  this_color <- ifelse(x$bySample, "red2", "royalblue")
  
  tmp <- 1- as.matrix(x$consensus.diss)[hhc$order,]
  
  pheatmap(tmp, 
           cluster_rows = FALSE, show_rownames = FALSE,
           cluster_cols = hhc, cutree_cols = G, show_colnames = FALSE,
           annotation_col=annotation_col, 
           annotation_row=annotation_row, 
           annotation_colors = ann_colors,
           color = colorRampPalette(c("white", this_color))(50),
           xlab = "ste") 
  
  annotation_row$consensus.col <- clust
  levels(annotation_row$consensus.col) <- clust.col[annotation_row$consensus.col]
  
  ans <- list(annotation = annotation_row, ann_colors = ann_colors, hclust = hhc,
              statistics = print.yaConsensus(x, verbose = FALSE)
  )
  class(ans) <- "yaConsensus_plot"
  invisible(ans)
}
