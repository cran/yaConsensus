reduce_dist <-
function(ddist, K = 500, given = NULL) {
  
  if(is.null(given)) {
    hhc <- hclust(ddist, method = "complete")
    clusters <- as.factor(cutree(hhc, k = K)) # get the clusters
    namesOfClusters <- levels(clusters) <- paste0("c", levels(clusters))
  } else {
    clusters <- as.factor(given) # get the clusters
    namesOfClusters <- levels(clusters)
    K <- length(unique(given))
  }
  
  ddist <- as.matrix(ddist)
  rdist<- matrix(0, ncol = K, nrow = K)
  colnames(rdist) <- rownames(rdist) <- levels(clusters)
  
  for(j in 1:K) for(i in j:K) {
    rrow <- names(clusters)[which(clusters == namesOfClusters[i])]
    ccol <- names(clusters)[which(clusters == namesOfClusters[j])]
    rdist[i, j] <- mean(ddist[rrow, ccol])
  }
  rdist <- as.dist(tmp <- rdist)
  rdist <- (rdist - min(rdist))/(max(rdist)-min(rdist))
  attr(rdist, "clustering") <- clusters
  invisible(rdist)
}
