get_consensus_dissimilarity <-
function (obj, G = 2) {
  
  ans_template <- rep(0, length(obj$labels))
  names(ans_template) <- obj$labels
  
  clusters <- sapply(obj$hclust, function(x) {
    ans_tmp <- ans_template
    tmp <- cutree(x, k = G)
    ans_tmp[names(tmp)] <- tmp
    return(ans_tmp)
  })
  
  obj$consensus.diss <- consensus.diss(t(clusters))
  obj$G <- G
  obj$hclust <- NULL
  invisible(obj)
}
