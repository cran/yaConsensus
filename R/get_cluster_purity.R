get_cluster_purity <-
function(theoretical, estimated)
  return(1-mean(unlist(lapply(unique(theoretical), function(x){cal_entropy(estimated[theoretical == x])})), na.rm = TRUE))
