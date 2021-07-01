get_cluster_accuracy <-
function(theoretical, estimated) 
  return(1-mean(unlist(lapply(unique(estimated), function(x){cal_entropy(theoretical[estimated == x])})), na.rm = TRUE))
