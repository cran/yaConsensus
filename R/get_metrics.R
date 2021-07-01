get_metrics <-
function(theoretical, estimated, verbose = TRUE) {
  
  eca <- get_cluster_accuracy(theoretical, estimated)
  ecp <- get_cluster_purity(theoretical, estimated)
  avg <- sqrt(eca * ecp)
  thNofC <- length(unique(theoretical))
  estNofC <- length(unique(estimated))
  
  ans <- c(round(c(eca, ecp, avg), 4), thNofC, estNofC)
  names(ans) <- c("eca", "ecp", "average", "thNofC", "estNofC")
  
  if(verbose) {
    message("cluster accuracy (eca): ", ans[1])
    message("cluster purity (ecp): ", ans[2])
    message("G index (geometric average of eca, and ecp): ", ans[3])
    message("no. of clusters in theoretical partition: ", ans[4])
    message("no. of clusters in estimated partition: ", ans[5])
  }
  
  return(ans)
}
