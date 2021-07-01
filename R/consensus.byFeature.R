consensus.byFeature <-
function(ddata, runs = 1000, epsilon = 0.1, # ddata is a sample X gene matrix
                                distMethod = "euclidean", hcMethod = "ward.D2", prefix = NULL) {
  #  if(class(ddata) != "matrix") stop("The input data is not a data matrix.")
  
#  suppressPackageStartupMessages(require(foreach))
  mm <- floor(epsilon * ncol(ddata))
  llabels <- rownames(ddata)
  
  
  ttime <- system.time({
    hhclust <- foreach(k = 1:runs) %dopar% {
      these_features <- sample.int(ncol(ddata), mm, replace = FALSE)
      ddist <- dist(ddata[, these_features], method = distMethod)
      return(hclust(ddist, method = hcMethod))}
  })[3]
  
  
  ans <- list()
  ans$distMethod <- distMethod
  ans$hcMethod <- hcMethod
  ans$labels <- llabels
  ans$bySample <- FALSE
  ans$epsilon <- epsilon
  ans$subsetDimension <- mm
  ans$runs <- runs
  ans$hclust <- hhclust
  names(ttime) <- NULL
  ans$elapsed_time <- ttime
  ans$ncores <- getDoParWorkers()
  class(ans) <- "yaConsensus"
  
  if(!is.null(prefix)) {
    ans$fname <- paste0(prefix, "_yaConsensus_eps", ans$epsilon * 100, 
                        "pct_runs", ans$runs, "_byFeatures.RData") 
    aConsensus <- ans
    save(aConsensus, file = ans$fname)
    message("yaConsensus data structure saved in ", ans$fname)
  } 
  invisible(ans)
}
