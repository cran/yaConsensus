consensus.bySample <-
function (ddata, runs = 1000, epsilon = .65, 
                                distMethod = "euclidean", hcMethod = "ward.D2", prefix = NULL) {
  
  #  if (class(ddata) != "dist") stop("A 'distance' object is required.")
  
#  require(foreach)
  
  if ("dist" %in% class(ddata)) {
    n <- attr(ddata, "Size") 
    mm <- floor(epsilon * n)
    distMethod <- attr(ddata, "method")
    llabels <- attr(ddata, "Labels")
    
    ttime <- system.time({
      hhclust <- foreach(k = 1:runs) %dopar% {
        these_samples <- sample.int(n, mm, replace = FALSE)
        tmp <- hclust(as.dist(as.matrix(ddata)[these_samples, these_samples]), method = hcMethod)
        return(tmp)}
    })[3]
  } else {
    n <- nrow(ddata) 
    mm <- floor(epsilon * n)
    llabels <- rownames(ddata)
    
    ttime <- system.time({
      hhclust <- foreach(k = 1:runs) %dopar% {
        these_samples <- sample.int(nrow(ddata), mm, replace = FALSE)
        ddist <- dist(ddata[these_samples,], method = distMethod)
        return(hclust(ddist, method = hcMethod))}
    })[3]
    
  }
  
  ans <- list()
  ans$distMethod <- distMethod
  ans$hcMethod <- hcMethod
  ans$labels <- llabels
  ans$bySample <- TRUE
  ans$epsilon <- epsilon
  ans$subsetDimension <- mm
  ans$runs <- runs
  ans$samples <- NULL
  ans$hclust <- hhclust
  names(ttime) <- NULL
  ans$elapsed_time <- ttime
  ans$ncores <- getDoParWorkers()
  class(ans) <- "yaConsensus"
  if(!is.null(prefix)) {
    ans$fname <- paste0(prefix, "_yaConsensus_eps", ans$epsilon * 100, 
                        "pct_runs", ans$runs, "_bySamples.RData")
    aConsensus <- ans
    save(aConsensus, file = ans$fname)
    message("yaConsensus data structure saved in ", ans$fname)
  } 
  #attr(ans,"package")
  invisible(ans)
}
