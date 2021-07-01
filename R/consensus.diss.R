consensus.diss <-
function (cclusters, similarity = FALSE)   {
    getLoadedDLLs()
    n <- ncol(cclusters)
    H <- nrow(cclusters)
    m <- n * (n - 1)/2
    tmp <- rep(0, m)
    storage.mode(cclusters) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(H) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(tmp) <- "integer"
    ans <- .Fortran("cconsensus", cclusters, n, H, tmp, m)[[4]]
    ans <- as.double(ans)
    
    if(sum(cclusters < 1) > 0) {
      cclusters  <-  0 + (cclusters > 0)
      tmp <- rep(0, m)
      storage.mode(cclusters) <- "integer"
      storage.mode(n) <- "integer"
      storage.mode(H) <- "integer"
      storage.mode(m) <- "integer"
      storage.mode(tmp) <- "integer"
      tmp <- .Fortran("cconsensus", cclusters, n, H, tmp, m)[[4]]
      tmp <- as.double(tmp)
    } else tmp  <- H
    
    ans <- ans/tmp
    
    if(!similarity)
      ans <- 1 - ans
    
    attr(ans, "Size") <- n
    attr(ans, "Labels") <- colnames(cclusters)
    attr(ans, "Diag") <- FALSE
    attr(ans, "method") <- "consensus"
    class(ans) <- "dist"
    return(ans)
  }
