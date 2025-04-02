monti <-
function(obj, gMax = nclass.Sturges(obj$labels), just_compute = FALSE) {
  
  if(is.null(obj$fname) & just_compute) {
    message("No file to updated with monti statistics.")
    invisible(obj)
  }
  
  if(is.null(obj$monti)) {
    obj$monti$x <- seq(from = 0, to = 1, length.out = 500)
    obj$monti$y <- matrix(NA, ncol = gMax-1, nrow = 500)
    obj$monti$area <- rep(NA, gMax-1)
    
    for(G in 2:gMax) {
      tmp <- get_consensus_dissimilarity(obj, G = G)
      #      obj <- get_consensus_dissimilarity(obj, G = G)
      tmp <- get_ecdf(tmp)
      obj$monti$y[, G-1] <- tmp$y
      obj$monti$area[G-1] <- tmp$area
    }
  }
  
  if(just_compute) {
    
    aConsensus <- obj # controllare che il file sia yet salvato
    save(aConsensus, file = aConsensus$fname) 
    message(aConsensus$fname, " updated with monti statistics.")
    
  } else {
    delta <- obj$monti$area
    m <- gMax-1
    y <- c(0, (delta[2:m] - delta[1:(m-1)])/delta[1:(m-1)])
    #    y <- c(0, (delta[2:m] - delta[1:(m-1)]))
    x <- 2:(m+1)
    
    ccol <- rep("black", m)
    for(g in 2:(m-1)) if (y[g - 1] < y[g] & y[g + 1] < y[g]) ccol[g] <- "red"
    
    col <- colorRampPalette(c("cyan2", "royalblue"))(gMax)
    col[1+which(ccol == "red")] <- "red"
    
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    layout(matrix(c(1,2), 1, 2, byrow = TRUE), respect = TRUE)
    plot(obj$monti$x, obj$monti$y[, 1], type = "l", ylim = c(0, 1), col = col[2], 
         xlab = "consensus values", ylab = "distribution functions")
    for(G in 2:(gMax-1)) lines(obj$monti$x, obj$monti$y[, G], col = col[G+1])
    
    plot(x, y, type = "l", col = "grey", xaxt = "n", bty = "n",
         xlab = "", ylab = "Delta(K)", xlim = c(1, gMax), 
         omi = c(0.2, 0.15, 0.2, 0.15))
    text(x, y, paste(2:gMax), col = ccol)
    par(def.par)  #- reset to default
  }
  invisible(obj)
}
