get_ecdf <-
function(obj, x = obj$monti$x)  {
  
  Fn <- ecdf(1 - unclass(obj$consensus.diss))
  y <- Fn(x)
  n <- length(y)
  area <- sum((y[2:n] + y[1:(n-1)])*(x[2:n] - x[1:(n-1)])/2)
  ans <- list(y = y, area = area)
  invisible(ans)
  
}
