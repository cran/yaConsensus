print.yaConsensus_plot <-
function(x, verbose = T, ...) {
  results <- x$statistics
  if(verbose) {
    message("yaConsensus statistics")
    for(k in 1:length(results)) message("\t", names(results)[k], results[[k]])
  }
  invisible(results)
}
