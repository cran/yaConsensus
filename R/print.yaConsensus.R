print.yaConsensus <-
function(x, verbose = T, ...) {
    results <- list('consensus type:' = ifelse(x$bySample, " sampling of samples", " sampling of features"),
                    'number of items: ' = length(x$labels),
                    'sampling rate: ' = paste0(round(100*x$epsilon, 2), "%"),
                    'number of samplings: ' = x$runs,
                    'inner distance: ' = x$distMethod,
                    'inner hc method: ' = x$hcMethod,
                    'outer distance: ' = "consensus",
                    'outer hc method: ' = "complete",
                    'running time (seconds): ' = round(x$elapsed_time, 4),
                    'number of cores: ' = x$ncores)
    if(verbose) {
      message("yaConsensus statistics")
      for(k in 1:length(results)) message("\t", names(results)[k], results[[k]])
    }
    invisible(results)
  }
