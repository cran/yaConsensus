summary.yaConsensus_plot <-
function(object, given = NULL, ...) {
  
  if(is.null(given)) results <- print.yaConsensus_plot(object) else {
    results <- print.yaConsensus_plot(object, verbose = FALSE)
    
    j <- 1
    
    if(!(given[j] %in% colnames(object$annotation))) stop(given[j], " not found as a-priori clustering.")
    
    wwhich <- which(colnames(object$annotation) == given[j])
    tmp <- get_metrics(object$annotation[, wwhich], object$annotation$consensus, verbose = FALSE)
    results$'entropy consensus accuracy' <- paste0(round(100 * tmp[1], 2), "%")
    results$'entropy consensus precision' <- paste0(round(100 * tmp[2], 2), "%")
    results$'entropy consensus (average)' <- paste0(round(100 * tmp[3], 2), "%")
    
    names(results)[length(results)-2] <- paste0(names(results)[length(results)-2], ", given ", given[j], ": ")
    names(results)[length(results)-1] <- paste0(names(results)[length(results)-1], ", given ", given[j], ": ")
    names(results)[length(results)] <- paste0(names(results)[length(results)], ", given ", given[j], ": ")
    message("yaConsensus statistics")
    for(k in 1:length(results)) message("\t", names(results)[k], results[[k]])
  }
  results <- object$statistics
  #results
  invisible(results)
}
