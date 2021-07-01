match_colors <-
function(refClust, consClust) {
  (ttable <- table(consClust, refClust))
  (ttable <- ttable/rowSums(ttable))
  (whichMaxTable <- apply(ttable, 1, which.max))
  ans <- colnames(ttable)[whichMaxTable]
  names(ans) <- rownames(ttable)
  return(ans)
}
