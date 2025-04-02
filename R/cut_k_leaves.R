cut_k_leaves <-
function(x, k) {
  
  cut_first_2leaves <- function(x) {
    
    x_merge_1 <- x$merge[1,]
    x$merge <- x$merge[-1,]
    
    x$merge[x$merge == 1] <- max(x_merge_1)
    x$labels[-max(x_merge_1)] <- paste(x$labels[-x_merge_1], collapse = ";")
    x$labels <- x$labels[min(x_merge_1)]
    x$order <- x$order[-which(x$order == -min(x_merge_1))]
    tmp <- x$order > -min(x_merge_1)
    x$order[tmp] <- x$order[tmp] - 1
    
    tmp <- x$merge > 1 # tutti i positivi
    x$merge[tmp] <- x$merge[tmp] - 1
    
    tmp <- x$merge < min(x_merge_1)
    x$merge[tmp] <- x$merge[tmp] + 1
    
    x$height <- x$height[-1]
    return(x)
  }
  
  if(length(x$labels) < 3) return(NULL)
  x$original_labels <- x$labels
  for(j in 1:k) x <- cut_first_2leaves(x)
  return(x)
}
