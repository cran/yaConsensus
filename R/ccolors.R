ccolors <-
function (k) {
  if (k > 7) {
    tmp <- grep("gray", colors())
    tmp <- c(tmp, grep("grey", colors()))
    tmp <- c(tmp, grep("white", colors()))
    tmp <- colors()[-tmp]
    tmp <- tmp[(1:k) * floor(length(tmp)/k)]
  }
  else tmp <- palette()[1:k]
  return(tmp)
}
