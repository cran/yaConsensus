cal_entropy <-
function(x) {
  freqs <- table(x)/length(x)
  freqs = freqs[freqs>0]
  return(-sum(freqs * log(freqs))/log(length(x)))
}
