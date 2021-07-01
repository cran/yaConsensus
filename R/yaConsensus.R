yaConsensus <-
function(ddata, runs = 1000, epsilon = .65, is_by_sample = TRUE, 
                        distMethod = "euclidean", hcMethod = "ward.D2", prefix = NULL) {
  
  if(is_by_sample) 
    ans <- consensus.bySample(ddata, runs = runs, epsilon = epsilon, hcMethod = hcMethod, distMethod = distMethod, prefix = prefix) else
      ans <- consensus.byFeature(ddata, runs = runs, epsilon = epsilon, distMethod = distMethod, hcMethod = hcMethod, prefix = prefix) 
    invisible(ans)
}
