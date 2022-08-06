`permutatedSamples` <-
function(data, B, cl) {
  samples <- matrix(nrow=B, ncol=ncol(data))
  for(i in seq_len(B)){
    samples[i,] <- sample(seq_len(ncol(data)))
  }
  return(samples)
}
