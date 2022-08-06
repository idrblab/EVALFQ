`calculateP` <- function(observed, permuted) {
  # Store order for later use
  observed_order <- order(abs(observed), decreasing=TRUE)
  
  # Sort observed and permuted values
  observed <- sort(abs(observed), decreasing=TRUE)
  permuted <- sort(abs(as.vector(permuted)), decreasing=TRUE)
  
  # Get p-values from C++ code
  # (expects ordered vectors)
  p <- pvalue(observed, permuted)
  
  # Revert to original ordering
  results <- vector(mode="numeric", length=length(p))
  for(i in seq_along(results)) {
    results[observed_order[i]] <-  p[i]
  }
  
  return(results)
}
