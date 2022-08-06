`calculateFDR` <-
function(observed, permuted, progress) {
  
   observed <- abs(observed)
   permuted <- abs(permuted)
   ord <- order(observed, decreasing=TRUE, na.last=TRUE)
   a <- observed[ord]
   
   A <- matrix(NA, nrow=length(a), ncol=ncol(permuted))
   if (progress) pb <- txtProgressBar(min=0, max=ncol(A), style=3)
   for(i in seq_len(ncol(A))) {
      a.rand <- sort(permuted[,i], decreasing=TRUE, na.last=TRUE)
      n.bigger <- biggerN(a, a.rand)
      A[ord,i] <- n.bigger/(seq_along(a))
      if (progress) setTxtProgressBar(pb, i)
   }
   if (progress) close(pb)
   
   FDR <- apply(A, 1, median)
   FDR[FDR>1] <- 1
	
   return(FDR)
}
