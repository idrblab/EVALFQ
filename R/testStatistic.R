`testStatistic` <-
 function(X, Y, paired){
   ## Calculates the test statistic for each row.
   ## X and Y are the data matrices of the two groups.
   ## Each row of these two matrices must contain at least TWO not NA values.
   ## Thus the "variance" always exists.  
   
   ## Row means
   mX <- rowMeans(X, na.rm=TRUE)
   mY <- rowMeans(Y, na.rm=TRUE)

   ## Pooled standard deviations for each row
   sX <- rowSums((X - mX)^2, na.rm=TRUE)
   sY <- rowSums((Y - mY)^2, na.rm=TRUE)

   if(!paired) {
      ## Number of not NA values in each row
      nX <- rowSums(!is.na(X))
      nY <- rowSums(!is.na(Y))
    
      ## d == difference between the group means for each row (==gene)
      ## s == pooled standard deviation for each row (==gene)        
      d <- mY - mX
      s <- sqrt(((sX + sY) / (nX + nY - 2)) * (1 / nX + 1 / nY))
    
      ## Cases with less than two non-missing values.
      ## Set d = 0, s = 1
      ind <- which( nY < 2 | nX < 2 )
      d[ind] <- 0
      s[ind] <- 1
   }
 
   if(paired) {
      ## Add for paired
      sXY <- rowSums((X - mX)*(Y - mY), na.rm=TRUE)

      ## Number of not NA values in each row
      n <- rowSums(!is.na(X*Y))
    
      ## d == difference between the group means for each row (==gene)
      ## s == pooled standard deviation for each row (==gene)        
      d <- mY - mX
      s <- sqrt(((sX + sY) / (n + n - 2)) * (2 / n) - 2/(n*n-n)*sXY)
    
      ## Cases with less than two non-missing values.
      ## Set d = 0, s = 1
      ind <- which( n < 2 )
      d[ind] <- 0
      s[ind] <- 1
   }
   
   return(list(d=d, s=s))
}
