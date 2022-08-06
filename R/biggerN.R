`biggerN` <-
function(x, y) {
 
   x <- sort(x, decreasing=TRUE, na.last=TRUE)		# sort x in decreasing order
   y <- sort(y, decreasing=TRUE, na.last=TRUE)		# sort y in decreasing order
   a <- match(x, x)				 # vector of the positions of (first) matches of the first argument in the second
   b <- x %in% y				   # a logical vector indicating if there is a match or not for its left operand
   z <- sort(c(x, y), decreasing=TRUE, na.last=TRUE)		# sort c(x,y) in decreasing order
   d <- match(x, z)				 # vector of the positions of (first) matches of the first argument in the second
 
   return(d - a + b)
}
