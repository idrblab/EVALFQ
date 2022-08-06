# x The first argument onesss
# percent The second argument onesss
quan.norm <- function(x,percent=50) {
  low <- 0.5*(100 - percent)/100
  high <- 0.5*(100 + percent)/100
  difference <- as.vector(diff(quantile(x, probs=c(low,high), na.rm = TRUE)))
  return(difference)
}


# x The first argument ones
# percent The second argument ones
# quan.norm
quartile.normalize <- function(x,percent=50) {
  quartile <- apply(x,2,quan.norm,percent=percent)
  max.quartile <- max(quartile)
  ratio <- (quartile/max.quartile)
  ratio.vect <- rep(ratio,nrow(x))
  adjusted <- matrix(ratio.vect, nrow=nrow(x),byrow=TRUE)
  normalized <- data.frame(x/adjusted)
  return(normalized)
}



# x The first argument one
# y The second argument one
lowess.normalize <- function(x,y)
{
  # x = log(cy3 or chip1) and y = log(cy5 or chip2)
  na.point <- (1:length(x))[!is.na(x) & !is.na(y)]
  x <- x[na.point]; y <- y[na.point]
  fit <- lowess(x+y, y-x)

  diff.fit <- approx(fit$x,fit$y,x+y,ties=mean)
  out <- y - diff.fit$y
  return(out)
}


# x The first argument one
# data.type The second argument one
# threshold The third argument one
# LOWESS The fourth argument one
# quartile.normalize
# lowess.normalize

preprocess <- function(x, data.type = "MAS5", threshold=1,LOWESS=FALSE) {

  # Removing NA values
  x <- as.matrix(na.exclude(x))

  # IQR normalization
  if (data.type =="MAS4" || data.type == "MAS5") {
    x <- quartile.normalize(x, percent=50)
  }

  # Thresholding to 'threshold' (default = 1)
  if (data.type == "MAS4" || data.type =="MAS5"|| data.type == "dChip") {
    if (length(x[x<threshold]) !=0) {
      x[x<threshold] <- threshold
    }
  }

  # Log based 2 transformation
  x <- logb(x,2)

  # Loess normalization of all the chips w.r.t. first one
  if (LOWESS) {
    y <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
    y[,1] <- x[,1]
    for (i in 2:ncol(x)) {
      y[,i] <- lowess.normalize(x[,1],x[,i])
    }
    x <- y
  }
  return(x)
}

# Above function preprocesses the data from MAS4/5 and dchip.
# First, IQR (inter-quartile normalization) is applied to the data
# from MAS 4/5. Then for MAS4/dChip data thresholding is applied
# at 1 and for MAS5 data, thresholding is applied at 0.1
# Finally, the data is log transformed to the base 2.

#lowess normalization
# datasim The first argument oness
# preprocess
LOWESS<-function(datasim){
  res<-preprocess(datasim,data.type="MAS5",LOWESS=TRUE)
  return(res)
}
