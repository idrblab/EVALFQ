#################Transformation Methods################################
Box_Cox <- function(data, lambda) {
  lambda<-0.3
  if (length(lambda)>1 || lambda!=0) data <- (sign(data)*abs(data)^lambda-1)/lambda else data <- log(data)
  data
}

#################Centering Methods###################################
MEC <- function(data) {
  centered.data <- data - apply(data, 1, function(x) mean(x, na.rm = TRUE))
  return(centered.data)
}

MDC <- function(data) {
  centered.data <- data - apply(data, 1, function(x) median(x, na.rm = TRUE))
  return(centered.data)
}
################Scaling Methods####################################
AUTO1 <- function(data) {
  scaling.auto <- apply(data, 1, function(x) sd(x, na.rm = TRUE))
  return(scaling.auto)
}
PARETO1 <- function(data) {
  scaling.pareto <- sqrt(apply(data, 1, function(x) sd(x, na.rm = TRUE)))
  return(scaling.pareto)
}
VAST1 <- function(data) {
  scaling.vast <- apply(data, 1, function(x) var(x, na.rm = TRUE))/apply(data, 1, function(x) mean(x, na.rm = TRUE))
  return(scaling.vast)
}
RANGE1 <- function(data) {
  scaling.range <- apply(data, 1, function(x) max(x, na.rm = TRUE))-apply(data, 1, function(x) min(x, na.rm = TRUE))
  return(scaling.range)
}

################Normalization Methods################
SMAD<-function(data2log2){
  mediandata<-apply(data2log2,2,"median",na.rm=T)
  maddata<-apply(data2log2,2,function(x) mad(x,na.rm=T))
  data2mad<-t(apply(data2log2,1,function(x) ((x-mediandata)/maddata)))
  data2mad<-data2mad+mean(mediandata)
  return(data2mad)
}

Basicfilter<-function(x,label,g1,g2){
  label_c<-as.factor(label)
  data<-data.frame(x,label=label_c)
  data_1<-data[data$label==levels(label_c)[1],]
  data_2<-data[data$label==levels(label_c)[2],]

  pos<-c()
  for(k in 1:ncol(x)){
    data_c1<-data_1[,k]
    data_c2<-data_2[,k]
    if(length(data_c1[is.na(data_c1)])>=g1||length(data_c2[is.na(data_c2)])>=g2)
      c(pos,k)->pos
  }

  if(is.null(pos)){
    filterdata<-x
  }else{
    filterdata<-x[,-pos]
  }

  return(filterdata)

}

MEAN <- function(data) {
  #a <- log(data)
  a<-data
  inputdata<-data.frame(as.factor(rep("sample",ncol(a))),t(a))
  norm_mean <- metabolomics::Normalise(inputdata, method = "mean")
  mean<-t(norm_mean$output[,-1])
  return(mean)
}

LINEAR <- function(data) {
  #data<-log(data)##########################################################
  linear.baseline <- apply(data,1,median)
  baseline.mean <- mean(linear.baseline)
  sample.means <- apply(data,2,mean)
  linear.scaling <- baseline.mean/sample.means
  linear.baseline.data <- t(t(data)*linear.scaling)
  return(linear.baseline.data)
}

MEDIAN <- function(data) {
  #data <- log(data)
  inputdata<-data.frame(as.factor(rep("sample",ncol(data))),t(data))
  norm_med <- metabolomics::Normalise(inputdata, method = "median")
  median<-t(norm_med$output[,-1])
  return(median)
}

PQN <- function(data) {
  #data<-log(data) #
  reference <- apply(data, 1, median, na.rm=T)
  quotient <- data/reference
  quotient.median <- apply(quotient, 2, median, na.rm=T)
  pqn.data <- t(t(data)/quotient.median)
  return(pqn.data)
}

QUANTILE <- function(data) {
  normalize.quantile <- get("normalize.quantiles",
                            envir = asNamespace("affy"))
  quantile.data <- normalize.quantile(data)
  rownames(quantile.data) <- rownames(data)
  colnames(quantile.data) <- colnames(data)
  return(quantile.data)
}

RLR1<-function(data2log2){
  mediandata<-apply(data2log2,1,"median",na.rm=T)
  flag1=1
  for(j in 1:ncol(data2log2))
  {
    LRfit <- MASS::rlm(as.matrix(data2log2[,j])~mediandata,na.action=na.exclude)
    Coeffs<-LRfit$coefficients
    a<-Coeffs[2]
    b<-Coeffs[1]
    if(flag1==1)
    {
      globalfittedRLR<-(data2log2[,j]-b)/a
      flag1=2
    }
    else
    {
      globalfittedRLR<-cbind(globalfittedRLR,(data2log2[,j]-b)/a)
    }
  }
  colnames(globalfittedRLR)<-colnames(data2log2)
  rownames(globalfittedRLR)<-rownames(data2log2)
  return(globalfittedRLR)
}

MSTUS <- function(data) {
  data_sum <- matrix(colSums(data, na.rm = TRUE), nrow = 1)
  uni <- matrix(rep(1, nrow(data)), ncol = 1)
  area.uni <- uni %*% data_sum
  MSTUS <- data/area.uni
  return(MSTUS)
}

TMM<-function(datasim){
  datasim[is.na(datasim)]<-0
  res<-tmm(datasim, long = 1000, lc = 0, k = 0)
  return(res)
}

VSN <- function(data) {
  # load package unless it is already loaded
  vsn.model <- vsn::vsn2(data)
  vsn.data <- predict(vsn.model, data)
  return(vsn.data)
}

################Imputation Methods################
back<-function(x){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata<-x
  filterdata[is.na(filterdata)]<-min(filterdata,na.rm=TRUE)
  cObs <- filterdata
  return(cObs)
}

bpca<-function(x,nPcs){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata<-x
  result <- pca(filterdata, method="bpca", nPcs=nPcs, center = FALSE)
  cObs <- completeObs(result)
  return(cObs)
}

censor<-function(x){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata<-x
  filterdata[is.na(filterdata)]<-min(filterdata,na.rm=TRUE)
  cObs <- filterdata
  return(cObs)
}

knn<-function(x,k){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  result<-impute::impute.knn(filterdata, k, rowmax = 0.5, colmax = 0.8, maxp = 1500)
  cObs <- result$data
  return(cObs)
}

svdm<-function(x,nPcs){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata<-x
  result <- pca(filterdata, method="svdImpute", nPcs=nPcs, center = FALSE)
  cObs <- completeObs(result)
  return(cObs)
}

zero<-function(x){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata[is.na(filterdata)]<-0
  cObs <- filterdata
  return(cObs)
}

#' @importFrom pcaMethods llsImpute
lls<-function(x,k){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata<-x
  result <- llsImpute(filterdata, k , correlation="pearson", allVariables=TRUE)
  cObs <- completeObs(result)
  return(cObs)
}

################End imputation################
