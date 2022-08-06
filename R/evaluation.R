
###==========================================================================###
### Other evaluation indexes, such as PEV, PMAD, PCV, and so on.
###==========================================================================###

#calculate PEV
PEV <- function(data0){
  data0<-as.matrix(data0)
  data0<-data0[order(data0[,1]),]
  label <- as.factor(data0[,1])
  data <- data0[,-1]
  data <- t(as.matrix(data))
  x<-levels(label)[1]
  z<-1
  y<-1
  flag<-1
  count<-0
  varmem<-vector()
  tempvar<-vector()
  nonmissingmat<-vector()
  for(i in 1:length(label))
  {
    if(x!=label[i] || i==length(label))
    {
      y<-i-1
      if(i==length(label))
      {
        y<-i
      }
      if(flag==1)
      {
        count<-count+1
        nonmissingmat<-(apply(data[,z:y],1,function(x) {((sum(!is.na(x))))}))-1
        tempvar<-nonmissingmat*apply(data[,z:y],1,function(x) {var(x,na.rm=TRUE)})
      }
      if(flag==2)
      {
        count<-count+1
        nonmissingmat<-(apply(data[,z:y],1,function(x) {((sum(!is.na(x))))}))-1
        tempvar<-nonmissingmat*apply(data[,z:y],1,function(x) {var(x,na.rm=TRUE)})
      }
      varmem<-c(varmem,((sum(tempvar,na.rm=T))/(sum(nonmissingmat,na.rm=T))))
      z<-i
      x<-label[i]
      flag=2;
    }
  }
  avgvarmem<-varmem
  names(avgvarmem)<-levels(label)
  return(avgvarmem)
}
#calculate PMAD
PMAD <- function(data0){
  data0<-as.matrix(data0)
  data0<-data0[order(data0[,1]),]
  label <- as.factor(data0[,1])
  data <- data0[,-1]
  data <- t(as.matrix(data))
  data <- apply(data, 1:2, as.numeric)
  x<-levels(label)[1]
  z<-1
  y<-1
  flag<-1
  count<-0
  madmem<-matrix(nrow=nrow(data),ncol=length(levels(as.factor(unlist(label)))),byrow=T)
  for(i in 1:length(label))
  {
    if(x!=label[i] || i==length(label))
    {
      y<-i-1
      if(i==length(label))
      {
        y<-i
      }
      if(flag==1)
      {
        count<-count+1
        madmem[,count]<-apply(data[,z:y],1,function(x) {mad(x,na.rm=T)})
      }
      if(flag==2)
      {
        count<-count+1
        madmem[,count]<-apply(data[,z:y],1,function(x) {mad(x,na.rm=T)})
      }
      z<-i
      x<-label[i]
      flag=2;
    }
  }
  avgmadmem<-apply(madmem,2,mean,na.rm=T)
  names(avgmadmem)<-levels(label)
  return(avgmadmem)
}

PMAD1 <- function(data0){
  data0<-as.matrix(data0)
  data0<-data0[order(data0[,1]),]
  label <- as.factor(data0[,1])
  data <- data0[,-1]
  data <- t(as.matrix(data))
  data <- apply(data, 1:2, as.numeric)
  x<-levels(label)[1]
  z<-1
  y<-1
  flag<-1
  count<-0
  madmem<-matrix(nrow=nrow(data),ncol=length(levels(as.factor(unlist(label)))),byrow=T)
  for(i in 1:length(label))
  {
    if(x!=label[i] || i==length(label))
    {
      y<-i-1
      if(i==length(label))
      {
        y<-i
      }
      if(flag==1)
      {
        count<-count+1
        madmem[,count]<-apply(data[,z:y],1,function(x) {mad(x,na.rm=T)})
      }
      if(flag==2)
      {
        count<-count+1
        madmem[,count]<-apply(data[,z:y],1,function(x) {mad(x,na.rm=T)})
      }
      z<-i
      x<-label[i]
      flag=2;
    }
  }
  avgmadmem<-apply(madmem,2,mean,na.rm=T)
  names(avgmadmem)<-levels(label)
  return(madmem)
}


#calculate PCV
#' @importFrom RcmdrMisc numSummary
PCV <- function(data0){
  data0<-as.matrix(data0)
  data0<-data0[order(data0[,1]),]
  label <- as.factor(data0[,1])
  data <- data0[,-1]
  data <- t(as.matrix(data))
  data <- apply(data, 1:2, as.numeric)

  tempcvmat<-matrix(nrow=nrow(data),ncol=length(levels(as.factor(unlist(label)))),byrow=T)
  for(i in 1:nrow(data))
  {
    tempcv<-RcmdrMisc::numSummary(data[i,],statistics=c("cv"),groups=unlist(label))
    tempcvmat[i,]<-tempcv$table
  }

  colnames(tempcvmat)[1:2]<-levels(label)
  temcvmatsum<-apply(tempcvmat,1,sum,na.rm=T)
  result<-cbind(tempcvmat,total=temcvmatsum)
  result<-na.omit(result)
  colnames(result)[1:2]<-levels(label)
  return(result)
}



PCV1 <- function(data0){
  data0<-as.matrix(data0)
  data0<-data0[order(data0[,1]),]
  label <- as.factor(data0[,1])
  data <- data0[,-1]
  data <- t(as.matrix(data))
  data <- apply(data, 1:2, as.numeric)

  tempcvmat<-matrix(nrow=nrow(data),ncol=length(levels(as.factor(unlist(label)))),byrow=T)
  for(i in 1:nrow(data))
  {
    tempcv<-RcmdrMisc::numSummary(data[i,],statistics=c("cv"),groups=unlist(label))
    tempcvmat[i,]<-tempcv$table
  }

  colnames(tempcvmat)[1:2]<-levels(label)
  temcvmatsum<-apply(tempcvmat,1,sum,na.rm=T)
  result<-cbind(tempcvmat,total=temcvmatsum)
  #result<-na.omit(result)
  result[is.na(result)]<-0
  colnames(result)[1:2]<-levels(label)
  return(result)
}

plotrCV<-function(k, estvar){

  #only signals without negativ variances
  #pos<-(estvar[,2:k]>0)
  #VarNeg<-unique(unlist(sapply(1:(k-2), function(i){which(pos[,i]==FALSE)})))
  #SiganzTrue<- nrow(estvar)-length(VarNeg)


  values<-unlist(as.vector(estvar))
  names<-rep(1:(k+1),each=nrow(estvar))
  colpl<-rainbow(k+1)
  ymax<-max(estvar)

  #plot in file

  par(mar=c(2,2,0.1,0.2))
  stripchart(values~names, method="jitter" ,jitter=0.3,vertical=TRUE,
             pch=8,cex=0.2, ylim=c(0,ymax),
             col= colpl,
             at=1:(k+1),
             group.names=colnames(estvar),
             xlab="", ylab="",cex.lab=2,
             main="", cex.main=1.2)

  sapply(1:(k+1), function(i){segments(i-0.4,mean(as.numeric(estvar[,i])), i+0.4, mean(as.numeric(estvar[,i])), lwd=1.8)})
  sapply(1:(k+1), function(i){text(i,ymax,round(1000*mean(as.numeric(estvar[,i])))/1000, cex=1.2, col= colpl[i]

  )})

  abline(v=k+0.5)

}

### End
