#' @title lfqvisualize
#' @details Draw heatmap and save as ProteoLFQ-OUTPUT.Figure-Top.XXX.workflows.pdf. For function definitions and descriptions please use "??ProteoLFQ" command in R.
#' @param object The input is the output file of the "lfq_access" or "lfq_spiked".
#' @param top The default 'top' value is 100. You can view the top ranking heatmap you want.
#' @rawNamespace import(gplots, except=lowess)
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @usage lfqvisualize(object, top = 100)
#' @export lfqvisualize
#' @examples
#' library(EVALFQ)
#' \donttest{data_q <- PrepareInuputFiles(acquisitionmethods=2,
#' rawdataset = "MaxQuant_proteinGroups_LFQ.txt", lable = "MaxQuant_LFQ_Label.txt")}
#' \donttest{lfqevalueall(data_q = data_q,
#' assum_a="Y", assum_b="Y", assum_c="Y", Ca="1", Cb="1", Cc="1", Cd="1")}
#' \donttest{lfqvisualize(object = "EVALFQ-OUTPUT.Data-Overall.Ranking.csv", top = 100)}

lfqvisualize <- function(object, top=100){

  message("The default 'top' value is 100.")
  allranks1 <- read.csv(file = object, header = TRUE)
  x <- allranks1
  if(length(grep("Rank.OverallRank",colnames(x)))==1)
  {
    ID<-which(as.numeric(as.character(x[,"Rank.OverallRank"]))>top)
    x<-x[-ID,]
  }else{
    x<-x
  }

  hmcols<-colorRampPalette(c("#116A39","white","#F6FBBA"))
  groups <- factor(x[, 1], levels = unique(x[, 1]))
  unique.groups <- levels(groups)
  cols<-colorRampPalette(c("#116A39","white","#F6FBBA"))(length(unique.groups))
  ColSideColors <- c(rep(NA, length(rownames(x))))
  for (ii in 1:length(x[, 1])) {
    selected <- which(unique.groups == x[, 1][ii])
    ColSideColors[ii] <- cols[selected]
  }

  x[,1]<-rank(x[,1],ties.method="first")
  x1 <- x[,c(2:5)]
  x1 <- as.matrix(x1)
  options(show.error.messages=FALSE)
  options(warn=-1)

  if(top<500){

    pdf(file=paste("./ProteoLFQ-OUTPUT.Figure-Top",top,"workflows.pdf",sep="."),width=20,height=56)
    heatmap.2(x1,Colv=NA,Rowv=NA,trace="none",col=hmcols,scale= "none",dendrogram="none",main=paste("ProteoLFQ-OUTPUT.Figure-Top",top,"workflows",sep="."),
              RowSideColors =ColSideColors,key.title=NA,symkey=FALSE,symbreaks=FALSE,keysize=0.1,
              labRow=paste(x[,1]," ",rownames(x),sep=""),labCol=colnames(x1),
              lwid=c(0.5, 4, 0.2), offsetRow=-151,margins = c(21,0.2),cexRow = 1.3,colsep=c(1:ncol(x)),rowsep=c(1:nrow(x)),sepwidth = c(0.01,0))
    dev.off()
    message("The 'figure' has been successfully saved in the current path.")
    message("Please use 'getwd()' to find the current path!")
  }else{

    pdf(file=paste("./ProteoLFQ-OUTPUT.Figure-Top",top,"workflows.pdf",sep=" "),width=20,height=56)
    heatmap.2(x1,Colv=NA,Rowv=NA,trace="none",col=hmcols,scale= "none",dendrogram="none",main=paste("ProteoLFQ-OUTPUT.Figure-Top",top,"workflows",sep="."),
              RowSideColors =ColSideColors,key.title=NA,symkey=FALSE,symbreaks=FALSE,keysize=0.1,
              labRow=paste(x[,1]," ",rownames(x),sep=""),labCol=colnames(x1),
              lwid=c(0.5, 4, 0.2), offsetRow=-151,margins = c(21,0.2),cexRow = 1,colsep=c(1:ncol(x)),rowsep=c(1:nrow(x)),sepwidth = c(0.01,0))
    dev.off()

    message("The 'figure' has been successfully saved in the current path.")
    message("Please use 'getwd()' to find the current path!")
  }}
