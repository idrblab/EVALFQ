#' @title Conduct LFQ and assess performance of all possible LFQ workflows.
#' @description The EVALFQ enables the label-free quantification of proteomic
#' data and the performance assessment of each LFQ workflow from multiple perspectives.  Moreover, it
#' provides the unique function of ranking all possible LFQ workflows (>3,000 random combinations of
#' transformation, normalization and imputation methods) based on their performances.
#' All in all, this tool makes the performance assessment of whole LFQ workflow possible
#' (collectively assessed by five well-established criteria with distinct underlying theories) and
#' gives the ranking results of all possible workflows based on the criteria preferred and selected by the users.
#' For function definitions and descriptions please use "??EVALFQ" command in R.
#' @param data_q This input file should be numeric type except the first and second column containing the names and label (control or case) of the studied samples, respectively. The intensity data should be provided in this input file with the following order: samples in row and proteins/peptides in column. Missing value (NA) of protein intensity are allowed.
#' @param assum_a all proteins were assumed to be equally important.The authors will be asked to input a letter “Y” to indicate the corresponding assumption is held for the studied dataset and a letter “N” to denote the opposite.
#' @param assum_b The level of protein abundance was assumed to be constant among all samples. The authors will be asked to input a letter “Y” to indicate the corresponding assumption is held for the studied dataset and a letter “N” to denote the opposite.
#' @param assum_c The intensities of the vast majority of the proteins were assumed to be unchanged under the studied conditions. The authors will be asked to input a letter “Y” to indicate the corresponding assumption is held for the studied dataset and a letter “N” to denote the opposite.
#' @param Ca Criterion (a): precision of LFQ based on the proteomes among replicates (Proteomics. 15:3140-51, 2015). If set 1, the user chooses to assess LFQ workflows using Criterion (a). If set 0, the user excludes Criterion (a) from performance assessment. The default setting of this value is “1”.
#' @param Cb Criterion (b): classification ability of LFQ between distinct sample groups (Nat Biotechnol. 28:83-9, 2010). If set 1, the user chooses to assess LFQ workflows using Criterion (b). If set 0, the user excludes Criterion (b) from performance assessment. The default setting of this value is “1”.
#' @param Cc Criterion (c): differential expression analysis by reproducibility-optimization (Nat Biotechnol. 32:896-902, 2014). If set 1, the user chooses to assess LFQ workflows using Criterion (c). If set 0, the user excludes Criterion (c) from performance assessment. The default setting of this value is “1”.
#' @param Cd Criterion (d): reproducibility of the identified protein markers among different datasets (Mol Biosyst. 11:1235-40, 2015). If set 1, the user chooses to assess LFQ workflows using Criterion (d). If set 0, the user excludes Criterion (d) from performance assessment. The default setting of this value is “1”.
#' @return preprocessed matrix
#' @import utils stats
#' @import metabolomics
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid grid.draw
#' @import futile.logger
#' @import affy vsn
#' @import MASS limma
#' @import ProteoMM ROTS
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @useDynLib EVALFQ
#' @importFrom Rcpp sourceCpp
#' @rawNamespace import(dplyr, except=c(filter,lag,select,combine))
#' @rawNamespace import(gplots, except=lowess)
#' @importFrom pcaMethods pca
#' @importFrom pcaMethods completeObs
#' @import impute
#' @usage lfqevalueall(data_q, assum_a="Y", assum_b="Y",
#' assum_c="Y", Ca="1", Cb="1", Cc="1", Cd="1")
#' @export lfqevalueall
#' @examples
#' library(EVALFQ)
#' \donttest{data_q <- PrepareInuputFiles(acquisitionmethods=2,
#' rawdataset = "MaxQuant_proteinGroups_LFQ.txt", lable = "MaxQuant_LFQ_Label.txt")}
#' \donttest{lfqevalueall(data_q = data_q,
#' assum_a="Y", assum_b="Y", assum_c="Y", Ca="1", Cb="1", Cc="1", Cd="1")}


lfqevalueall <- function(data_q, assum_a="Y", assum_b="Y", assum_c="Y", Ca="1", Cb="1", Cc="1", Cd="1"){

  cat("\n")
  cat("EVALFQ is Running ...","\n")
  cat("\n")

  cat("*************************************************************************","\n")
  cat("Depending on the size of your input dataset","\n")
  cat("Several mintues or hours may be needed for this assessment","\n")
  cat("*************************************************************************","\n")
  cat("\n")

  cat("STEP 1: Prepare input file in standard formats of EVALFQ", "\n")
  cat("\n")

  cat("STEP 2: The assumption(s) held as indicated by users","\n")
  cat("Assumption A: all proteins were equally important (Y/N): ", assum_a, "\n")
  cat("Assumption B: the level of protein abundance was constant among all samples (Y/N): ", assum_b, "\n")
  cat("Assumption C: the intensity of the majority of proteins were unchanged (Y/N): ", assum_c, "\n")

  #YNN
  if(any(assum_a == "Y") && any(assum_b == "N") && any(assum_c == "N")){
    cat("According to your selection, there are 248 methods are included and 2880 methods are excluded", "\n")
  }
  #NYN
  if(any(assum_a == "N") && any(assum_b == "Y") && any(assum_c == "N")){
    cat("According to your selection, there are 72 methods are included and 3056 methods are excluded", "\n")
  }
  #NNY
  if(any(assum_a == "N") && any(assum_b == "N") && any(assum_c == "Y")){
    cat("According to your selection, there are 136 methods are included and 2992 methods are excluded", "\n")
  }
  #YYN
  if(any(assum_a == "Y") && any(assum_b == "Y") && any(assum_c == "N")){
    cat("According to your selection, there are 968 methods are included and 2160 methods are excluded", "\n")
  }
  #YNY
  if(any(assum_a == "Y") && any(assum_b == "N") && any(assum_c == "Y")){
    cat("According to your selection, there are 1928 methods are included and 1200 methods are excluded", "\n")
  }
  #NYY
  if(any(assum_a == "N") && any(assum_b == "Y") && any(assum_c == "Y")){
    cat("According to your selection, there are 216 methods are included and 2912 methods are excluded", "\n")
  }
  #YYY
  if(any(assum_a == "Y") && any(assum_b == "Y") && any(assum_c == "Y")){
    cat("According to your selection, there are 3128 methods are included and 0 methods are excluded", "\n")
  }

  cat("\n")
  cat("STEP 3: The criteira selected by users for this assessment","\n")
  cat("Criterion A: Precision (1/0): ", Ca, "\n")
  cat("Criterion B: Classification Ability (1/0): ", Cb, "\n")
  cat("Criterion C: Differential Expression (1/0): ", Cc, "\n")
  cat("Criterion D: Reproducibility (1/0): ", Cd, "\n")
  cat("\n")

  cat("EVALFQ is Running ...","\n")
  cat("\n")

  #### 样本在行，特征在列
  trans <- function(data,n){
    matrix <- switch(
      n,
      Box_Cox(data),
      log2(data),
      data
    )
    return(matrix)
  }

  cents <- function(data,n){
    matrix <- switch(
      n,
      MEC(data),
      MDC(data),
      data
    )
    return(matrix)
  }

  scals<- function(data,n){
    matrix <- switch(
      n,
      1,
      AUTO1(data),
      PARETO1(data),
      VAST1(data),
      RANGE1(data)
    )
    return(matrix)
  }

  norm <- function(data,n){
    matrix <- switch(
      n,
      t(fastlo(as.matrix(data))),
      t(EIGENMS(data, label)),
      t(LINEAR(data)),
      t(LOWESS(data)),
      t(SMAD(data)),
      t(MEAN(data)),
      t(MEDIAN(data)),
      t(data),
      t(PQN(data)),
      t(QUANTILE(as.matrix(data))),
      t(RLR1(data)),
      t(MSTUS(data)),
      t(TMM(data)),
      t(VSN(as.matrix(data)))
    )
    return(matrix)
  }

  impute <- function(data,n){
    matrix <- switch(
      n,
      filter_train_data,
      t(back(filter_train_data)),
      t(bpca(filter_train_data,nPcs=3)),
      t(censor(filter_train_data)),
      t(knn(filter_train_data,k=10)),
      t(lls(filter_train_data,k=10)),
      t(svdm(filter_train_data,nPcs=3)),
      t(zero(filter_train_data))
    )
    return(matrix)
  }

  consistency <-  function(fold = 5, top = 20) {
    folds <- fold
    control.label <- control.y # variable-1
    test.fold1 <- split(sample(1:length(control.label)), 1:folds) #ignore warning
    case.label <- case.y # variable-2
    test.fold2 <- split(sample(1:length(case.label)), 1:folds) #ignore warning

    DEG <- list()
    for (i in 1:folds) {
      com.x <- cbind(control.x[, test.fold1[[i]]], case.x[, test.fold2[[i]]]) # variable-3 & 4.
      lab.ct <- test.fold1[[i]]
      lab.ca <- test.fold2[[i]]
      design <- cbind(Grp1 = 1, Grp2vs1 = c(rep(0, length(lab.ct)), rep(1, length(lab.ca))))
      fit <- limma::lmFit(com.x, design)
      fit <- limma::eBayes(fit)
      DEG[[i]] <- rownames(limma::topTable(fit, coef = 2, number = nrow(com.x)))
    }
    names(DEG) <- LETTERS[1:folds]

    top.n <-top # Extracting the top n genes.
    DEG.list <- DEG
    for (g in 1:length(DEG.list)) {
      DEG.list[[g]] <- DEG.list[[g]][1:top.n]
    }

    # Calculating consistency score:

    setlist <- DEG.list
    OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
    con.score <- 0
    VennList <- OLlist$Venn_List
    for (i in 1:length(VennList)) {
      insect.n <- nchar(names(VennList[i]))
      if (insect.n < 2) next
      num.i <- 2^(insect.n - 2) * length(VennList[[i]])
      con.score <- con.score + num.i
    }

    return(con.score) # consistense score

  }

  # Stable consistense score, with 20 repeats

  stabel.score <- function(repeats = 20, fold = 5, top = 10) {
    score <- 0
    for (r in 1:repeats) {
      score <- score + consistency(fold, top)
    }
    return(score/repeats)
  }

  #####################################################################

  iName<-c("BOX","LOG","NON")

  oName<-c("MEC","MDC","NON")

  pName<-c("NON","ATO","PAR","RAN","VAS")

  jName<-c("CYC","EIG","LIN","LOW","MAD","MEA","MED","NON","PQN","QUA","RLR","TIC","TMM","VSN")

  gName<-c("NON","BAK","BPC","CEN","KNN","LLS","SVD","ZER")

  dataa<-data_q
  rownames(dataa)<-dataa[,1]
  label <- dataa[,2]
  frame <- dataa[, -(1:2)]
  frame <- t(frame)
  frame <- data.frame(frame)
  train_data <- data.matrix(frame, rownames.force = NA)

  train_data[train_data == 0] <- NA
  #train_data <- as.matrix(train_data)
  data2 <- train_data
  col_f <- apply(data2, 2, function(x) length(which(is.na(x)))/length(x))

  if (length(which(col_f >0.2))==0){
    data2_f <- data2
  }else {
    data2_f <- data2[, -which(col_f > 0.2)]
  }

  train_data_t <- data2_f

  Fpcv <- list()
  Fscore <- list()
  Faccuracy <- list()
  Fpmad <- list()

  Fbar1num <- list()
  Fbaronummean <- list()
  Fbarosd <- list()
  Fbarorsd <- list()
  Fbaro_rsdr_to_bar1num <- list()
  spike <- list()
  backgound <- list()

  #time <- 0

  center <- c(1,2,3)
  scaling <- c(1,2,3,4,5)
  normalization <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)

  if(is.na(match("Y",assum_a))){
    center<-3
    scaling<-1
  }

  if(is.na(match("Y",assum_b)) & is.na(match("Y",assum_c))) {
    normalization <- c(8)
  }

  if(is.na(match("Y",assum_b)) & !is.na(match("Y",assum_c))){
    normalization <- c(1,2,3,4,9,10,11,13)
  }

  if(!is.na(match("Y",assum_b)) & is.na(match("Y",assum_c))){
    normalization <- c(5,6,7,12)
  }

  newname=c()
  newname1=c()
  for(ii in 1:3){
    if(ii!=3){
      for(o in  center){
        for(p in  scaling){
          for(j in normalization){
            for(g in 1:8){
              namebind <- paste(iName[ii],oName[o],pName[p],jName[j],gName[g],sep="+")
              namebind1 <- paste(ii,o,p,j,g,sep="")
              newname1=rbind(newname1,namebind1)
              newname=rbind(newname,namebind)
            }
          }
        }
      }
    }else{
      for(o in  3){
        for(p in  1){
          for(j in 14){
            for(g in 1:8){
              namebind <- paste(iName[ii],oName[o],pName[p],jName[j],gName[g],sep="+")
              namebind1 <- paste(ii,o,p,j,g,sep="")
              newname1=rbind(newname1,namebind1)
              newname=rbind(newname,namebind)
            }
          }
        }
      }
    }
  }
  newnameall <- cbind(newname,newname1)
  rownames(newnameall) <- NULL
  colnames(newnameall) <- c("NAME","NUMBER")
  write.csv(newnameall,file = "InlusionExclusionList.csv")

  dir.create(paste0("Criteria.A-Precision"))
  dir.create(paste0("Criteria.B-Classification.Ability"))
  dir.create(paste0("Criteria.C-Differential.Expression"))
  dir.create(paste0("Criteria.D-Reproducibility"))

  for(ii in 1:3){

    train_data_tran <- try(trans(train_data_t,ii))
    if(class(train_data_tran) == "try-error")
    { next }

    if(ii!=3){

      for(o in  center){

        tran_train_data <- train_data_tran
        tran_train_data[is.infinite(data.matrix(tran_train_data))] <- NA
        cen_train_data <- try(cents(tran_train_data,o))

        if(class(cen_train_data) == "try-error")
        { next }

        for(p in  scaling){

          scal_factor <- try(scals(tran_train_data,p))

          if(class(scal_factor) == "try-error")
          { next }

          scal_train_data <- cen_train_data/scal_factor

          if(class(scal_train_data)=="try-error")
          { next }

          for(j in normalization){

            scal_train_data[is.nan(scal_train_data)] <- NA
            scal_train_data[is.infinite(scal_train_data)] <- NA
            normalized_data <- try(norm(scal_train_data,j))

            if(class(normalized_data)=="try-error")
            { next }

            label_c <- as.factor(label)
            g1 <- table(label_c)[levels(label_c)[1]]*0.8
            g2 <- table(label_c)[levels(label_c)[2]]*0.8

            train_data_filtering <- try(Basicfilter(normalized_data, label, g1=2, g2=2))

            if(class(train_data_filtering) == "try-error")
            { next }

            filter_train_data <- train_data_filtering

            for(g in 1:8){

              if(g==1){

                imputed_data <- filter_train_data
              }

              if(g!=1){

                imputed_data <- try(impute(filter_train_data,g))

                if(class(imputed_data)=="try-error")
                { next }
              }

              ##### count number
              namebind11 <- paste(iName[ii],oName[o],pName[p],jName[j],gName[g],sep="+")
              id <- which(newname==namebind11, arr.ind = TRUE)
              id <- as.data.frame(id)
              id1 <- id$row

              ##### Feature Selection
              #time = time+1
              dataa <- data_q
              frame <- imputed_data
              label <- dataa[,2]

              rots.out <-try(ROTS(data = t(frame), groups = as.character(label), B = 200, K = 500 , seed = 1234, log = FALSE))

              if(class(rots.out)=="try-error")
              { next }

              frame <- imputed_data
              label <- as.factor(as.character(label))
              im.data <- data.frame(label=label,frame)

              #(a) Precision of LFQ Based on the Proteomes among Replicates

              if( Ca == 1 ){

                #dir.create(paste0("Criteria.A-Precision"))

                cat(paste("Assessing" , paste(id1,"/",length(newname),":",sep=""), paste(iName[ii],"-",oName[o],"-",pName[p],"-", jName[j], "-",gName[g],sep=""),"Under Criteria A: Precision"),"\n")

                pdf(file=paste("./Criteria.A-Precision/Criteria.A_",ii,o,p,j,g,".pdf",sep=""),width = 25,height = 12,onefile = T)

                data <- im.data
                result <- PCV1(data)
                try(plotrCV(k=2,result))
                pcv <- sapply(1:3, function(i){round(1000*mean(as.numeric(result[,i])))/1000})[3]

                pmad <- try(mean(PMAD(data)))
                if(class(pmad) == "try-error")
                { next }

                Fpmad[paste(ii,o,p,j,g,sep="")] <- pmad
                dev.off()
              }else{
                message("'Criteria A: Precision' cannot be evaluated, Please Check!")
              }

              #(b) Classification Ability of LFQ between Distinct Sample Groups
              if( Cb == 1 ){

                #dir.create(paste0("Criteria.B-Classification.Ability"))

                cat(paste("Assessing" , paste(id1,"/",length(newname),":",sep=""), paste(iName[ii],"-",oName[o],"-",pName[p],"-", jName[j], "-",gName[g],sep=""),"Under Criteria B: Classification.Ability"),"\n")

                pdf(file=paste("./Criteria.B-Classification.Ability/Criteria.B_",ii,o,p,j,g,".pdf",sep=""),width = 25,height = 12,onefile = T)

                data <- im.data
                rots.out <- rots.out
                col_pos <- which(rots.out$FDR<0.05)
                if(length(col_pos) <= 10) {
                  markerid <-order(rots.out$pvalue)[1:20]
                }else {
                  markerid <- which(rots.out$FDR<0.05)+1
                }

                #markerid<-which(rots.out$pvalue<0.05)+1
                heat<-try( HeatMap(data[,c(1,markerid)], margins = c(6, 6),
                                   dendrogram="both",main = "",
                                   key=FALSE))
                if(class(heat)=="try-error")
                { next }
                dev.off()

                clusters <- hclust(dist(data[,markerid]))
                clusterCut <- cutree(clusters, 2)
                dataa <- data_q
                label <- dataa[,2]
                tmatrix <- table(clusterCut, label)
                tru <- as.numeric(data[,1])
                accuracy <- (tmatrix[1,1]+ tmatrix[2,2])/length(label)
                Faccuracy[paste(ii,o,p,j,g,sep="")] <- accuracy

              }else{
                message("'Criteria B: Classification.Ability' cannot be evaluated, Please Check!")
              }

              #(c) Differential Expression Analysis Based on Reproducibility-optimization

              if( Cc == 1 ){

                #dir.create(paste0("Criteria.C-Differential.Expression"))

                cat(paste("Assessing" ,paste(id1,"/",length(newname),":",sep=""), paste(iName[ii],"-", oName[o],"-",pName[p],"-",jName[j], "-",gName[g],sep=""),"Under Criteria C: Differential.Expression"),"\n")

                pdf(file=paste("./Criteria.C-Differential.Expression/Criteria.C_",ii,o,p,j,g,".pdf",sep=""),width = 22,height = 12,onefile = T)

                rots.out <- rots.out
                try(plot(rots.out, type="pvalue"))
                dev.off()

                breaks <- seq(0,1,0.05)
                sdres <- affy::hist(rots.out$pvalue,breaks=breaks)
                bar1num <- sdres$counts[1]
                baronummean <- mean(sdres$counts[-1])
                barosd <- sd(sdres$counts[-1])
                barorsd <- barosd/baronummean

                Fbar1num[paste(ii,o,p,j,g,sep="")] <- bar1num
                Fbaronummean[paste(ii,o,p,j,g,sep="")] <- baronummean
                Fbarosd[paste(ii,o,p,j,g,sep="")] <- barosd
                Fbarorsd[paste(ii,o,p,j,g,sep="")] <- barorsd
                Fbaro_rsdr_to_bar1num[paste(ii,o,p,j,g,sep="")] <- barorsd/bar1num

              }else{
                message("'Criteria C: Differential.Expression' cannot be evaluated, Please Check!")
              }

              #(d) Reproducibility of the Identified Protein Markers among Different Datasets
              if( Cd == 1 && length(label) >= 20 ){

                #dir.create(paste0("Criteria.D-Reproducibility"))

                cat(paste("Assessing" , paste(id1,"/",length(newname),":",sep=""), paste(iName[ii],"-", oName[o],"-",pName[p],"-",jName[j], "-",gName[g],sep=""),"Under Criteria D: Reproducibility"),"\n")

                test_data <- imputed_data
                label.vector <- names(table(label))
                control.x <- as.data.frame(t(test_data[label == label.vector[1], -1]))
                case.x <- as.data.frame(t(test_data[label == label.vector[2], -1]))
                control.y <- rep(0, table(label)[1])
                case.y <- rep(1, table(label)[2])

                score <- try(stabel.score(repeats = 200, fold = 5, top = 20))
                pdf(file=paste("./Criteria.D-Reproducibility/Criteria.D_",ii,o,p,j,g,".pdf",sep=""),width = 22,height = 12,onefile = T)

                sink(file=paste("./Criteria.D-Reproducibility/Criteria.D_",ii,o,p,j,g,".txt",sep=""))
                print(score)
                sink()

                Fscore[paste(ii,o,p,j,g,sep="")] <- score

                fold = 5
                top = 20
                folds <- fold
                control.label <- control.y # variable-1
                test.fold1 <- split(sample(1:length(control.label)), 1:folds) #ignore warning
                case.label <- case.y # variable-2
                test.fold2 <- split(sample(1:length(case.label)), 1:folds) #ignore warning

                DEG <- list()
                for (i in 1:folds) {
                  com.x <- cbind(control.x[, test.fold1[[i]]], case.x[, test.fold2[[i]]]) # variable-3 & 4.
                  lab.ct <- test.fold1[[i]]
                  lab.ca <- test.fold2[[i]]
                  design <- cbind(Grp1 = 1, Grp2vs1 = c(rep(0, length(lab.ct)), rep(1, length(lab.ca))))
                  fit <- lmFit(com.x, design)
                  fit <- eBayes(fit)
                  DEG[[i]] <- rownames(topTable(fit, coef = 2, number = nrow(com.x)))
                }
                names(DEG) <- LETTERS[1:folds]

                top.n <-top # Extracting the top n genes.
                DEG.list <- DEG
                for (g in 1:length(DEG.list)) {
                  DEG.list[[g]] <- DEG.list[[g]][1:top.n]
                }

                #setwd(tempdir())
                DEG3.list <- DEG.list
                #save(DEG3.list, file = "DEG3.list.RData")
                cols <- rainbow(length(DEG.list))
                if (length(DEG.list) == 3) {

                  setlist3 <- DEG.list;
                  #OLlist3 <- overLapper(setlist = setlist3, sep="_", type = "vennsets")
                  #counts <- list(sapply(OLlist3$Venn_List, length), sapply(OLlist3$Venn_List, length))
                  #par(mar=c(2,2,3,0.2))
                  #vennPlot(counts = counts[[1]], yoffset = c(0.3, -0.2))
                  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
                  venn.plot <- try(VennDiagram::venn.diagram(x = DEG.list, filename = NULL, col=cols))
                  try(grid::grid.draw(venn.plot))


                } else if (length(DEG.list) == 4) {

                  setlist4 <- DEG.list
                  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
                  venn.plot <- try(VennDiagram::venn.diagram(x = DEG.list, filename = NULL, col=cols))

                  #pdf(file=paste("./Criteria.D-Reproducibility/Criteria.D_",i,o,p,j,g,".pdf",sep=""),width = 22,height = 12,onefile = T)
                  try(grid::grid.draw(venn.plot))
                  #dev.off()
                  #OLlist4 <- overLapper(setlist=setlist4, sep="_", type="vennsets")
                  #counts <- list(sapply(OLlist4$Venn_List, length), sapply(OLlist4$Venn_List, length))
                  #par(mar=c(2,2,3,0.2))
                  #vennPlot(counts=counts[[1]], yoffset=c(0.3, -0.2))

                } else {

                  setlist5 <- DEG.list
                  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
                  venn.plot <- try(VennDiagram::venn.diagram(x = DEG.list, filename = NULL, col=cols))

                  #pdf(file=paste("./Criteria.D-Reproducibility/Criteria.D_",i,o,p,j,g,".pdf",sep=""),width = 22,height = 12,onefile = T)
                  try(grid::grid.draw(venn.plot))

                  #OLlist5 <- overLapper(setlist=setlist5, sep="_", type="vennsets")
                  #counts <- sapply(OLlist5$Venn_List, length)
                  #par(mar=c(2,2,3,0.2))
                  #vennPlot(counts=counts, ccol=c(rep(1,30),2), lcex=1.5, ccex=c(rep(1.5,5), rep(0.6,25),1.5)) # Plots a non-proportional 5-way Venn diagram.

                }
                dev.off()

              }else{
                message("'Criteria.D-Reproducibility' cannot be evaluated, Please Check!")
              }

            }

          }

        }
      }

    }else{

      for(o in  3){

        tran_train_data <- train_data_tran
        tran_train_data[is.infinite(data.matrix(tran_train_data))] <- NA
        cen_train_data <- try(cents(tran_train_data,o))

        if(class(cen_train_data) == "try-error")
        { next }

        for(p in  1){

          scal_factor <- try(scals(tran_train_data,p))

          if(class(scal_factor)=="try-error")
          { next }

          scal_train_data <- cen_train_data/scal_factor

          if(class(scal_train_data) == "try-error")
          { next }

          for(j in 14){

            scal_train_data[is.nan(scal_train_data)] <- NA
            scal_train_data[is.infinite(scal_train_data)] <- NA
            normalized_data <- try(norm(scal_train_data,j))

            if(class(normalized_data) == "try-error")
            { next }

            label_c <- as.factor(label)
            g1 <- table(label_c)[levels(label_c)[1]]*0.8
            g2 <- table(label_c)[levels(label_c)[2]]*0.8

            train_data_filtering <- try(Basicfilter(normalized_data,label,g1=2,g2=2))

            if(class(train_data_filtering)=="try-error")
            { next }

            filter_train_data <- train_data_filtering

            for(g in 1:8){

              if(g==1){

                imputed_data <- filter_train_data
              }

              if(g!=1){

                imputed_data <- try(impute(filter_train_data,g))

                if(class(imputed_data) == "try-error")
                { next }
              }
              ##### count number
              namebind11 <- paste(iName[ii],oName[o],pName[p],jName[j],gName[g],sep="+")
              id <- which(newname==namebind11, arr.ind = TRUE)
              id <- as.data.frame(id)
              id1 <- id$row
              ##### Feature Selection
              #time = time+1
              dataa <- data_q
              frame <- imputed_data
              label <- dataa[,2]

              rots.out <- try(ROTS(data = t(frame), groups = as.character(label), B = 200, K = 500 , seed = 1234,log = FALSE))

              if(class(rots.out) == "try-error")
              { next }

              frame <- imputed_data
              label <- as.factor(as.character(label))
              im.data <- data.frame(label=label,frame)

              #(a) Precision of LFQ Based on the Proteomes among Replicates

              if( Ca == 1 ){
                #dir.create(paste0("Criteria.A-Precision"))

                cat(paste("Assessing" , paste(id1,"/",length(newname),":",sep=""), paste(iName[ii],"-",oName[o],"-",pName[p],"-", jName[j], "-",gName[g],sep=""),"Under Criteria A: Precision"),"\n")

                pdf(file=paste("./Criteria.A-Precision/Criteria.A_",ii,o,p,j,g,".pdf",sep=""),width = 25,height = 12,onefile = T)

                data <- im.data
                result <- PCV1(data)
                try(plotrCV(k=2,result))

                pcv <- sapply(1:3, function(i){round(1000*mean(as.numeric(result[,i])))/1000})[3]

                pmad <- try(mean(PMAD(data)))
                if(class(pmad) == "try-error")
                { next }

                Fpmad[paste(ii,o,p,j,g,sep="")] <- pmad
                dev.off()
              }else{
                message("'Criteria A: Precision' cannot be evaluated, Please Check!")
              }

              #(b) Classification Ability of LFQ between Distinct Sample Groups
              if( Cb == 1 ){
                #dir.create(paste0("Criteria.B-Classification.Ability"))

                cat(paste("Assessing" , paste(id1,"/",length(newname),":",sep=""), paste(iName[ii],"-",oName[o],"-",pName[p],"-", jName[j], "-",gName[g],sep=""),"Under Criteria B: Classification.Ability"),"\n")

                pdf(file=paste("./Criteria.B-Classification.Ability/Criteria.B_",ii,o,p,j,g,".pdf",sep=""),width = 25,height = 12,onefile = T)

                data <- im.data
                rots.out <- rots.out
                col_pos <- which(rots.out$FDR<0.05)
                if(length(col_pos) <= 10) {
                  markerid <- order(rots.out$pvalue)[1:20]
                }else {
                  markerid <- which(rots.out$FDR < 0.05)+1
                }

                #markerid<-which(rots.out$pvalue<0.05)+1
                heat<-try( HeatMap(data[,c(1,markerid)], margins = c(6, 6),
                                   dendrogram="both",main = "",
                                   key=FALSE))
                if(class(heat)=="try-error")
                { next }
                dev.off()

                clusters <- hclust(dist(data[,markerid]))
                clusterCut <- cutree(clusters, 2)
                dataa <- data_q
                label <- dataa[,2]
                tmatrix <- table(clusterCut, label)
                tru <- as.numeric(data[,1])
                accuracy <- (tmatrix[1,1]+ tmatrix[2,2])/length(label)
                Faccuracy[paste(ii,o,p,j,g,sep="")] <- accuracy

              }else{
                message("'Criteria B: Classification.Ability' cannot be evaluated, Please Check!")
              }

              #(c) Differential Expression Analysis Based on Reproducibility-optimization

              if( Cc == 1 ){

                #dir.create(paste0("Criteria.C-Differential.Expression"))

                cat(paste("Assessing" ,paste(id1,"/",length(newname),":",sep=""), paste(iName[ii],"-", oName[o],"-",pName[p],"-",jName[j], "-",gName[g],sep=""),"Under Criteria C: Differential.Expression"),"\n")

                pdf(file=paste("./Criteria.C-Differential.Expression/Criteria.C_",ii,o,p,j,g,".pdf",sep=""),width = 22,height = 12,onefile = T)

                rots.out <- rots.out
                try(plot(rots.out, type="pvalue"))
                dev.off()

                breaks <- seq(0,1,0.05)
                sdres <- affy::hist(rots.out$pvalue,breaks=breaks)
                bar1num <- sdres$counts[1]
                baronummean <- mean(sdres$counts[-1])
                barosd <- sd(sdres$counts[-1])
                barorsd <- barosd/baronummean

                Fbar1num[paste(ii,o,p,j,g,sep="")] <- bar1num
                Fbaronummean[paste(ii,o,p,j,g,sep="")] <- baronummean
                Fbarosd[paste(ii,o,p,j,g,sep="")] <- barosd
                Fbarorsd[paste(ii,o,p,j,g,sep="")] <- barorsd
                Fbaro_rsdr_to_bar1num[paste(ii,o,p,j,g,sep="")] <- barorsd/bar1num

              }else{
                message("'Criteria C: Differential.Expression' cannot be evaluated, Please Check!")
              }

              #(d) Reproducibility of the Identified Protein Markers among Different Datasets
              if( Cd == 1 && length(label) >= 20){

                #dir.create(paste0("Criteria.D-Reproducibility"))

                cat(paste("Assessing" , paste(id1,"/",length(newname),":",sep=""), paste(iName[ii],"-", oName[o],"-",pName[p],"-",jName[j], "-",gName[g],sep=""),"Under Criteria D: Reproducibility"),"\n")

                test_data <- imputed_data
                label.vector <- names(table(label))
                control.x <- as.data.frame(t(test_data[label == label.vector[1], -1]))
                case.x <- as.data.frame(t(test_data[label == label.vector[2], -1]))
                control.y <- rep(0, table(label)[1])
                case.y <- rep(1, table(label)[2])

                score <- try(stabel.score(repeats = 200, fold = 5, top = 20))

                pdf(file=paste("./Criteria.D-Reproducibility/Criteria.D_",ii,o,p,j,g,".pdf",sep=""),width = 22,height = 12,onefile = T)

                sink(file=paste("./Criteria.D-Reproducibility/Criteria.D_",ii,o,p,j,g,".txt",sep=""))
                print(score)
                sink()

                Fscore[paste(ii,o,p,j,g,sep="")] <- score

                fold = 5
                top = 20
                folds <- fold
                control.label <- control.y # variable-1
                test.fold1 <- split(sample(1:length(control.label)), 1:folds) #ignore warning
                case.label <- case.y # variable-2
                test.fold2 <- split(sample(1:length(case.label)), 1:folds) #ignore warning

                DEG <- list()
                for (i in 1:folds) {
                  com.x <- cbind(control.x[, test.fold1[[i]]], case.x[, test.fold2[[i]]]) # variable-3 & 4.
                  lab.ct <- test.fold1[[i]]
                  lab.ca <- test.fold2[[i]]
                  design <- cbind(Grp1 = 1, Grp2vs1 = c(rep(0, length(lab.ct)), rep(1, length(lab.ca))))
                  fit <- lmFit(com.x, design)
                  fit <- eBayes(fit)
                  DEG[[i]] <- rownames(topTable(fit, coef = 2, number = nrow(com.x)))
                }
                names(DEG) <- LETTERS[1:folds]

                top.n <-top # Extracting the top n genes.
                DEG.list <- DEG
                for (g in 1:length(DEG.list)) {
                  DEG.list[[g]] <- DEG.list[[g]][1:top.n]
                }

                #setwd(tempdir())
                DEG3.list <- DEG.list
                #save(DEG3.list, file = "DEG3.list.RData")
                cols <- rainbow(length(DEG.list))
                if (length(DEG.list) == 3) {

                  setlist3 <- DEG.list;
                  #OLlist3 <- overLapper(setlist = setlist3, sep="_", type = "vennsets")
                  #counts <- list(sapply(OLlist3$Venn_List, length), sapply(OLlist3$Venn_List, length))
                  #par(mar=c(2,2,3,0.2))
                  #vennPlot(counts = counts[[1]], yoffset = c(0.3, -0.2))
                  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
                  venn.plot <- try(venn.diagram(x = DEG.list, filename = NULL, col=cols))
                  try(grid::grid.draw(venn.plot))


                } else if (length(DEG.list) == 4) {

                  setlist4 <- DEG.list
                  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
                  venn.plot <- try(venn.diagram(x = DEG.list, filename = NULL, col=cols))

                  #pdf(file=paste("./Criteria.D-Reproducibility/Criteria.D_",i,o,p,j,g,".pdf",sep=""),width = 22,height = 12,onefile = T)
                  try(grid::grid.draw(venn.plot))
                  #dev.off()
                  #OLlist4 <- overLapper(setlist=setlist4, sep="_", type="vennsets")
                  #counts <- list(sapply(OLlist4$Venn_List, length), sapply(OLlist4$Venn_List, length))
                  #par(mar=c(2,2,3,0.2))
                  #vennPlot(counts=counts[[1]], yoffset=c(0.3, -0.2))

                } else {

                  setlist5 <- DEG.list
                  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
                  venn.plot <- try(venn.diagram(x = DEG.list, filename = NULL, col=cols))

                  #pdf(file=paste("./Criteria.D-Reproducibility/Criteria.D_",i,o,p,j,g,".pdf",sep=""),width = 22,height = 12,onefile = T)
                  try(grid::grid.draw(venn.plot))

                  #OLlist5 <- overLapper(setlist=setlist5, sep="_", type="vennsets")
                  #counts <- sapply(OLlist5$Venn_List, length)
                  #par(mar=c(2,2,3,0.2))
                  #vennPlot(counts=counts, ccol=c(rep(1,30),2), lcex=1.5, ccex=c(rep(1.5,5), rep(0.6,25),1.5)) # Plots a non-proportional 5-way Venn diagram.

                }
                dev.off()

              }else{
                message("'Criteria.D-Reproducibility' cannot be evaluated, Please Check!")
              }

            }
          }
        }
      }
    }
  }

  save(Fpmad,file="./Criteria.A-Precision.Rdata")
  save(Faccuracy,file="./Criteria.B-Classification.Ability.Rdata")
  save(Fbaro_rsdr_to_bar1num,file="./Criteria.C-Differential.Expression.Rdata")
  save(Fscore,file="./Criteria.D-Reproducibility.Rdata")

  ###################################################Step-2
  if(!is.null(unlist(Faccuracy))){
    Acc_revison <- unlist(Faccuracy)
    Acc_revison_ID <- which(Acc_revison<0.5)
    Acc_revison[Acc_revison_ID] <- 1-Acc_revison[Acc_revison_ID]
  }

  result <- dplyr::bind_rows("Precision"=unlist(Fpmad),"Classcification.Ability"=Acc_revison,"Differential.Expression"=unlist(Fbaro_rsdr_to_bar1num),"Reproducibility"=unlist(Fscore),.id = "id")
  result1 <- t(result)
  colnames(result1) <- result1["id",]
  result2 <- result1[-1,]
  result2 <- data.frame(result2, check.names=FALSE)
  rownames(result2) <- names[match(rownames(result2), names[,1]), 2]

  Rank <- apply(result2, 2, function(x){rank(as.numeric(as.character(x)), ties.method="min", na.last = "keep")})
  rownames(Rank) <- rownames(result2)

  if(length(grep("Reproducibility", colnames(Rank))) == 1){
    Rank[,"Reproducibility"] <- rank(-as.numeric(as.character(result2[,"Reproducibility"])), ties.method="min", na.last = "keep")
  }else{
    Rank <- Rank
  }

  if(length(grep("Classcification.Ability", colnames(Rank))) == 1){
    Rank[,"Classcification.Ability"] <- rank(-as.numeric(as.character(result2[,"Classcification.Ability"])),ties.method="min",na.last = "keep")
  }else{
    Rank <- Rank
  }

  Rank_revision <- apply(Rank, 2, function(x){x[is.na(x)] <- nrow(Rank);return(x)})

  Ranksum0 <- apply(Rank_revision, 1, sum)
  Rankres0 <- cbind("OverallRank" = Ranksum0,Rank_revision)

  zuihou0 <- cbind("Rank"=Rankres0,"Value"=result2)

  if(length(grep("Precision",colnames(result2))) == 1){
    ID <- which(as.numeric(as.character(result2[,"Precision"]))>0.7)
    AA <- zuihou0[-ID,]
    AA0 <- AA[order(AA[,1], decreasing = FALSE),]
    AA0[,1] <- rank(AA0[,1], ties.method = "min")
    BB <- zuihou0[ID,]
    BB0 <- BB[order(BB[,1], decreasing = FALSE),]
    BB0[,1] <- rank(BB0[,1], ties.method = "min") + max(AA0[,1])
    zuihou2 <- rbind(AA0,BB0)
  }else{
    zuihou1 <- zuihou0[order(Rankres0[,"OverallRank"], decreasing = FALSE),]
    zuihou1[,1] <- rank(zuihou1[,1], ties.method = "min")
    zuihou2 <- zuihou1
  }

  write.csv(zuihou2,file = "./EVALFQ-OUTPUT.Data-Overall.Ranking.csv")
  cat("\n")
  cat("*************************************************************************","\n")
  cat("Congratulations! Assessment Successfully Completed!","\n")
  cat("Thanks for Using EVALFQ. Wish to See You Soon ...","\n")
  cat("*************************************************************************","\n")
  cat("\n")
  #return(zuihou2)
  #return(zuihou2)
  ###############################################End##########################################################################
}
