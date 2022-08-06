`summary.ROTS` <-
function(object, fdr=NULL, num.genes=NULL, verbose=TRUE, ...){
  
  if( !is.null(fdr) || !is.null(num.genes)){
    ## Sort ROTS-statistic values (abs)
    sorted.rots <- abs(object$d)
    names(sorted.rots) <- rownames(object$data)
    sorted.rots <- sort(sorted.rots, decreasing=TRUE)

    ## Add feature-names to fdr and ROTS values
    names(object$FDR) <- rownames(object$data)
    names(object$d) <- rownames(object$data)
	names(object$pvalue) <- rownames(object$data)

    ## Result matrix, columns are: "Row", "ROTS-statistic", "pvalue" and "FDR" 
    result <- numeric(0)
    result <- cbind(match(names(sorted.rots), rownames(object$data)),
                    object$d[names(sorted.rots)],
					object$pvalue[names(sorted.rots)],
                    object$FDR[names(sorted.rots)])
    colnames( result ) <- c("Row", "ROTS-statistic", "pvalue", "FDR")

    ## Show only num.gene top rows or rows whose false discovery rate <= fdr
    if(!is.null(fdr))
      result <- result[ result[,4] <= fdr, ,drop=FALSE]
    else
      result <- result[ 1 : min(num.genes, nrow(result)), ,drop=FALSE]

    if( verbose ){
      cat("ROTS results:", "\n\n")
      cat("Number of resamplings: ", object$B, "\n\n")
      cat("a1:                    ", object$a1,"\n")
      cat("a2:                    ", object$a2,"\n")
      cat("Top list size:         ", object$k, "\n")
      cat("Reproducibility value: ", object$R, "\n")
      cat("Z-score:               ", object$Z, "\n\n")
      
      cat(nrow(result), "rows satisfy the condition.")
      ## Print oly 10 first rows
      if(nrow(result) > 10){
        cat(" Only ten first rows are" ,"\n",
            "displayed, see the return value for the whole output.\n")
        print(result[1:10, ])
        cat("...", "\n")
      }
      else{
        cat("\n")
        print(result)
      }
    }
    
    return(invisible(result))
  }
}
