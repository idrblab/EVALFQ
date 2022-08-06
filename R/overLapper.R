##########################################
## Intersect and Venn Diagram Functions ##
##########################################
## Author: Thomas Girke
## Last update: March 24, 2012
## Please note: an improved version of the functionalities provided by this script is available
## in the systemPipeR package: http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html

## Utilities:
## (1) Venn Intersects
##     Computation of Venn intersects among 2-20 or more sample sets using the typical
##     'only in' intersect logic of Venn comparisons, such as: objects present only in
##     set A, objects present only in the intersect of A & B, etc. Due to this restrictive
##     intersect logic, the combined Venn sets contain no duplicates.
## (2) Regular Intersects
##     Computation of regular intersects among 2-20 or more sample sets using the
##     following intersect logic: objects present in the intersect of A & B, objects present
##     in the intersect of A & B & C, etc. The approach results usually in many duplications
##     of objects among the intersect sets.
## (3) Graphical Utilities
##     - Venn diagrams of 2-5 sample sets.
##     - Bar plots for the results of Venn intersect and all intersect approaches derived
##       from many samples sets.
##
## Detailed instructions for using the functions of this script are available on this page:
##     http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn
##
## Revision history:
##     March 24, 2012: fixed substring problem in plotVenn function

#######################################
## Define Generic Intersect Function ##
#######################################
## Computation of (1) Venn Intersects and (2) Regular Intersects
overLapper <- function(setlist=setlist, complexity=1:length(setlist), sep="-", cleanup=FALSE, keepdups=FALSE, type) {
  ## Clean up of sample sets to minimize formatting issues
  if(cleanup==TRUE) {
    ## Set all characters to upper case
    setlist <- sapply(setlist, function(x) gsub("([A-Z])", "\\U\\1", x, perl=T, ignore.case=T))
    ## Remove leading and trailing spaces
    setlist <- sapply(setlist, function(x) gsub("^ {1,}| {1,}$", "", x, perl=T, ignore.case=T))
  }

  ## Append object counter to retain duplicates
  if(keepdups==TRUE) {
    dupCount <- function(setlist=setlist) {
      count <- table(setlist)
      paste(rep(names(count), count), unlist(sapply(count, function(x) seq(1, x))), sep=".")
    }
    mynames <- names(setlist)
    setlist <- lapply(setlist, function(x) dupCount(x)) # lapply necessary for numeric data!
    names(setlist) <- mynames
  }

  ## Create intersect matrix (removes duplicates!)
  setunion <- sort(unique(unlist(setlist)))
  setmatrix <- sapply(names(setlist), function(x) setunion %in% unique(setlist[[x]]))
  rownames(setmatrix) <- setunion
  storage.mode(setmatrix) <- "numeric"

  ## Create all possible sample combinations within requested complexity levels
  labels <- names(setlist)
  allcombl <- lapply(complexity, function(x) combn(labels, m=x, simplify=FALSE))
  allcombl <- unlist(allcombl, recursive=FALSE)
  complevels <- sapply(allcombl, length)

  ## Return intersect list for generated sample combinations
  if(type=="intersects") {
    OLlist <- sapply(seq(along=allcombl), function(x) setunion[rowSums(setmatrix[, rep(allcombl[[x]], 2)]) == 2 * length(allcombl[[x]])])
    names(OLlist) <- sapply(allcombl, paste, collapse=sep)
    return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Intersect_List=OLlist))
  }

  ## Return Venn intersect list for generated sample combinations
  if(type=="vennsets") {
    vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
      mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
      mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
      cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
      cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
      return(setunion[cond1 & cond2])
    }
    vennOLlist <- sapply(seq(along=allcombl), function(x) vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x))
    names(vennOLlist) <- sapply(allcombl, paste, collapse=sep)
    return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Venn_List=vennOLlist))
  }
}

###########################################
