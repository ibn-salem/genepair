#'
#' This is a collection of functions that work on pairs of genes. Thse gene
#' pairs consists usually of a data.frame that hold indices of genes in an
#' associate GRanges object (or strings correponding to names in a GRanges
#' object)
#'


#'get subset of gene pairs that are located on the same chromosome
#'
#'@param gp A data.frame like object holding gene pairs
#'@param genesGR A \code{\}link{GRanges}} object
#'@return A data.frame objec with subset of rows in \code{gp} with pair on the
#'  same chromosome.
#'@export
getCisPairs <- function(gp, genesGR){

  # get chromosomes of gene pairs
  c1 = as.character(seqnames(genesGR[gp[,1]]))
  c2 = as.character(seqnames(genesGR[gp[,2]]))

  # subset of gene pairs that are located on the same chromosome
  return(gp[c1==c2,])
}


#' Remove duplicated entries of gene pairs of the form A-B and B-A.
#' @export
uniquePair <- function(gp){

  # get string of sorted IDs as unique pair ID
  pairID = apply(apply(gp[,1:2], 1, sort), 2, paste, collapse="_")

  gp[!duplicated(pairID),]
}



#' Get only one pair per unique gene.
#' @export
uniquePairPerGene <- function(gp){

  seen = c()  # set off seen IDs
  uniqPairs = rbind() # initialize output pairs

  for (i in seq(nrow(gp))){
    g1 = as.character(gp[i,1])
    g2 = as.character(gp[i,2])

    # check if one of the pairs was seen before
    if (! ( (g1 %in% seen) | (g2 %in% seen) )){
      uniqPairs = rbind(uniqPairs, gp[i,])
    }
    # update seen ID set
    seen = c(seen, g1, g2)
  }

  return(uniqPairs)
}


#' Creates an ID for each gene pair by concatenating both gene names.
#'
getPairID <- function(gp) paste( gp[,1], gp[,2], sep="_" )


#' Get only one pair per unique gene by choosing the gene pair with highest
#' sequence similarity first.
#'
#' @param gp data.frame like object with gene pairs
#' @param similarty numeric vector of the same lenght as rows in \code{gp}
#' @return a data frame that is filtered in a way that each gene occures only in
#' on pair. The pairs with samller similarity values are keeped preferentially.
#' @export
uniquePairPerGeneBySim <- function(gp, similarity){

  # get maximal weight matching of the graph G induced by the paris of genes
  # with similarity as weight.
  # This command call the function form inside the python script
  matching = python.call( "getMaxWeightMatchingAsDict", gp[,1], gp[,2], similarity)

  # convert the matching to a data frame of gene pairs
  uniqPairs <- data.frame(names(matching), matching, stringsAsFactors=FALSE)

  # for all unique pairs get indices in the input set of pairs
  orgIDs = match(getPairID(uniqPairs), getPairID(gp))

  # add annotation columns from original data.frame
  uniqPairs = cbind(uniqPairs, gp[orgIDs,3:ncol(gp)])

  # set colom names to those of the input data.frame and delete column names
  names(uniqPairs) <- names(gp)
  row.names(uniqPairs) = NULL

  return(uniqPairs)
}


