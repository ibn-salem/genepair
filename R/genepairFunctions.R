
#' Get all possible gene pairs within a given distance range as data.frame
#'
#' @param genesGR A \code{\link{GRanges}} object.
#' @param maxDist single numeric value for maximal allowed distance between
#'   start position.
#' @param minDist single numeric value for minimal allowed distance between
#'   start postions.
#' @return A data.frame with
#' @export
#' @author Jonas Ibn-Salem and Dimitris Polychronopoulos
getAllCisPairs <- function(genesGR, maxDist=10^6, minDist=0){

  # calculate overlap all possible gene pairs within maxDist bp
  # posGR <- resize(genesGR, width=1, fix="start")
  # windowGR <- resize(posGR, 2*maxDist+1, fix="center")
  #
  # hits = findOverlaps(posGR, windowGR, ignore.strand=TRUE)
  #
  posGR <- resize(genesGR, width=1, fix="start")
  hits = findOverlaps(posGR,
                      maxgap=maxDist,
                      drop.redundant=TRUE,
                      drop.self=TRUE,
                      ignore.strand=TRUE)

  # convert hits objec to data.frame
  gp = as.data.frame(hits)
  names(gp) <- c("g1", "g2")

  # if no gene pairs left, return empty data.frame
  if (nrow(gp) == 0) return(gp)

  # sort columns according to index in genesGR
  firstSmaller <- gp[,1] <= gp[,2]

  first <- gp[cbind(1:nrow(gp), ifelse(firstSmaller, 1, 2))]
  second <- gp[cbind(1:nrow(gp), ifelse(firstSmaller, 2, 1))]

  gp[,1] <- first
  gp[,2] <- second

  # sort rows according to index in genesGR
  gp <- gp[order(gp[,1], gp[,2]), ]

  # add distance
  gp <- addPairDist(gp, genesGR)

  # remove pars with distance smaller than minDist or larger than maxDist
  gp <- gp[abs(gp$dist) >= minDist,]

  return(gp)

}


#' Get subset of gene pairs that are located on the same chromosome.
#'
#' @param gp A data.frame like object holding gene pairs.
#' @param genesGR A \code{\link{GRanges}} object.
#' @return A data.frame objec with subset of rows in \code{gp} with pair on the
#'  same chromosome.
#' @import GenomicRanges
#' @export
filterForCisPairs <- function(gp, genesGR){

  # get chromosomes of gene pairs
  c1 = as.character(seqnames(genesGR[gp[,1]]))
  c2 = as.character(seqnames(genesGR[gp[,2]]))

  # subset of gene pairs that are located on the same chromosome
  return(gp[c1 == c2,])
}


#' Add linear distance between gene paires or NA if on differnet chromosomes.
#'
#' Distance is caluclated based on start positon of genes in basepair and can be
#' negative, if second gene has a lower start postion than the first gene in
#' pairs.
#'
#' @import GenomicRanges
#' @export
#' @param gp data.frame like object with gene pairs. Assumes indexes or names of
#'   ranges in 'genesGR' in the first two columns.
#' @param genesGR A \code{\link{GRanges}} object.
addPairDist <- function(gp, genesGR, colname="dist", ignore.strand=FALSE){

  # get chromosomes of gene pairs
  c1 <- as.character(seqnames(genesGR[gp[,1]]))
  c2 <- as.character(seqnames(genesGR[gp[,2]]))
  sameChrom <- c1 == c2

  # get coordinate of start positon of genes
  s1 = start(resize(genesGR[gp[,1]], 1, ignore.strand = ignore.strand))
  s2 = start(resize(genesGR[gp[,2]], 1, ignore.strand = ignore.strand))

  # add a new column "dist" to the data.frame
  gp[, colname] = ifelse(sameChrom, s2 - s1, NA)
  return(gp)
}


#' Add same starnd information.
#'
#' @import GenomicRanges
#' @export
addSameStrand <- function(gp, genesGR, colname="sameStrand"){

  # check if both genes have strand information
  s1 = as.character(strand(genesGR[gp[,1]]))
  s2 = as.character(strand(genesGR[gp[,2]]))
  hasStrandInfo = s1 != "*" & s2 != "*"

  #check if they are equal
  sameStrand = s1==s2

  gp[, colname] = as.vector(ifelse(hasStrandInfo, sameStrand, NA))

  return(gp)
}


#' Remove duplicated entries of gene pairs of the form A-B and B-A.
#' @export
uniquePair <- function(gp){

  if(nrow(gp) == 0) return(gp)

  # get string of sorted IDs as unique pair ID
  pairID = getPairIDsorted(gp)

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
#' @export
getPairID <- function(gp) paste( gp[,1], gp[,2], sep="_" )


#' Returns a unique gene ID for each pair that is independent of pair
#' order (e.g. sorted)
#' @export
getPairIDsorted <- function(gp){

  firstSmaler <- gp[,1] <= gp[,2]

  first <- gp[cbind(1:nrow(gp), ifelse(firstSmaler, 1, 2))]
  second <- gp[cbind(1:nrow(gp), ifelse(firstSmaler, 2, 1))]

  id <- paste(
    first,
    second,
    sep="_"
  )

  return(id)
}


#' Get only one pair per unique gene by choosing the gene pair with highest
#' sequence similarity first.
#'
#' @param gp data.frame like object with gene pairs
#' @param similarty numeric vector of the same lenght as rows in \code{gp}
#' @return a data frame that is filtered in a way that each gene occures only in
#' on pair. The pairs with samller similarity values are keeped preferentially.
#' @export
uniquePairPerGeneBySim <- function(gp, similarity){

  # (1) sort by similarity while keeping track of the original order

  orgOrder <- 1:nrow(gp)
  sortedOrder <- order(similarity)
  sortedGP <- gp[sortedOrder,]

  # (2) filter with "seen" set of geens
  dupMat <- t(duplicated(t(as.matrix(sortedGP[,1:2])), MARGIN=0))

  uniq <- apply(dupMat, 1, function(x) !any(x))

  # debug print
  # cbind(sortedGP, dupMat, uniq)

  # (3) reorder to original order
  uniqPairs <- gp[uniq[sortedOrder],]


  return(uniqPairs)
}



#' Returns the percentage of gene pairs that are contained in another set of
#' gene pairs.
#'
#' @export
percentIncluded <- function(gPa, gPb){
  a <- getPairIDsorted(gPa)
  b <- getPairIDsorted(gPb)
  return( 100 * sum(a %in% b) / length(a))
}


#' Test if gene pairs are contained in another set of gene pairs.
#'
#' Compaires gene pair by combintion of IDs by using the
#' \code{\link{getPairIDsorted}} function.
#'
#' @param gp a data.frame with gene pairs.
#' @param other another data.frame with gene pairs.
#' @return a logical vector for each pair in \code{gp} indicating whether it is
#'   contained in \code{other} or not.
#' @export
containsGenePairs <- function(gp, other){

  # get unique ID for pairs in gp and other

  gpID <- getPairIDsorted(gp)
  otherID <- getPairIDsorted(other)

  inOther <- gpID %in% otherID
}


#' Test if gene pairs are non-overlapping. That is that two paird genes do not
#' overlap each other. in the genome (on the same strand).
#'
#' @import GenomicRanges
#' @export
nonOverlappingGenePairs <- function(gp, genesGR, useIDs=FALSE){

  genesHitDF <- as.data.frame(findOverlaps(genesGR, genesGR))
  ovl <- containsGenePairs(gp, negPairs=genesHitDF, gPidx=useIDs, nPidx=TRUE, gr=genesGR)
  return( ! ovl )

}


#' Add annotatino column of each gene to gene pairs
#'
#'@export
addGeneAnnotation <- function(gp, genesGR, colname){

  gp[, paste0(colname, "_g1")] <- mcols(genesGR[gp[,1]])[, colname]
  gp[, paste0(colname, "_g2")] <- mcols(genesGR[gp[,2]])[, colname]

  return(gp)

}


#' Make GRange object of range between start position of paired genes.
#'
#' It assumes that genes in pair are on the same chromosome.

#' @param gp A data.frame like object holding gene pairs.
#' @param genesGR A \code{\link{GRanges}} object.
#' @return A \code{\link{GRanges}} object with the range of start positions of
#'   paired ranges.
#'
#' @import GenomicRanges
#' @seealso \code{\link{getPairsAsGRL}}
#' @export
getPairAsGR <- function(gp, genesGR){

  # get chromosomes of gene pairs
  c1 <- seqnames(genesGR[gp[,1]])
  c2 <- seqnames(genesGR[gp[,2]])
  if(!all(c1 == c2)) stop("Gene pairs are not on the same chromosome.")

  chrom = c1
  s1 = start(genesGR[gp[,1]])
  s2 = start(genesGR[gp[,2]])

  up = apply(cbind(s1, s2), 1, min)
  down = apply(cbind(s1, s2), 1, max)

  outGR = GRanges(chrom, IRanges(up, down), seqinfo=seqinfo(genesGR))

  # add all columns of gp as annotation columns
  mcols(outGR) = gp

  return(outGR)
}

#' Make GRangesList object of range between start position of paired genes.
#'
#' Does NOT assume that genes in pair are on the same chromosome.

#' @param gp A data.frame like object holding gene pairs.
#' @param genesGR A \code{\link{GRanges}} object.
#' @return A \code{\link{GRangesList}} object with the range of start positions of
#'   paired ranges. GRanges between different chromosomes are empty.
#'
#' @import GenomicRanges
#' @seealso \code{\link{getPairsAsGRL}}
#' @export
getPairAsGRL <- function(gp, genesGR){

  # get chromosomes and starts of gene pairs
  c1 <- seqnames(genesGR[gp[,1]])
  c2 <- seqnames(genesGR[gp[,2]])

  s1 = start(genesGR[gp[,1]])
  s2 = start(genesGR[gp[,2]])

  # assign start / end coordinates depending on gene order
  up = apply(cbind(s1, s2), 1, min)
  down = apply(cbind(s1, s2), 1, max)

  pairGRL <- GRangesList(
    mapply(function(g1_chr, g1_start, g2_chr, g2_start)
    {
      # when on same chr combine
      if(g1_chr == g2_chr){
        return(GRanges(seqnames=c(g1_chr),
                       ranges=IRanges(start=g1_start, end=g2_start),
                       seqinfo=seqinfo(genesGR)))
      } else {
        return(GRanges())
      }
    },
    as.character(c1), up, as.character(c2), down)
  )
  return(pairGRL)
}



#' Add column to indicate that two query region overlap the same subset of
#' subject regions. If both query regions do not overlap any subject, FALSE is
#' returned.
#'
#' @import GenomicRanges
#' @export
addSubTADmode <- function(gp, tadGR, genesGR, colName="subTAD"){

  # compute overlap of all genes with TADs
  hitsDF <- as.data.frame(findOverlaps(genesGR, tadGR, type="within"))

  # get indices of genes in tssGR
  if ( !is.numeric(gp[,1]) ){
    idx1 <- match(gp[,1], names(genesGR))
    idx2 <- match(gp[,2], names(genesGR))
  }else{
    idx1 <- gp[,1]
    idx2 <- gp[,2]
  }

  # get for each gene the set of TADs that overlap
  set1 <- lapply(idx1, function(i) hitsDF[hitsDF[,1]==i, 2])
  set2 <- lapply(idx2, function(i) hitsDF[hitsDF[,1]==i, 2])

  # test if sets are equal for each gene in pair
  commonSubset <- mapply(setequal, set1, set2)

  # compute overlap with any common TAD
  genePairGR <- getPairAsGR(gp, genesGR)
  commonTAD <- countOverlaps(genePairGR, tadGR, type="within") >= 1

  # get combinations
  subTADmode <- NA
  subTADmode[!commonTAD & commonSubset] <- "no TAD"
  subTADmode[!commonTAD & !commonSubset] <- "diff TAD"
  subTADmode[commonTAD & !commonSubset] <- "diff sub TAD"
  subTADmode[commonTAD & commonSubset] <- "same sub TAD"
  subTADmode <- factor(subTADmode, levels=c("no TAD", "diff TAD", "diff sub TAD", "same sub TAD"))

  # add column to gene pairs and return
  gp[, colName] <- subTADmode
  return(gp)
}


#' Inter-chromosomal gene pairs counts matrix
#'
#' @import GenomicRanges
#' @export
interChromPairMatrix <- function(gp, genesGR, symmetric=FALSE){

  # get vector of all unique chromosome names
  chroms = seqnames(seqinfo(genesGR))

  # initialize matrix with zero counts
  n = length(chroms)
  mat = matrix(rep(0, n*n), n, dimnames=list(chroms, chroms))

  # get chromsome names of gene pairs
  c1 = seqnames(genesGR[gp[,1]])
  c2 = seqnames(genesGR[gp[,2]])

  # count pairwise occurrences
  counts = count(data.frame(c1, c2, stringsAsFactors=FALSE))

  # iterate over all found pairs and increase counter
  for (i in 1:nrow(counts)){
    mat[counts[i,1], counts[i,2]] = mat[counts[i,1], counts[i,2]] + counts[i,3]

    # if option symmetric is FALSE, count pair as cA-cB and cB-cA
    if (!symmetric){
      mat[counts[i,2], counts[i,1]] =  mat[counts[i,2], counts[i,1]] + counts[i,3]
    }
  }

  # make single letter chromosome names
  rownames(mat) = gsub("chr", "", rownames(mat))
  colnames(mat) = gsub("chr", "", colnames(mat))

  return(mat)
}


#' Adds the expression values for both genes
#'
#' @param expDF data.frame that can be accesd by rows using the indecies or
#'   names in the first two columns of \code{gp}.
#' @export
addPairExp <- function(gp, expDF, expCol, colname="exp"){

  g1exp <- expDF[gp[,1], expCol]
  g2exp <- expDF[gp[,2], expCol]

  gp[,paste0("g1_", colname)] <- g1exp
  gp[,paste0("g2_", colname)] <- g2exp

  return(gp)
}



#' Returns the Pearson correlation coefficient for expression of two input
#' genes.
#'
#' If more than one gene pair is provied in \code{gp} a pairwise matirx of
#' correlation coeficients is returned.
#'
#'
getCor <- function(gp, expDF){

  # correct intput to matrix (in case of only two element vector)
  gp = matrix(gp, ncol=2)

  # get correlation values of all cells/conditions
  # this will make a vector of NA's if the gene is not contained in the expression data set
  # furthermore, cbind(c(.)) guarantees that cor() will deal with column-vectors
  x = t(as.vector(expDF[gp[,1],]))
  y = t(as.vector(expDF[gp[,2],]))

  # return pearson correlation coefficient
  cor(x, y, method="pearson")

}


#' Adds Pearson correlation coefficient for all gene pairs
#'
#' @seealso \code{\link{applyToCisPair}}s, \code{\link{applyToCisPairs}}
#' @export
addCor <- function(gp, expDF, colName="expCor"){

  pairsAsChars = sapply(gp[,1:2], as.character)

  gp[,colName] = apply(pairsAsChars, 1, getCor, expDF=expDF)

  return(gp)
}

#' Checks that none of the genes in pairs is NA.
#'
#' @return a logical vector that has TRUE for each gene pair that has no NA and
#'   FALSE otherwise.
#' @export
pairNotNA <- function(gp){
  return( !is.na(gp[,1]) & !is.na(gp[,2]) )
}

