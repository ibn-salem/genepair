

#' Apply a function to a pair of genomic ranges.
#'
#' This function returns a vector with the resulting functin call for each pair..
#'
#' @param gp data.frame of pairs of genomic ranges. First two columns hold
#'   indices of the range in \code{rangesGR}.
#' @param rangesGR GenomicRanges object with all ranges used in \code{gp}.
#' @param datamat a matrix of with data values associated to each rang in
#'   \code{rangesGR}. Assumes the ranges in rows while \code{row.names(datamat)}
#'   are the \code{names(rangesGR)}.
#' @param fun A function that takes two numeric vectors as imput to compute a
#'   summary statsitic. Default is \code{\link{cor()}}.
#' @return A Numeric vector
#' @examples
#' rangesGR <- GRanges(
#'  rep(c("chr1", "chr2"), c(3,2)),
#'  IRanges(
#'       c(100, 200, 300, 100, 200),
#'       c(150, 250, 350, 150, 250)
#'    ))
#' gp <- data.frame(
#'  g1=c(1,4,2,1,4),
#'  g2=c(2,4,3,3,5)
#' )
#'
#' datamat <- rbind(
#'     c(10, 20, 30),
#'     c(15, 26, 40),
#'     c(100, 2, 0),
#'     c(10, 20, 30),
#'     c(15, 26, 40)
#'   )
#' applyToCisPairs(gp, rangesGR, datamat)
#'
#' @export
applyToCisPairs <- function(gp, rangesGR, datamat, fun=cor){

  # Algorithm
  # (1) group pairs by chromosome
  # (2) get mapping from single chr index to rangesIdx
  # (3) M <- empty matrix with NAs of dim nxn
  # (4) iterate over all chroms
  #     - select ranges on this chrom
  #     - comp corelation matrix of all pairs
  #     - Fill
  # (5) query matrix with all pairs

  n <- length(rangesGR)
  idx1 <- gp[,1]
  idx2 <- gp[,2]
  regChr <- as.vector(GenomeInfoDb::seqnames(rangesGR))

  # to check that paris are from same chromosome
  if (!all(regChr[idx1] == regChr[idx2])) {
    strop("Pairs of ranges are not all from the same chromosome.")
  }
  # check that ranges GR is annotated with a seqinof object

  chromosomes <- unique(regChr)

  countChr <- plyr::count(regChr)

  #-----------------------------------------------------------------------------
  # (3) initialize empty matrix
  #-----------------------------------------------------------------------------
  m <- Matrix::Matrix(0,nrow=n, ncol=n)

  #-----------------------------------------------------------------------------
  # (4) iteratre over all chromosomes
  #-----------------------------------------------------------------------------
  subCorList <- BiocParallel::bplapply(chromosomes, function(chr){

    # to get indexes of regions on the current chrom
    subIdx <- which(chr == regChr)

    # select coresponding data
    subDat <- datamat[subIdx,]

    # apply function to subset of data
    fun(t(subDat))

  })

  names(subCorList) <- chromosomes

  # update Matrix with submatrixes
  for (chr in chromosomes){

    subIdx <- which(chr == regChr)

    subCor <- subCorList[[chr]]

    m[subIdx, subIdx] <- subCor

  }
  #-----------------------------------------------------------------------------
  # (5) query matrix with all pairs.
  #-----------------------------------------------------------------------------
  values <- m[cbind(idx1, idx2)]

  return(values)

}


#' Apply a function to pairs of close genomic regions.
#'
#' This function returns a vector with the resulting functin call for each input pair.
#'
#'
#' @param gp data.frame of pairs of genomic ranges. First two columns hold
#'   indices of the range in \code{rangesGR}.
#' @param rangesGR GenomicRanges object with all ranges used in \code{gp}.
#' @param datamat a matrix of with data values associated to each rang in
#'   \code{rangesGR}. Assumes the ranges in rows while \code{row.names(datamat)}
#'   are the \code{names(rangesGR)}.
#' @param fun A function that takes two numeric vectors as imput to compute a
#'   summary statsitic. Default is \code{\link{cor()}}.
#' @return A Numeric vector
#' @examples
#' rangesGR <- GRanges(
#'  rep(c("chr1", "chr2"), c(3,2)),
#'  IRanges(
#'       c(100, 200, 300, 100, 200),
#'       c(150, 250, 350, 150, 250)
#'    ))
#' gp <- data.frame(
#'  g1=c(1,4,2,1,4),
#'  g2=c(2,4,3,3,5)
#' )
#'
#' datamat <- rbind(
#'     c(10, 20, 30),
#'     c(15, 26, 40),
#'     c(100, 2, 0),
#'     c(10, 20, 30),
#'     c(15, 26, 40)
#'   )
#' applyToCisPairs(gp, rangesGR, datamat)
#'
#' @export
applyToClosePairs <- function(gp, rangesGR, datamat, fun=cor, maxDist=10^6){

  # Algorithm
  # (1) Group genome in overlapping bins of size 2*maxDist
  # (2) Run pairwise correlation for all ranges in each bin
  # (3) Combine correlations to data.frame with proper id1 and id2 in first columns
  # (4) Query data frame with input pairs
  #   /uses inner_join() from dplyr like in
  #  http://stackoverflow.com/questions/26596305/match-two-data-frames-based-on-multiple-columns

  #-----------------------------------------------------------------------------
  # (1) group ranges in by bins
  #-----------------------------------------------------------------------------

  # create GRanges object for entire genome
  genomeGR <- GenomicRanges::GRanges(GenomeInfoDb::seqinfo(rangesGR))

  # tile genoe in overlapping bins of with 2*maxDist
  binGR <- unlist(GenomicRanges::slidingWindows(genomeGR, 2*maxDist, maxDist))

  hits <- GenomicRanges::findOverlaps(binGR, rangesGR)
  #-----------------------------------------------------------------------------
  # (2) compute pairwise correlatin for all ranges in each bin
  #-----------------------------------------------------------------------------
  corMatList <- lapply(1:length(binGR), function(i){

    # get regions in this bin
    idxs <- S4Vectors::subjectHits(hits)[S4Vectors::queryHits(hits) == i]

    n = length(idxs)

    # compute pairwise correlations for all regions in this bin
    m <- cor(t(datamat[idxs,]))

    corMat <- cbind(
      rep(idxs, n),
      rep(idxs, each=n),
      array(m)
      )
  })

  #-----------------------------------------------------------------------------
  # (3) combine all data.frames
  #-----------------------------------------------------------------------------
  corDF <- data.frame(do.call("rbind", corMatList))

  #-----------------------------------------------------------------------------
  # (4) Query with input pairs
  #-----------------------------------------------------------------------------
  names(corDF) <- c("id1", "id2", "val")
  names(gp)[1:2] <- c("id1", "id2")

  matches <- dplyr::semi_join(corDF, gp, by=c("id1", "id2"))


  return(matches$val)

}
