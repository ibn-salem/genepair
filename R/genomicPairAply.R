

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
# rangesGR <- GRanges(
#  rep(c("chr1", "chr2"), c(3,2)),
#  IRanges(
#       c(100, 200, 300, 100, 200),
#       c(150, 250, 350, 150, 250)
#    ))
#
# datamat <- rbind(
#     c(10, 20, 30),
#     c(15, 26, 40),
#     c(100, 2, 0),
#     c(10, 20, 30),
#     c(15, 26, 40)
#   )
#
# applyToCisPairs(gp, rangesGR, datamat)
#
#' @keywords internal
#'
applyToCisPairs <- function(gp, rangesGR, datamat, fun=cor){

  # Algorithm
  # (1) group pairs by chromosome
  # (2) get mapping from single chr index to rangesIdx
  #     M <- empty matrix with NAs of dim nxn
  # (3) iterate over all chroms
  #     - select ranges on this chrom
  #     - comp corelation matrix of all pairs
  #     - Fill

  idx1 <- gp[,1]
  idx2 <- gp[,2]

  # to check that paris are from same chromosome
  chr <- as.vector(seqnames(rangesGR))
  if(!all(chr[idx1] == chr[idx2])){strop("Pairs of ranges are not all from the same chromosome.")}
  # check that ranges GR is annotated with a seqinof object

  chrPair <- chr[idx1]
  # just for testing:
  #chrPair <- c("chr1", "chr2", "chr1", "chr1", "chr2")

  chromosomes <- unique(chrPair)

  countChr <- plyr::count(chrPair)

  chrID <- rep(NA, length(chrPair))

  for(chr in countChr$x){
    chrSub <- chr == chrPair
    chrID[chrSub] <- seq(sum(chrSub))
  }

  chrIdxList <- lapply(chromosomes, function(chr) which(chr == chrPair))

  vecList <- lapply(chromosomes, function(chr){

    # get subset of matrix
    m <- datamat[chrIdxList[[chr]],]

  })

}
