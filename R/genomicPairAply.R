

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
  for (chr in chromosomes){

    subIdx <- which(chr == regChr)

    subDat <- datamat[subIdx,]

    subCor <- fun(t(subDat))

    m[subIdx, subIdx] <- subCor


  }
  #-----------------------------------------------------------------------------
  # (5) query matrix with all pairs.
  #-----------------------------------------------------------------------------
  values <- m[cbind(idx1, idx2)]

  return(values)

}
