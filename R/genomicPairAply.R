

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
#' gp <- data.frame(g1=c(1,2), g2=c(2,3))
#' rangesGR <- GRanges(chr1, IRanges(c(100, 200, 300), c(150, 250, 350)))
#' datamat <- rbind(c(10, 20, 30), c(15, 26, 40), c(100, 2, 0))
#'
#' applyToCisPairs(gp, rangesGR, datamat)
#'
#' @keywords internal
#'
applyToCisPairs <- function(gp, rangesGR, datamat, fun=cor){

  # to check that paris are from same chromosome
  chr <- seqnames(rangesGR)

}
