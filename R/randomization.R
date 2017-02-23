



#' Sample randomly pairs from set of genes with equal probabilities
#' @param genesGR a \code{\link[GenomicRanges]{GRanges}} with genes.
#' @param n number of pairs to sample.
#' @return A \code{\link{data.frame} with indexes of ranges in \code{genesGR}
#' in the frist two columns.
#' @export
getRandomPairs <- function(genesGR, n){

    gp = data.frame(
      t(replicate(
          n,
          sample(1:length(genesGR), size=2, replace=FALSE)
        )
      )
    )

  # sort ids that the first <= second
  gp <- t(apply(gp, 1, sort))

  # sort pairs according to first column
  gp <- gp[order(gp[,1], gp[,2]),]

  gp <- as.data.frame(gp)
  names(gp) = c("g1", "g2")

  return(gp)
}



#' Compute sampling weights.
#'
#' Compute sampling weights to sample from a population with osberved
#' probabilities by binning.
#'
#' @param obs vector of observed values according to which the sampling should
#'   be done (e.g distances observed for paralog gene pairs).
#' @param population vector of all values in the total population from which one
#'   wants to sample (e.g distances of all gene pairs)
#' @param breaks breaks used for sampling resolution (see \code{breaks} argument
#'   in \code{\link[graphics]{hist}} function).
#'
#' @return numeric vector with sampling weights to sample from \code{samp} with
#'   probabilities observed in \code{obs}.
#'
#' @export
weightsByBin <- function(obs, population, breaks=50){

  breaksAll <- hist(c(obs, population), breaks=breaks, plot=FALSE)$breaks

  # calculate the number of observation for nBin equal sized bins
  hObs <- hist(obs, breaks=breaksAll, plot=FALSE)

  # get for each individual in the population the bin index
  binPop <- .bincode(population, breaksAll, include.lowest=TRUE)

  # get counts per bin in the population
  hPop <- hist(population, breaks=breaksAll, plot=FALSE)

  # get the number of observed counts normalized to one as weight
  # normalize the number of observed counts by the bias observed in the population
  weight <- hObs$counts[binPop] / hPop$counts[binPop]

  # remove NA's, e,g, bis not observed in obs but in population. Set their probability to zero
  weight[is.na(weight)] = 0

  # normalize the weights to sum up to 1
  weightNormed <- weight / sum(weight)

  return(weightNormed)
}


