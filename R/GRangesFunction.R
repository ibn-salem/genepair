



#' Add column to indicate that query lies within at least one subject object.
#'
#' @param query A GRanges object
#' @param subject A GRanges object
#' @param colName A single string that is used as columnname for the new
#'   annotation column in \code{query}.
#' @return A GRanges object like \code{query} with an additioanl column of
#'   logical values indicating whehter each range is within one range in
#'   \code{subject} or not.
#' @export
addWithinSubject <- function(query, subject, colName="inRegion"){

  mcols(query)[, colName] = countOverlaps(query, subject, type="within") >= 1
  return(query)

}

