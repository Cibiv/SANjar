#' Fill empty (NA) cells in rates table
#' 
#' Empty (NA) cells are set to the next non-empty value in the same column,
#' or to zero if no such value exists.
#'
#' @param ratestable a ratestable (a data.frame or data.table)
#' @return a modified table without empty (NA) cells
#' @export
ratestable.fill.na <- function(ratestable) {
  for(c in SAN.RATENAMES) {
    # Fetch row corresponding to rate x
    x <- ratestable[[c]]
    x.set <- !is.na(x)
    # Set empty entries to the next non-empty entry or zero if no such entry exists
    ratestable[[c]] <- c(x[x.set], 0)[c(1, head(cumsum(x.set) + 1, n=-1))]
  }
  return(ratestable)
}
