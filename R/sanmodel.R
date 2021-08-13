#' Create a SANModel from a rates table and an initial number of cells
#'
#' @export
san_model <- function(rates, s0) {
  if (!is.numeric(s0) || (length(s0) != 1) || !is.finite(s0) || (s0 %% 1 != 0) || (s0 < 0))
    stop("s0 must be a single finite and non-negative integer")
  
  if (!is.data.frame(rates) || (nrow(rates) == 0))
    stop("rates must be a non-empty data.table or.data.frame)")
  rates <- as.data.table(rates)
  if (!all(SAN.RATENAMES %in% colnames(rates)))
    stop("rates ", paste0(SAN.RATENAMES[!(SAN.RATENAMES %in% colname(rates))], collapse=", "), " are missing")
  for (c in SAN.RATENAMES) {
    if (!is.numeric(rates[[c]]) || any(is.infinite(rates[[c]])) || !all(is.na(rates[[c]]) || (rates[[c]] >= 0)))
      stop(paste(n, " must contain non-negative and finite numeric values or NA"))
  }
  if (!is.numeric(rates$Tmax) || !all(is.finite(rates$Tmax)) || is.unsorted(rates$Tmax))
    stop(paste("Tmax must contain non-negative and finite numeric values, and must increase monotonically"))

  return(structure(list(rates=rates, s0=s0),
                   class="SANModel")) 
}

#' Fill empty (NA) cells in rates table
#' 
#' @export
ratestable.fill.na <- function(...) UseMethod("ratestable.fill.na")

#' Fill empty (NA) cells in rates table
#' 
#' Empty (NA) cells are set to the next non-empty value in the same column,
#' or to zero if no such value exists.
#'
#' @param ratestable a ratestable (a data.frame or data.table)
#' @return a modified table without empty (NA) cells
#' @export
ratestable.fill.na.data.frame <- function(ratestable) {
  for(c in SAN.RATENAMES) {
    # Fetch row corresponding to rate x
    x <- ratestable[[c]]
    x.set <- !is.na(x)
    # Set empty entries to the next non-empty entry or zero if no such entry exists
    ratestable[[c]] <- c(x[x.set], 0)[c(1, head(cumsum(x.set) + 1, n=-1))]
  }
  return(ratestable)
}

#' Fill empty (NA) cells in the model's rates table
#' 
#' Empty (NA) cells are set to the next non-empty value in the same column,
#' or to zero if no such value exists.
#'
#' @param model a SAN model created with `san_model()`
#' @return a model whose rates table does not contain empty (NA) cells
#' @export
ratestable.fill.na.SANModel <- function(model) {
  model$rates <- ratestable.fill.na(model$rates)
  return(model)
}
