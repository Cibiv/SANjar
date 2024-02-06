#' Create a SANModel
#' 
#' @export
san_model <- function(x, ...) UseMethod("san_model")

#' Simulate a SANModel
#' 
#' @export
simulate <- function(x, ...) UseMethod("simulate")



#' Create a SANModel from a rates table and an initial number of cells
#'
#' @export
san_model.default <- function(rates, s0) {
  if (!is.numeric(s0) || (length(s0) != 1) || !is.finite(s0) || (s0 %% 1 != 0) || (s0 < 0))
    stop("s0 must be a single finite and non-negative integer")
  
  if (!is.data.frame(rates) || (nrow(rates) == 0))
    stop("rates must be a non-empty data.table or.data.frame)")
  rates <- as.data.table(rates)
  if (!all(SAN.RATENAMES %in% colnames(rates)))
    stop("rates ", paste0(SAN.RATENAMES[!(SAN.RATENAMES %in% colname(rates))], collapse=", "), " are missing")
  for (c in SAN.RATENAMES) {
    if (!is.numeric(rates[[c]]) || any(is.infinite(rates[[c]])) || any(rates[[c]] < 0))
      stop(paste(c, " must contain non-negative and finite numeric values or NA"))
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

#' Return the final time, i.e. the time up to which the model defines rates
#'
#' @export
Tmax <- function(model) UseMethod("Tmax")

#' Return the final time, i.e. the time up to which the model defines rates
#'
#' @export
Tmax.SANModel <- function(model) {
  return(model$rates[nrow(model$rates), Tmax])
}

#' Update the model's final time, i.e. the time up to which it defines rates
#' 
#' If the final time exceeds the model's current final time, the validity interval
#' of the last rate set will be extended until the new final time. Otherwise, the
#' validity period of the rates set active at the new final time will be set to the
#' new final time and any later rate sets will be removed.
#' 
#' @export
`Tmax<-` <- function(model, value) set_Tmax(model, value)
set_Tmax <- function(model, value) UseMethod("set_Tmax")

set_Tmax.SANModel <- function(model, value) {
  # Fill in empty cell
  model$rates <- ratestable.fill.na(model$rates)
  # Keep all rows whose validity period ends before the new Tmax, plus
  # one more row because that row's validity period then covers Tmax
  # (if no such row exists, )
  keep <- model$rates$Tmax < value
  keep[which.min(keep)] <- TRUE
  model$rates <- model$rates[keep]
  # Limit the validity period of the last rate set to the new Tmax
  model$rates[nrow(model$rates), Tmax := value]
  return(model)
}

#' Simulate lineage tracing data 
#'
#' @export
simulate.SANModel <- function(model, template, steps=c("aliases", "seqsim", "threshold"), method.seqsim=NULL, p_cutoff=1e-3) {
  lineagesizes <- template$sequencing[, {
    # Simulate stochastic model
    r <- san_stochastic(L=model$s0, s0=1, rates=ratestable.fill.na(model$rates), Tmax=max(day), p_cutoff=p_cutoff)
    # Rename "t" column to "day"
    colnames(r) <- ifelse(colnames(r) == "t", "day", colnames(r))
    # Compute total lineage sizes C
    r[, C := S + A + N ]
    # Add "sid" column
    r[, sid := 0 ]
    # Simulate sequencing
    r <- san_sequencing(r, .SD, steps=steps, method.seqsim=method.seqsim)
    # Generate lineage size table
    if ("seqsim" %in% steps)
      r[, list(day, reads=R, cells=C.obs, S, A, N)]
    else
      r[, list(day, cells=C, S, A, N)]
  }, by=c(template$groups, "sid")]
  
  organoidsizes <- template$sequencing[, {
    san_deterministic(s0=model$s0, rates=ratestable.fill.na(model$rates), samples_per_day=1)[, list(
      day=t, cells=S+A+N, S, A, N
    )]
  }, by=template$groups]
  
  return(LTData(lineagesizes, organoidsizes, template$sequencing, unit=template$unit, groups=template$groups))
}
