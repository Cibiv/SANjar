#' The estimated SAN model parameters for the LT47 dataset
"lt47.ps"

#' The LT47 cerebral organoid lineage tracing dataset
"lt47"

#' Create a lineage tracing data object
#' 
#' @export
LTData <- function(lineagesizes=NULL, organoidsizes=NULL, sequencing=NULL, unit="reads", groups=character()) {
  # Sanitize arguments
  if (!is.null(organoidsizes)) {
    organoidsizes <- as.data.table(organoidsizes)
    if (!all(c(groups, "day", "sid", "cells") %in% colnames(organoidsizes)))
      stop("organoidsizes table must contain the columns ", paste0(c(groups, "day", "sid", "cells"), collapse=", "))
    setkeyv(organoidsizes, c(groups, "day", "sid"))
  }
  
  unit <- match.arg(unit, c("cells", "reads"))
  
  if (!is.null(lineagesizes)) {
    lineagesizes <- as.data.table(lineagesizes)
    if (!all(c(groups, "day", "sid", unit) %in% colnames(lineagesizes)))
      stop("lineagesizes table must contain the columns ", paste0(c(groups, "day", "sid", unit), collapse=", "))
    setkeyv(lineagesizes, c(groups, "day", "sid"))
  }

  if (!is.null(sequencing)) {
    sequencing <- as.data.table(sequencing)
    if (!all(c(groups, "day", "sid", "library_size", "pcr_efficiency", "phantom_threshold") %in%
             colnames(sequencing)))
      stop("sequencing parameter table must contain the columns ",
           paste0(c(groups, "day", "sid", "library_size", "pcr_efficiency", "phantom_threshold"),
                  collapse=", "))
    setkeyv(sequencing, c(groups, "day", "sid"))
  }
  
  # Return LTData object
  return(structure(list(organoidsizes=organoidsizes,
                         lineagesizes=lineagesizes,
                         sequencing=sequencing,
                         unit=as.name(unit),
                         groups=groups),
                    class="LTData"))
}

#' Estimate sequencing parameters (library size, PCR efficiency, phantom threshold)
#' 
#' @export
estimate_sequencing_parameters <- function(lt, ...) UseMethod("estimate_sequencing_parameters")

#' Infer the total number of cells each sub-library of an organoid corresponds to
#'
#' @export
partial_organoid_sizes <- function(partial, ...) UseMethod("partial_organoid_sizes")

#' Normalize lineage sizes within each sample group to the same sequencing library size
#' 
#' @export
normalize_library_sizes <- function(lt, ...) UseMethod("normalize_library_sizes")

#' Subset of lineage tracing data
#' 
#' @export
`[.LTData`  <- function(lt, ...) {
  # Convert arguments to index list and then to filter expression
  idx <- list(...)
  f <- Reduce(function(e, p) substitute(e & p, list2env(list(e=e, p=p))),
              Map(function(k) bquote((.(as.name(k))==.(idx[[k]]))), names(idx)))
  
  return(structure(list(organoidsizes=lt$organoidsizes[eval(f)],
                        lineagesizes=lt$lineagesizes[eval(f)],
                        sequencing=lt$sequencing[eval(f)],
                        unit=as.name(lt$unit),
                        groups=lt$groups),
                   class="LTData"))
}

#' Estimate sequencing parameters (library size, PCR efficiency, phantom threshold)
#' 
#' @import gwpcR
#' @export
estimate_sequencing_parameters.LTData <- function(lt, pcr_efficiency=NA, molecules=1) {
  if (lt$unit == "cells")
    stop("sequencing parameter estimation requires read counts, not cell counts")
  
  # Estimate PCR efficiency from day-0 data if possible, fall back to hard-coded default of 33% otherwise
  if (is.na(pcr_efficiency)) {
    day0 <- lt$lineagesizes[day==0]
    pcr_efficiency <- if (nrow(day0) > 0) {
      day0.gwpcr <- pcr_efficiency <- day0[, {
        gwpcrpois.est(reads, threshold=min(reads), molecules=molecules)
      }, by=c(lt$groups, "day", "sid")]
      day0.gwpcr[, signif(median(efficiency), 3) ]
    } else {
      warning("PCR efficiency cannot be estimated due to missing day 0 data, setting to 33%")
      0.33
    }
  }
  
  # Determine each sample's library size and phantom threshopld
  lt$sequencing <- lt$lineagesizes[, list(
    library_size=sum(reads),
    pcr_efficiency=pcr_efficiency,
    phantom_threshold=min(reads)
  ), keyby=c(lt$groups, "day", "sid")]
  
  return(lt)
}

#' Infer the total number of cells each sub-library of an organoid corresponds to
#'
#' @export
partial_organoid_sizes.LTData <- function(partial, whole, sublib) {
  if (!all(whole$groups %in% partial$groups))
    stop("whole-organoid group identifiers must be a subset of partial-organoid group identifiers")
  if (missing(sublib)) {
    sublib <- setdiff(partial$groups, whole$groups)
    message("no sub-library identifiers specified, will use ", paste0(sublib, collapse=","))
  }
  if (!all(sublib %in% partial$groups))
    stop("sub-library identifiers must be a subset of partial-organoid group identifiers")
  if (any(sublib %in% whole$groups))
    stop("sub-library identifiers cannot also be whole-organoid group identifiers")
  groups=setdiff(partial$groups, sublib)

  # Handle each sample group separately
  partial_organoidsizes <- partial$sequencing[, {
    # Compute average (in log-space) of fractional sizes of the sub-libraries on each day.
    # First, group by day and sample to compute the fractional sub-library sizes separately
    # for each sample, then group by sub-library and day to average the per-sample fractions.
    f <- .SD[, c(
      .SD[, sublib, with=FALSE],
      list(logf=log(library_size) - log(sum(library_size)),
           f=library_size / sum(library_size))
    ), by=c("day", "sid")][, list(
      logf.mean=mean(logf),
      logf.var=var(logf)
    ), keyby=c(sublib, "day")]

    # Fetch whole-organoid sizes for the correct group
    wholesizes <- if (length(whole$groups) > 0) {
      stopifnot(all(whole$groups %in% names(.BY)))
      whole$organoidsizes[.BY[whole$groups]]
    } else
      whole$organoidsizes
    
    # Compute per-sub-library organoid sizes by scaling the whole-organoid sizes accordingly
    f[, {
      # Create interpolation functions for the mean and variance of the sub-library fractions
      logf.mean.fun <- approxfun(day, logf.mean)
      logf.var.fun <- if (sum(!is.na(logf.var)) >= 2) 
        approxfun(day, logf.var)
      else
        function(day) { 0 }
      # Scale organoid sizes for each day separately
      wholesizes[, {
        # Scale the organoid sizes by shifting the log-mean according to the sub-library
        # fractions. To account for uncertainty in these fractions, the log-residuals (i.e
        # the deviations of the individual log-measurements from the per-day log-mean) are
        # scaled. Specifically, scaled with s such that the log-variance of the final cell
        # counts is increased by the log-variance of the fractional sub-library sizes.
        logc <- log(cells)
        logc.mean <- mean(logc)
        s <- sqrt(1 + logf.var.fun(day) / var(logc))
        if (!is.finite(s)) s <- 1
        copy(.SD)[, cells := exp(s*(logc - logc.mean) + logc.mean + logf.mean.fun(day))]
      }, by=.(day)]
    }, by=sublib]
  }, by=groups]
  setkeyv(partial_organoidsizes, c(partial$groups, "day", "sid"))
  
  partial$organoidsizes <- partial_organoidsizes
  return(partial)
}

#' Normalize lineage sizes within each sample group to the same sequencing library size
#' 
#' @export
normalize_library_sizes.LTData <- function(lt, method="scale") {
  switch(match.arg(method, choices=c("scale")),
         `scale`=normalize_library_sizes.LTData.scale(lt))
}

normalize_library_sizes.LTData.scale <- function(lt) {
  if (is.null(lt$sequencing))
    stop("LTData must contain sequencing parameter table to normalize library sizes, ",
         "call estimate_sequencing_parameters first")
  
  # Determine library-size-normalized sequencing parameters 
  sequencing.norm <- lt$sequencing[, {
    # The target library size for the current group is the smallest library size among replicates
    library_size.min <- min(library_size)
    # The assumed PCR efficiency is the median PCR efficiency across replicates
    pcr_efficiency.med <- median(pcr_efficiency)
    # The phantom threshold is the *maximum* across replicates, after taking the
    # library size scaling that happens below into account. This ensures that we
    # only move the threshold *upwards* when we re-apply the new threshold below.
    phantom_threshold.max <- max(ceiling(phantom_threshold * library_size.min / library_size))
    # Return new per-sample sequencing parameters
    list(sid=sid,
         library_size=library_size.min,
         pcr_efficiency=pcr_efficiency.med,
         phantom_threshold=phantom_threshold.max)
  }, by=c(lt$groups, "day")]
  setkeyv(sequencing.norm, c(lt$groups, "day", "sid"))

  # Transform read counts
  lineagesizes.norm <- copy(lt$lineagesizes)
  lineagesizes.norm[, reads := {
    # Fetch original and normalized sequencing parameters for the current sample
    seq <- lt$sequencing[.BY]
    seq.norm <- sequencing.norm[.BY]
    if ((nrow(seq) < 1) || (nrow(seq.norm) < 1))
      stop("missing sequencing parameters for sample ", paste(names(.BY), "=", as.character(.BY), collapse=","))
    
    # Scale read counts
    reads.norm <- round(reads * seq.norm$library_size / seq$library_size)
    # Re-apply phantom threshold
    ifelse(reads.norm >= seq.norm$phantom_threshold, reads.norm, 0)
  }, keyby=c(lt$groups, "day", "sid")]
  
  # Return updated LTData object
  lt$lineagesizes <- lineagesizes.norm
  lt$sequencing <- sequencing.norm
  return(lt)
}
