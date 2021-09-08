#' The estimated SAN model parameters for the LT47 dataset
"lt47.model"

#' The LT47 cerebral organoid lineage tracing dataset
"lt47"

#' Create a lineage tracing data object
#' 
#' @export
LTData <- function(lineagesizes=NULL, organoidsizes=NULL, sequencing=NULL, unit="reads", groups=character()) {
  # Sanitize arguments
  unit <- match.arg(unit, c("cells", "reads"))
  
  # Validate lineagesizes table
  if (!is.null(lineagesizes)) {
    # Check table layout
    lineagesizes <- as.data.table(lineagesizes)
    if (!all(c(groups, "day", "sid", unit) %in% colnames(lineagesizes)))
      stop("lineagesizes table must contain the columns ", paste0(c(groups, "day", "sid", unit), collapse=", "))
    setkeyv(lineagesizes, c(groups, "day", "sid"))
    # Extract list of samples
    samples <- lineagesizes[, list(dummy=0), keyby=c(groups, "day", "sid")]
  }
  
  # Validate organoidsizes table
  if (!is.null(organoidsizes)) {
    # Check table layout
    organoidsizes <- as.data.table(organoidsizes)
    if (!all(c(groups, "day", "cells") %in% colnames(organoidsizes)))
      stop("organoidsizes table must contain the columns ", paste0(c(groups, "day", "cells"), collapse=", "))
    setkeyv(organoidsizes, c(groups, "day"))
    
    # Check contents
    if (length(groups) > 0) {
      organoidsizes[samples, {
        if (.N < 1)
          stop("organoidsizes table is missing an entry for ",
               paste0(names(.BY), "=", .BY, collapse=", "))
        NULL
      }, by=.EACHI, nomatch=NA, on=groups]
    }
    
    # Add missing columns to organoidsizes table
    if (!("fraction" %in% names(organoidsizes))) {
      organoidsizes[, fraction := 1.0]
    }
  }
  
  # Validate sequencing table
  if (!is.null(sequencing)) {
    # Check table layout
    sequencing <- as.data.table(sequencing)
    if (!all(c(groups, "day", "sid") %in%
             colnames(sequencing)))
      stop("sequencing parameter table must contain the columns ",
           paste0(c(groups, "day", "sid"),
                  collapse=", "))
    setkeyv(sequencing, c(groups, "day", "sid"))
    
    # Check contents
    sequencing[samples, {
      if (.N < 1)
        stop("sequencing table is missing an entry for ",
             paste0(names(.BY), "=", .BY, collapse=", "))
      if (.N > 1)
        stop("sequencing table contains multiple entries for ",
             paste0(names(.BY), "=", .BY, collapse=", "))
      NULL
    }, by=.EACHI, nomatch=NA]
  } else {
    # Create sequencing parameter table
    sequencing <- lineagesizes[, list(0), keyby=c(groups, "day", "sid")][, -length(groups)-3, with=FALSE]
  }

  # Add missing columns to sequencing table
  for(col in c("library_size", "pcr_efficiency", "phantom_threshold", "reads_per_cell", "lineage_aliases")) {
    if (!(col %in% names(sequencing)))
      sequencing[, eval(col) := NA_real_]
  }
  
  # Return LTData object
  return(structure(list(organoidsizes=organoidsizes,
                         lineagesizes=lineagesizes,
                         sequencing=sequencing,
                         unit=as.name(unit),
                         groups=groups),
                    class="LTData"))
}

#' Estimate sequencing parameters (library size, PCR efficiency, phantom threshold, reads per cell)
#' 
#' @export
estimate_sequencing_parameters <- function(lt, ...) UseMethod("estimate_sequencing_parameters")

#' Re-estimate the number of reads per cell without relying on the total organoid sizes
#' 
#' @export
estimate_reads_per_cell <- function(lt, ...) UseMethod("estimate_reads_per_cell")

#' Infer the total number of cells each sub-library of an organoid corresponds to
#'
#' @export
partial_organoid_sizes <- function(partial, ...) UseMethod("partial_organoid_sizes")

#' Normalize library sizes within each sample group to the same sequencing library size
#' 
#' @export
normalize_library_sizes <- function(lt, ...) UseMethod("normalize_library_sizes")

#' Transform lineage sizes from (relative) read counts into (absolute) cells counts
#' 
#' @export
absolute_lineage_sizes <- function(lt, ...) UseMethod("absolute_lineage_sizes")

#' Subset of lineage tracing data
#' 
#' @export
subset.LTData  <- function(lt, ...) {
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

#' Subset of lineage tracing data
#' 
#' @export
`[.LTData` <- subset.LTData

#' Estimate sequencing parameters (library size, PCR efficiency, phantom threshold, read_per_cell)
#' 
#' @import gwpcR
#' @import matrixStats
#' @export
estimate_sequencing_parameters.LTData <- function(lt, pcr_efficiency=NA, lineage_aliases=NA, molecules=1, replace=FALSE) {
  if (lt$unit == "cells")
    stop("sequencing parameter estimation requires read counts, not cell counts")
  
  # Copy sequencing parameter table before modifying it below
  lt$sequencing <- copy(lt$sequencing)
  
  # Estimate PCR efficiency and number of lineage aliases from day-0 data if possible
  # fall back to hard-coded default of 33% and no aliasing otherwise
  if (is.na(pcr_efficiency) || is.na(lineage_aliases)) {
    seq.day0 <- lt$lineagesizes[, {
      # Fetch organoidsizes for current group
      os <- if (length(lt$groups) > 0)
        lt$organoidsizes[.BY]
      else
        lt$organoidsizes
      # Compute median number of cells on day 0
      cells0 <- os[day==0, median(cells, na.rm=TRUE)]
      # Fetch lineagesizes for day 0
      ls0 <- .SD[day==0]
      r <- if (nrow(ls0) > 0) {
        # We have day 0 data for the current group, estimate PCR efficiency and
        # total number of lineages (including unobserved).
        # When computing the median number of lineages present (including unobserved!),
        # we take the accuracy of the estimates into account by weighting them with
        # (1-loss)^4. The exponent 4 stems from applying the delta method to the
        # factor 1 / (1-loss) which is used to compute n.tot from n.obs. 
        pcr0 <- ls0[, gwpcrpois.est(reads, threshold=min(reads), molecules=molecules), by=.(sid)]
        pcr0[, list(pcr_efficiency=median(efficiency, na.rm=TRUE),
                    lineage_aliases=weightedMedian(n.tot/cells0, (1 - loss)^4, na.rm=TRUE))]
      } else {
        # No day-0 data for the current group
        list(pcr_efficiency=NA_real_, lineage_aliases=NA_real_)
      }
      # Replace NAs by hardcoded defaults
      if (is.na(r$pcr_efficiency)) {
        warning("Unable to estimate PCR efficiency for ",
                paste0(names(.BY), "=", lapply(.BY, as.character), collapse=", "),
                ", setting to hard-coded default of 33%")
        r$pcr_efficiency <- 0.33
      }
      if (is.na(r$lineage_aliases)) {
        message("Unable to estimate number of lineage aliases for ",
                paste0(names(.BY), "=", lapply(.BY, as.character), collapse=", "),
                ", setting to 1")
        r$lineage_aliases <- 1
      }
      r
    }, by=eval(lt$groups)]
    if (length(lt$groups) > 0)
      seq.day0 <- seq.day0[lt$sequencing[, c(lt$groups, "day", "sid"), with=FALSE], on=c(lt$groups)]
    else
      seq.day0 <- seq.day0[, c(lt$sequencing[, c("day", "sid"), with=FALSE], .SD) ]
    setkeyv(seq.day0, c(lt$groups, "day", "sid"))
  }

  # Update PCR efficiency in sequencing table
  if (is.na(pcr_efficiency)) {
    lt$sequencing[seq.day0, pcr_efficiency :=
                    ifelse(is.na(pcr_efficiency) | replace, i.pcr_efficiency, pcr_efficiency) ]
  } else {
    lt$sequencing[is.na(is.na(pcr_efficiency) | replace), pcr_efficiency := ..pcr_efficiency]
  }

  # Update number of lineage aliases in sequencing table
  if (is.na(lineage_aliases)) {
    lt$sequencing[seq.day0, lineage_aliases :=
                    ifelse(is.na(lineage_aliases) | replace, i.lineage_aliases, lineage_aliases) ]
  } else {
    lt$sequencing[is.na(is.na(lineage_aliases) | replace), lineage_aliases := ..lineage_aliases]
  }
  
  # Determine each sample's library size and phantom threshold and update sequencing table
  seq <- lt$lineagesizes[, list(
    library_size=sum(reads),
    phantom_threshold=min(reads)
  ), keyby=c(lt$groups, "day", "sid")]
  lt$sequencing[seq, library_size :=
                  ifelse(is.na(library_size) | replace, i.library_size, library_size)      ]
  lt$sequencing[seq, phantom_threshold :=
                  ifelse(is.na(phantom_threshold) | replace, i.phantom_threshold, phantom_threshold) ]

  # Estimate the number of reads per cell from the total organoid size and library size
  if (!is.null(lt$organoidsizes)) {
    lt$sequencing[, reads_per_cell := {
      # Fetch organoidsizes for current sample group
      organoidsizes <- if(length(lt$groups) > 0)
        lt$organoidsizes[.BY]
      else
        lt$organoidsizes
      # Interpolate organoid sizes for current sample group
      t <- organoidsizes[, list(logc=mean(log(cells))), by=.(day)]
      logf <- approxfun(t$day, t$logc, rule=2)
      f <- function(...) exp(logf(...))
      # Return new reads_per_cell column
      ifelse(is.na(reads_per_cell) | replace, library_size / f(day), reads_per_cell)
    }, by=eval(lt$groups)]
  }
  
  return(lt)
}

#' Re-estimate the number of reads per cell without relying on the total organoid sizes
#' 
#' @export
estimate_reads_per_cell.LTData <- function(lt, method="singleton_mode") {
  switch(match.arg(method, choices=c("singleton_mode")),
         `singleton_mode`=estimate_reads_per_cell.LTData.singleton_mode(lt))
}

estimate_reads_per_cell.LTData.singleton_mode <- function(lt, replace=TRUE) {
  if (lt$unit == "cells")
    stop("reads per cell estimation requires read counts, not cell counts")
  
  # Copy sequencing parameter table before modifying it below
  lt$sequencing <- copy(lt$sequencing)
  
  # Determine mode of log-lineagesizes, and assume these lineages represent
  # single-cell lineages.
  seq <- lt$lineagesizes[, {
    f <- density(log10(reads))
    list(reads_per_cell=10^(f$x[which.max(f$y)]))
  }, keyby=c(lt$groups, "day", "sid")]
  lt$sequencing[seq, reads_per_cell :=
                  ifelse(is.na(reads_per_cell) | replace, i.reads_per_cell, reads_per_cell) ]
  
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
      list(logf=log10(library_size) - log10(sum(library_size)),
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
      logf.mean.fun <- approxfun(day, logf.mean, rule=2)
      logf.var.fun <- if (sum(!is.na(logf.var)) >= 2) 
        approxfun(day, logf.var, rule=2)
      else
        function(day) { 0 }
      # Scale organoid sizes for each day separately
      wholesizes[, {
        # Scale the organoid sizes by shifting the log-mean according to the sub-library
        # fractions, accounting for the uncertainty in these fractions.
        logc <- log10(cells)
        logc.mean <- mean(logc)
        logc.var <- var(logc)
        if (logc.var > 0) {
          # To account for uncertainty in the fractions f, the log-residuals (i.e
          # the deviations of the individual log-measurements from the per-day log-mean) are
          # scaled. Specifically, scaled with s such that the log-variance of the final cell
          # counts is increased by the log-variance of the fractional sub-library sizes.
          s <- sqrt(1 + logf.var.fun(day) / var(logc))
        } else if (logf.var.fun(day) > 0) {
          # If the cell counts are all identical, scaling the deviations from the mean
          # does not work. Instead, we spread them according to the quantiles of a
          # normal distribution around the mean, selected to (approximately) yield
          # the standard deviation of the fractions.
          # TODO: Figure out how to compute quantiles which yield the correct std.
          # dev *exactly*. The current formula is a result of trial-and-error, and
          # the "fudge factor" d was found by experimentation. Ugh.
          k <- length(logc)
          d <- 0.171010
          logc <- qnorm(seq(from=1/(2*k+d), to=(2*k+d-1)/(2*k+d), length.out=k),
                        mean=logc.mean, sd=sqrt(logf.var.fun(day)))
          s <- 1
        } else {
          warning("Unable to compute mean-deviation scaling factor for some samples")
          s <- 1
        }
        # Scale cell counts, and also output the fraction we found
        # Note that "cells" and "fraction" are required columns for LTData,
        # but fraction.logsd is not. We just output the latter for validation
        # purposes.
        sd <- copy(.SD)
        sd[, cells := 10^(s*(logc - logc.mean) + logc.mean + logf.mean.fun(day))]
        sd[, fraction := 10^(logf.mean.fun(day)) ]
        sd[, fraction.logsd := sqrt(logf.var.fun(day)) ]
        sd
      }, by=.(day)]
    }, by=sublib]
  }, by=groups]
  setkeyv(partial_organoidsizes, c(partial$groups, "day"))
  
  partial$organoidsizes <- partial_organoidsizes
  return(partial)
}

#' Normalize lineage sizes within each sample group to the same sequencing library size
#' 
#' @export
normalize_library_sizes.LTData <- function(lt, method="scale") {
  if (any(is.na(lt$sequencing[, library_size])))
    stop("LTData must contain valid library sizes in sequencing parameter table to normalize library sizes, ",
         "call estimate_sequencing_parameters first")
  
  method <- match.arg(method, c("scale", "subsample"))
  
  # Determine library-size-normalized sequencing parameters 
  sequencing.norm <- lt$sequencing[, {
    # The target library size for the current group is the smallest library size among replicates
    library_size.min <- min(library_size)
    # The assumed PCR efficiency is the median PCR efficiency across replicates
    pcr_efficiency.med <- median(pcr_efficiency)
    # The phantom threshold is the *maximum* across replicates, after taking the
    # library size scaling that happens below into account. This ensures that we
    # only move the threshold *upwards* when we re-apply the new threshold below.
    phantom_threshold.max <- max(ceiling(phantom_threshold * (library_size.min / library_size)))
    # Return new per-sample sequencing parameters
    list(sid=sid,
         library_size=library_size.min,
         pcr_efficiency=pcr_efficiency.med,
         phantom_threshold=phantom_threshold.max,
         reads_per_cell=reads_per_cell * (library_size.min / library_size),
         lineage_aliases=lineage_aliases)
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
    reads.norm <- if (method == "scale")
      round(as.numeric(reads) * (seq.norm$library_size / seq$library_size))
    else if (method == "subsample")
      rbinom(.N, size=round(as.numeric(reads)), prob=seq.norm$library_size / seq$library_size)
    else
      stop("unknown normalization method ", method)
    # Re-apply phantom threshold
    ifelse(reads.norm >= seq.norm$phantom_threshold, reads.norm, 0)
  }, keyby=c(lt$groups, "day", "sid")]
  
  # Return updated LTData object
  lt$lineagesizes <- lineagesizes.norm
  lt$sequencing <- sequencing.norm
  return(lt)
}

#' Transform lineage sizes from (relative) read counts into (absolute) cells counts
#' 
#' @export
absolute_lineage_sizes.LTData <- function(lt) {
  if (lt$unit == "cells")
    stop("lineage lineages are already expressed in cells")
  
  # Copy lineagesizes before modifying it below
  lt$lineagesizes <- copy(lt$lineagesizes)
  
  # Convert read count to cell count using the estimated number of reads per cell
  lt$lineagesizes[lt$sequencing, cells := reads / reads_per_cell]
  
  # Lineage size units are now cells, not reads
  lt$unit <- "cells"
  
  return(lt)
}
