#' Create SAN model parametrization
#'
#' @export
san_parametrization <- function(ranges, basemodel, s0="autodetect", cc_days="autodetect",
                                rs_days="autodetect", rs_ranks="autodetect", sc_days="autodetect") {
  # Validate parameter names (I)
  pnames <- names(ranges)
  pvalid <- grepl("^([0-9]+)([A-Z]+)$", pnames)
  if (!all(pvalid))
    stop("invalid parameter names ", paste0(pnames[!pvalid], collapse=", "))
  
  s0 <- match.arg(s0, c("autodetect", "basemodel"))
  
  # Create map that translates parameter names to rate table columns
  # rateslist.pmap is a that contains as many entries as the rate table has rows, i.e.
  # it mimic the format of the rate table itself as stored in SANParametrization
  # object (see split() below).
  basemodel <- ratestable.fill.na(basemodel)
  rateslist.pmap <- lapply(basemodel$rates$Tmax, function(tmax) {
    # Find rates which affect the current time interval (up to tmax)
    pattern <- paste0("^", tmax, "([A-Z]+)$")
    for.tmax <- grepl(pattern, pnames)
    # Extract rate column name(s) from the parameter names (e.g. 40SA translates to columns S, A)
    p <- pnames[for.tmax]
    r <- gsub("X", "0", sub(pattern, "\\1", p))
    # Create map, i.e. a vector with rate names as values and target columns as names
    k <- nchar(r)
    map <- rep(p, k)
    if (length(map) > 0)
      names(map) <- paste0("r_", unlist(strsplit(r, split="")))
    # Check that target columns exist (II)
    if (!(all(names(map) %in% colnames(basemodel$rates))))
      stop("invalid parameter names ", paste0(map[!(names(map) %in% colnames(basemodel$rates))], collapse=", "),
           " (no such rates)")
    return(map)
  })
  # Check that all parameter names could be translated (III)
  if (!all(pnames %in% unlist(rateslist.pmap)))
    stop("invalid parameter names ", paste0(pnames[!(pnames %in% unlist(rateslist.pmap))], collapse=", "),
         " (no such time interval endpoints)")

  # Convert default ratetable into list, and replace the last Tmax by infinity
  # to allow simulations to run arbitrarily long
  rateslist.base <- lapply(split(ratestable.fill.na(basemodel$rates), by="Tmax"), as.list)
  rateslist.base[[length(rateslist.base)]]$Tmax <- Inf
  
  # Create SANParametrization object
  return(structure(list(
    names=as.character(pnames),
    ranges=as.list(ranges),
    rateslist.base=rateslist.base,
    s0.base=if (is.character(s0) && (s0 == "autodetect")) NA else as.integer(basemodel$s0),
    rateslist.pmap=rateslist.pmap,
    cc_days=if (is.character(cc_days) && (cc_days == "autodetect")) NULL else as.integer(cc_days),
    rs_days=if (is.character(rs_days) && (rs_days == "autodetect")) NULL else as.integer(rs_days),
    rs_ranks=if (is.character(rs_ranks) && (rs_ranks == "autodetect")) NULL else as.integer(rs_ranks),
    sc_days=if (is.character(sc_days) && (sc_days == "autodetect")) NULL else as.integer(sc_days)
  ), class="SANParametrization"))
}

#' Convert a parameter vector to a rate list and inital cell count (s0)
#' 
#' @export
rateslist_s0 <- function(parameterization, params) UseMethod("rateslist_s0")

#' Adapt a SAN model parametrization (SANP instance) to a specific dataset
#' 
#' @export
autodetect.metaparameters <- function(parameterization, ...) UseMethod("autodetect.metaparameters")

#' Convert a parameter vector to a rate list and initial cell count (s0)
#' 
#' @export
rateslist_s0.SANParametrization <- function(parameterization, params) {
  rl <- parameterization$rateslist.base
  rateslist.pmap <- parameterization$rateslist.pmap
  for(i in  seq_along(rl)) {
    m <- rateslist.pmap[[i]]
    if (length(m) > 0)
      rl[[i]][names(m)] <- params[m]
  }
  
  s0 <- parameterization$s0.base
  if ("s0" %in% parameterization$names)
    s0 <- params["s0"]
  
  return(list(s0=s0, rateslist=rl))
}

#' Convert a parameter vector to a SAN model instance
#' 
#' @export
san_model.SANParametrization <- function(parameterization, params, lt=NULL, s0=NA, Tmax=NA) {
  # If a dataset was provided, use it to fill in missing meta-parameters
  if (!is.null(lt))
    parameterization <- autodetect.metaparameters(parameterization, lt)
  
  # Extract the parameters from the parameter vector or single-row table
  if (is.data.table(params)) {
    if (nrow(params) != 1)
      stop("cannot create a SAN model instance from multiple parameter combinations")
    params <- unlist(params[, parameterization$names, with=FALSE])
  } else if (is.data.frame(params)) {
    if (nrow(params) != 1)
      stop("cannot create a SAN model instance from multiple parameter combinations")
    params <- unlist(params[, parameterization$names])
  } else if (is.vector(params)) {
    params <- params[parameterization$names]
  } 
    
  # Convert parameter vector to ratelist, and set s0 unless it was overriden manually
  rl_s0 <- rateslist_s0(parameterization, params)
  if (is.na(s0))
    s0 <- rl_s0$s0
  # Determine a sensible Tmax unless one was provided
  # (The original Tmax of the basemodel that the parametrizations was created with
  # is lost by san_parametrization)
  rl_s0$rateslist[[length(rl_s0$rateslist)]]$Tmax <- if (!is.na(Tmax))
    Tmax
  else max(
    parameterization$cc_days,
    parameterization$rs_days,
    parameterization$sc_days
  )
  
  # Return SANModel instance
  return(san_model(rbindlist(rl_s0$rateslist), s0))
}

#' Adapt a SAN model parametrization (SANP instance) to a specific dataset
#' 
#' @export
autodetect.metaparameters.SANParametrization <- function(parameterization, lt) {
  # Auto-detect 
  if (is.na(parameterization$s0.base)) {
    day0 <- lt$organoidsizes[day==0]
    if (nrow(day0) < 1)
      stop("No organoidsizes for day 0, unable to auto-detect s0.base")
    parameterization$s0.base <- round(median(day0$cells))
    message("Auto-detected s0 to be ", parameterization$s0.base)
  }
  
  if (is.null(parameterization$cc_days))
    parameterization$cc_days <- sort(unique(lt$organoidsizes$day))

  if (is.null(parameterization$rs_days))
    parameterization$rs_days <- sort(unique(lt$lineagesizes$day))

  if (is.null(parameterization$rs_ranks))
    parameterization$rs_ranks <- c(1, 2, 5, 10, 15, 25, 40, 60, 100, 150, 250, 400, 600, 1000, 1500, 2500, 4000, 6000, 10000, 15000, 25000)

  if (is.null(parameterization$sc_days))
    parameterization$sc_days <- parameterization$rs_days
  
  return(parameterization)
}

#' Create a posterior distribution instance from a SAN model parametrization and a dadtaset
#' 
#' @export
san_posterior<- function(parametrization, lt, cc.cutoff=1e7, p.cutoff=1e-2, ll.site.min=-Inf,
                         min.size=0.1, min.logsd=0.1, seqsim.coarse=FALSE)
{
  # Fill in dataset-dependent parameters
  parametrization <- autodetect.metaparameters(parametrization, lt)

  # Compute log-CIs for cell counts
  logcis.cc <- lt$organoidsizes[list(parametrization$cc_days), list(
    logmean=pmax(mean(log10(cells), na.rm=TRUE), log10(min.size), na.rm=TRUE),
    logsd=pmax(sd(log10(cells), na.rm=TRUE), min.logsd, na.rm=TRUE),
    samples=sum(is.finite(cells))
  ), on=.(day), keyby=.EACHI]

  # Compute log-CIs for lineage rank-sizes
  rs <- rank_size(lt$lineagesizes[list(parametrization$rs_days), list(day, sid, size=eval(as.name(lt$unit))), on=.(day)])
  logcis.rs <- rs[CJ(day=parametrization$rs_days, rank=parametrization$rs_ranks), list(
    logmean=pmax(mean(log10(size), na.rm=TRUE), log10(min.size), na.rm=TRUE),
    logsd=pmax(sd(log10(size), na.rm=TRUE), min.logsd, na.rm=TRUE),
    samples=sum(is.finite(size))
  ), on=.(day, rank), keyby=.EACHI]

  # Extract days at which the rates change
  params.na <- rep(NA_real_, length(parametrization$names))
  names(params.na) <- parametrization$names
  days_ratechanges <- lapply(rateslist_s0(parametrization, params.na)$rateslist, function(r) r$Tmax)
  # Extract list of days at which to compute cell-counts
  day_prev <- -1
  cc_days <- lapply(days_ratechanges, function(day) {
    r <- parametrization$cc_days[(day_prev < parametrization$cc_days) & (parametrization$cc_days <= day)]
    day_prev <<- day
    r
  })
  # Prepare vectors of cell-count log-means and log-sds to compute likelihoods
  stopifnot(logcis.cc[, day] == unlist(cc_days))
  cc_days_logmean <- logcis.cc$logmean
  cc_days_logsd <- logcis.cc$logsd
  if (!all(is.finite(cc_days_logmean)) || !all(is.finite(cc_days_logsd)))
    stop("LTData does not contain sufficient organoid size data for days ",
         paste0(parametrization$cc_days[!is.finite(cc_days_logmean) | !is.finite(cc_days_logsd)],
                collapse=", "))

  # Extract list of days at which to compute S-cell-counts
  day_prev <- -1
  sc_days <- lapply(days_ratechanges, function(day) {
    r <- parametrization$sc_days[(day_prev < parametrization$sc_days) & (parametrization$sc_days <= day)]
    day_prev <<- day
    r
  })
  
  # Extract list of days at which to compute rank-sizes
  # Prepare vectors of size log-means and log-sds to compute likelihoods
  rs_days <- parametrization$rs_days
  rs_ranks <- parametrization$rs_ranks
  rs_days_logmean <- lapply(rs_days, function(rs_day) {
    if (logcis.rs[day==rs_day, max(samples)] == 0)
      stop("LTData does not contain sufficient lineagesize data for day ", rs_day)
    stopifnot(logcis.rs[day==rs_day, rank] == rs_ranks)
    logcis.rs[day==rs_day, logmean]
  })
  rs_days_logsd <- lapply(rs_days, function(rs_day) {
    stopifnot(logcis.rs[day==rs_day, rank] == rs_ranks)
    logcis.rs[day==rs_day, logsd]
  })

  # PCR and sequencing parameters
  rs_libsize <- as.list(lt$sequencing[list(rs_days), ceiling(median(library_size)),
                                      on=.(day), by=.EACHI, nomatch=NULL]$V1)
  rs_pcreff <- as.list(lt$sequencing[list(rs_days), median(pcr_efficiency),
                                     on=.(day), by=.EACHI, nomatch=NULL]$V1)
  rs_th <- as.list(lt$sequencing[list(rs_days), ceiling(median(phantom_threshold)),
                                 on=.(day), by=.EACHI, nomatch=NULL]$V1)
  rs_aliaslambda <- as.list(lt$sequencing[list(rs_days), tpois.lambda(median(lineage_aliases)),
                                          on=.(day), by=.EACHI, nomatch=NULL]$V1)

  # Result vector template
  res_ll_rs <- rep(list(NA_real_), length(rs_days))
  if (length(parametrization$rs_days) > 0)
    names(res_ll_rs) <- paste0("ll_rs_", rs_days)
  res_cc <- rep(list(NA_real_), length(parametrization$cc_days))
  if (length(parametrization$cc_days) > 0)
    names(res_cc) <- paste0("CC", as.character(parametrization$cc_days))
  res_sc <- rep(list(NA_real_), length(parametrization$sc_days))
  if (length(parametrization$sc_days) > 0)
    names(res_sc) <- paste0("SC", as.character(parametrization$sc_days))
  res_rs <- rep(list(NA_real_), length(rs_days)*length(rs_ranks))
  if (length(parametrization$rs_days) > 0)
    names(res_rs) <- as.vector(outer(rs_ranks, rs_days,
                                     function(r, d) { paste0("RS", d, "_", r) }))
  res <- c(list(ll_tot=0, ll_cc=NA_real_), res_ll_rs, res_cc, res_sc, res_rs)
  
  # Create environment to evaludate functions in
  env <- new.env(parent=.SANjar.env)
  env$parametrization <- parametrization
  env$cc.cutoff <- cc.cutoff
  env$p.cutoff <- p.cutoff
  env$ll.site.min <- ll.site.min
  env$min.size <- min.size
  env$min.logsd <- min.logsd
  env$unit <- lt$unit
  env$cc_days <- cc_days
  env$cc_days_logmean <- cc_days_logmean
  env$cc_days_logsd <- cc_days_logsd
  env$rs_days <- rs_days
  env$rs_ranks <- rs_ranks
  env$rs_days_logmean <- rs_days_logmean
  env$rs_days_logsd <- rs_days_logsd
  env$sc_days <- sc_days
  env$res <- res
  env$res_ll_rs <- res_ll_rs
  env$res_cc <- res_cc
  env$res_sc <- res_sc
  env$res_rs <- res_rs
  env$rs_libsize <- rs_libsize
  env$rs_pcreff <- rs_pcreff
  env$rs_th <- rs_th
  env$rs_aliaslambda <- rs_aliaslambda
  env$method_seqsim <- if (seqsim.coarse) "gamma" else NULL

  # Compute cellcounts
  cellcounts <- function(s0, rl) {
    x0 <- c(S=s0, A=0, N=0)
    t0 <- 0
    cc_list <- rep(list(NULL), length(rl))
    sc_list <- rep(list(NULL), length(rl))
    for(ri in seq_along(rl)) {
      # Times (relative to the previous stopping time t0) at which to evaluate the model
      ts <- c(cc_days[[ri]] - t0, sc_days[[ri]] - t0, rl[[ri]]$Tmax - t0)
      ti_cc <- seq_along(cc_days[[ri]])
      ti_sc <- length(cc_days[[ri]]) + seq_along(sc_days[[ri]])
      ti_tmax <- length(ts)
      # Evaluate model
      san.out <- san_deterministic_eval_fixedrates(x0=x0, times=ts, rates=rl[[ri]])
      # Extract results (sum S, A, N to get total cellcount)
      cc_list[[ri]] <- san.out[ti_cc, c("S", "A", "N")] %*% c(1,1,1)
      sc_list[[ri]] <- unname(san.out[ti_sc, "S"])
      # Update final time and state
      x0 <- san.out[ti_tmax,]
      t0 <- rl[[ri]]$Tmax
    }

    return(list(cc=unlist(cc_list), sc=unlist(sc_list)))
  }
  environment(cellcounts) <- env
  env$cellcounts <- cellcounts
  
  # Simulate lineage sizes and compute rank-sizes
  ranksizes_aliases_pcrseq <- function(san.out, day, param_i) {
    # Extract total size sizes
    ls <- san.out[list(day), S+A+N, on=.(t)]
    # Simulate lineage aliasing (multiple labels per lineage)
    ls.labs <- if (is.finite(rs_aliaslambda[[param_i]]) && (rs_aliaslambda[[param_i]] > 0))
      rep(ls, extraDistr::rtpois(length(ls), lambda=rs_aliaslambda[[param_i]], a=0))
    else
      ls
    # Simulate sequencing using gwpcR
    ls.reads <- seqsim(ls.labs, reads.target=ceiling(as.numeric(rs_libsize[[param_i]])),
                       efficiency=rs_pcreff[[param_i]], method=method_seqsim)
    # Apply threshold, set read count to zero for filtered lineages
    ls.reads[ls.reads < rs_th[[param_i]]] <- 0
    if (unit == "cells") {
      # Normalize reads to cells to produce observed lineage sizes
      ls.obs <- ls.reads * sum(ls.labs) / (rs_libsize[[param_i]])
    } else if (unit == "reads") {
      # Use raw number of reads as "observed lineage size"
      ls.obs <- ls.reads
    }
    # Rank-sizes for selected ranks
    ifelse(rs_ranks <= length(ls.obs), sort(ls.obs, decreasing=TRUE)[rs_ranks], 0)
  }
  environment(ranksizes_aliases_pcrseq) <- env
  env$ranksizes_aliases_pcrseq <- ranksizes_aliases_pcrseq
  
  # Evaluate partial log-likelihood for cellcounts
  ll_cc <- function(cc) {
    ll.site <- -0.5 * ((log10(pmax(cc, min.size)) - cc_days_logmean) / cc_days_logsd)^2
    sum(pmax(ll.site, ll.site.min))
  }
  environment(ll_cc) <- env
  env$ll_cc <- ll_cc
  
  # Evaluate partial log-likelihood for rank-sizes
  ll_rs <- function(rs, rs_i) {
    ll.site <- -0.5 * ((log10(pmax(rs, min.size)) - rs_days_logmean[[rs_i]]) / rs_days_logsd[[rs_i]])^2
    sum(pmax(ll.site, ll.site.min))
  }
  environment(ll_rs) <- env
  env$ll_rs <- ll_rs
  
  # Evaluate total log-likelihood for a single rate vector
  ll_evaluate <- function(params, ll_cutoff) {
    # Get initial cell count and rates table
    m <- rateslist_s0(parametrization, params)
    rl <- m$rateslist
    s0 <- m$s0

    # Evaluate deterministic SAN model
    # Exit early if the number of rcells grows unreasonably large
    cc.out <- cellcounts(s0, rl)
    stopifnot(length(names(res_cc)) == length(cc.out$cc))
    res[names(res_cc)] <- cc.out$cc
    stopifnot(length(names(res_sc)) == length(cc.out$sc))
    res[names(res_sc)] <- cc.out$sc
    if (max(cc.out$cc) > cc.cutoff) {
      res[["ll_tot"]] <- (-Inf)
      return(res)
    }
    
    # Compute partial likelihood L_cc for cell counts
    # Exit early if the likelihood is unreasonably small
    ll_part <- ll_cc(cc.out$cc)
    res[["ll_tot"]] <- res[["ll_tot"]] + ll_part
    res[["ll_cc"]] <- ll_part
    # Early abort
    if (res[["ll_tot"]] < ll_cutoff) {
      res[["ll_tot"]] <- (-Inf)
      return(res)
    }

    # Evaluate stochastic SAN model and compute partial likelihoods L_rs_<day> for rank-sizes
    # Exit early if combined partial likelihoods become unreasonably small
    rt <- rbindlist(rl)
    san.out <- NULL
    for(rs_i in seq_along(rs_days)) {
      rs_day <- rs_days[rs_i]
      # Evaluate model
      san.out <- tryCatch(
        san_stochastic(L=s0, s0=1, previous=san.out, rates=rt,
                       Tmax=rs_day, samples_per_day=1, p_cutoff=p.cutoff),
        error=function(e) {
          # If cell counts overflow set likelihood to -infinity (see alow below)
          if (grepl("^S, A or N cell count overflowed", conditionMessage(e)))
            return(NULL)
          signalCondition(e)
        }
      )
      if (is.null(san.out)) {
        res[["ll_tot"]] <- (-Inf)
        return(res)
      }
      # Simulate sequencing and compute rank-sizes
      rs <- tryCatch(
        ranksizes_aliases_pcrseq(san.out, day=rs_day, param_i=rs_i),
        error=function(e) {
          # This shouldn't fail, but if it does, don't abort the whole MCMC run,
          # but instead log the error and return -Inf as the likelhood.
          f <- paste0(tempfile(pattern = "ranksizes_aliases_pcrseq.", tmpdir = "."), ".rd")
          warning("Caught an unexpected error while simulating sequencing!\n",
                  "  Error: ", conditionMessage(e), "\n",
                  "  Parameters: ", paste0(names(params), "=", signif(params, 6), collapse=", "), "\n",
                  "  Dump: ", f)
          tryCatch(
            save(san.out, params, rs_day, rs_i, rs_libsize, rs_pcreff, rs_th, unit, rs_ranks,
                 file=f),
            error=function(e) {
              warning("Unable to dump to file ", f)
            }
          )
          return(NULL)
        }
      )
      if (is.null(rs)) {
        res[["ll_tot"]] <- (-Inf)
        return(res)
      }
      res[names(res_rs)[seq(to=rs_i*length(rs_ranks), along.with=rs_ranks, by=1)]] <- rs
      # Compute partial likelihood
      ll_part <- ll_rs(rs, rs_i)
      res[["ll_tot"]] <- res[["ll_tot"]] + ll_part
      res[[names(res_ll_rs)[rs_i]]] <- ll_part
      # Early abort
      if (res[["ll_tot"]] < ll_cutoff) {
        res[["ll_tot"]] <- (-Inf)
        return(res)
      }
    }
    
    return(res)
  }
  environment(ll_evaluate) <- env
  env$ll_evaluate <- ll_evaluate
  
  # Evaluate total log-likelihood for each row of the matrix `params`
  loglikelihood <- function(params, cutoffs=-Inf) {
    stopifnot(is.parameter.vector(params) || is.parameter.matrix(params))
    stopifnot((length(cutoffs) == 1) || (is.parameter.matrix(params) && (length(cutoffs) == nrow(params))))
    if (is.parameter.vector(params))
      ll_evaluate(params, cutoffs)
    else
      rbindlist(mapply(ll_evaluate, asplit(params, MARGIN=1), as.list(cutoffs), SIMPLIFY=FALSE),
                use.names=FALSE)
  }
  environment(loglikelihood) <- env
  env$loglikelihood <- loglikelihood
  
  return(structure(list(
    loglikelihood=loglikelihood,
    parametrization=parametrization,
    parameters=parametrization$ranges,
    auxiliary=list(ll=c("ll_cc", names(res_ll_rs)), cc=names(res_cc), sc=names(res_sc), rs=names(res_rs)),
    data=list(cc=logcis.cc, rs=logcis.rs, unit=lt$unit),
    arguments=list(cc.cutoff=cc.cutoff, p.cutoff=p.cutoff, ll.site.min=ll.site.min,
                   min.size=min.size, min.logsd=min.logsd, seqsim.coarse=seqsim.coarse)
  ), class="SANPosterior"))
}

#' Combines multiple posteriors by summing up the likelihoods of the components
#'
#' @export
san_posterior_combine <- function(..., components=list(), parametrization=NULL, args=list()) {
  # Handle components of combined posterior distribution
  components <- c(list(...), components)
  if ((length(components) > 1) && (is.null(names(components)) || any(names(components) == "")))
    stop("san_posterior_combine expects datasets and posterior distribution arguments to be named")
  if (is.null(parametrization)) {
    # Expect posterior distribution arguments
    if (!all(sapply(components, function(p) "SANPosterior" %in% class(p))))
      stop("if no parametrization is specified, san_posterior_combine() expects instances of SANPosterior (e.g. created with san_posterior())")
  } else {
    # Expect datasets distribution arguments
    if (!all(sapply(components, function(p) "LTData" %in% class(p))))
      stop("if a parametrization is specified, san_posterior_combine() expects instances of LTData (such as lt74)")
    # Construct posterior distribution objects from datasets
    components <- lapply(components, function(lt) do.call(san_posterior, c(list(parametrization, lt), args)))
  }

  # Check that the parameters are the same
  parameters <- NULL
  for(p in components) {
    if (is.null(parameters))
      parameters <- p$parameters
    
    if ((length(p$parameters) != length(parameters)) || any(names(p$parameters) != names(parameters)))
      stop("posterior distributions to be combined must have the same parameters")
  }
  
  # Create environment to evaluate functions in
  env <- new.env(parent=.SANjar.env)
  env$loglikelihoods <- lapply(components, function(p) p$loglikelihood)
    
  # Evaluate total log-likelihood for each row of the matrix `params`
  loglikelihood <- function(params, cutoffs=-Inf) {
    stopifnot(is.parameter.vector(params) || is.parameter.matrix(params))
    stopifnot((length(cutoffs) == 1) || (is.parameter.matrix(params) && (length(cutoffs) == nrow(params))))
    if (is.parameter.matrix(params) && (length(cutoffs) == 1))
      cutoffs <- rep(cutoffs, nrow(params))
    # Evaluate individual posterior likelihoods
    lls <- lapply(loglikelihoods, function(ll) {
      r <- ll(params, cutoffs)
      cutoffs <<- ifelse(is.finite(r$ll_tot), cutoffs - r$ll_tot, Inf)
      r
    })
    # Compute total likelihood
    ll_tot <- Reduce(function(ll_tot, ll) ll_tot + ll$ll_tot, lls, 0)
    # Return either list or table containing the total likelihood and all auxiliary data
    # from all individual likelihoods, with column names prefixed with the likelihood label
    result <- c(list(ll_tot=ll_tot), lls)
    if (is.parameter.vector(params))
      do.call(c, result)
    else
      do.call(cbind, result)
  }
  environment(loglikelihood) <- env
  env$loglikelihood <- loglikelihood
  
  return(structure(list(
    loglikelihood=loglikelihood,
    parameters=parameters,
    auxiliary=list(ll=unlist(lapply(names(components), function(n) paste0(n, ".", c("ll_tot", components[[n]]$auxiliary$ll)))),
                   cc=unlist(lapply(names(components), function(n) paste0(n, ".", components[[n]]$auxiliary$cc))),
                   rs=unlist(lapply(names(components), function(n) paste0(n, ".", components[[n]]$auxiliary$sc))),
                   sc=unlist(lapply(names(components), function(n) paste0(n, ".", components[[n]]$auxiliary$rs)))),
    components=components
  ), class=c("SANCombinedPosterior", "SANPosterior")))
}

#' Parallelizes likelihood evaluations for different parameter combinations
#' 
#' @export
san_posterior_parallel <- function(posterior, cluster) {
  if (!("SANPosterior" %in% class(posterior)))
    stop("posterior argument must be an instance of SANPosterior (e.g. created with san_posterior())")

  if (is.null(cluster))
    return(posterior)
  
  # Send likelihood evaluation function to cluster
  id <- paste0(lapply(1:4, function(i) sprintf("%04x", as.integer(sample.int(size=1, n=2**16)-1))), collapse="")
  loglikelihood.worker.name <- paste(".san_posterior_parallel", id, "loglikelihood.worker", sep=".")
  setup.worker <- function(i) {
    # Load libraries on worker
    library(data.table)
    library(SANjar)
    library(gwpcR)

    # Assign name in the global environment to the likelihood function 
    if (exists(loglikelihood.worker.name, envir=globalenv()))
      stop("name ", loglikelihood.worker.name, " is already taken on node ", i)
    assign(loglikelihood.worker.name, loglikelihood.worker.impl, envir=globalenv())
  }
  environment(setup.worker) <- list2env(list(
    loglikelihood.worker.name=loglikelihood.worker.name,
    loglikelihood.worker.impl=posterior$loglikelihood
  ), parent=baseenv())
  clusterApply(cluster, 1:length(cluster), setup.worker)
  
  # Create wrapper function to evaluate likelihoods on workers
  loglikelihood.worker <- eval(bquote(function(...) .(as.name(loglikelihood.worker.name))(...)))
  environment(loglikelihood.worker) <- globalenv()

  # Create parallel-evaluation log-likelihood function
  env <- new.env(parent=.SANjar.env)
  env$cluster <- cluster
  env$loglikelihood.worker <- loglikelihood.worker
  env$loglikelihood.serial <- posterior$loglikelihood
  posterior$loglikelihood <- function(params, cutoffs=-Inf) {
    stopifnot(is.parameter.vector(params) || is.parameter.matrix(params))
    stopifnot((length(cutoffs) == 1) || (is.parameter.matrix(params) && (length(cutoffs) == nrow(params))))
    if (is.parameter.vector(params))
      loglikelihood.serial(params, cutoffs)
    else
      rbindlist(parMLapply(cl=cluster, loglikelihood.worker, asplit(params, MARGIN=1), as.list(cutoffs)),
                use.names=FALSE)
  }
  environment(posterior$loglikelihood) <- env

  return(posterior)
}

is.parameter.vector <- function(x) {
  is.vector(x) || (is.array(x) && (length(dim(x)) == 1))
}

is.parameter.matrix <- function(x) {
  is.data.table(x) || (is.array(x) && (length(dim(x)) == 2))
}

.SANjar.env <- parent.frame()$envir
