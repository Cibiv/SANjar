#' Create SAN model parametrization
#'
#' @export
san_parametrization <- function(ranges, basemodel, match.s0.to.data=TRUE, cc_days=NULL, rs_days=NULL, rs_ranks=NULL, sc_days=NULL) {
  # Validate parameter names (I)
  pnames <- names(ranges)
  pvalid <- grepl("^([0-9]+)([A-Z]+)$", pnames)
  if (!all(pvalid))
    stop("invalid parameter names ", paste0(pnames[!pvalid], collapse=", "))
  
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
    r <- gsub("0", "X", sub(pattern, "\\1", p))
    # Create map, i.e. a vector with rate names as values and target columns as names
    k <- nchar(r)
    map <- rep(p, k)
    if (length(map) > 0)
      names(map) <- paste0("r_", unlist(strsplit(r, split="")))
    # Check that target columns exist (II)
    if (!(all(names(map) %in% colnames(basemodel$rates))))
      stop("invalid parameter names ", paste0(map[!(names(map) %in% colnames(basemodel$rates))], collapse=", "))
    return(map)
  })
  # Check that all parameter names could be translated (III)
  if (!all(pnames %in% unlist(rateslist.pmap)))
    stop("invalid parameter names ", paste0(pnames[!(pnames %in% unlist(rateslist.pmap))], collapse=", "))

  # Convert default ratetable into list, and replace the last Tmax by infinity
  # to allow simulations to run arbitrarily long
  rateslist.base <- lapply(split(ratestable.fill.na(basemodel$rates), by="Tmax"), as.list)
  rateslist.base[[length(rateslist.base)]]$Tmax <- Inf
  
  # Create SANParametrization object
  return(structure(list(
    names=as.character(pnames),
    ranges=as.list(ranges),
    rateslist.base=rateslist.base,
    s0.base=if (match.s0.to.data) NA else as.integer(basemodel$s0),
    rateslist.pmap=rateslist.pmap,
    cc_days=if (!is.null(cc_days)) as.integer(cc_days) else NULL,
    rs_days=if (!is.null(rs_days)) as.integer(rs_days) else NULL,
    rs_ranks=if (!is.null(rs_ranks)) as.integer(rs_ranks) else NULL,
    sc_days=if (!is.null(sc_days)) as.integer(sc_days) else NULL
  ), class="SANParametrization"))
}

#' Convert a parameter vector to a rate list and inital cell count (s0)
#' 
#' @export
rateslist_s0 <- function(parameterization, params) UseMethod("rateslist_s0")

#' Adapt a SAN model parametrization (SANP instance) to a specific dataset
#' 
#' @export
adapt <- function(parameterization, ...) UseMethod("adapt")

#' Convert a parameter vector to a model instance
#' 
#' @export
rateslist_s0.SANParametrization <- function(parameterization, params) {
  rl <- parameterization$rateslist.base
  rateslist.pmap <- parameterization$rateslist.pmap
  for(i in 1:length(rl)) {
    m <- rateslist.pmap[[i]]
    if (length(m) > 0)
      rl[[i]][names(m)] <- params[m]
  }
  
  s0 <- parameterization$s0.base
  if ("s0" %in% parameterization$names)
    s0 <- params["s0"]
  
  return(list(s0=s0, rateslist=rl))
}

#' Adapt a SAN model parametrization (SANP instance) to a specific dataset
#' 
#' @export
adapt.SANParametrization <- function(parameterization, lt) {
  # Auto-detect 
  if (is.na(parameterization$s0.base)) {
    day0 <- lt$organoidsizes[day==0]
    if (nrow(day0) < 1)
      stop("No organoidsizes for day 0, unable to auto-detect s0.base")
    parameterization$s0.base <- round(median(day0$cells))
    message("Auto-detected s0 to be ", parameterization$s0.base)
  }
  
  if (is.null(parameterization$cc_days)) {
    parameterization$cc_days <- sort(unique(lt$organoidsizes$day))
  }

  if (is.null(parameterization$rs_days))
    parameterization$rs_days <- sort(unique(lt$lineagesizes$day))

  if (is.null(parameterization$rs_ranks))
    parameterization$rs_ranks <- c(1, 2, 5, 10, 15, 25, 40, 60, 100, 150, 250, 400, 600, 1000, 1500, 2500, 4000, 6000, 10000, 15000, 25000)
  
  return(parameterization)
}

#' Create a posterior distribution instance from a SAN model parametrization and a dadtaset
#' 
#' @export
san_posterior<- function(parametrization, lt, cc.cutoff=1e7, p.cutoff=1e-2, ll.site.min=-Inf,
                         min.size=0.1, min.logsd=0.1, cluster=NULL)
{
  # Setup cluster
  if (!is.null(cluster)) {
    clusterEvalQ(cluster, {
      library(OpenMPController)
      library(data.table)
      library(SANjar)
      library(gwpcR)
      
      # No nested parallelism
      data.table::setDTthreads(1)
      OpenMPController::omp_set_num_threads(1)
    })
  }
  MLapply <- if (is.null(cluster)) {
    function(FUN, ...) mapply(FUN, ..., SIMPLIFY=FALSE)
  } else
    function(FUN, ...) parMLapply(cl=cluster, FUN, ...)
  
  # Fill in dataset-dependent parameters
  parametrization <- adapt(parametrization, lt)

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
                                      on=.(day), by=.EACHI]$V1)
  rs_pcreff <- as.list(lt$sequencing[list(rs_days), median(pcr_efficiency),
                                     on=.(day), by=.EACHI]$V1)
  rs_th <- as.list(lt$sequencing[list(rs_days), ceiling(median(phantom_threshold)),
                                 on=.(day), by=.EACHI]$V1)
  rs_aliaslambda <- as.list(lt$sequencing[list(rs_days), tpois.lambda(median(lineage_aliases)),
                                          on=.(day), by=.EACHI]$V1)

  # Result vector template
  res_ll_rs <- rep(list(NA_real_), length(rs_days))
  names(res_ll_rs) <- paste0("ll_rs_", rs_days)
  res_cc <- rep(list(NA_real_), length(parametrization$cc_days))
  if (length(parametrization$cc_days) > 0)
  names(res_cc) <- paste0("CC", as.character(parametrization$cc_days))
  res_sc <- rep(list(NA_real_), length(parametrization$sc_days))
  if (length(parametrization$sc_days) > 0)
    names(res_sc) <- paste0("SC", as.character(parametrization$sc_days))
  res_rs <- rep(list(NA_real_), length(rs_days)*length(rs_ranks))
  names(res_rs) <- as.vector(outer(rs_ranks, rs_days,
                                   function(r, d) { paste0("RS", d, "_", r) }))
  res <- c(list(ll_tot=0, ll_cc=NA_real_), res_ll_rs, res_cc, res_sc, res_rs)
  
  # Compute cellcounts
  cellcounts <- function(s0, rl) {
    x0 <- c(S=s0, A=0, N=0)
    t0 <- 0
    cc_list <- rep(list(NULL), length(rl))
    sc_list <- rep(list(NULL), length(rl))
    for(ri in 1:length(rl)) {
      # Times (relative to the previous stopping time t0) at which to evaluate the model
      ts <- c(cc_days[[ri]] - t0, sc_days[[ri]] - t0, rl[[ri]]$Tmax - t0)
      ti_cc <- 1:length(cc_days[[ri]])
      ti_sc <- length(cc_days[[ri]]) + 1:length(sc_days[[ri]])
      ti_tmax <- length(ts)
      # Evaluate model
      san.out <- san_deterministic_eval_fixedrates(x0=x0, times=ts, rates=rl[[ri]])
      # Extract results (sum S, A, N to get total cellcount)
      cc_list[[ri]] <- san.out[ti_cc, c("S", "A", "N")] %*% c(1,1,1)
      if (length(sc_days[[ri]]) > 0)
        sc_list[[ri]] <- san.out[ti_sc, "S"]
      # Update final time and state
      x0 <- san.out[ti_tmax,]
      t0 <- rl[[ri]]$Tmax
    }

    return(list(cc=unlist(cc_list), sc=unlist(sc_list)))
  }
  
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
                       efficiency=rs_pcreff[[param_i]])
    # Apply threshold, set read count to zero for filtered lineages
    ls.reads[ls.reads < rs_th[[param_i]]] <- 0
    if (lt$unit == "cells") {
      # Normalize reads to cells to produce observed lineage sizes
      ls.obs <- ls.reads * sum(ls.labs) / (rs_libsize[[param_i]])
    } else if (lt$unit == "reads") {
      # Use raw number of reads as "observed lineage size"
      ls.obs <- ls.reads
    }
    # Rank-sizes for selected ranks
    ifelse(rs_ranks <= length(ls.obs), sort(ls.obs, decreasing=TRUE)[rs_ranks], 0)
  }
  
  # Evaluate partial log-likelihood for cellcounts
  ll_cc <- function(cc) {
    ll.site <- -0.5 * ((log10(pmax(cc, min.size)) - cc_days_logmean) / cc_days_logsd)^2
    sum(pmax(ll.site, ll.site.min))
  }
  
  # Evaluate partial log-likelihood for rank-sizes
  ll_rs <- function(rs, rs_i) {
    ll.site <- -0.5 * ((log10(pmax(rs, min.size)) - rs_days_logmean[[rs_i]]) / rs_days_logsd[[rs_i]])^2
    sum(pmax(ll.site, ll.site.min))
  }
  
  # Evaluate total log-likelihood for a single rate vector
  ll_params <- function(params, ll_cutoff) {
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
    for(rs_i in 1:length(rs_days)) {
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
            save(san.out, params, rs_day, rs_i, rs_libsize, rs_pcreff, rs_th, lt$unit, rs_ranks,
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
  
  # Evaluate total log-likelihood for each row of the matrix `params`
  loglikelihood <- function(params, cutoffs=-Inf) {
    stopifnot((length(cutoffs) == 1) || (length(cutoffs) == nrow(params)))
    rbindlist(MLapply(ll_params, asplit(params, MARGIN=1), as.list(cutoffs)))
  }
  
  return(structure(list(
    loglikelihood=loglikelihood, ll_params=ll_params,
    cellcounts=cellcounts, ranksizes_aliases_pcrseq=ranksizes_aliases_pcrseq,
    ll_cc=ll_cc, ll_rs=ll_rs,
    cc_days=parametrization$cc_days, sc_days=sc_days,
    rs_days=rs_days, rs_ranks=rs_ranks,
    parameters=parametrization$ranges,
    parametrization=parametrization,
    unit=lt$unit,
    logcis_cc=logcis.cc, logcis_rs=logcis.rs
  ), class="SANPosterior"))
}
