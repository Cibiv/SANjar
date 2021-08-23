#' The MCMC results for the LT47 dataset
"lt47.mcmc"

#' Find a suitable bandwidth matrix for kernel density estimation 
#' 
#' @export
bandwidth.matrix <- function(x, ...) UseMethod("bandwidth.matrix")

#' Locates the distribution's mode(s) using the mean-shift algorithm
#'
#' @export
locate.modes <- function(x, ...) UseMethod("locate.modes")

#' Computes various summary statistics
#' 
#' @export
summarystats <- function(x, ...) UseMethod("summarystats")

#' Samples from a distribution using the Metropolis-Hastings MCMC algorithm  
#' 
#' @export
mcmc <- function(llfun, variables, fixed=character(),
                 steps=NA, accepts=NA, chains, candidate.samples=chains,
                 candidates=NULL, initial.samplewithreplacement=FALSE, initial=NULL, 
                 keep.history=FALSE, keep.initial=FALSE, keep.candidates=FALSE,
                 reevaluate.likelihood=FALSE, acceptance.target=0.234,
                 initial.Rproposal=1.0, minimal.Rproposal=0.2,
                 initial.Vproposal=NULL, minimal.Vproposal.CV=0.1, maximal.Vproposal.corr=0.95,
                 proposal.update.rate=0.2, delta.ll.cutoff=15, initial.ll.cutoff=NA,
                 verbose=FALSE)
{
  # Negative values would reject any proposal that doesn't *decrease* the likelihood sufficiently
  stopifnot(delta.ll.cutoff >= 0)
  
  # Either steps or accepts must be specified
  stopifnot(!is.na(steps) || !is.na(accepts))
  
  # Prepare a list of variable names, ensuring that lapply(varnames) will return a *named* list
  varnames <- names(variables)
  names(varnames) <- varnames
  
  # Initial state distribution
  if (is.null(initial)) {
    if (is.null(candidates)) {
      if (length(fixed) > 0)
        stop("If some variables are held fixed, candidate samples must be provided")
      # Sample parameter space
      if (verbose)
        message("Generating ", candidate.samples, " candidate samples")
      candidates <- sample.variables(candidate.samples, variables)
    }
    # Determine initial.ll.cutoff
    if (is.na(initial.ll.cutoff)) {
      candidates.cutoffest <- candidates[1:ceiling(sqrt(.N))]
      if (verbose)
        message("Evaluating posterior likelihood for ", nrow(candidates.cutoffest), " candidate samples to estimate initial.ll.cutoff")
      r <- as.data.table(llfun(candidates.cutoffest, cutoffs=-Inf))
      initial.ll.cutoff <- max(r[, 1])
    }
    # Evaluate likelihood
    if (verbose)
      message("Evaluating posterior likelihood for ", nrow(candidates), " candidate samples ",
              "with initial.ll.cutoff ", signif(initial.ll.cutoff, 3))
    r <- as.data.table(llfun(candidates, cutoffs=initial.ll.cutoff))
    r.ll <- r[[ colnames(r)[1] ]]
    if (ncol(r) > 1)
      candidates.meta <- r[, 2:ncol(r)]
    else
      candidates.meta <- NULL
    # Remove samples with log-likelihood -infinity
    r.ll.finite <- is.finite(r.ll)
    candidates <- candidates[r.ll.finite]
    candidates[, ll := r.ll[r.ll.finite] ]
    if (!is.null(candidates.meta))
      candidates.meta <- candidates.meta[r.ll.finite]
    if (nrow(candidates) == 0)
      stop("No candidate sample had finite likelihood, aborting")
    # Remove samples whose log-likelihood is too small to yield a positive probability when
    # exponentiated Since the latter could, in principle, affect *all* likelihoods,
    # log-likelihoods are normalized to make the largest log-likelihood 1 before exponentiation.
    # Since likelihoods are only defined up to a constant scaling factor, this additional#
    # scaling does not affect the results.
    ll.max <- candidates[, max(ll)]
    candidates[, ll.norm.exp :=  exp(ll-ll.max) ]
    ll.norm.exp.positive <- candidates[, (ll.norm.exp > 0)]
    candidates <- candidates[ll.norm.exp.positive]
    if (!is.null(candidates.meta))
      candidates.meta <- candidates.meta[ll.norm.exp.positive]
    if (verbose) {
      message("The following ", nrow(candidates), " candidate samples are usable:")
      message(paste0("  ", capture.output(print(signif(candidates, 3))), collapse="\n"))
    }
    # Sample initial states accordingly to likelihood from the uniformly spaced samples
    if (initial.samplewithreplacement) {
      # Sample the initial states
      if (verbose)
        message("Sampling initial states for ", chains, " chains from ", nrow(candidates), " candidates with replacement")
      i <- candidates[, sample.int(.N, size=chains, prob=ll.norm.exp, replace=TRUE) ]
    } else {
      # Repeat existing sample often enough to not run out of samples when sampling the initial states
      k <- ceiling(chains / nrow(candidates))
      if (verbose)
        message("Sampling initial states for ", chains, " chains from ", nrow(candidates), " candidates, ",
                "using each candidate at most ", k, " times")
      candidates <- rbindlist(rep(list(candidates), k))
      if (!is.null(candidates.meta))
        candidates.meta <- rbindlist(rep(list(candidates.meta), k))
      # Sample the initial states *without* replacement
      # The idea is that this ensures a reasonable dispersion of initial values
      i <- candidates[, sample.int(.N, size=chains, prob=ll.norm.exp, replace=FALSE) ]
    } 
    initial <- candidates[i, c(list(chain=1:.N, naccepts=0), .SD[, c("ll", fixed, varnames), with=FALSE]) ]
    if (!is.null(candidates.meta))
      initial.meta <- candidates.meta[i]
    else
      initial.meta <- NULL
  } else {
    # Initial states were specified
    if (!is.null(candidates))
      warning("Initial states were specified, provided candidates not being used")
    if (!missing(chains) && (chains != nrow(initial)))
      warning("Number of initial states (", nrow(initial), ") overrides specified number of chains (", chains, ")")
    chains <- nrow(initial)
    if (ncol(initial) > (length(fixed) + length(varnames) + 1))
      initial.meta <- initial[, .SD, .SDcols=!c("ll", fixed, varnames)]
    else
      initial.meta <- NULL
    initial <- initial[, c(list(chain=1:.N, naccepts=0), .SD[, c("ll", fixed, varnames), with=FALSE]) ]
  }
  if (verbose) {
    message("Initial states")
    initial.signif <- initial[, lapply(.SD, function(c) { if (is.numeric(c)) signif(c, 3) else c } )]
    message(paste0("  ", capture.output(print(initial.signif)), collapse="\n"))
  }
  
  # Initialize proposal distribution
  # If we started with ABC-style uniform sampling of the parameter space, we compute
  # the initial covariance estimate using all samples and their likelihoods, since
  # this gives a more precise estimate than using just the initial states.
  Vproposal <- if (is.null(initial.Vproposal)) {
    if (!is.null(candidates[, varnames, with=FALSE]))
      cov.wt(candidates[, varnames, with=FALSE], wt=candidates$ll.norm.exp, method="ML")$cov
    else
      cov(initial[, varnames, with=FALSE])
  } else initial.Vproposal
  Rproposal <- initial.Rproposal
  if (verbose) {
    message("Initial proposal distribution")
    message("  covariance matrix V:")
    message(paste0("    ", capture.output(print(signif(Vproposal, 3))), collapse="\n"))
    if (!is.null(acceptance.target))
      message("  radius R: ", signif(Rproposal, 3))
  }
  
  # Setup MCMC variables
  states <- initial
  states.meta <- initial.meta
  if (keep.history) {
    history <- list(states[, c(list(step=0), .SD) ])
  }

  # Prevent the variances of the proposal distribution from becoming too small
  Vproposal.clamp <- function() {
    # Compute vector of adjustment factors which make the i-th CV equal to minimal.Vproposal.CV
    sd.min <- sapply(varnames, function(v) { minimal.Vproposal.CV * mean(states[[v]]) })
    var <- diag(Vproposal)
    f <- sd.min / sqrt(var)
    if (any(f > 1)) {
      overshoot <- 2
      adj <- ifelse(f>1, overshoot*f, 1)
      if (verbose)
        message("Variances of proposal distribution clamped to ",
                paste0(varnames[f>1], "=", signif(overshoot*sd.min[f>1]^2, 3), collapse=", "))
      # Adjust for i=1,..,k, adjust i-th row and i-th column of V by the i-th adjustment factor
      M <- outer(adj, adj)
      # The adjustment may be infinite if V contained zeros. To handle that, first set any
      # non-finite entries in M to zero to avoid computing "0 * inf" which yields NaN
      M[!is.finite(M)] <- 0
      Vproposal <<- Vproposal * M
      # Then set any zeros in the diagonal of V to the desired minimal variance
      adj.inf <- !is.finite(adj)
      diag(Vproposal)[adj.inf] <<- overshoot*sd.min[adj.inf]^2
    }

    # Restrict correlations to at most maximal.Vproposal.corr
    if (maximal.Vproposal.corr == 0) {
      # Remove all covariances
      Vproposal <- diag(diag(Vproposal))
    } else if (maximal.Vproposal.corr == 1.0) {
      # Nothing to do
    } else {
      # Compute adjustment factors which ensure that no correlation exceeds maximal.Vproposal.corr
      Vcorr <- Vproposal / sqrt(diag(Vproposal) %o% diag(Vproposal))
      f <- Vcorr / maximal.Vproposal.corr
      diag(f) <- 1
      if (any(f > 1)) {
        # Find maximal adjustment necessary for every variable
        f.rowmax = apply(f, MARGIN=1, FUN=max)
        f.colmax = apply(f, MARGIN=2, FUN=max)
        s <- sqrt(pmax(f.rowmax, f.colmax, 1))
        # Conceptually, we now increase the variances of all variables according to the factors
        # in s^2 to reduce the correlations, and then scale the whole variable (i.e. its
        # row *and* column in V) down by s This restores the original variances, but keeps
        # the reduced correlations.
        M <- s %o% s
        diag(M) <- 1
        Vproposal <<- Vproposal / M
        if (verbose) {
          p <- (f>1) & upper.tri(Vproposal)
          message("Covariances of proposal distribution reduced to ",
                  paste0(outer(varnames, varnames, function(a,b) paste0("cov(",a,",",b,")"))[p],
                         "=", signif(Vproposal[p], digits=3), collapse=", "))
        }
      }
    }
  }
  Vproposal.clamp()

  # Run k Metropolis-Hasting MCMC steps
  step <- 0
  while(  ( is.na(steps) || (step < steps) ) && ( is.na(accepts) || (sum(states$naccepts < accepts) > 0) )  ) {
    step <- step + 1

    # Copy states, since data.table does not copy on update
    states <- copy(states)
    
    # Determine which chains need extension
    extend <- if (!is.na(accepts))
      states$naccepts < accepts
    else
      rep(TRUE, nrow(states))

    # Generate proposals by sampling displacements using current covariances and radius,
    # and displacing the previous samples accordingly. We only really generate proposals
    # for chains that need extension.
    displacements <- mvtnorm::rmvnorm(n=sum(extend), sigma=Vproposal) * Rproposal
    colnames(displacements) <- names(varnames)
    proposals <- data.table(valid=extend)
    if (length(fixed) > 0)
      proposals[valid==TRUE, (fixed) := states[extend, fixed, with=FALSE] ]
    proposals[valid==TRUE, (varnames) := lapply(varnames, function(v) {
      states[extend, eval(as.name(v))] + displacements[, v]
    })]

    # Accept or reject proposals
    # Proposals outside the parameter's domain are skipped.
    accept <- rep(FALSE, nrow(states))
    for(v in varnames)
      proposals[valid==TRUE, valid := valid & (variables[[v]][1] <= eval(as.name(v))) & (eval(as.name(v)) <= variables[[v]][2]) ]
    if (any(proposals$valid)) {
      # Re-evaluate likelihoods of states with valid proposals (if requested)
      if (reevaluate.likelihood) {
        r <- as.data.table(llfun(states[proposals$valid, c(fixed, varnames), with=FALSE],
                                 cutoffs=-Inf))
        stopifnot(is.finite(unlist(r[, 1])))
        if (!is.null(states.meta) && (ncol(r) > 1))
          states.meta[proposals$valid, colnames(r)[2:ncol(r)] := r[, 2:ncol(r)] ]
        else
          states.meta <- NULL
        states[proposals$valid, ll := r[, 1] ]
      }
      
      # Evaluate likelihood of valid proposals
      r <- as.data.table(llfun(proposals[valid==TRUE, c(fixed, varnames), with=FALSE],
                               cutoffs=states[proposals$valid, ll] - delta.ll.cutoff))
      # Separate meta information
      if (ncol(r) > 1)
        proposals.valid.meta <- r[, 2:ncol(r)]
      else
        proposals.valid.meta <- NULL
      # Set proposal likelihoods
      proposals[, ll := -Inf]
      proposals[valid==TRUE, ll := r[, 1]]
      
      # Accept or reject proposals using the Metropolis-Hastings rule
      mh <- (proposals$valid & is.finite(proposals$ll))
      accept[mh] <- (runif(sum(mh)) <= exp(proposals$ll[mh] - states$ll[mh]))
      # Update parameters if a proposal was accepted
      states[accept, c("ll", varnames) := proposals[accept, c("ll", varnames), with=FALSE] ]
      states[accept, naccepts := naccepts + 1]
      # Replace meta information if a proposal was accepted
      if (!is.null(states.meta) && !is.null(proposals.valid.meta)) {
        stopifnot(!accept[!proposals$valid])
        states.meta[accept, colnames(proposals.valid.meta) :=
                      proposals.valid.meta[accept[proposals$valid]] ]
      }
      else
        states.meta <- NULL
    }

    if (verbose) {
      message("States after step ", step)
      states.signif <- states[, lapply(.SD, function(c) { if (is.numeric(c)) signif(c, 3) else c } )]
      message(paste0("  ", capture.output(print(states.signif)), collapse="\n"))
      message("Chain extension attempts: ", sum(extend))
      message("Valid proposal ratio over extension attempts: ", signif(proposals[extend, mean(valid)], 3))
      message("Acceptance ratio over extension attempts: ", signif(mean(accept[extend]), 3))
      message("Average log-likelihood over all chains: ", signif(states[, mean(ll)], 3))
    }

    # Update the proposal distribution, smoothly over about proposal.update.rate steps.
    if (proposal.update.rate > 0) {
      Vproposal <- Vproposal * (1 - proposal.update.rate) + cov(states[, varnames, with=FALSE]) * proposal.update.rate
      Vproposal.clamp()
    }
    
    # Update the proposal radius to meet the acceptance target, smoothly over about proposal.update.rate steps.
    # Widening the proposal radius reduces the acceptance rate, we thus adjust the radius based on the ratio
    # of actual vs. targetted acceptance rate. To avoid the radius becoming zero (where it would then remain)
    # updates are limited to at most a factor of two, regardless of the update rate.
    if (!is.null(acceptance.target) && (proposal.update.rate > 0))
      Rproposal <- max(Rproposal * max(0.5, min((mean(accept[extend]) / acceptance.target)^proposal.update.rate, 2)), minimal.Rproposal)
    
    if (verbose && proposal.update.rate) {
      message("Updated proposal distribution after step ", step)
      message("  covariance matrix V:")
      message(paste0("    ", capture.output(print(signif(Vproposal, 3))), collapse="\n"))
      if (!is.null(acceptance.target))
        message("  radius R: ", signif(Rproposal, 3))
    }
    
    # Store samples
    if (keep.history)
      history[[step+1]] <- states[accept, c(list(step=..step), .SD) ]
  }
  
  # Return result
  result.data <- list(
    variables=variables,
    llfun=llfun,
    arguments=list(
      steps=steps, accepts=accepts, chains=chains, candidate.samples=chains,
      initial.samplewithreplacement=initial.samplewithreplacement,
      reevaluate.likelihood=reevaluate.likelihood, acceptance.target=acceptance.target,
      initial.Rproposal=initial.Rproposal, minimal.Rproposal=minimal.Rproposal,
      initial.Vproposal=initial.Vproposal, minimal.Vproposal.CV=minimal.Vproposal.CV, maximal.Vproposal.corr=maximal.Vproposal.corr,
      proposal.update.rate=proposal.update.rate, delta.ll.cutoff=delta.ll.cutoff, initial.ll.cutoff=initial.ll.cutoff))

  if (keep.candidates)
    result.data[["candidates"]] <- candidates
  
  if (keep.initial)
    result.data[["initial"]] <- cbind(initial, initial.meta)
  
  if (keep.history)
    result.data[["history"]] <- rbindlist(history)
  
  result.data[["final"]] <- cbind(states, states.meta)

  return(structure(result.data, class="SANMCMC"))
}

#' Sample from a multivariate uniform distribution
#' 
#' @export
sample.variables <- function(n, variables) {
  vs <- names(variables)
  names(vs) <- names(variables)
  return(as.data.table(lapply(vs, function(v) {
    runif(n, min=variables[[v]][1], max=variables[[v]][2])
  })))
}

#' Sample from the posterior distribution
#'
#' @export
sample.posterior <- function(n, mcmc, H="auto") {
  # Estimate bandwidth matrix
  H <- bandwidth.matrix(mcmc, H)
  # Sample randomly from the posterior samples
  k <- nrow(mcmc$final)
  s <- mcmc$final[sample.int(k, size=n, replace=TRUE), names(mcmc$variables), with=FALSE]
  # Add displacements samples from a multivariate normal with covariance matrix H
  # We retry until the displaced value lies within the variable's domain.
  vs <- names(mcmc$variables)
  names(vs) <- names(mcmc$variables)
  r <- list()
  while(nrow(s) > 0) {
    # Sample displacements
    m <- nrow(s)
    dp <- mvtnorm::rmvnorm(n=m, sigma=H)
    colnames(dp) <- names(mcmc$variables)
    # Displace samples, check if result lies within domain otherwise set to NA
    sp <- as.data.frame(lapply(vs, function(v) {
      # Fetch variable's domain [lb, ub]
      lb <- mcmc$variables[[v]][1]
      ub <- mcmc$variables[[v]][2]
      x <- s[[v]] + dp[, v]
      ifelse((lb <= x) & (x <= ub), x, NA_real_)
    }), optional=TRUE)
    # Append accepted samples to result r and remove from s
    a <- !Reduce(`|`, Map(is.na, sp), rep(FALSE, m))
    r <- c(r, list(sp[a,]))
    s <- s[!a]
  }
  
  # Merge results from all lopp iterations and return
  return(rbindlist(r))
}

#' Find a suitable bandwidth matrix for kernel density estimation 
#' 
#' @export
bandwidth.matrix.default <- function(x, H="auto") {
  x <- as.matrix(x)
  d <- ncol(x)
  H <- if (is.numeric(H)) {
    # Bandwidth H was specified numerically
    if (is.vector(H)) {
      if (length(H) %in% c(1, d)) diag(H, nrow=d, ncol=d)
      else stop("If H is a numeric vector, it must either specify a single ",
                "bandwith for all variables, or one bandwith per variable.")
    } else if (is.matrix(H)) {
      if ((nrow(H) == d) && (ncol(H) == d)) H
      else stop("If H is a numeric matrix, if must have as many rows ",
                "and columns as the posterior has variables")
    } else stop("If H is numeric, it must either be a vector or a matrix")
  } else {
    if (as.character(H) == "auto")
      H <- if (d <= 6) "Hpi" else "silverman.cov"
    # Bandwith H must be estimated
    if (as.character(H) == "silverman.cov") {
      # Determine bandwith matrix H by scaling the empirical covariance matrix
      # by a factor R. The choice R=n^(-1/5) is loosely based on the rule-of-thumb
      # univariate bandwidth estimate H = 0.9 * min(sigma, IQR/1.34) * n^(-1/5).
      # (https://en.wikipedia.org/wiki/Kernel_density_estimation)
      # Note: A more precise estimate for the bandwith matrix H would be desirable,
      # but typical algorithms like the ones implemented in the package ks seems to
      # fail already for a moderate number of dimensions (say, 6).
      V <- cov(x)
      R <- nrow(x)^(-1/5)
      V * (R^2)
    } else if (as.character(H) == "Hpi") {
      ks::Hpi(x, deriv.order=1, nstage=2-(d > 2))
    } else if (as.character(H) == "Hpi.diag") {
      ks::Hpi.diag(x, deriv.order=1, nstage=2-(d > 2))
    } else stop("unknown bandwith estimation method ", as.character(H))
  }
  colnames(H) <- colnames(x)
  rownames(H) <- colnames(x)
  return(H)
}

#' Find a suitable bandwith matrix for kernel density estimation and sampling of the posterior
#' 
#' @export
bandwidth.matrix.SANMCMC <- function(sanmcmc, H="auto") {
  x <- as.matrix(sanmcmc$final[, names(sanmcmc$variables), with=FALSE])
  return(bandwidth.matrix(x, H))
}

#' Locates the posterior distribution's mode(s) using the mean-shift algorithm
#'
#' @export
locate.modes.default <- function(x, tolerance=0.1, adjust=1.0, H="Hpi") {
  x <- as.matrix(x)
  H <- bandwidth.matrix(x, H)
  ks::kms(x, H=H*(adjust^2), min.clust.size=0.1*nrow(x), tol.clust=tolerance)
}

#' Locates the posterior distribution's mode(s) using the mean-shift algorithm
#'
#' @export
locate.modes.SANMCMC <- function(sanmcmc, tolerance=0.1, adjust=1.0, H="Hpi") {
  locate.modes(as.matrix(sanmcmc$final[, names(sanmcmc$variables), with=FALSE]),
               tolerance, adjust=adjust, H=H)
}

#' Computes various summary statistics of the posterior distribution
#' 
#' @export
summarystats.SANMCMC <- function(sanmcmc, modes, expressions=names(sanmcmc$variables)) {
  stats <- c("mean", "std. dev.", "median", "mad*1.48", "mode", "ML")
  
  # Run mean-shift algorithm if necessary and find mode
  if (missing(modes))
    modes <- locate.modes(sanmcmc)
  i.mode <- which.max(modes$nclust.table)
  mode <- modes$mode[i.mode,]
  
  sanmcmc$final[, {
    stats <- lapply(expressions, function(expr) {
      # Translate strings into expressions. For added convenience, column
      # names are quoted automatically, this makes e.g. writing "40S-40A" work.
      if (is.character(expr)) {
        expr <- gsub(paste0("(", paste0("(?:", colnames(.SD), ")", collapse="|"), ")"), "`\\1`", expr)
        expr <- parse(text=expr)
      }
      if (!is.expression(expr) && !is.name(expr))
        stop("invalid expression of type ", class(expr))
      
      # Evaluate expression
      v <- eval(expr)
      
      # For MCMC variables, the mode is the multi-variate mode as determined
      # by the mean-shift algorithm. For other expressions, the mode is determined
      # using univariate KDE.
      m <- eval(expr, envir=as.list(mode))
      if (is.null(m)) {
        d <- density(v)
        m <- d$x[which.max(d$y)]
      }
      
      # Compute summary statistics 
      c(`mean`=mean(v),
        `std. dev.`=sd(v),
        `median`=median(v),
        `mad*1.48`=mad(v),
        `mode`=m,
        `ML`=v[which.max(ll)])
    })
    
    # Prepend "stat" column which indicates which row represents
    # which statistic and name the columns
    names(stats) <- as.character(expressions)
    c(list(statistic=names(stats[[1]])), stats)
  }]
}