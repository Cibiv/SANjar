#' The MCMC results for the LT47 dataset
"lt47.mcmc"

#' Find a suitable bandwidth matrix for kernel density estimation 
#' 
#' @export
bandwidth.matrix <- function(x, ...) UseMethod("bandwidth.matrix")

#' Computes the MAP (maximum a-posteriori) estimate from MCMC results
#' 
#' @export
map.estimate <- function(x, ...) UseMethod("map.estimate")

#' Computes various summary statistics
#' 
#' @export
summarystats <- function(x, ...) UseMethod("summarystats")

#' Samples from a distribution using the Metropolis-Hastings MCMC algorithm  
#' 
#' @export
mcmc <- function(llfun, variables, fixed=character(), llfun.average=1,
                 steps=NA, accepts=NA, chains, candidate.samples=chains,
                 candidates=NULL, initial.samplewithreplacement=FALSE, initial=NULL, 
                 keep.history=FALSE, keep.initial=FALSE, keep.candidates=FALSE,
                 acceptance.target=if ((length(variables) == 1) || directional.metropolis.gibbs) 0.44 else 0.234,
                 directional.metropolis.gibbs=FALSE, initial.Rproposal=1.0, minimal.Rproposal=0.2,
                 initial.Vproposal=NULL, minimal.Vproposal.CV=0.1, maximal.Vproposal.corr=0.95,
                 proposal.update.rate=0.2, delta.ll.cutoff=15, initial.ll.cutoff=NA,
                 verbose=FALSE, extremely.verbose=FALSE)
{
  # Negative values would reject any proposal that doesn't *decrease* the likelihood sufficiently
  stopifnot(delta.ll.cutoff >= 0)
  
  # Either steps or accepts must be specified
  stopifnot(!is.na(steps) || !is.na(accepts))
  
  # Prepare a list of variable names, ensuring that lapply(varnames) will return a *named* list
  varnames <- names(variables)
  names(varnames) <- varnames
  
  # Average likelihoods if requested
  average.likelihoods <- function(parameters, cutoffs) {
    if (llfun.average > 1) {
      # Evaluate likelihoods multiple times for each proposal
      i <- rep(1:nrow(parameters), each=llfun.average)
      parameters <- parameters[i,]
      if (length(cutoffs) > 1)
        cutoffs = rep(cutoffs, each=llfun.average)[i]
      rbindlist(lapply(split(as.data.table(llfun(parameters, cutoffs=cutoffs)), i), function(g) {
        # Average likelihoods (NOT log-likelihoods!) across proposals
        # Since we can't easily average the meta-information, we take the one
        # that corresponds to the highest likelihood.
        ll <- g[, 1][[1]]
        j.max <- which.max(ll)
        ll.max <- ll[j.max]
        ll.avg <- log(mean(exp(ll - ll.max))) + ll.max
        cbind(ll.avg, g[j.max, -1])
      }))
    } else {
      as.data.table(llfun(parameters, cutoffs=cutoffs))
    }
  }

  # Proposal generator for multivariate normal proposals
  fill.proposals.mvnorm <- function(proposals, states) {
    # Generate proposals by sampling displacements using current covariances and radius,
    # and displacing the previous samples accordingly. We only really generate proposals
    # for chains that need extension.
    displacements <- mvtnorm::rmvnorm(n=sum(proposals$valid), sigma=Vproposal) * Rproposal
    colnames(displacements) <- names(varnames)
    proposals[valid==TRUE, (varnames) := lapply(varnames, function(v) {
      states[proposals$valid, eval(as.name(v))] + displacements[, v]
    })]
  }
  
  # Proposal generator for Metropolis-within-Gibbs sampling
  fill.proposals.gibbs <- function(proposals, direction.index, states) {
    if (direction.index == 1) {
      # Re-compute the factorization of the covariance matrix before starting a new iteration
      ev <- eigen(Vproposal, symmetric=TRUE)
      A <- t(ev$vectors) * sqrt(pmax(ev$values, 0))
      colnames(A) <- varnames
      Vproposal.directions <<- asplit(A, MARGIN=1)
      message("Updated Gibbs sampling directions:")
      message(paste0("  ", capture.output(print(signif(do.call(rbind, Vproposal.directions), 3))),
                     collapse="\n"))
    }
    
    # Compute displacements along the direction.index-ith direction
    displacements <- (matrix(rnorm(sum(proposals$valid)), ncol=1)
                      %*% (Vproposal.directions[[direction.index]]
                           * Rproposal[direction.index]))
    proposals[valid==TRUE, (varnames) := lapply(varnames, function(v) {
      states[proposals$valid, eval(as.name(v))] + displacements[, v]
    })]
  }

  # Split table returned by the llfun into the likelihood column and the meta-data table
  split.llfun.out <- function(llfun.out) {
    ll <- llfun.out[[ colnames(llfun.out)[1] ]]
    if (ncol(llfun.out) > 1)
      meta <- llfun.out[, 2:ncol(llfun.out)]
    else
      meta <- NULL
    return(list(ll=ll, meta=meta))
  }
  
  # Likelihood evaluation
  evaluate.loglikelihoods <- function(proposals, states) {
    split.llfun.out(average.likelihoods(parameters=proposals[valid==TRUE, c(fixed, varnames), with=FALSE],
                                        cutoffs=states[proposals$valid, ll] - delta.ll.cutoff))
  }
  
  # Perform MH step
  metropolis.hastings <- function(proposals, results, states, states.meta, varnames) {
    # Check that we got the correct number of results
    stopifnot(nrow(proposals) == nrow(states))
    stopifnot(is.null(states.meta) || (nrow(states) == nrow(states.meta)))
    stopifnot(nrow(results) == sum(proposals$valid))
    
    # Set proposal likelihoods
    proposals[, ll := -Inf]
    proposals[valid==TRUE, ll := results$ll ]
    
    # Accept or reject proposals using the Metropolis-Hastings rule
    accept <- rep(FALSE, nrow(states))
    mh <- (proposals$valid & is.finite(proposals$ll))
    accept[mh] <- (runif(sum(mh)) <= exp(proposals$ll[mh] - states$ll[mh]))
    # Update parameters if a proposal was accepted
    states[accept, c("ll", varnames) := proposals[accept, c("ll", varnames), with=FALSE] ]
    states[accept, naccepts := naccepts + 1]
    # Replace meta information if a proposal was accepted
    if (!is.null(states.meta) && !is.null(results$meta)) {
      stopifnot(!accept[!proposals$valid])
      states.meta[accept, colnames(results$meta) := results$meta[accept[proposals$valid]] ]
    }
    else
      states.meta <- NULL
    
    return(accept)
  }
  
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
    # Check that we have enough candidates
    if (nrow(candidates) < chains)
      stop("need at least as many candidate samples as there are chains, ideally much more")
    # Run MCMC with independent proposals to initialize the chains
    initial <- data.table(chain=1:chains, naccepts=0, ll=-Inf)
    initial[, c(fixed, varnames) := rep(list(NA_real_), length(fixed) + length(varnames)) ]
    initial.meta <- NULL
    i <- 1
    while (i + chains - 1 <= nrow(candidates)) {
      # Fetch next <chains> proposals from the candidates table
      proposals <- data.table(valid=rep(TRUE, chains), ll=-Inf)
      proposals[, c(fixed, varnames) := candidates[i:(i+chains-1), c(fixed, varnames), with=FALSE] ]
      # Evaluate proposal likelihoods
      results <- evaluate.loglikelihoods(proposals, initial)
      if ((i == 1) && !is.null(results$meta))
        initial.meta <- results$meta
      # Perform MH step
      accept <- metropolis.hastings(proposals, results, initial, initial.meta, c(fixed, varnames))
      # Log
      if (extremely.verbose) {
        message("States after initialization step ", 1 + (i-1)/chains)
        initial.signif <- initial[, lapply(.SD, function(c) { if (is.numeric(c)) signif(c, 3) else c } )]
        message(paste0("  ", capture.output(print(initial.signif)), collapse="\n"))
      } else
        message("Completed initialization step ", 1 + (i-1)/chains)
      message("Valid proposal ratios: ", signif(proposals[, mean(valid)], 3))
      message("Finite-LL proposal ratios: ", signif(proposals[, mean(is.finite(ll))], 3))
      message("Acceptance ratio: ", signif(mean(accept), 3))
      message("Accepted proposals in each chain: min=", min(initial$naccepts), " median=", median(initial$naccepts), " max=", max(initial$naccepts))
      message("Average log-likelihood over all chains: ", signif(initial[, mean(ll)], 3))
      # Continue until condidates are exhausted
      i <- i + chains
    }
    if (all(is.infinite(initial$ll)))
      stop("Unable to find initial parameters with finite likelihood, aborting")
    if (any(is.infinite(initial$ll))) {
      message("Some chains have likelihood -inf after initialization, re-sampling their parameters from other chains")
      # Select replacement parameters from the initial parameters with finite likelihood
      ll.inf <- is.infinite(initial$ll)
      ll.fin.idx <- (1:nrow(initial))[!ll.inf]
      ll.inf.repidx <- ll.fin.idx[sample.int(length(ll.fin.idx), size=sum(ll.inf), replace=TRUE)]
      # Replace parameters
      initial[ll.inf, c("ll", varnames) := initial[ll.inf.repidx, c("ll", varnames), with=FALSE] ]
      # Replace meta information
      if (!is.null(initial.meta))
        initial.meta[ll.inf, colnames(initial.meta) := initial.meta[ll.inf.repidx] ]
    }
    # Reset number of accepted proposals before main MCMC part begins
    reset.naccepts <- TRUE
  } else {
    # Initial states were specified
    if (!is.null(candidates))
      warning("Initial states were specified, provided candidates not being used")
    if (!missing(chains) && (chains != nrow(initial)))
      warning("Number of initial states (", nrow(initial), ") overrides specified number of chains (", chains, ")")
    chains <- nrow(initial)
    # Split off meta-data columns, which are the columns after the last variable (or fixed parameter)
    ci.lastvar <- max((1:ncol(initial))[colnames(initial) %in% c(fixed, varnames)])
    if (ncol(initial) > ci.lastvar)
      initial.meta <- initial[, -(1:ci.lastvar)]
    else
      initial.meta <- NULL
    # Add chain and naccepts columns to initial states
    if (!("chain" %in% colnames(initial)))
      initial[, chain := 1:.N]
    if (!("naccepts" %in% colnames(initial)))
      initial[, naccepts := 0]
    # Reorder columns to match the expected order
    initial <- initial[, c("chain", "naccepts", "ll", fixed, varnames), with=FALSE]
    # Don't reset number of accepted proposals before main MCMC part begins
    reset.naccepts <- FALSE
  }
  if (llfun.average > 1) {
    message("Will average the likelihood over ", llfun.average, " evaluations, re-evaluating initial likelihoods now")
    # Evaluate likelihood of valid proposals
    r <- average.likelihoods(initial[, c(fixed, varnames), with=FALSE], cutoffs=-Inf)
    initial[, ll := r[, 1]]
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
  Vproposal <- if (is.null(initial.Vproposal))
    cov(initial[, varnames, with=FALSE])
  else
    initial.Vproposal
  # For Metropolis-within-Gibbs, we scale each component separately, but
  # for multivariate normal proposals, there is only a single scaling factor
  Rproposal <- if (directional.metropolis.gibbs)
    sapply(varnames, function(v) initial.Rproposal)
  else
    initial.Rproposal
  if (verbose) {
    message("Initial proposal distribution")
    message("  covariance matrix V:")
    message(paste0("    ", capture.output(print(signif(Vproposal, 3))), collapse="\n"))
    if (!is.null(acceptance.target))
      message("  radius R: ", paste(signif(Rproposal, 3), collapse=" "))
  }
  Vproposal.directions <- NULL
  
  # Setup MCMC variables
  states <- copy(initial)
  if (reset.naccepts)
    states[, naccepts := 0 ]
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
    
    # If we're doing metropolis-within-gibbs, we iterate through the possible directions individually
    direction <- if (directional.metropolis.gibbs)
      direction <- 1 + (step - 1) %% length(variables)
    else 1

    # Copy states, since data.table does not copy on update
    states <- copy(states)
    
    # Determine which chains need extension
    extend <- if (!is.na(accepts))
      states$naccepts < accepts
    else
      rep(TRUE, nrow(states))

    # Setup proposals table, and fill in fixed variables
    proposals <- data.table(valid=extend, ll=-Inf)
    if (length(fixed) > 0)
      proposals[valid==TRUE, (fixed) := states[extend, fixed, with=FALSE] ]
    
    # Generate proposals
    if (directional.metropolis.gibbs)
      fill.proposals.gibbs(proposals, direction, states)
    else
      fill.proposals.mvnorm(proposals, states)
    
    # Accept or reject proposals
    # Proposals outside the parameter's domain are skipped.
    accept <- rep(FALSE, nrow(states))
    for(v in varnames)
      proposals[valid==TRUE, valid := valid & (variables[[v]][1] <= eval(as.name(v))) & (eval(as.name(v)) <= variables[[v]][2]) ]
    if (any(proposals$valid)) {
      # Evaluate likelihood of valid proposals
      results <- evaluate.loglikelihoods(proposals, states)
      # Perform MH updates
      accept <- metropolis.hastings(proposals, results, states, states.meta, varnames)
    }

    if (verbose) {
      if (extremely.verbose) {
        message("States after step ", step, if (directional.metropolis.gibbs) paste0(" (direction ", direction, ")") else "")
        states.signif <- states[, lapply(.SD, function(c) { if (is.numeric(c)) signif(c, 3) else c } )]
        message(paste0("  ", capture.output(print(states.signif)), collapse="\n"))
      } else
        message("Completed step ", step, if (directional.metropolis.gibbs) paste0(" (direction ", direction, ")") else "")
      message("Chain extension attempts: ", sum(extend))
      message("Valid proposal ratio over extension attempts: ", signif(proposals[extend, mean(valid)], 3))
      message("Finite-LL proposal ratio over extension attempts: ", signif(proposals[extend, mean(is.finite(ll))], 3))
      message("Acceptance ratio over extension attempts: ", signif(mean(accept[extend]), 3))
      message("Accepted proposals in each chain: min=", min(states$naccepts), " median=", median(states$naccepts), " max=", max(states$naccepts))
      message("Average log-likelihood over all chains: ", signif(states[, mean(ll)], 3))
    }

    # Update the proposal distribution, smoothly over about proposal.update.rate steps.
    Vproposal.update <- ((proposal.update.rate > 0) &&
                         (!directional.metropolis.gibbs || (direction == length(variables))))
    if (Vproposal.update) {
      Vproposal <- Vproposal * (1 - proposal.update.rate) + cov(states[, varnames, with=FALSE]) * proposal.update.rate
      Vproposal.clamp()
    }
    
    # Update the proposal radius to meet the acceptance target, smoothly over about proposal.update.rate steps.
    # Widening the proposal radius reduces the acceptance rate, we thus adjust the radius based on the ratio
    # of actual vs. targeted acceptance rate. To avoid the radius becoming zero (where it would then remain)
    # updates are limited to at most a factor of two, regardless of the update rate. To avoid large fluctuations
    # the update factor is dampened according to proposal.update.rate, and according to how many chains we
    # attempted to extend.
    if (!is.null(acceptance.target) && (proposal.update.rate > 0)) {
      f.raw <- mean(accept[extend]) / acceptance.target
      gamma <- 0.1
      r <- proposal.update.rate * gamma * sum(extend) / (1 + gamma * sum(extend))
      f <- max(0.5, min(f.raw ^ r, 2))
      Rproposal[direction] <- max(minimal.Rproposal, Rproposal[direction] * f)
    }
    
    if (verbose && Vproposal.update) {
      message("Updated proposal distribution after step ", step)
      message("  covariance matrix V:")
      message(paste0("    ", capture.output(print(signif(Vproposal, 3))), collapse="\n"))
      if (!is.null(acceptance.target))
        message("  radius R: ", paste(signif(Rproposal, 3), collapse=" "))
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
      llfun.average=llfun.average, steps=steps, accepts=accepts, chains=chains, candidate.samples=chains,
      initial.samplewithreplacement=initial.samplewithreplacement,
      acceptance.target=acceptance.target, directional.metropolis.gibbs=directional.metropolis.gibbs,
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

#' Computes the MAP (maximum a-posteriori) estimate from MCMC results
#' 
#' @export
map.estimate.SANMCMC <- function(sanmcmc, ms.tolerance=0.1, H.adjust=1.0, H="Hpi", method="max.local.lh.avg", ms=NULL) {
  # Get posterior samples as a matrix with one column per variable, one row per sample 
  final.ll <- sanmcmc$final[, ll]
  final <- as.matrix(sanmcmc$final[, names(sanmcmc$variables), with=FALSE])
  # Run mean-shift algorithm unless mean-shift results are provided (see meanshift()).
  # This algorithm finds a set of modes (i.e. local density maxima), and assigns each
  # sample to one of these modes. The mean-shift algorithm thus is not only a mode-finding
  # algorithm, but also an (unsupervised) clustering algorithm.
  if (is.null(ms)) {
    x <- as.matrix(sanmcmc$final[, names(sanmcmc$variables), with=FALSE])
    H <- bandwidth.matrix(x, H)
    ms <- ks::kms(x, H=H*(H.adjust^2), min.clust.size=0.1*nrow(x), tol.clust=ms.tolerance)
  }
  # Pick the "maximal" mode
  # Which mode constitutes the "maximum" depends on the metric used to compare the modes
  method <- match.arg(method, c("max.local.lh.avg", "max.cluster.lh.sum", "largest.cluster"))
  i.mode <- if (method == "max.local.lh.avg") {
    # Metric is the average likelihoods of the MCMC samples in the neighbourhood of each node.
    # The neighborhoods are defined by the bandwidth matrix;  each sample's LH is weighted
    # by its contribution to the density at the mode. The idea behind this approach is that
    # it should be robust against partially converged MCMC runs. If some MCMC chains got stuck
    # in local optima and hence have low likelihoods, this approach should not pick these
    # clusters.
    H <- bandwidth.matrix(sanmcmc, H=H)
    ll.scale <- max(final.ll)
    modes.ll <- sapply(1:ms$nclust, function(i) {
      w <- mvtnorm::dmvnorm(final, mean=ms$mode[i,], sigma=H*(H.adjust^2))
      log(weighted.mean(exp(final.ll - ll.scale), w=w)) + ll.scale
    })
    which.max(modes.ll)
  } else if (method == "max.cluster.lh.sum") {
    # Metric is the sum of likelihoods within each cluster. This is a version of "largest.cluster"
    # this is (somehwat) robust against partially converged MCMC runs because the influence of
    # low-likelihood samples is reduced. However, larger sub-optimal clusters may still dominate
    # over smaller clusters containing better likelihoods.
    ll.scale <- max(final.ll)
    modes.ll <- sapply(1:ms$nclust, function(i) {
      log(sum(exp(final.ll[ms$label==i] - ll.scale))) + ll.scale
    })
    print(modes.ll)
    which.max(modes.ll)
  } else if (method == "largest.cluster") {
    # Metric is simply the cluster size. This is robust against noisy likelihood evaluations,
    # but not robust against partially converged MCMC results because low-likelihood samples can
    # dominate the mode selection
    which.max(ms$nclust.table)
  }
  return(ms$mode[i.mode,])
}

sanmcmc.evaluate <- function(sanmcmc, which="final", expressions=names(sanmcmc$variables), extra=character()) {
  varnames <- names(sanmcmc$variables)
  data <- sanmcmc[[which]]
  varidx <- which(colnames(data) %in% varnames)
  ididx <- 1:(min(varidx)-1)
  
  labels <- list()
  r <- do.call(cbind, c(list(data[, ididx, with=FALSE]),
                   lapply(c(as.list(expressions), as.list(extra)), function(expr) {
    # Translate strings into expressions. For added convenience, column
    # names are quoted automatically, this makes e.g. writing "40S-40A" work.
    if (is.character(expr)) {
      label <- expr
      expr <- gsub(paste0("(", paste0("(?:", varnames, ")", collapse="|"), ")"), "`\\1`", expr)
      expr <- parse(text=expr)
    } else {
      label <- as.character(expr)
    }
    if (!is.expression(expr) && !is.name(expr))
      stop("invalid expression of type ", class(expr))
    
    labels <<- c(labels, list(label))
                     
    # Evaluate expression and return single-column data.table
    r <- data[, list(eval(expr)) ]
    colnames(r) <- label
    r
  })))
  attr(r, "expressions") <- as.character(labels)
  r
}

#' Computes various summary statistics of the posterior distribution
#' 
#' @export
summarystats.SANMCMC <- function(sanmcmc, map, expressions=names(sanmcmc$variables)) {
  stats <- c("mean", "std. dev.", "median", "mad*1.48", "MAP", "ML")
  
  # Run mean-shift algorithm if necessary and find mode
  if (missing(map))
    map <- map.estimate(sanmcmc)

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
      
      # For MCMC variables, the MAP (maximum á posteriori) is the maximal multi-variate mode
      # For other expressions, the mode is determined using univariate KDE.
      m <- if (is.name(expr[[1]]) && (length(expr) == 1)) {
        # Expression is a single variable name, use multi-variate MAP estimates computed above
        eval(expr, envir=as.list(map))
      } else {
        # Expression is a complex expression, compute univariate mode
        d <- density(v)
        d$x[which.max(d$y)]
      }

      # Compute summary statistics 
      c(`mean`=mean(v),
        `std. dev.`=sd(v),
        `median`=median(v),
        `mad*1.48`=mad(v),
        `MAP`=m,
        `ML`=v[which.max(ll)])
    })
    
    # Prepend "stat" column which indicates which row represents
    # which statistic and name the columns
    names(stats) <- as.character(expressions)
    c(list(statistic=names(stats[[1]])), stats)
  }]
}