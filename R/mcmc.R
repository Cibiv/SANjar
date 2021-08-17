library(data.table)
library(mvtnorm)

#' Samples from a distribution using the Metropolis-Hastings MCMC algorithm  
#' 
#' @export
mcmc <- function(llfun, variables, steps=NA, accepts=NA, chains,
                 initial.samples=chains, initial.withreplacement=FALSE, initial.states=NULL,
                 keep.history=FALSE, keep.initial=FALSE, keep.uniform=FALSE,
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
  
  parameter.samples <- function(n) {
    as.data.table(c(list(ll=NA_real_), lapply(varnames, function(v) {
      runif(n, min=variables[[v]][1], max=variables[[v]][2])
    })))
  }
  
  # Initial state distribution
  if (is.null(initial.states)) {
    # Initial approximation (ABC-like).
    # To improve performance, we try to determine a reasonable initial.ll.cutoff value
    if (is.na(initial.ll.cutoff)) {
      preinitial.samples <- ceiling(sqrt(initial.samples))
      if (verbose)
        message("Generating ", preinitial.samples, " parameter samples to estimate initial.ll.cutoff")
      preinitial <- parameter.samples(preinitial.samples)
      if (verbose)
        message("Evaluating posterior likelihood for ", preinitial.samples, " parameter samples to estimate initial.ll.cutoff")
      r <- as.data.table(llfun(preinitial[, varnames, with=FALSE], cutoffs=-Inf))
      initial.ll.cutoff <- max(r[, 1])
    }
    # Sample parameter space uniformly, evaluate posterior likelihood, and re-sample according to this likelihood.
    if (verbose)
      message("Generating ", initial.samples, " initial parameter samples")
    uniform <- parameter.samples(initial.samples)
    if (verbose)
      message("Evaluating posterior likelihood for ", initial.samples, " initial parameter samples ",
              "with initial.ll.cutoff ", signif(initial.ll.cutoff, 3))
    r <- as.data.table(llfun(uniform[, varnames, with=FALSE], cutoffs=initial.ll.cutoff))
    if (ncol(r) > 1)
      meta <- r[, 2:ncol(r)]
    else
      meta <- NULL
    uniform[, ll := r[, 1]]
    # Remove sample with likelihood -infinity
    ll.finite <- is.finite(uniform$ll)
    uniform <- uniform[ll.finite]
    if (!is.null(meta))
      meta <- meta[ll.finite]

    if (nrow(uniform) == 0)
      stop("No parameter sample had finite likelihood, aborting")
    if (verbose) {
      message("The following ", nrow(uniform), " samples are usable:")
      message(paste0("  ", capture.output(print(signif(uniform, 3))), collapse="\n"))
    }

    # Sample initial states accordingly to likelihood from the uniformly spaced samples
    if (verbose)
      message("Selecting initial states for ", chains, " chains from ", nrow(uniform), " initial parameter samples")
    if (initial.withreplacement) {
      # Sample the initial states
      i <- uniform[, sample.int(.N, size=chains, prob=exp(ll), replace=TRUE) ]
    } else {
      # Repeat existing sample often enough to not run out of samples when sampling the initial states
      k <- ceiling(chains / sum(exp(uniform$ll) > 0))
      uniform <- rbindlist(rep(list(uniform), k))
      if (!is.null(meta))
        meta <- rbindlist(rep(list(meta), k))
      # Sample the initial states *without* replacement
      # The idea is that this ensures a reasonable dispersion of initial values
      i <- uniform[, sample.int(.N, size=chains, prob=exp(ll), replace=FALSE) ]
    } 
    initial.states <- uniform[i, c(list(chain=1:.N, naccepts=0), .SD[, c("ll", varnames), with=FALSE]) ]
    if (!is.null(meta))
      initial.meta <- meta[i]
    else
      initial.meta <- NULL
  } else {
    uniform <- NULL
    if (!missing(chains) && (chains != nrow(initial.states)))
      warning("Number of initial states (", nrow(initial.states), ") overrides specified number of chains (", chains, ")")
    chains <- nrow(initial.states)
    if (ncol(initial.states) > (length(varnames) + 1))
      initial.meta <- initial.states[, .SD, .SDcols=!c("ll", varnames)]
    else
      initial.meta <- NULL
    initial.states <- initial.states[, c(list(chain=1:.N, naccepts=0), .SD[, c("ll", varnames), with=FALSE]) ]
  }
  if (verbose) {
    message("Initial states")
    initial.states.signif <- initial.states[, lapply(.SD, function(c) { if (is.numeric(c)) signif(c, 3) else c } )]
    message(paste0("  ", capture.output(print(initial.states.signif)), collapse="\n"))
  }
  
  # Initialize proposal distribution
  # If we started with ABC-style uniform sampling of the parameter space, we compute
  # the initial covariance estimate using all samples and their likelihoods, since
  # this gives a more precise estimate than using just the initial states.
  Vproposal <- if (is.null(initial.Vproposal)) {
    if (is.null(initial.states))
      cov.wt(uniform[, varnames, with=FALSE], wt=exp(uniform$ll))$cov
    else
      cov(initial.states[, varnames, with=FALSE])
  } else
    initial.Vproposal
  Rproposal <- initial.Rproposal
  if (verbose) {
    message("Initial proposal distribution")
    message("  covariance matrix V:")
    message(paste0("    ", capture.output(print(signif(Vproposal, 3))), collapse="\n"))
    if (!is.null(acceptance.target))
      message("  radius R: ", signif(Rproposal, 3))
  }
  
  # Setup MCMC variables
  states <- initial.states
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
      # And force the radius to one, otherwise the actual proposal radius may still be tiny
      Rproposal <<- 1
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
    displacements <- rmvnorm(n=sum(extend), sigma=Vproposal) * Rproposal
    colnames(displacements) <- names(varnames)
    proposals <- data.table(valid=extend)
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
        r <- as.data.table(llfun(states[proposals$valid, varnames, with=FALSE], cutoffs=-Inf))
        stopifnot(is.finite(unlist(r[, 1])))
        if (!is.null(states.meta) && (ncol(r) > 1))
          states.meta[proposals$valid, colnames(r)[2:ncol(r)] := r[, 2:ncol(r)] ]
        else
          states.meta <- NULL
        states[proposals$valid, ll := r[, 1] ]
      }
      
      # Evaluate likelihood of valid proposals
      r <- as.data.table(llfun(proposals[valid==TRUE, varnames, with=FALSE], cutoffs=states[proposals$valid, ll] - delta.ll.cutoff))
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
    arguments=list(
      steps=steps, accepts=accepts, chains=chains,
      initial.withreplacement=initial.withreplacement, reevaluate.likelihood=reevaluate.likelihood,
      acceptance.target=acceptance.target,
      initial.Rproposal=initial.Rproposal, minimal.Rproposal=minimal.Rproposal,
      initial.Vproposal=initial.Vproposal, minimal.Vproposal.CV=minimal.Vproposal.CV, maximal.Vproposal.corr=maximal.Vproposal.corr,
      proposal.update.rate=proposal.update.rate, delta.ll.cutoff=delta.ll.cutoff, initial.ll.cutoff=initial.ll.cutoff))

  if (keep.uniform)
    result.data[["uniform"]] <- uniform
  
  if (keep.initial)
    result.data[["initial"]] <- cbind(initial.states, initial.meta)
  
  if (keep.history)
    result.data[["history"]] <- rbindlist(history)
  
  result.data[["final"]] <- cbind(states, states.meta)

  return(structure(result.data, class="SANMCMC"))
}
