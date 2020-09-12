library(inline)
library(lamW)
library(data.table)

#' Internal pure S-model simulation function
#'
#' Simulates the pure S-model starting from the initial counts in \code{S}
#' for the number(s) of steps specified in \code{steps}.
discrete_s_c <- rcpp(
  signature(C_init="list", p_S="numeric", p_R="numeric", steps="integer"),
  includes=c("#include <R_ext/Utils.h>", "#include <Rmath.h>"),
  Rcpp=TRUE,
  cxxargs=c("-O3", "-fast-math"),
  body='
using namespace Rcpp;
  
/* Extract & validate parameters */
List C_init_(C_init);
IntegerVector C_curr_S(clone(as<IntegerVector>(C_init_["S"])));
const int L_ = IntegerVector(C_curr_S).size();
if ((DoubleVector(p_S).size() != 1) || (DoubleVector(p_R).size() != 1))
  stop("p_S, p_R must have length 1");
const double p_S_ = DoubleVector(p_S)[0];
const double p_R_ = DoubleVector(p_R)[0];
const IntegerVector stepsv = steps;

/* Prepare result */
List result;

/* Run simulation over all requested time intervals */
RNGScope rngScope;
for(int i=0; i < stepsv.size(); ++i) {
  int* const C_curr_S_ = INTEGER(C_curr_S);

  /* Proceed until the end of the current interval */
  for(int l=0, l_end=stepsv[i]; l < l_end; ++l) {
    /* Process lineages one by one */
    for(int j=0; j < L_; ++j) {
      /* Prevent overflows */
      if ((C_curr_S_[j] >= 1000000000))
        stop("S cell count reached one billion, stopping to avoid integer overflow");
        
      /* Determine the number of cells affected by each type of event */
      const int dS = (p_S_ > 0) ? R::rbinom(C_curr_S_[j], p_S_) : 0;
      const int dR = (p_R_ > 0) ? R::rbinom(C_curr_S_[j], p_R_) : 0;

      /* Update cell counts */
      C_curr_S_[j] = std::max(C_curr_S_[j] + dS - dR, 0);
    }
  }
  
  /* Append current per-lineage cell counts to the result list */
  List C_curr = List::create(
    Named("i") = IntegerVector(1, i+1),
    Named("S") = clone(C_curr_S)
  );
  result.push_back(C_curr);
  
  /* Allow interruptions */
  Rcpp::checkUserInterrupt();
}

/* Return result */
return result;
')

#' Internal SAN model simulation function
#'
#' Simulates the SAN model starting from the initial counts in \code{C_init}
#' for the number(s) of steps specified in \code{steps}.
discrete_san_c <- rcpp(
  signature(C_init="list",
            p_S="numeric", p_0="numeric", p_R="numeric",
            p_A="numeric", p_N="numeric", p_D="numeric",
            steps="integer"),
  includes=c("#include <R_ext/Utils.h>", "#include <Rmath.h>"),
  Rcpp=TRUE,
  cxxargs=c("-O3", "-fast-math"),
  body='
using namespace Rcpp;
  
/* Extract & validate parameters */
List C_init_(C_init);
IntegerVector C_curr_S(clone(as<IntegerVector>(C_init_["S"])));
IntegerVector C_curr_A(clone(as<IntegerVector>(C_init_["A"])));
IntegerVector C_curr_N(clone(as<IntegerVector>(C_init_["N"])));
const int L_ = IntegerVector(C_curr_S).size();
if ((C_curr_A.size() != L_) || (C_curr_N.size() != L_))
  stop("list S_init must contain initial count vectors S, A, N of the same length");
if ((DoubleVector(p_S).size() != 1) || (DoubleVector(p_0).size() != 1) || (DoubleVector(p_R).size() != 1) ||
    (DoubleVector(p_A).size() != 1) || (DoubleVector(p_N).size() != 1) ||
    (DoubleVector(p_D).size() != 1))
  stop("p_S, p_0, p_A, p_N, p_D must have length 1");
const double p_S_ = DoubleVector(p_S)[0];
const double p_0_ = DoubleVector(p_0)[0];
const double p_R_ = DoubleVector(p_R)[0];
const double p_A_ = DoubleVector(p_A)[0];
const double p_N_ = DoubleVector(p_N)[0];
const double p_D_ = DoubleVector(p_D)[0];
const IntegerVector stepsv = steps;

/* Prepare result */
List result;

/* Run simulation over all requested time intervals */
RNGScope rngScope;
for(int i=0; i < stepsv.size(); ++i) {
  int* const C_curr_S_ = INTEGER(C_curr_S);
  int* const C_curr_A_ = INTEGER(C_curr_A);
  int* const C_curr_N_ = INTEGER(C_curr_N);
  
  /* Proceed until the end of the current interval */
  for(int l=0, l_end=stepsv[i]; l < l_end; ++l) {
    /* Process lineages one by one */
    for(int j=0; j < L_; ++j) {
      /* Prevent overflows */
      if ((C_curr_S_[j] >= 1000000000) || (C_curr_A_[j] >= 1000000000) || (C_curr_N_[j] >= 1000000000))
        stop("S, A or N cell count overflowed (>= 1,000,000,000)");
        
      /* Determine the number of cells affected by each type of event */
      const int dS = (p_S_ > 0) ? R::rbinom(C_curr_S_[j], p_S_) : 0;
      const int d0 = (p_0_ > 0) ? R::rbinom(C_curr_S_[j], p_0_) : 0;
      const int dR = (p_R_ > 0) ? R::rbinom(C_curr_S_[j], p_R_) : 0;
      const int dA = (p_A_ > 0) ? R::rbinom(C_curr_S_[j], p_A_) : 0;
      const int dN = (p_N_ > 0) ? R::rbinom(C_curr_A_[j], p_N_) : 0;
      const int dD = (p_D_ > 0) ? R::rbinom(C_curr_A_[j], p_D_) : 0;
      
      /* Update cell counts */
      C_curr_S_[j] = std::max(C_curr_S_[j] + dS - d0 - dR - dA          , 0);
      C_curr_A_[j] = std::max(C_curr_A_[j]                + dA     - dD , 0);
      C_curr_N_[j] = std::max(C_curr_N_[j]           + dR      +dN + dD , 0);
    }
  }
  
  /* Append current per-lineage cell counts to the result list */
  List C_curr = List::create(
    Named("i") = IntegerVector(1, i+1),
    Named("S") = clone(C_curr_S),
    Named("A") = clone(C_curr_A),
    Named("N") = clone(C_curr_N)
  );
  result.push_back(C_curr);
}

/* Return result */
return result;
')

#' Determine a sensible time step
#' 
#' The time step is chosen to reduce the probability of two sequential
#' events (i.e. where the second event affects a cell created or destroyed
#' by the first event) below \code{p_cutoff}
#' 
#' @param rate the maximal rate in the model
#'
find_time_step <- function(rate, p_cutoff) {
  # Compute a time step dt such that for events occuring with the given
  # rate per time unit the probability of more than one event per time unit
  # is at most p_cutoff. In other words, determine lambda such that
  #   P(0 or 1 event) = e^-lambda + lambda*e^-lambda >= 1 - p_cutoff,
  # and then set
  #   dt = lambda / rate
  (-lambertWm1((p_cutoff-1)/exp(1))-1) / rate
}

#' Simulate the S model
#' 
simulate_s <- function(L, s0=1, rates, samples_per_day=1, p_cutoff=1e-3) {
  # The names of the rate parameters
  RATES <- c("r_S", "r_R")
  
  # rates can change over times, and are therefore passed as a list of lists. Each
  # list must specify rates r_S and r_R plus a Tmax. The Tmax must increase
  # monotonically, and each parameter set is used for times between the previous set's
  # Tmax and the set's own Tmax. If the list contains no nested list, or contains a value
  # Tmax, the code assumes that only a single set or parameters is to be used for all times t.
  if (!is.list(rates[[1]]) || !is.null(rates$Tmax))
    rates <- list(rates)
  
  # Check that all necessary parameters were specified correctly
  if (!is.numeric(L) || (length(L) != 1) || (L %% 1 != 0))
    stop("L must be a single integer")
  if (!is.numeric(s0) || (length(s0) != 1) || !is.finite(s0) || (s0 < 0) || (s0 %% 1 != 0))
    stop("s0 must be be a single non-negative and finite integral value")
  if (!is.numeric(samples_per_day) || (length(samples_per_day) != 1) || !is.finite(samples_per_day) || (samples_per_day <= 0) || (samples_per_day %% 1 != 0))
    stop("samples_per_day must be be a single positive and finite integral value")
  if (!is.numeric(p_cutoff) || (length(p_cutoff) != 1) || !is.finite(p_cutoff) || (p_cutoff < 0) || (p_cutoff >= 1))
    stop("p_cutoff must be be a single numeric value within [0,1)")
  if (!is.list(rates))
    stop("rates must be a list")
  Tmax <- 0
  for(r in rates) {
    # Check rates
    for(n in RATES) {
      if (is.null(r[[n]]) || !is.numeric(r[[n]]) || (length(r[[n]]) != 1) || !is.finite(r[[n]]) || (r[[n]] < 0))
        stop(paste(n, " must be a a single non-negative and finite numeric value"))
    }
    # Check Tmax
    if (!is.numeric(r$Tmax) || (length(r$Tmax) != 1) || !is.finite(r$Tmax) || (r$Tmax < 0) || (r$Tmax %% 1 != 0) || (r$Tmax < Tmax))
      stop(paste("Tmax must be a a single non-negative and finite integral value, and must increase monotonically"))
    Tmax <- r$Tmax
  }
  
  # Count vectors for the S-cell counts
  state <- data.table(t=0, dt=NA_real_, lid=1L:L, S=rep(as.integer(s0), L))
  
  # Simulate
  ri <- 0
  r <- list(Tmax=0)
  rows <- list(state)
  t <- 0
  while (t < Tmax) {
    # Switch to the next set of rates if necessary
    while (r$Tmax <= t) {
      # Next set of rates
      ri <- ri + 1
      r <- rates[[ri]]
      
      # Determine time step.
      # First, determine maximum time step that makes the probability of more than one event
      # per cell and time step less than p_cutoff
      dt.max <- find_time_step(do.call(max, r[RATES]), p_cutoff=p_cutoff)
      # Now find the smallest integral value steps_per_sample such that
      #       1/(steps_per_sample*samples_per_day) <= dt.max
      #   <=> steps_per_sample*samples_per_day >= 1/dt.max
      #   <=> steps_per_sample >= 1/(samples_per_day*dt.max)
      steps_per_sample = max(1, ceiling(1 / (samples_per_day*dt.max)))
      # Finally, compute the actual time step dt
      dt <- 1/(steps_per_sample*samples_per_day)
      stopifnot(dt <= dt.max)
    }
    
    # Simulate as far as the current set of rates is valid, output the state at the end of each day
    #
    # See the discussion in simulate_san for a discussion of how p_S and p_R are computed
    res <- discrete_s_c(state,
                        p_S=expm1(r$r_S*dt), p_R=expm1(r$r_R*dt),
                        steps=rep(steps_per_sample, (r$Tmax - t)*samples_per_day))
    # Append to results
    rows.res <- rbindlist(res)[, list(t=t+i/samples_per_day, dt=dt, lid=1L:L, S)]
    rows <- c(rows, list(rows.res))
    # Update state
    state <- res[[length(res)]]
    stopifnot(rows.res[, max(t)] == r$Tmax)
    t <- r$Tmax
  }
  
  # Return per-lineage cell-count table
  return(rbindlist(rows))
}

#' Simulate the SAN model
#' 
simulate_san <- function(L, rates, samples_per_day=1, p_cutoff=1e-3) {
  # The names of the rate parameters
  RATES <- c("r_S", "r_0", "r_R", "r_A", "r_N", "r_D")

  # rates can change over times, and are therefore passed as a list of lists. Each
  # list must specify rates r_S, r_0, r_R, r_A, r_N, r_D plus a Tmax. The Tmax must increase
  # monotonically, and each parameter set is used for times between the previous set's
  # Tmax and the set's own Tmax. If the list contains no nested list, or contains a value
  # Tmax, the code assumes that only a single set or parameters is to be used for all times t.
  if (!is.list(rates[1]) || !is.null(rates$Tmax))
    rates <- list(rates)

  # Check that all necessary parameters were specified correctly
  if (!is.numeric(L) || (length(L) != 1) || (L %% 1 != 0))
    stop("L must be a single integer")
  if (!is.numeric(samples_per_day) || (length(samples_per_day) != 1) || !is.finite(samples_per_day) || (samples_per_day <= 0) || (samples_per_day %% 1 != 0))
    stop("samples_per_day must be be a single positive and finite integral value")
  if (!is.numeric(p_cutoff) || (length(p_cutoff) != 1) || !is.finite(p_cutoff) || (p_cutoff < 0) || (p_cutoff >= 1))
    stop("p_cutoff must be be a single numeric value within [0,1)")
  if (!is.list(rates))
    stop("rates must be a list")
  Tmax <- 0
  for(r in rates) {
    # Check rates
    for(n in RATES) {
      if (!is.numeric(r[[n]]) || (length(r[[n]]) != 1) || !is.finite(r[[n]]) || (r[[n]] < 0))
        stop(paste(n, " must be a a single non-negative and finite numeric value"))
    }
    # Check Tmax
    if (!is.numeric(r$Tmax) || (length(r$Tmax) != 1) || !is.finite(r$Tmax) || (r$Tmax < 0) || (r$Tmax %% 1 != 0) || (r$Tmax < Tmax))
      stop(paste("Tmax must be a a single non-negative and finite numeric value, and must increase monotonically"))
    Tmax <- r$Tmax
  }
  
  # Count vectors for the three cell types S, A, N
  state <- data.table(t=0, dt=NA_real_, lid=1L:L, S=rep(1L, L), A=rep(0L, L), N=rep(0L, L))
  
  # Simulate
  ri <- 0
  r <- list(Tmax=0)
  rows <- list(state)
  t <- 0
  while (t < Tmax) {
    # Switch to the next set of rates if necessary
    while (r$Tmax <= t) {
      # Next set of rates
      ri <- ri + 1
      r <- rates[[ri]]
      
      # Determine time step.
      # First, determine maximum time step that makes the probability of more than one event
      # per cell and time step less than p_cutoff
      dt.max <- find_time_step(do.call(max, r[RATES]), p_cutoff=p_cutoff)
      # Now find the smallest integral value steps_per_sample such that
      #       1/(steps_per_sample*samples_per_day) <= dt.max
      #   <=> steps_per_sample*samples_per_day >= 1/dt.max
      #   <=> steps_per_sample >= 1/(samples_per_day*dt.max)
      steps_per_sample = max(1, ceiling(1 / (samples_per_day*dt.max)))
      # Finally, compute the actual time step dt
      dt <- 1/(steps_per_sample*samples_per_day)
      stopifnot(dt <= dt.max)
    }

    # Simulate as far as the current set of rates is valid, output the state at the end of each day
    #
    # discrete_san_c simulates the system in discrete time steps, and requires the probabilities
    # p_S, p_0, ... of an event to occurs within one time step to be specified. We chose those
    # probabilities such that the *expected* number of events matches the time-continuous process.
    #
    # In a simple one-type birth-only time-continuous process with birth rate r, the expected number
    # of individuals at time t+dt is C(t) * e^(r*dt) if C(t) is the number of individuals at time t.
    # The expected number of birth events per individual within a time interval dt is thus
    #    e^(r*dt) - 1,
    # and generalizing this, we use the probability
    #    p_x=e^(r_x*dt) - 1
    # for the per-individual occurence of an event of type x within a time interval of length dt.
    #
    # Note if |z| << 1, e^z is close to 1, and using `exp(z)-1` to evaluate e^z - 1 thus suffers
    # from loss of precision. We thus use the function `expm1` which directly evaluates e^z - 1 to
    # circumvent this cancellation issue.
    res <- discrete_san_c(state,
                          p_S=expm1(r$r_S*dt), p_0=expm1(r$r_0*dt), p_R=expm1(r$r_R*dt),
                          p_A=expm1(r$r_A*dt), p_N=expm1(r$r_N*dt), p_D=expm1(r$r_D*dt),
                          steps=rep(steps_per_sample, (r$Tmax - t)*samples_per_day))
    
    # Append to results
    rows.res <- rbindlist(res)[, list(t=t+i/samples_per_day, dt=dt, lid=1L:L, S, A, N)]
    rows <- c(rows, list(rows.res))
    # Update state
    state <- res[[length(res)]]
    stopifnot(rows.res[, max(t)] == r$Tmax)
    t <- r$Tmax
  }
  
  # Return per-lineage cell-count table
  return(rbindlist(rows))
}

RATES_SIMON <- list(
  list(Tmax=3,  r_S=0,     r_0=0,   r_R=0, r_A=0,   r_N=0,   r_D=0),
  list(Tmax=9,  r_S=0.349, r_0=0,   r_R=0, r_A=0,   r_N=0,   r_D=0),
  list(Tmax=41, r_S=10,    r_0=9.6, r_R=0, r_A=0.4, r_N=0.2, r_D=2.6)
)
