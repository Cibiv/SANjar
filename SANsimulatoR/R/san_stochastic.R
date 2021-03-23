#' Determine a sensible time step
#' 
#' The time step is chosen to reduce the probability of two sequential
#' events (i.e. where the second event affects a cell created or destroyed
#' by the first event) below \code{p_cutoff}
#' 
#' @param rate the maximal rate in the model
#' @param p_cutoff the maximal allowed probability of two sequential events
#'                 within one time step
#'
#' @import lamW
#' @export
find_time_step <- function(rate, p_cutoff) {
  # Compute a time step dt such that for events occuring with the given
  # rate per time unit the probability of more than one event per time unit
  # is at most p_cutoff. In other words, determine lambda such that
  #   P(0 or 1 event) = e^-lambda + lambda*e^-lambda >= 1 - p_cutoff,
  # and then set
  #   dt = lambda / rate
  (-lambertWm1((p_cutoff-1)/exp(1))-1) / rate
}

#' Simulate the stochastic SAN model
#' 
#' The stochastic SAN model describes the stochastic behavior of lineages
#' consisting of S, A and N-cells which can undergo the following conversions
#' in a memoryless fashion, making the model a Markov process.
#' 
#' | event                         | rate      |
#' | \eqn{S \to S S}{S -> S S}     | \eqn{r_S} |
#' | \eqn{S \to \emptyset}{S -> 0} | \eqn{r_0} |
#' | \eqn{S \to N}{S -> N}         | \eqn{r_R} |
#' | \eqn{S \to A}{S -> A}         | \eqn{r_A} |
#' | \eqn{A \to A N}{S -> A N}     | \eqn{r_N} |
#' | \eqn{A \to N}{A -> N}         | \eqn{r_D} |
#' 
#' The picewise constant rates are specified through the table `rates` with
#' columns `Tmax`, `r_S`, `r_0`, `r_R`, `r_A`, `r_N`. The parameters in each row
#' are in effect between the previous row's `Tmax` and the row's own `Tmax`.
#' 
#' *Note*: Typically, the deterministic SAN model is applied to a system consisting
#' of \eqn{s_0} (\eqn{\gg 1}) initial S-cells. The point of the *stochastic* SAN model
#' is typically to simulate the lineages arising from these cells separately, meaning
#' that in the stochastic model, one will typically set \eqn{s_0} to 1, and \eqn{L} to
#' the initial number of S-cells. In a way, it is thus the parameter \eqn{L} here that
#' corresponds to the parameter \eqn{s_0} of the deterministic model. 
#'
#' @param L the number of lineages
#' @param s0 the initial number of S-cells *in each lineage* (default: `1`)
#' @param previous the result of a previously san_stochastic invocation to be continued from
#' @param rates a `data.table` with columns `Tmax`, `r_S`, `r_0`, `r_R`, `r_A`, `r_N`
#'              and `r_D`  and monotonically increasing values in the column `Tmax`.
#' @param Tmax the stopping time (defaults to the largest Tmax in `rates`)
#' @param samples_per_day the number of (equally spaced) times point per day at
#'                        which to output the cell counts
#' @param p_cutoff the maximal allowed probability of two sequential events
#'                 within one time step
#'                        
#' @return a `data.table` with columns `t`, `S`, `A`, `N` containing the cell counts
#'         at each day from \eqn{t=0} to \eqn{t=T}, where \eqn{T} is the value of
#'         `Tmax` in the last row of the rates table. For each day, the table contains
#'         `samples_per_day` rows with equally spaced evaluation times within that day. 
#'
#' While the model as defined above is time-continuous, this function simulates
#' the model in discrete time steps. The parameter `p_cutoff` controls the accuracy
#' of this discretization, see \code{\link[=func]{find_time_step()}}
#'
#' @import data.table
#' @useDynLib SANsimulatoR, .registration = TRUE
#' @export 
san_stochastic <- function(L=NA, s0=1, rates, previous=NULL, Tmax=max(rates$Tmax), samples_per_day=1, p_cutoff=1e-3) {
  # Check that all necessary parameters were specified correctly
  if (is.null(previous)) {
    if (!is.numeric(L) || (length(L) != 1) || !is.finite(L) || (L %% 1 != 0))
      stop("L must be a single non-negative integer")
    if (!is.numeric(s0) || (length(s0) != 1) || !is.finite(s0) || (s0 %% 1 != 0) || (s0 < 0))
      stop("s0 must be a single finite and non-negative integer")
  } else {
    if (!is.data.table(previous) || any(colnames(previous) != c("t", "dt", "lid", "S", "A", "N")))
      stop("previous must be the result of a previous san_stochastic call")
  }
  if (!is.numeric(samples_per_day) || (length(samples_per_day) != 1) || !is.finite(samples_per_day) || (samples_per_day <= 0) || (samples_per_day %% 1 != 0))
    stop("samples_per_day must be be a single positive and finite integral value")
  if (!is.numeric(p_cutoff) || (length(p_cutoff) != 1) || !is.finite(p_cutoff) || (p_cutoff < 0) || (p_cutoff >= 1))
    stop("p_cutoff must be be a single numeric value within [0,1)")
  if (!is.data.frame(rates) || (nrow(rates) == 0))
    stop("rates must be a non-empty data.table (or.data.frame")
  rates <- as.data.table(rates)
  for (c in SAN.RATENAMES)
    if (!is.numeric(rates[[c]]) || !all(is.finite(rates[[c]])) || any(rates[[c]] < 0))
      stop(paste(n, " must contain non-negative and finite numeric values"))
  if (is.unsorted(rates$Tmax))
    stop(paste("Tmax must contain non-negative and finite numeric values, and must increase monotonically"))
  if (Tmax > max(rates$Tmax))
    stop("specified stopping time Tmax exceeds largest Tmax in rates table")

  # Initialize state, either for t=0 or to continue from where we previously left off
  if (is.null(previous)) {
    # Start fresh
    t <- 0
    state <- list(S=rep(as.integer(s0), L), A=rep(0L, L), N=rep(0L, L))
    rows <- list(data.table(t=0, dt=NA_real_, lid=1L:L, S=state$S, A=state$A, N=state$N))
  } else {
    # Continue from previous output
    t <- previous[, max(t)]
    t_ <- t
    p <- previous[t == t_][order(lid)]
    L <- p$lid[nrow(p)]
    if (any(p$lid != 1:L))
      stop("column lid of argument previous must contain LIDs 1 through L")
    state <- p[, list(S, A, N)]
    rows <- list(previous)
  }
  
  # Simulate
  ri <- 0
  r <- list(Tmax=0)
  while (t < Tmax) {
    # Switch to the next set of rates if necessary
    while (r$Tmax <= t) {
      # Next set of rates
      ri <- ri + 1
      r <- rates[ri]
    }

    # Determine stopping time for the current set of rates
    r_Tmax <- min(r$Tmax, Tmax)

    # Determine time step.
    # First, determine maximum time step that makes the probability of more than one event
    # per cell and time step less than p_cutoff
    dt.max <- find_time_step(do.call(max, r[, ..SAN.RATENAMES]), p_cutoff=p_cutoff)
    # Now find the smallest integral value steps_per_sample such that
    #       1/(steps_per_sample*samples_per_day) <= dt.max
    #   <=> steps_per_sample*samples_per_day >= 1/dt.max
    #   <=> steps_per_sample >= 1/(samples_per_day*dt.max)
    steps_per_sample = max(1, ceiling(1 / (samples_per_day*dt.max)))
    # Finally, compute the actual time step dt
    dt <- 1/(steps_per_sample*samples_per_day)
    stopifnot(dt <= dt.max)

    # Simulate as far as the current set of rates is valid, output the state at the end of each day
    #
    # discrete_san_c simulates the system in discrete time steps, and requires the probabilities
    # p_S, p_0, ... of an event to occurs within one time step to be specified.
    res <- san_timediscrete_c(state,
                              p_S=r$r_S*dt, p_0=r$r_0*dt, p_R=r$r_R*dt,
                              p_A=r$r_A*dt, p_N=r$r_N*dt, p_D=r$r_D*dt,
                              steps=rep(steps_per_sample, (r_Tmax - t)*samples_per_day))
    

    # Append to results
    rows.res <- rbindlist(res)[, list(t=t+i/samples_per_day, dt=dt, lid=1L:L, S, A, N)]
    rows <- c(rows, list(rows.res))
    # Update state
    state <- res[[length(res)]]
    stopifnot(rows.res[, max(t)] == r_Tmax)
    t <- r_Tmax
  }
  
  # Return per-lineage cell-count table
  return(rbindlist(rows))
}
