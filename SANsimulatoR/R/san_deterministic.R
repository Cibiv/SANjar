# The names of the rate parameters
SAN.INITIALVALUES <- c("s0", "a0", "n0")
SAN.RATENAMES <- c("r_S", "r_0", "r_R", "r_A", "r_N", "r_D")

# Abbreviations used in the analytical solutions below
SAN.DETERMINISTIC.SOLUTIONS.ABBREVIATIONS <- list(
  gamma       =quote(r_S - r_0 - r_R - r_A),
  theta       =quote(gamma + r_D),
  rho         =quote(r_R / gamma),
  alpha       =quote(r_A / gamma),
  alpha.tilde =quote(r_A / theta),
  nu          =quote(r_N / gamma),
  delta       =quote(r_D / gamma),
  beta        =quote(r_A / r_D),
  omega       =quote((r_D + r_N) / r_D)
)

# Substitute abbreviations and additional substitutions to create a solution expression
make_san_solution <- function(expr, subst) {
  # Get unenvaluated expression passed for the "expr" parameter
  expr <- substitute(expr)
  # Substitute for the times parameter to make sure the evaluation times are a row vector
  # Note: This is necessary to ensure that %*% gets interpreted as an outer product
  # in the solution expressions. Using %o% would be preferrable, but it doesn't seem to
  # be supported by the Deriv package for symbolic differentations, so we stick with %*%
  expr <- do.call(substitute, list(expr, list(times=quote(matrix(times, nrow=1)),
                                              const=quote(matrix(1, nrow=1, ncol=length(times))))))
  # Substitute the substitutions in "subst"
  expr <- do.call(substitute, list(expr, subst))
  # Substitute the abbreviations.
  # Must substitute thrice because the abbreviations are defined in three tiers,
  # i.e. alpha.tilde depends on theta depends on gamma.
  expr <- do.call(substitute, list(expr, SAN.DETERMINISTIC.SOLUTIONS.ABBREVIATIONS))
  expr <- do.call(substitute, list(expr, SAN.DETERMINISTIC.SOLUTIONS.ABBREVIATIONS))
  expr <- do.call(substitute, list(expr, SAN.DETERMINISTIC.SOLUTIONS.ABBREVIATIONS))
}

# Creates the 3x1 column vector (expr_s, expr_a, expr_n)^T from the *unevaluated* arguments
make_vector <- function(expr_s, expr_a, expr_n) {
  substitute(matrix(c(expr_s, expr_a, expr_n)))
}

# The analytical solutions for the 5 different parameter regimes.
# All abbreviations are substituted away, the expressions only refer to the
# model parameters, the inital values s0, a0, n0 and the evaluation times t.
# The expressions in the solutions below consist of sums of terms of the form
#   (b * v) %*% F
# or
#   (b * v) %*% C
# where b is a scalar, v = (v_S, v_A, v_N)^T a 3x1 column vector,
# F = ( f(t_1) ... f(t_k) ) a 1xk row vector containing transformed versions
# of the k evaluation times t1_1 ... t_k, and C the constant 1xk row vector
# (1 ... 1). Since the LHS of the matrix product %*% is a column- and the RHS
# a row-vector, %*% is an outer product and will return the matrix
#   ( d * v_S * f(t_1)   d * v_S * f(t_2)   ...   d * v_S * f(t_k) )
#   ( d * v_A * f(t_1)   d * v_A * f(t_2)   ...   d * v_A * f(t_k) )
#   ( d * v_N * f(t_1)   d * v_N * f(t_2)   ...   d * v_N * f(t_k) )
# respectively
#   ( d * v_S            d * v_S            ...   d * v_S          )
#   ( d * v_A            d * v_A            ...   d * v_A          )
#   ( d * v_N            d * v_N            ...   d * v_N          ).
# The result of evaluation these expressions is a matrix of the form
#   ( s(t_1)             s(t_2)             ...   s(t_k)           )
#   ( a(t_1)             a(t_2)             ...   a(t_k)           )
#   ( n(t_1)             n(t_2)             ...   n(t_k)           ).
# Note that the substitutions done by make_san_solution guarantee that F
# and C have the desired shapes.
SAN.DETERMINISTIC.SOLUTION.REGIMES <- list(
  s_noneq.a_fin.general=make_san_solution(
    (  (s0 * v.g.s0)                             %*% exp(gamma * times)
     + (s0 * v.rD.s0 + a0 * v.rD.a0)             %*% exp(-r_D * times)
     + (s0 * v.c.s0 + a0 * v.c.a0 + n0 * v.c.n0) %*% const),
    list(
      v.g.s0 =make_vector(1,  alpha.tilde, alpha.tilde * delta * omega + rho),
      v.rD.s0=make_vector(0, -alpha.tilde,               alpha.tilde * omega),
      v.rD.a0=make_vector(0,            1,                           - omega),
      v.c.s0 =make_vector(0,            0,             - alpha * omega - rho),
      v.c.a0 =make_vector(0,            0,                             omega),
      v.c.n0 =make_vector(0,            0,                                 1)
    )),
  
  s_noneq.a_fin.g_eq_mrD=make_san_solution(
    (  (s0 * v.g.s0 + a0 * v.g.a0)               %*% exp(gamma * times)
     + (s0 * v.tg.s0)                            %*% (times * exp(gamma * times))
     + (s0 * v.c.s0 + a0 * v.c.a0 + n0 * v.c.n0) %*% const),
    list(
      v.g.s0 =make_vector(1,            0,               alpha * omega + rho),
      v.g.a0 =make_vector(0,            1,                            -omega),
      v.tg.s0=make_vector(0,          r_A,                     - r_A * omega),
      v.c.s0 =make_vector(0,            0,             - alpha * omega - rho),
      v.c.a0 =make_vector(0,            0,                             omega),
      v.c.n0 =make_vector(0,            0,                                 1)
    )),
  
  s_noneq.a_inf=make_san_solution(
    (  (s0 * v.g.s0)                             %*% exp(gamma * times)
     + (s0 * v.t.s0 + a0 * v.t.a0)               %*% times
     + (s0 * v.c.s0 + a0 * v.c.a0 + n0 * v.c.n0) %*% const),
    list(
      v.g.s0=make_vector(1,        alpha,                   alpha * nu + rho),
      v.t.s0=make_vector(0,            0,                      - alpha * r_N),
      v.t.a0=make_vector(0,            0,                                r_N),
      v.c.s0=make_vector(0,       -alpha,                 - alpha * nu - rho),
      v.c.a0=make_vector(0,            1,                                  0),
      v.c.n0=make_vector(0,            0,                                  1)
    )),
  
  s_equil.a_fin=make_san_solution(
    (  (s0 * v.rD.s0 + a0 * v.rD.a0)             %*% exp(-r_D * times)
     + (s0 * v.t.s0)                             %*% times
     + (s0 * v.c.s0 + a0 * v.c.a0 + n0 * v.c.n0) %*% const),
    list(
      v.rD.s0=make_vector(0,        -beta,                      beta * omega),
      v.rD.a0=make_vector(0,            1,                           - omega),
      v.t.s0 =make_vector(0,            0,                 r_A * omega + r_R),
      v.c.s0 =make_vector(1,         beta,                    - beta * omega),
      v.c.a0 =make_vector(0,            0,                             omega),
      v.c.n0 =make_vector(0,            0,                                 1)
    )),
  
  s_equil.a_inf=make_san_solution(
    (  (s0 * v.t2.s0)                            %*% (times^2)
     + (s0 * v.t.s0 + a0 * v.t.a0)               %*% times
     + (s0 * v.c.s0 + a0 * v.c.a0 + n0 * v.c.n0) %*% const),
    list(
      v.t2.s0=make_vector(0,            0,                     r_A * r_N / 2),
      v.t.s0 =make_vector(0,          r_A,                               r_R),
      v.t.a0 =make_vector(0,            0,                               r_N),
      v.c.s0 =make_vector(1,            0,                                 0),
      v.c.a0 =make_vector(0,            1,                                 0),
      v.c.n0 =make_vector(0,            0,                                 1)
    ))
)

if (FALSE) {
# Disabled as long as we dont use it to save loading time
SAN.DETERMINISTIC.SOLUTION.REGIMES.DIFF <- lapply(SAN.DETERMINISTIC.SOLUTION.REGIMES, function(expr) {
  params <- c(SAN.INITIALVALUES, SAN.RATENAMES)
  d <- lapply(params, function(x) {
    Deriv(expr, x=x)
  })
  names(d) <- as.vector(params)
  return(d)
})
}

#' Evaluates the analytical solution of the deterministic SAN model
#' 
#' Evaluates the anyltical solution at the specified times, starting
#' from the cell counts specified in `x0` at time `t=0`, and uses
#' the rates specified via `r`.
#' 
#' @param x0 initial cell counts at time `t=0`
#' @param times the times to evaluate the solution at
#' @param rates a names list with non-negative finite values for r_S, r_0 r_R, r_A, r_N, r_D.
#' @param eps the epsilon value used when determining the parameter regime
#' @return a matrix with columns `t`, `S`, `A`, `N` and one row per time point
san_deterministic_eval_fixedrates <- function(x0, times, rates, eps=1e-3) {
  # Check parameters
  if (!is.numeric(x0) || (length(x0) != 3) || !all(is.finite(x0)) || !all(x0 >= 0) ||
      is.null(names(x0)) || !all(names(x0) == c("S", "A", "N")))
    stop("initial values in x0 must be numeric, non-negative, finite, and named 'S', 'A', 'N'")
  if (!is.numeric(times) || !all(is.finite(times)) || !all(times >= 0))
    stop("evaluation times must be numeric, non-negative and finite")
  if (!is.list(rates))
    stop("rates must be a named list")
  for(n in SAN.RATENAMES) {
    if (!is.numeric(rates[[n]]) || (length(rates[[n]]) != 1) || !is.finite(rates[[n]]) || (rates[[n]] < 0))
      stop(paste("rate ", n, " must be a a single non-negative and finite numeric value"))
  }

  # Select regime
  regime <- with(rates, {
    gamma <- eval(SAN.DETERMINISTIC.SOLUTIONS.ABBREVIATIONS$gamma)
    theta <- eval(SAN.DETERMINISTIC.SOLUTIONS.ABBREVIATIONS$theta)
    if ((abs(gamma) >= eps) && (abs(r_D) >= eps) && (abs(theta) >= eps))
      # S-cell non-equilibrium, finite A-cell lifetime, general case
      "s_noneq.a_fin.general"
    else if ((abs(gamma) >= eps) && (abs(r_D) >= eps) && (abs(theta) < eps))
      # S-cell non-equilibrium, finite A-cell lifetime, special case r_D == -gamma
      "s_noneq.a_fin.g_eq_mrD"
    else if ((abs(gamma) >= eps) && (abs(r_D) < eps))
      # S-cell non-equilibrium, infinite A-cell lifetime
      "s_noneq.a_inf"
    else if ((abs(gamma) < eps) && (abs(r_D) >= eps))
      # S-cell equilibrium, finite A-cell lifetime,
      "s_equil.a_fin"
    else if ((abs(gamma) < eps) && (abs(r_D) < eps))
      # S-cell equilibrium, infinite A-cell lifetime
      "s_equil.a_inf"
  })

  # Evaluate solution and return results.
  expr <- SAN.DETERMINISTIC.SOLUTION.REGIMES[[regime]]
  x <- eval(expr, enclos=baseenv(), envir=c(
    list(s0=x0['S'], a0=x0['A'], n0=x0['N'], times=times),
    rates
  ))
  rownames(x) <- c("S", "A", "N")
  r <- cbind(t=as.vector(times), t(x))
  attr(r, "regime") <- regime
  return(r)
}

#' Computes the number of S-, A- and N-cells under the deterministic SAN model
#' 
#' The deterministic SAN model describes the *average* behavior of the stochastic
#' SAN model, and is the defined by the ODE
#' \deqn{S' = (r_S - r_0 - r_R - r_A)\cdot S}{S' = (r_S - r_0 - r_R - r_A) * S}
#' \deqn{A' = r_A\cdot S - r_D\cdot A}{A' = r_A * S - r_D * A}
#' \deqn{N' = r_R\cdot S + r_N\cdot A}{N' = r_R * S + r_N * A}
#' This function computes \eqn{(S(t), A(t), N(t))} from \eqn{t=0} onwards
#' given \eqn{s_0=S(0)} and time-dependent but piece-wise constant rates
#' \eqn{r_S}, \eqn{r_0}, \eqn{r_R}, \eqn{r_A}, \eqn{r_N}, \eqn{r_D}. 
#' 
#' The picewise constant rates are specified through the table `rates` with
#' columns `Tmax`, `r_S`, `r_0`, `r_R`, `r_A`, `r_N`. The parameters in each row
#' are in effect between the previous row's `Tmax` and the row's own `Tmax`.
#' 
#' @param s0 the initial number of S-cells
#' @param rates a `data.table` with columns `Tmax`, `r_S`, `r_0`, `r_R`, `r_A`, `r_N`
#'              and `r_D`  and monotonically increasing values in the column `Tmax`.
#' @param samples_per_day the number of (equally spaced) times point per day at
#'                        which to output the cell counts
#'                        
#' @return a `data.table` with columns `t`, `S`, `A`, `N` containing the cell counts
#'         at each day from \eqn{t=0} to \eqn{t=T}, where \eqn{T} is the value of
#'         `Tmax` in the last row of the rates table. For each day, the table contains
#'         `samples_per_day` rows with equally spaced evaluation times within that day. 
#'
#' @import data.table
#' @export
san_deterministic <- function(s0, rates, samples_per_day=1) {
  # Check that all necessary parameters were specified correctly
  if (!is.numeric(samples_per_day) || (length(samples_per_day) != 1) || !is.finite(samples_per_day) || (samples_per_day <= 0) || (samples_per_day %% 1 != 0))
    stop("samples_per_day must be be a single positive and finite integral value")
  if (!is.list(rates))
    stop("rates must be a list")
  if (!is.data.frame(rates) || (nrow(rates) == 0))
    stop("rates must be a non-empty data.table (or.data.frame")
  rates <- as.data.table(rates)
  for (c in SAN.RATENAMES)
    if (!is.numeric(rates[[c]]) || !all(is.finite(rates[[c]])) || any(rates[[c]] < 0))
      stop(paste(n, " must contain non-negative and finite numeric values"))
  if (is.unsorted(rates$Tmax))
    stop(paste("Tmax must contain non-negative and finite numeric values, and must increase monotonically"))
  Tmax <- max(rates$Tmax)
  
  # Count vectors for the three cell types S, A, N
  state <- data.table(t=0, S=s0, A=0, N=0, regime=NA_character_)
  
  # Simulate
  ri <- 0
  r <- list(Tmax=0)
  rows <- list(state)
  tstart <- 0
  while (tstart < Tmax) {
    # Switch to the next set of rates if necessary
    while (r$Tmax <= tstart) {
      # Next set of rates
      ri <- ri + 1
      r <- rates[ri]
    }
    
    # Evaluate as far as the current set of rates is valid, output the state at the end of each day
    res <- san_deterministic_eval_fixedrates(
      x0=c(S=state$S, A=state$A, N=state$N),
      times=tail(seq(from=0, to=r$Tmax - tstart, length.out=(r$Tmax - tstart + 1)*samples_per_day), n=-1),
      rates=r
    )

    # Append to results
    rows.res <- as.data.table(res)[, list(t=tstart+t, S, A, N, regime=attr(res, "regime"))]
    rows <- c(rows, list(rows.res))
    # Update state
    state <- rows.res[nrow(res),]
    stopifnot(rows.res[, max(t)] == r$Tmax)
    tstart <- r$Tmax
  }
  
  # Return per-lineage cell-count table
  return(rbindlist(rows))
}
