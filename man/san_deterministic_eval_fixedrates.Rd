% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/san_deterministic.R
\name{san_deterministic_eval_fixedrates}
\alias{san_deterministic_eval_fixedrates}
\title{Evaluates the analytical solution of the deterministic SAN model}
\usage{
san_deterministic_eval_fixedrates(x0, times, rates, eps = 0.001)
}
\arguments{
\item{x0}{initial cell counts at time `t=0`}

\item{times}{the times to evaluate the solution at}

\item{rates}{a names list with non-negative finite values for r_S, r_0 r_R, r_A, r_N, r_D.}

\item{eps}{the epsilon value used when determining the parameter regime}
}
\value{
a matrix with columns `t`, `S`, `A`, `N` and one row per time point
}
\description{
If the goal is simply to evaluate the deterministic SAN model, `san_deterministic` should
be used which allows rates to vary over time, provides a user-friedly API,
and checks parameters for validity.
}
\details{
Evaluates the analytical solution at the specified times, starting
from the cell counts specified in `x0` at time `t=0`, and uses
the rates specified via `r`.
}
