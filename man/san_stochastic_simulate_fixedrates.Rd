% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/san_stochastic.R
\name{san_stochastic_simulate_fixedrates}
\alias{san_stochastic_simulate_fixedrates}
\title{Simulate the stochastic SAN model for with constant rates in discrete steps}
\usage{
san_stochastic_simulate_fixedrates(x0, rates, steps, dt)
}
\arguments{
\item{x0}{a list containing vectors `S`, `A`, `N` with the initial cell-counts per lineage}

\item{rates}{a names list with non-negative finite values for r_S, r_0 r_R, r_A, r_N, r_D.}

\item{steps}{a vector listing the number of steps between each sample}

\item{dt}{the time each step corresponds to}
}
\description{
See @seealso san_stochastic for details about the stochastic SAN model
}
