#' Simulate sequencing
#' 
#' To pass the results of `san_stochastic` as the `lineagesizes` parameters, the
#' `t` column has to be renamed to `day`, and columns `C` containing the total number
#' of cells (i.e., S+A+N) and `sid` (containing an arbitrary sample id) must be added.
#' 
#' @export
san_sequencing <- function(lineagesizes, parameters, steps=c("aliases", "seqsim", "threshold"), method.seqsim=NULL) {
  # Join sequencing parameters to lineage sizes
  ls <- lineagesizes[parameters, on=.(day), nomatch=NULL]

  if ("aliases" %in% steps) {
    # Compute lambda parameter of a zero-truncated Poisson distribution which yields
    # the requested mean number of aliases per lineage
    ls[, lineage_aliases_lambda := tpois.lambda(lineage_aliases[1]), by=.(lineage_aliases) ]
  
    # Add zero-truncated Poisson-distributed number of aliases of each lineage such
    # that the average number of aliases is lineage_aliases
    ls <- ls[, {
      stopifnot(all(lineage_aliases_lambda == lineage_aliases_lambda[1]))
      if (!is.na(lineage_aliases_lambda[1]) && (lineage_aliases_lambda[1] > 0))
        .SD[rep(1:.N, extraDistr::rtpois(.N, lambda=lineage_aliases_lambda[1], a=0))]
      else 
        SD
    }, keyby=.(day, sid) ]
  }

  if ("seqsim" %in% steps) {
    ls[, R := seqsim(C, reads.target=library_size[1], efficiency=pcr_efficiency[1], method=method.seqsim)
       , keyby=.(day, sid)]
    ls[, reads_per_cell := library_size[1] / sum(C)
       , keyby=.(day, sid) ]
    ls[, C.obs := R / reads_per_cell ]
    
    if ("threshold" %in% steps)
      ls <- ls[R >= phantom_threshold]
  }
  
  return(ls)
}
