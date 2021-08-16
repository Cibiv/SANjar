# Equation to solve to determine the lambda parameter
# for a zero-truncated Poisson distribution with mean m
tpois.lambda.eq <- deriv(expression((m - l/(1-exp(-l)))^2), namevec=c("l"),
                         function.arg=c("l", "m"))

#' Compute lambda parameter for a zero-truncated Poisson distribution with mean m
#' 
#' @param mean the desired mean of the zero-truncated Poisson distribution
#' @return the lambda Parameter of the truncated Poisson distribution
tpois.lambda <- function(mean) {
  if (is.na(mean) || (mean <= 1))
    return(NA_real_)
  if (is.infinite(mean))
    return(Inf)
  r <- nlm(tpois.lambda.eq, p=mean, m=mean)
  return(r$estimate)
}

