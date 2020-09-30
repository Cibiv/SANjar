library(data.table)

#' Compute rank-size table from a table with observed lineage sizes
#' 
#' @param subset a `data.table` containing columns `sid` (sample id) and `lsize` (lineage size)
#' @return a `data.table`  with columns `sid` (sample id), `rank` (lineage size rank) and `size` (lineage size)
rank_size <- function(subset, by=c("sample", "sid", "day")) {
  subset[, {
    .SD[order(lsize, decreasing=TRUE)][, list(rank=1:.N, size=lsize)]
  }, by=mget(by)]
}

#' Align rank-size curves
align_rank_size <- function(ranked) {
  rank.max.all <- max(ranked$rank)
  size.min.all <- min(ranked$size)
  ranked[, {
    rank.scale <- rank.max.all / (max(rank)-1)
    size.scale <- size.min.all / min(size)
    sm <- min(size)
    list(rank, size,
         rank.aligned=(rank-1) * rank.scale + 1,
         size.aligned=size * size.scale,
         rank.scale, size.scale)
  }, keyby="sid"]
}

fit_pareto <- function(sizes) {
  # https://en.wikipedia.org/wiki/Pareto_distribution#Estimation_of_parameters
  xm <- min(sizes)
  alpha <- length(sizes) / sum(log(sizes/xm))
  return(list(xm=xm, alpha=alpha))
}

#' Fit the powerlaw model to rank-size data
#' 
#' @param ranked  a `data.table`  with columns `sid` (sample id), `rank` (lineage size rank) and `size` (lineage size)
#' @param alpha the Pareto index alpha to use. By default, alpha is estimated for each sample separately.
#' @param r.large a two-component vector which specifies which range of ranks corresponds to "large" lineages 
#' @return a `data.table` with columns `sid` (sample id), `alpha` (Pareto index), `k` (Zipf exponent),
#'         `d` (Zipf intersect), `r` (smallest Zipf-goverened rank), `s` (largest Zipf-goverend lineage size).
fit_powerlaw_model <- function(ranked, alpha=NA, r.large=c(15, 100), by=c("sid", "sample", "day")) {
  r <- ranked[, .SD[order(rank), {
    stopifnot(all(rank == 1:length(rank)))
    stopifnot(!is.unsorted(rev(size)))
    # 1. Fit powerlaw
    #
    # If no global alpha was specified, compute a sample-specific Pareto index alpha
    a <- if (!is.na(alpha)) alpha else fit_pareto(size)$alpha
    # In log-log-space, the powerlaw model for the rank-size relationship is
    #   log(size) = k * log(rank) + d
    # where k = -1/alpha if the sizes follow a Pareto distribution with index alpha,
    # and d is found such that the powerlaw intersects the rank-size curve in its
    # rightmost (and thus minimal) point, i.e. in the point
    #   p = (max(rank), min(size))
    k <- -1/a
    r.max <- length(rank)
    s.min <- size[r.max]
    d <- log10(s.min) + log10(r.max)*(1/a)
    #
    # 2. Fit another powerlaw for large lineages (i.e. small ranks)
    m.small <- lm(log10(size) ~ log10(rank), .SD[(r.large[1] <= rank) & (rank <= r.large[2])])
    #
    # 3. Intersect the two powerlaws (lines in log-log-space) to find the knee (r,s),
    # i.e. the rank r from which onward the curve shows a powerlaw tail and the correspoinding
    # size s.
    r <- round(10^((d - m.small$coefficients[1]) / (m.small$coefficients[2] - k)))
    s <- size[r]
    #
    # 4. Return model fit
    list(pareto.alpha=a, zipf.k=k, zipf.d=d, zipf.rank.min=r, zipf.size.max=s)
  }], by=mget(by)]
  setkeyv(r, by[1])
  r
}

# Return true if x and y have the same value, including both being NA
is.samevalue <- function(x, y) {
  if (is.atomic(x) && is.atomic(y)) {
    xv <- !is.na(x)
    yv <- !is.na(y)
    return(!((length(x) != length(y)) || any(xv != yv) || any(x[xv] != y[xv])))
  }
  else
    stop("non-atomic values are currently unsupported")
}
