#' Compute rank-size table from a table with observed lineage sizes
#' 
#' @param data a `data.table`
#' @param groups the name(s) of a set of columns of `data` which define the group within which sizes are ranked
#' @param size the name of a column of `data` containing the lineage sizes
#' @return a `data.table` containing the grouping columns, `rank` (lineage size rank) and `size` (lineage size)
rank_size <- function(data, groups=c("day", "sid")) {
  r <- data[, {
    .SD[order(size, decreasing=TRUE)][, list(rank=1:.N, size)]
  }, by=groups]
  setkeyv(r, c(groups, "rank"))
  r
}

fit_pareto <- function(sizes) {
  sizes[, {
    # https://en.wikipedia.org/wiki/Pareto_distribution#Estimation_of_parameters
    xm <- min(size)
    alpha <- length(size) / sum(log(size/xm))
    list(xm=xm, alpha=alpha)
  }, keyby=.(sid, day)]
}

#' Fit the powerlaw model to rank-size data
#' 
#' @param ranked  a `data.table`  with columns `sid` (sample id), `rank` (lineage size rank) and `size` (lineage size)
#' @param alpha the Pareto index alpha to use. By default, alpha is estimated for each sample separately.
#' @param r.nonzipf a two-component vector which specifies which range of ranks corresponds to "non-Zipfian" lineages
#' @return a `data.table` with columns `sid` (sample id), `alpha` (Pareto index), `k` (Zipf exponent),
#'         `d` (Zipf intersect), `r` (smallest Zipf-goverened rank), `s` (largest Zipf-goverend lineage size).
fit_powerlaw_model <- function(ranked, alpha=NA, r.nonzipf=c(NA_integer_, NA_integer_)) {
  r <- ranked[, .SD[order(rank), {
    n <- length(rank)
    stopifnot(all(rank == 1:n))
    stopifnot(!is.unsorted(rev(size)))
    stopifnot(length(r.nonzipf)==2)
    log_rank=log10(rank)
    log_size=log10(size)
    # 1. Fit powerlaw
    #
    # If no global alpha was specified, compute a sample-specific Pareto index alpha
    a <- if (!is.na(alpha)) alpha else fit_pareto(.SD[, list(sid, day, size)])$alpha
    # In log-log-space, the powerlaw model for the rank-size relationship is
    #   log(size) = k * log(rank) + d
    # where k = -1/alpha if the sizes follow a Pareto distribution with index alpha,
    # and d is found such that the powerlaw intersects the rank-size curve in its
    # rightmost (and thus minimal) point, i.e. in the point
    #   p = (max(rank), min(size))
    k <- -1/a
    d <- log_size[n] - k * log_rank[n]
    #
    # 3. Auto-detect the rank range corresponding to "large" lineages
    #
    # By default, the range includes all lineages up to (but not including)
    # the first lineages whose size exceeds the powerlaw prediction, but at
    # most the largest sqrt(n) lineages.
    r.nonzipf.min <- if (is.na(r.nonzipf[1])) 1 else r.nonzipf[1]
    r.nonzipf.max <- if (is.na(r.nonzipf[2])) {
      log_size.pred <- k*log_rank + d
      r.below.pred <- min(match(FALSE, log_size.pred > log_size) - 1, n, na.rm=T)
      round(r.below.pred^0.5)
    } else r.nonzipf[2]
    r.nonzip.max <- max(r.nonzipf.min, r.nonzipf.max)
    #
    # 2. Fit another powerlaw for large lineages (i.e. small ranks)
    m.nonzipf <- if (r.nonzipf.min < r.nonzip.max) {
      lm(log10(size) ~ log10(rank), .SD[(r.nonzipf.min <= rank) & (rank <= r.nonzipf.max)])
    } else NULL
    #
    # 3. Intersect the two powerlaws (lines in log-log-space) to find the knee (r,s),
    # i.e. the rank r from which onward the curve shows a powerlaw tail and the correspoinding
    # size s.
    r <- if (!is.null(m.nonzipf)) {
      round(10^((d - m.nonzipf$coefficients[1]) / (m.nonzipf$coefficients[2] - k)))
    } else 1
    s <- size[r]
    #
    # 4. Return model fit
    list(pareto.alpha=a, zipf.k=k, zipf.d=d, zipf.rank.min=r, zipf.size.max=s,
         nonzipf.k=if (!is.null(m.nonzipf)) m.nonzipf$coefficients[2] else NA_real_,
         nonzipf.d=if (!is.null(m.nonzipf)) m.nonzipf$coefficients[1] else NA_real_)
  }], by=.(sid, day)]
  setkey(r, sid, day)
  r
}
