library(data.table)

#' Compute rank-size table from a table with observed lineage sizes
#' 
#' @param subset a `data.table` containing columns `sid` (sample id) and `lsize` (lineage size)
#' @return a `data.table`  with columns `sid` (sample id), `rank` (lineage size rank) and `size` (lineage size)
rank_size <- function(subset) {
  r <- subset[, {
    .SD[order(lsize, decreasing=TRUE)][, list(rank=1:.N, size=lsize)]
  }, by=.(sid, day)]
  setkey(r, sid, day, rank)
  r
}

#' Align rank-size curves
align_rank_size <- function(ranked, rank.max=NULL, size.min=NULL) {
  r <- ranked[, {
    rank.max.day <- if (is.null(rank.max)) max(rank) else rank.max
    size.min.day <- if (is.null(size.min)) min(size) else size.min
    .SD[, {
      rank.scale <- rank.max.day / (max(rank)-1)
      size.scale <- size.min.day / min(size)
      sm <- min(size)
      list(rank, size,
           rank.aligned=(rank-1) * rank.scale + 1,
           size.aligned=size * size.scale,
           rank.scale, size.scale)
    }, by=sid]
  }, by=day]
  setkey(r, sid, day, rank)
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
#' @param r.large a two-component vector which specifies which range of ranks corresponds to "large" lineages 
#' @return a `data.table` with columns `sid` (sample id), `alpha` (Pareto index), `k` (Zipf exponent),
#'         `d` (Zipf intersect), `r` (smallest Zipf-goverened rank), `s` (largest Zipf-goverend lineage size).
fit_powerlaw_model <- function(ranked, alpha=NA, r.large=c(NA_integer_, NA_integer_)) {
  r <- ranked[, .SD[order(rank), {
    n <- length(rank)
    stopifnot(all(rank == 1:n))
    stopifnot(!is.unsorted(rev(size)))
    stopifnot(length(r.large)==2)
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
    r.large.min <- if (is.na(r.large[1])) 1 else r.large[1]
    r.large.max <- if (is.na(r.large[2])) {
      log_size.pred <- k*log_rank + d
      r.below.pred <- min(match(FALSE, log_size.pred > log_size) - 1, n, na.rm=T)
      round(r.below.pred^0.5)
    } else r.large[2]
    if (r.large.min >= r.large.max) stop("r.large auto-detection failed")
    #
    # 2. Fit another powerlaw for large lineages (i.e. small ranks)
    m.small <- lm(log10(size) ~ log10(rank), .SD[(r.large.min <= rank) & (rank <= r.large.max)])
    #
    # 3. Intersect the two powerlaws (lines in log-log-space) to find the knee (r,s),
    # i.e. the rank r from which onward the curve shows a powerlaw tail and the correspoinding
    # size s.
    r <- round(10^((d - m.small$coefficients[1]) / (m.small$coefficients[2] - k)))
    s <- size[r]
    #
    # 4. Return model fit
    list(pareto.alpha=a, zipf.k=k, zipf.d=d, zipf.rank.min=r, zipf.size.max=s,
         large.k=m.small$coefficients[2], large.d=m.small$coefficients[1])
  }], by=.(sid, day)]
  setkey(r, sid, day)
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
