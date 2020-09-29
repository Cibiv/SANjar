library(data.table)
library(gwpcR)

#' Compute rank-size table from a table with observed lineage sizes
#' 
#' @param subset a `data.table` containing columns `sid` (sample id) and `lsize` (lineage size)
#' @return a `data.table`  with columns `sid` (sample id), `rank` (lineage size rank) and `size` (lineage size)
rank_size <- function(subset) {
  subset[, {
    .SD[order(lsize, decreasing=TRUE)][, list(rank=1:.N, size=lsize)]
  }, by=c("sample", "sid", "day")]
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
fit_powerlaw_model <- function(ranked, alpha=NA, r.large=c(15, 100)) {
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
  }], by=c("sample", "sid", "day")]
  setkey(r, "sid")
  r
}

# Load lineage size data
message("*** Loading LT47 from lt47.rd")
load("lt47.rd")

message("*** Removing low-quality samples")
LT47 <- LT47[!(sid %in% c(15, 16))]

# Compute number of observed lineages at each time point
message("*** Computing LT47.NLINEAGES")
LT47.NLINEAGES <- LT47[, list(nlineages=as.integer(sum(lsize > 0))), by=c("day", "sid")] 

message("*** Computing LT47.RANKSIZE")
LT47.RANKSIZE <- rank_size(LT47)

message("*** Computing LT47.POWERLAW")
# Compute powerlaw fits
LT47.LIMIT.ALPHA.STARTDAY <- 11
LT47.LIMIT.ALPHA <- signif(LT47[, fit_pareto(lsize), by=c("day", "sid")][
  day >= LT47.LIMIT.ALPHA.STARTDAY, mean(alpha)], digits=2)
LT47.POWERLAW <- fit_powerlaw_model(LT47.RANKSIZE, alpha=LT47.LIMIT.ALPHA)
LT47.POWERLAW[, zipf.k := signif(zipf.k, digits=3) ]
LT47.POWERLAW[, zipf.d := signif(zipf.d, digits=3) ]

# Estimate PCR and sequencing parameters
message("*** Estimating PCR and sequencing parameters")
LT47.GWPCR.PARAMETERS <- LT47[,
                         .SD[, list(
                           library_size=sum(lsize),
                           phantom_threshold=min(lsize)
                         ), by="sid"][, list(
                           library_size=signif(round(median(library_size)), 2),
                           phantom_threshold=round(median(phantom_threshold))
                         )]
                         , by="day"]
LT47.GWPCR.PARAMETERS[, pcr_efficiency := {
  LT47[day==0][, {
    gwpcrpois.est(lsize, threshold=min(lsize), molecules=1)
  }, by="sid"][, list(
    pcr_efficiency=signif(median(efficiency), 2)
  )]
}]

message("*** Saving LT47, LT47.NLINEAGES, LT47.GWPCR.PARAMETERS to lt47.processed.rd")
save(LT47, LT47.NLINEAGES, LT47.RANKSIZE, LT47.LIMIT.ALPHA, LT47.POWERLAW, LT47.GWPCR.PARAMETERS,
     file="lt47.processed.rd", version=2)
