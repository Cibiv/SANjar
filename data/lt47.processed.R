source("../utilities.R")
library(gwpcR)

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
LT47.POWERLAW <- fit_powerlaw_model(LT47.RANKSIZE, alpha=LT47.LIMIT.ALPHA, by=c("sid", "sample", "day"))
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

message("*** Saving LT47, LT47.NLINEAGES, LT47.RANKSIZE, LT47.LIMIT.ALPHA.STARTDAY, LT47.LIMIT.ALPHA, LT47.POWERLAW, LT47.GWPCR.PARAMETERS to lt47.processed.rd")
save(LT47, LT47.NLINEAGES, LT47.RANKSIZE, LT47.LIMIT.ALPHA.STARTDAY, LT47.LIMIT.ALPHA, LT47.POWERLAW, LT47.GWPCR.PARAMETERS,
     file="lt47.processed.rd", version=2)
