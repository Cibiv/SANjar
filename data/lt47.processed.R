library(data.table)
library(gwpcR)

# Load lineage size data
message("*** Loading LT47 from lt47.rd")
load("lt47.rd")

message("*** Removing low-quality samples")
LT47 <- LT47[!(sid %in% c(15, 16))]

# Compute number of observed lineages at each time point
message("*** Computing LT47.NLINEAGES")
LT47.NLINEAGES <- LT47[, list(nlineages=as.integer(sum(lsize > 0))), by=c("day", "sid")] 

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
save(LT47, LT47.NLINEAGES, LT47.GWPCR.PARAMETERS, file="lt47.processed.rd", version=2)
