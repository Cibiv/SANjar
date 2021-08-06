library(data.table)
library(gwpcR)
devtools::load_all(".")

message("*** Importing FACS-based organoid sizes from h9_organoidsizes_facs.tsv")
h9_organoidsizes_facs <- data.table(read.table("data-raw/h9_organoidsizes_facs.tsv", header=TRUE, sep="\t"))
h9_organoidsizes_facs <- h9_organoidsizes_facs[, list(
  day, source="facs", rep, cells=count
)]

message("*** Importing area-based organoid sizes from h9_organoidsizes_area.tsv")
h9_organoidsizes_area <- data.table(read.table("data-raw/h9_organoidsizes_area.tsv", header=TRUE, sep="\t",
                                               check.names=FALSE))
h9_organoidsizes_area <- h9_organoidsizes_area[, list(
  day=as.integer(substr(`Day`, 4, 6)),
  source="area",
  rep=paste0("b", trimws(substr(`Batch`, 6, 8))),
  volume=`Volume`
)]
h9_organoidsizes_area[, rep := paste0(rep, "-", 1:.N), by=.(rep)]

message("*** Combinding FACS- and area-based organoid sizes")
# Combine FACS- and area-based data based on the (interpolated) FACS count for day 10.
# The FACS-based measurements underestimate the organoid size from day 25 onwards because
# they fail to include dead cells, and are hence excluded.
f9 <- mean(h9_organoidsizes_facs[day==9, log(cells)])
f13 <- mean(h9_organoidsizes_facs[day==13, log(cells)])
a10 <- mean(h9_organoidsizes_area[day==10, log(volume)])
s <- exp((3/4)*f9+(1/4)*f13-a10)
h9_organoidsizes <- rbind(
  h9_organoidsizes_facs,
  h9_organoidsizes_area[, list(day, source, rep, cells=round(s*volume))]
)[((source=="facs") & (day < 25)) | (source=="area")]
setkey(h9_organoidsizes, day, rep)

message("*** Importing lineage sizes from data-raw/lt47_h9_only.tsv")
lt47_ls <- data.table(read.table("data-raw/lt47_h9_only.tsv", header=TRUE,
                                 sep="\t", quote="\"", stringsAsFactors=FALSE))

message("*** Converting lineage sizes into data.table")
lt47_ls <- lt47_ls[, {
  #    m <- regexec("^(\\d+)-H9-day(\\d\\d)-(\\d)+(-PCR\\d)?$", Sample)
  m <- regexec("^H9-day(\\d\\d)-(\\d+)$", Sample)
  stopifnot(all(sapply(m, `[`, 1) > 0))
  sid <- as.integer(sapply(regmatches(Sample, m), `[`, 3))
  day <- as.integer(sapply(regmatches(Sample, m), `[`, 2))
  list(sample=as.factor(Sample), sid=sid, day=day, reads=as.integer(nLineages))
}]
last_sid <- 0L
lt47_ls[, sid := { last_sid <<- last_sid + 1L }, by=c("day", "sid")]
setkey(lt47_ls, day, sid)

message("*** Creating LTData instance")
lt47 <- SANjar::LTData(organoidsizes=h9_organoidsizes,
                       lineagesizes=lt47_ls[!(sid %in% c(15, 16))],
                       unit="reads")

message("*** Estimating PCR and sequencing parameters")
lt47 <- SANjar::estimate_sequencing_parameters(lt47)

message("*** Translating lineage sizes from #reads to #cells")
lt47 <- SANjar::estimate_reads_per_cell(lt47, method="singleton_mode")
lt47 <- SANjar::absolute_lineage_sizes(lt47)

message("*** Saving lt47 to data/lt47.RData")
save(lt47, file="data/lt47.RData", version=2)
