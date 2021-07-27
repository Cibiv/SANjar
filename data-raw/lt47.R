library(data.table)
library(gwpcR)

message("*** Importing FACS-based organoid sizes from h9_organoidsizes_facs.tsv")
h9_organoidsizes_facs <- data.table(read.table("data-raw/h9_organoidsizes_facs.tsv", header=TRUE, sep="\t"))
h9_organoidsizes_facs <- h9_organoidsizes_facs[, list(
  day, source="facs", sid=rep, cells=count
)]

message("*** Importing area-based organoid sizes from h9_organoidsizes_area.tsv")
h9_organoidsizes_area <- data.table(read.table("data-raw/h9_organoidsizes_area.tsv", header=TRUE, sep="\t",
                                               check.names=FALSE))
h9_organoidsizes_area <- h9_organoidsizes_area[, list(
  day=as.integer(substr(`Day`, 4, 6)),
  source="area",
  sid=paste0("b", trimws(substr(`Batch`, 6, 8))),
  volume=`Volume`
)]
h9_organoidsizes_area[, sid := paste0(sid, "-", 1:.N), by=.(sid)]

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
  h9_organoidsizes_area[, list(day, source, sid, cells=round(s*volume))]
)[((source=="facs") & (day < 25)) | (source=="area")]
setkey(h9_organoidsizes, day, sid)

message("*** Saving organoid sizes to data/h9_organoidsizes.RData")
save(h9_organoidsizes, file="data/h9_organoidsizes.RData", version=2)

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

message("*** Removing low-quality samples")
lt47_ls <- lt47_ls[!(sid %in% c(15, 16))]

message("*** Translating lineage sizes from #reads to #cells")
reads_per_cell <- lt47_ls[, {
  f <- density(log10(reads))
  list(reads.per.cell=10^(f$x[which.max(f$y)]))
}, keyby=.(day, sid)]

lt47_ls_cells <- lt47_ls[reads_per_cell, {
  list(day=day,
       cells=reads/reads.per.cell)
}, keyby=.EACHI]

message("*** Estimating PCR and sequencing parameters")
sequencing <- lt47_ls[, list(
  library_size=sum(reads),
  phantom_threshold=min(reads)
), keyby=.(day, sid)]
# Compute PCR efficiency from day-0 samples
sequencing[, pcr_efficiency := {
  lt47_ls[day==0][, {
    gwpcrpois.est(reads, threshold=min(reads), molecules=1)
  }, by="sid"][, list(
    pcr_efficiency=signif(median(efficiency), 2)
  )]
}]

message("*** Constructing lt47 dataset object")
lt47 <- structure(list(organoidsizes=h9_organoidsizes,
                       lineagesizes=lt47_ls_cells,
                       sequencing=sequencing,
                       unit=as.name("cells"),
                       groups=character()),
                  class="LTData")

message("*** Saving lt47 to data/lt47.RData")
save(lt47, file="data/lt47.RData", version=2)
