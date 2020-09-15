library(data.table)

message("*** Importing lineagetracing_H9_only.tsv")
D <- data.table(read.table("lineagetracing_H9_only.tsv", header=TRUE,
                           sep="\t", quote="\"", stringsAsFactors=FALSE))

message("*** Converting into the data.table LT47")
LT47 <- D[, {
  #    m <- regexec("^(\\d+)-H9-day(\\d\\d)-(\\d)+(-PCR\\d)?$", Sample)
  m <- regexec("^H9-day(\\d\\d)-(\\d+)$", Sample)
  stopifnot(all(sapply(m, `[`, 1) > 0))
  sid <- as.integer(sapply(regmatches(Sample, m), `[`, 3))
  day <- as.integer(sapply(regmatches(Sample, m), `[`, 2))
  list(sample=as.factor(Sample), sid=sid, day=day, lsize=as.integer(nLineages))
}]
last_sid <- 0L
LT47[, sid := { last_sid <<- last_sid + 1L }, by=c("day", "sid")]
setkey(LT47, sid)

message("*** Saving as lt47.rd")
save(LT47, file="lt47.rd", version=2)
