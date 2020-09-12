library(data.table)

if (!exists("LT47") && !file.exists("lt47.rd")) {
  message("*** Importing lineagetracing_H9_only.tsv")
  D <- data.table(read.table("lineagetracing_H9_only.tsv", header=TRUE,
                             sep="\t", quote="\"", stringsAsFactors=FALSE))
  LT47 <- D[, {
    #    m <- regexec("^(\\d+)-H9-day(\\d\\d)-(\\d)+(-PCR\\d)?$", Sample)
    m <- regexec("^H9-day(\\d\\d)-(\\d+)$", Sample)
    stopifnot(all(sapply(m, `[`, 1) > 0))
    sid <- as.integer(sapply(regmatches(Sample, m), `[`, 3))
    day <- as.integer(sapply(regmatches(Sample, m), `[`, 2))
    list(sample=Sample, sid=sid, day=day, lsize=nLineages)
  }]
  last_sid <- 0
  LT47[, sid := { last_sid <<- last_sid + 1 }, by=c("day", "sid")]
  setkey(LT47, sid)
  save(LT47, file="lt47.rd", version=2)
} else if (!exists("LT47") && file.exists("lt47.rd")) {
  message("*** Loading lt47.rd")
  load("lt47.rd")
}
