library(data.table)

# Load cell count data
message("*** Loading celltypes.tab")
CELLTYPES <- data.table(read.table("celltypes.tab", header=TRUE, colClasses=c(
  day="integer", antibody="factor", rep="integer", percent="numeric")))[, list(
    day, source=as.factor("facs"), antibody, rep, percent
) ]

message("*** Saving as celltypes.rd")
save(CELLTYPES, file="celltypes.rd", version=2)
