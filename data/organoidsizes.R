library(data.table)

# Load cell count data
message("*** Loading organoidsizes.tab")
ORGANOIDSIZES <- data.table(read.table("organoidsizes.tab", header=TRUE, colClasses=c(
  day="integer", source="factor", rep="factor", count="integer", hq="logical")))

message("*** Saving as organoidsizes.rd")
save(ORGANOIDSIZES, file="organoidsizes.rd", version=2)
