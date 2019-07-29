
library(data.table)

# Target n
targetsize <- 100000

# File names
inputFile <- "Daten/Adrian Synthetic Perfect.txt"
inputFileB <- "Daten/Adrian Synthetic 10 percent.txt"

# Read input
clearTextA <- fread(inputFile, colClasses = "character", stringsAsFactors = FALSE)
clearTextB <- fread(inputFileB, colClasses = "character", stringsAsFactors = FALSE)


## Deduplicate by ID
clearTextA <- clearTextA[!duplicated(clearTextA$`Group Id`),]
clearTextB <- clearTextB[!duplicated(clearTextB$`Group Id`),]

# Sample for size
set.seed(42)
newfileA <- clearTextA[sample(nrow(clearTextA), targetsize,replace=FALSE),]

set.seed(42)
newfileB <- clearTextB[sample(nrow(clearTextA), targetsize,replace=FALSE),]

# Write new file
write.table(newfileA, "Daten/Adrian_Perfect_100k.csv", row.names = FALSE, sep="\t")
write.table(newfileB, "Daten/Adrian_10percent_100k.csv", row.names = FALSE, sep="\t")
