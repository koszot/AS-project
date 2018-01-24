library(ASpli)

setwd("~/Desktop/alternative_splicing/pchrysosporium/")

TxDb <- makeTxDbFromGFF(file="pchrysosporium_genome/pchrysosporium_AS_annotation.gtf", format="gtf")

features <- binGenome(TxDb)
