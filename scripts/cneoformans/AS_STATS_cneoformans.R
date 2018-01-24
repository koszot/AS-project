library(ASpli)

setwd("~/Desktop/alternative_splicing/cneoformans/")

TxDb <- makeTxDbFromGFF(file="cneoformans_genome/cneoformans_AS_annotation.gtf", format="gtf")

features <- binGenome(TxDb)
