library(ASpli)

setwd("~/Desktop/alternative_splicing/aampla/")

TxDb <- makeTxDbFromGFF(file="aampla_genome/aampla_AS_annotation.gtf", format="gtf")

features <- binGenome(TxDb)
