library(ASpli)

setwd("~/Desktop/alternative_splicing/rmellea/")

TxDb <- makeTxDbFromGFF(file="rmellea_genome/rmellea_AS_annotation.gtf", format="gtf")

features <- binGenome(TxDb)
