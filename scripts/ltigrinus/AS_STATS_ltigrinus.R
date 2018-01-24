library(ASpli)

setwd("~/Desktop/alternative_splicing/ltigrinus/")

TxDb <- makeTxDbFromGFF(file="ltigrinus_genome/ltigrinus_AS_annotation.gtf", format="gtf")

features <- binGenome(TxDb)
