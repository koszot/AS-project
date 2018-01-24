library(ASpli)

setwd("~/Desktop/alternative_splicing/ccinerea_AmutBmut/")

TxDb <- makeTxDbFromGFF(file="ccinerea_AmutBmut_genome/ccinerea_AmutBmut_AS_annotation.gtf", format="gtf")

features <- binGenome(TxDb)
