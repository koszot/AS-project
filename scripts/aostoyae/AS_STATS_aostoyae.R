library(ASpli)

setwd("~/Desktop/working/aostoyae")

TxDb <- makeTxDbFromGFF(file="genome/RRPM_transcripts.fixed.gtf", format="gtf")

features <- binGenome(TxDb)
