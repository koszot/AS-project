library(ASpli)

setwd("~/Desktop/working/scommune")

TxDb <- makeTxDbFromGFF(file="genome/RRPM_transcripts.fixed.gtf", format="gtf")

features <- binGenome(TxDb)
