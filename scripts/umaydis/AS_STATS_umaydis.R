library(ASpli)

setwd("~/Desktop/alternative_splicing/umaydis/")

TxDb <- makeTxDbFromGFF(file="umaydis_genome/umaydis_AS_annotation.gtf", format="gtf")

features <- binGenome(TxDb)
