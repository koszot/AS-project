library(tidyr)
library(readr)
library(dplyr)
library(stringr)

setwd("~/Desktop/alternative_splicing/umaydis/")

headers <- read_tsv("umaydis_transcripts.fasta.transdecoder_dir/headers.cds", col_names = F)

longest_orf <- function(headers_FUN) {
  headers.longest <- data.frame()
  ID <- unique(headers_FUN$X2)
  for (x in 1:length(ID)) {
    headers_temp <- headers_FUN[headers_FUN$X2 %in% ID[x],]
    headers.longest <- rbind(headers.longest, headers_temp[1,])
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(headers.longest)
}

headers_longest <- longest_orf(headers)

filter <- headers_longest$X1

write_lines(filter, "umaydis_transcripts.fasta.transdecoder_dir/filter")
