library(dplyr)
library(readr)
library(tidyr)

setwd("~/Desktop/alternative_splicing/pchrysosporium/")

transcript_lengths <- read_tsv("pchrysosporium_genome/gene_length", col_names = F)

transcripts <- read_tsv("pchrysosporium_genome/pchrysosporium_onlyexon.gtf", 
                  col_names = c("chr", "maker", "type", "start", "end", "att1", "strand", "att2", "attributes"))  

transcripts <- transcripts %>%
  separate(attributes, c("transcriptID_label", "transcriptID",
                         "geneID_label", "geneID"), sep = " ")   

TID <- unique(transcripts$transcriptID)
introns <- tibble()

for (x in 1:length(TID)) {
  transcripts_temp <- transcripts[transcripts$transcriptID %in% TID[x],] 
  introns_temp <- tibble()
  for (y in 1:length(transcripts_temp$chr)) {
    introns_temp[y,1] <- as.integer(transcripts_temp[(y + 1),4]) - as.integer(transcripts_temp[y, 5])
  }
  introns <- rbind(introns, introns_temp)
}

introns <- filter(introns, !is.na(V1))

message(c("A maximális transzkriptméret: ", max(transcript_lengths$X2)))
message(c("Intron min: ", min(introns$V1), " Intron max: ", max(introns$V1)))
