library(readr)
library(dplyr)
library(tidyr)
library(stringr)

setwd("~")

isoforms_tracking <- read_tsv("Desktop/alternative_splicing/ltigrinus/ltigrinus_RRPM_CUFFout/isoforms.fpkm_tracking")

annotation <- read_tsv("Desktop/alternative_splicing/ltigrinus/ltigrinus_genome/ltigrinus_AS_annotation.gtf", 
                       col_names = c("chr", "maker", 
                                     "type", "start", 
                                     "end", "att1", 
                                     "strand", "att2", 
                                     "attributes"))

annotation <- annotation %>%
  separate(attributes, c("transcriptID_label", "transcriptID",
                         "geneID_label", "geneID"), sep = " ") 

annotation$transcriptID <- annotation$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")    

annotation$geneID <- annotation$geneID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")    

transcriptID <- unique(annotation$transcriptID)

filtered_isoforms_tracking <- data.frame()

for (x in transcriptID) {
  temp_isoforms_tracking <- isoforms_tracking %>%
    filter(tracking_id == x)
  filtered_isoforms_tracking  <- rbind(filtered_isoforms_tracking , temp_isoforms_tracking)
}

filtered_isoforms_tracking <- select(filtered_isoforms_tracking, tracking_id, gene_id)

final_tracking <- data.frame()
transcriptID <- unique(filtered_isoforms_tracking$tracking_id)

for (x in transcriptID) {
  temp_annotation <- annotation[annotation$transcriptID %in% x,]
  temp_filtered_isoforms_tracking <- filtered_isoforms_tracking[filtered_isoforms_tracking$tracking_id %in% x,]
  temp_filtered_isoforms_tracking[1,3] <- temp_annotation[1,12]
  final_tracking <- rbind(final_tracking , temp_filtered_isoforms_tracking)
}

final_tracking <- select(final_tracking, gene_id, geneID)

final_tracking_2 <- distinct(final_tracking, gene_id, geneID)

candidates <- final_tracking_2 %>%
  group_by(geneID) %>%
  summarise(sum = length(geneID)) %>%
  filter(sum > 1) %>%
  select(geneID)

message(candidates)
