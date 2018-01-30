library(topGO)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

#######################################################
#     LENTINUS GO-kat tartalmazó fájl elkészítése     #
#######################################################

raw_interpro <- read_lines("~/Desktop/alternative_splicing/ltigrinus/Sisbr1_GeneCatalog_proteins_20130805.aa.fasta.tsv")

IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)

colnames(IPR) <- "all"

IPR <- IPR %>%
  separate(all, c("transcript_id", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t")

GO <- IPR %>%
  filter(GO != "") %>%
  select(transcript_id, GO) %>%
  distinct(GO, transcript_id, GO)

TID_GO <- unique(GO$transcript_id)
GO_final_merged <- tibble()

for (x in 1:length(TID_GO)) {
  GO_temp <- GO[GO$transcript_id %in% TID_GO[x],]
  GO_merged <- vector()
  for (y in 1:length(GO_temp$transcript_id)) {
    GO_merged[y] <- strsplit(as.character(GO_temp[y,2]), "\\|")
  }
  GO_final_merged[x,1] <- TID_GO[x]
  GO_final_merged[x,2] <- paste(unique(unlist(GO_merged)), collapse = ", ")
  pb <- txtProgressBar(min = 1, max = length(TID_GO), style = 3)                      # progress bar
  setTxtProgressBar(pb, x, title = NULL, label = NULL)
}

colnames(GO_final_merged) <- c("fasta_transcript_id", "GO_annotations")

dictionary <- read_tsv("~/Desktop/alternative_splicing/ltigrinus/ltigrinus_dictionary.tsv", col_names = T, cols(fasta_transcript_id = col_character(),gene_id = col_character(),annotation_transcript_id = col_character()))

GO_final_formated <- left_join(GO_final_merged, dictionary, by = "fasta_transcript_id")

GO_final_formated <- GO_final_formated[,c(4,1:3)] %>%
  select(annotation_transcript_id, GO_annotations)

colnames(GO_final_formated) <- c("transcript_id", "GO_annotations")

write_tsv(GO_final_formated, "~/Desktop/alternative_splicing/ltigrinus/ltigrinus_GO_dict.tsv", col_names = F)

#######################################################
#     LENTINUS GO-kat tartalmazó fájl elkészítése     #
#######################################################

raw_interpro <- read_lines("~/Desktop/alternative_splicing/ltigrinus/ltigrinus_proteins_all.fasta.tsv")

IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)

colnames(IPR) <- "all"

IPR <- IPR %>%
  separate(all, c("transcript_id", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t")

GO <- IPR %>%
  filter(GO != "") %>%
  select(transcript_id, GO) %>%
  distinct(GO, transcript_id, GO)

TID_GO <- unique(GO$transcript_id)
GO_final_merged <- tibble()

for (x in 1:length(TID_GO)) {
  GO_temp <- GO[GO$transcript_id %in% TID_GO[x],]
  GO_merged <- vector()
  for (y in 1:length(GO_temp$transcript_id)) {
    GO_merged[y] <- strsplit(as.character(GO_temp[y,2]), "\\|")
  }
  GO_final_merged[x,1] <- TID_GO[x]
  GO_final_merged[x,2] <- paste(unique(unlist(GO_merged)), collapse = ", ")
  pb <- txtProgressBar(min = 1, max = length(TID_GO), style = 3)                      # progress bar
  setTxtProgressBar(pb, x, title = NULL, label = NULL)
}

colnames(GO_final_merged) <- c("transcript_id", "GO_annotations")

write_tsv(GO_final_merged, "~/Desktop/alternative_splicing/ltigrinus/ltigrinus_transcript_GO_dict.tsv", col_names = F)
