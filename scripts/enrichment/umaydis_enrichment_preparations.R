library(dplyr)
library(readr)
library(stringr)
library(tidyr)

##############################################
#     GO-kat tartalmazó fájl elkészítése     #
##############################################

raw_interpro <- read_lines("~/Desktop/alternative_splicing/umaydis/p3_t237631_Ust_maydi_v2GB.prot.tsv")

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

write_tsv(GO_final_merged, "~/Desktop/alternative_splicing/umaydis/umaydis_GO_dict.tsv", col_names = F)
