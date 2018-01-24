#!/usr/bin/Rscript

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

#######################################
#    ADATOK BETÖLTÉSE, ÁTALAKÍTÁSA    #
#######################################

setwd("Desktop/working/aostoyae/")  # beállítjuk a working directoryt

genes <- read_tsv("genome/original_onlygene.gtf", 
                  col_names = c("chr", "maker", "type", "start", "end", "att1", "strand", "att2", "attributes"))     # beolvassuk az eredeti GTF fájlt

raw <- read_tsv("RRPM_CUFFout/transcripts.gtf", 
                col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes"))        # beolvassuk a Cufflinks output GTF fájlt

genes <- genes %>%
  separate(attributes, c("transcriptID_label", "transcriptID",
                         "geneID_label", "geneID"), sep = " ")                    # szétbontjuk megfelelő oszlopokra az eredeti csak a géneket tartalmazó GTF fájlt

raw <- raw %>%
  separate(attributes, c("geneID_label", "geneID", 
                         "transcriptID_label", "transcriptID", 
                         "FPKM_label", "FPKM", 
                         "frac_label", "frac",
                         "conf_lo_label", "conf_lo",
                         "conf_hi_label", "conf_hi",
                         "cov_label", "cov",
                         "full_read_support_label", "full_read_support"), sep = " ")      # szétbontjuk megfelelő oszlopokra a Cufflinks output GTF fájlt

stats <- list()                                                  # a stats list bevezetése
stats[[1]] <- paste0("Cufflinks output                  Gene count: ", 
                     nrow(unique(raw[,"geneID"])), 
                     "        Transcript count: ", 
                     nrow(unique(raw[,"transcriptID"])))         # raw adat statisztikái

raw$geneID <- raw$geneID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # geneID-k átalakítása

raw$transcriptID <- raw$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # transcriptID-k átalakítása

raw$FPKM <- raw$FPKM %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # FPKM-ek átalakítása    

raw$full_read_support <- raw$full_read_support %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # full_read_support-ok átalakítása    

################################
#     EXPRESSION FILTERING     #
################################

expression_filtering <- function(raw_FUN) 
{
  raw_IDfilter <- raw_FUN %>%                                                   # kiszedjük azon transcriptek ID-it amiknél az FPKM érték 0.0000000000
    filter(type == "transcript") %>%
    filter(FPKM == "0.0000000000") %>%
    select(transcriptID)
  IDs_zeroFPKM <- raw_IDfilter$transcriptID                                     # átalakítjuk karaktervektorrá
  for (x in 1:length(IDs_zeroFPKM)) {                                           # for ciklussal végigmegyünk az ID-kat tartalmazó listán
    temp_ID <- IDs_zeroFPKM[x]                                                  # kiemeljük a lista aktuális elemét
    raw_FUN <- raw_FUN %>%                                                      # kitöröljük az aktuális transcripthez tartozó sorokat a raw tibble-ből
      filter(transcriptID != temp_ID)
    pb <- txtProgressBar(min = 1, max = length(IDs_zeroFPKM), style = 3)        # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(raw_FUN)
}

expression_filtered <- expression_filtering(raw)

stats[[2]] <- paste0("Expression filtering output       Gene count: ", 
                     nrow(unique(expression_filtered[,"geneID"])), 
                     "        Transcript count: ", 
                     nrow(unique(expression_filtered[,"transcriptID"])))         # az exoression filtered statisztikái

#######################################
#     FULL READ SUPPORT FILTERING     #
#######################################

full_read_support_filtering <- function(expression_filtered_FUN) 
{
  expression_filtered_IDfilter <- expression_filtered_FUN %>%                   # kiszedjük azon transcriptek ID-it amiknél nincs full read coverage
    filter(type == "transcript") %>%
    filter(full_read_support == "no") %>%
    select(transcriptID)
  IDs_noFRS <- expression_filtered_IDfilter$transcriptID                        # átalakítjuk karaktervektorrá
  for (x in 1:length(IDs_noFRS)) {                                              # for ciklussal végigmegyünk az ID-kat tartalmazó listán
    temp_ID <- IDs_noFRS[x]                                                     # kiemeljük a lista aktuális elemét
    expression_filtered_FUN <- expression_filtered_FUN %>%                      # kitöröljük az aktuális transcripthez tartozó sorokat a raw tibble-ből
      filter(transcriptID != temp_ID)
    pb <- txtProgressBar(min = 1, max = length(IDs_noFRS), style = 3)           # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(expression_filtered_FUN)
}

read_support_filtered <- full_read_support_filtering(expression_filtered)

stats[[3]] <- paste0("Read support filtering output     Gene count: ", 
                     nrow(unique(read_support_filtered[,"geneID"])), 
                     "        Transcript count: ", 
                     nrow(unique(read_support_filtered[,"transcriptID"])))         # a full read support filtered statisztikái

###############################
#     CONTEXT RESTORATION     #
###############################

read_support_filtered <- read_support_filtered %>%
  select(chr, maker, type, start, end, att1, 
         strand, att2, transcriptID_label,
         transcriptID, geneID_label, geneID)                  # kitöröljük azokat a Cufflinks attributeokat amik már nemkellenek

context_restored <- read_support_filtered                     # létrehozzuk a context_restored tibble-t
context_restored$geneID <- context_restored$chr               # átmásoljuk a chr ID-kat a geneID-k helyére

genes$geneID <- genes$geneID %>%
  str_replace("\"", "") %>%
  str_replace("\";", "")                                      # geneID-k átalakítása a genes tibble-ben

# helyreállítjuk a Cufflinks output file kontextusát

context_restoration <- function(context_restored_FUN, genes_FUN) 
{
  context_restored_final <- data.frame()                                        # létrehozunk egy üres dataframe-et
  geneIDs <- unique(genes_FUN$geneID)                                           # kiszedjük a geneID-kat a genes tibbleből
  for (x in 1:length(geneIDs))                                                  # végiglépkedünk a geneIDs-en
  {
    temp_geneID <- geneIDs[x]                                                   # kiszedjük a geneID aktuális elemét
    raw_temp <- context_restored[context_restored$chr %in% temp_geneID,]        # egy temporary tibble-be kiemeljük az aktuális ID-hoz tartozó részt a context_restored-ből
    genes_temp <- genes[genes$geneID %in% temp_geneID,]                         # egy temporary tibble-be kiemeljük az aktuális ID-hoz tartozó részt a genes-ből
    if (length(raw_temp$chr) > 0)                                               # amennyiben a raw_temp tartalmaz sorokat csak akkor fusson le a restoration
    {
      raw_temp[,1] <- genes_temp[1,1]                                           # az eredeti kromoszómanév visszaállítása
      raw_temp[,4] <- raw_temp[,4] + as.integer(genes_temp[1,4]) - 1            # az eredeti start point visszaállítása
      raw_temp[,5] <- raw_temp[,5] + as.integer(genes_temp[1,4]) - 1            # az eredeti end point visszaállítása
      context_restored_final <- rbind(context_restored_final, raw_temp)         # összefűzzük a már készülő listával
    }
    pb <- txtProgressBar(min = 1, max = length(geneIDs), style = 3)             # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL) 
  }
  return(context_restored_final)
}

context_restored_final <- context_restoration(context_restored, genes)

stats[[4]] <- paste0("Context restoration output        Gene count: ", 
                     nrow(unique(context_restored_final[,"geneID"])), 
                     "        Transcript count: ", 
                     nrow(unique(context_restored_final[,"transcriptID"])))         # a context restoration statisztikái

############################
#     STRAND FILTERING     #
############################

strand_filtering <- function(context_restored_final_FUN, genes_FUN) 
{
  strand_filtered <- data.frame()                                                                               # létrehozunk egy üres dataframe-et
  geneIDs <- unique(genes$geneID)                                                                               # kiszedjük a geneID-kat a genes tibbleből
  for (x in 1:length(geneIDs))                                                                                  # végiglépkedünk a geneIDs-en
  {
    temp_geneID <- geneIDs[x]                                                                                   # kiszedjük a geneID aktuális elemét
    context_restored_final_temp <- context_restored_final[context_restored_final$geneID %in% temp_geneID,]      # egy temporary tibble-be kiemeljük az aktuális ID-hoz tartozó részt a context_restored_final-ből
    genes_temp <- genes[genes$geneID %in% temp_geneID,]                                                         # egy temporary tibble-be kiemeljük az aktuális ID-hoz tartozó részt a genes-ből
    context_restored_final_temp <- context_restored_final_temp %>%
      filter(context_restored_final_temp[,7] == as.character(genes_temp[1,7]))                                  # kiszedjük azokat a transzkripteket aminél megegyezik a strand a referenciával
    strand_filtered <- rbind(strand_filtered, context_restored_final_temp)                                      # összefűzzük a már készülő listával
    pb <- txtProgressBar(min = 1, max = length(geneIDs), style = 3)                                             # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL) 
  }
  return(strand_filtered)
}

strand_filtered <- strand_filtering(context_restored_final, genes)

stats[[5]] <- paste0("Strand filtering output           Gene count: ", 
                     nrow(unique(strand_filtered[,"geneID"])), 
                     "        Transcript count: ", 
                     nrow(unique(strand_filtered[,"transcriptID"])))         # a strand filtering statisztikái




strand_filtered_final <- strand_filtered

strand_filtered_final$geneID <- strand_filtered_final$geneID %>%
  str_replace("^", "\"") %>%
  str_replace("$", "\";")  

strand_filtered_final$transcriptID <- strand_filtered_final$transcriptID %>%
  str_replace("^", "\"") %>%
  str_replace("$", "\";")  

strand_filtered_final <- unite(strand_filtered_final , attributes, 9:12, sep = " ")

write_tsv(strand_filtered_final, "genome/RRPM_transcripts.gtf", col_names = F)

############################
#     ANNOTATION MERGE     #
############################

for_merge <- strand_filtered %>%
  filter(type != "transcript")

for_merge$geneID <- for_merge$geneID %>%
  str_replace("^", "\"") %>%
  str_replace("$", "\";")   

for_merge$transcriptID <- for_merge$transcriptID %>%
  str_replace("^", "\"") %>%
  str_replace("$", "\";")   

original <- read_tsv("genome/original_onlyexon.gtf", 
                  col_names = c("chr", "maker", "type", "start", "end", "att1", "strand", "att2", "attributes"))     # beolvassuk az eredeti GTF fájlt

original <- original %>%
  separate(attributes, c("transcriptID_label", "transcriptID",
                         "geneID_label", "geneID"), sep = " ")                    # szétbontjuk megfelelő oszlopokra az eredeti csak a géneket tartalmazó GTF fájlt

to_join <- anti_join(original, for_merge, by = "geneID") 

for_merge$att1 <- as.character(for_merge$att1)

merged <- bind_rows(for_merge, to_join)

merged <- arrange(merged, geneID, chr)

stats[[6]] <- paste0("Merged output                     Gene count: ", 
                     nrow(unique(merged[,"geneID"])), 
                     "        Transcript count: ", 
                     nrow(unique(merged[,"transcriptID"])))         # a merged statisztikái

merged <- unite(merged , attributes, 9:12, sep = " ")

write_tsv(merged, "genome/RRPM_transcripts.fixed.gtf", col_names = F)

write_lines(stats, "stats.log")

