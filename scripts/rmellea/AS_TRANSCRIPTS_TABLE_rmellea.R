library(tidyr)
library(readr)
library(dplyr)
library(stringr)

setwd("~/Desktop/alternative_splicing/rmellea/")

isoforms.fpkm <- read_tsv("rmellea_CUFFdiff/isoforms.fpkm_tracking", col_names = T, cols(.default = col_guess(), tracking_id = col_character(), gene_id = col_character()))  # betöltjük az isoform FPKM táblát
genes.fpkm <- read_tsv("rmellea_CUFFdiff/genes.fpkm_tracking", col_names = T, cols(.default = col_guess(), tracking_id = col_character(), gene_id = col_character()))           # betöltjük az genes FPKM táblát

isoforms.fpkm <- isoforms.fpkm %>%
  select(transcript_id = tracking_id, gene_id, VM_FPKM, P_FPKM, YFB_FPKM, FB_K_FPKM, FB_T_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat

genes.fpkm <- genes.fpkm %>%
  select(tracking_id, gene_id, VM_FPKM, P_FPKM, YFB_FPKM, FB_K_FPKM, FB_T_FPKM)              # kiválasztjuk a minket érdeklő oszlopokat

stats <- read_tsv("rmellea_stats.log")                                                 # a stats list bevezetése

##########################################
#     Alternaive Splicing filterezés     #
##########################################

alternative_splicing_associated_genes <- function(isoforms.fpkm_FUN) {
  geneids <- unique(isoforms.fpkm_FUN$gene_id)                                           # kimentjük az egyedi geneID-kat
  isoforms.fpkm.AS <- data.frame()                                                       # megcsináljuk az üres data framet
  for (x in 1:length(geneids)) {                                                         # egy ciklussal végigmegyünk az gén ID-ken
    geneID <- geneids[x]                                                                 # kiszedjük az aktuális geneID-t
    isoforms_temp <- isoforms.fpkm_FUN[isoforms.fpkm_FUN$gene_id %in% geneID,]           # aktuális geneID-hoz tartozó isoformák kiszedése
    isoforms_temp$isoformcount <- nrow(isoforms_temp)                                    # az adott génhez tartozó isoformák leszámolása
    isoforms.fpkm.AS <- rbind(isoforms.fpkm.AS, isoforms_temp)       
    pb <- txtProgressBar(min = 1, max = length(geneids), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL) 
  }
  isoforms.fpkm.AS <- filter(isoforms.fpkm.AS, isoformcount != 1)                # kiszedjük azokat az isoformákat amik részt vesznek AS-ben
  return(isoforms.fpkm.AS)
}

isoforms.fpkm.AS <- alternative_splicing_associated_genes(isoforms.fpkm)       # lefuttatjuk az AS szűrést

stats[7,1] <- "Isoforms associated /w Alternative Splicing"
stats[7,2] <- nrow(unique(isoforms.fpkm.AS[,"gene_id"]))
stats[7,3] <- nrow(unique(isoforms.fpkm.AS[,"transcript_id"]))

#################################################################################
#     Azon gének kiszedése ahol az FPKM nem éri el a 2-őt egyik fázisban se     #
#################################################################################

genes_FPKM_filtering <- function(genes.fpkm_FUN, isoforms.fpkm.AS_FUN) {
  genes.fpkm.FPKM <- genes.fpkm_FUN %>% filter(VM_FPKM > 2 | 
                                                 P_FPKM > 2 | 
                                                 YFB_FPKM > 2 | 
                                                 FB_K_FPKM > 2 | 
                                                 FB_T_FPKM > 2)                              # kiszedjük azokat a géneket amiknél az FPKM érték nagyobb mint 2 bármelyik fejlődési szakaszban
  geneids <- unique(genes.fpkm.FPKM$gene_id)                                             # kimentjük az egyedi geneID-kat
  isoforms.fpkm.AS.FPKM <- data.frame()                                                  # megcsináljuk az üres data framet
  for (x in 1:length(geneids)) {                                                         # egy ciklussal végigmegyünk az geneID-ken
    geneID <- geneids[x]                                                                 # kiszedjük az aktuális geneID-t
    isoforms_temp <- filter(isoforms.fpkm.AS_FUN, gene_id == geneID)                     # aktuális geneID-hoz tartozó isoformák kiszedése
    isoforms.fpkm.AS.FPKM <- rbind(isoforms.fpkm.AS.FPKM, isoforms_temp)       
    pb <- txtProgressBar(min = 1, max = length(geneids), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL) 
  }
  isoforms.fpkm.AS.FPKM <- isoforms.fpkm.AS.FPKM %>%
    select(gene_id, transcript_id, VM_FPKM, P_FPKM, YFB_FPKM, FB_K_FPKM, FB_T_FPKM) 
  return(isoforms.fpkm.AS.FPKM)
}

isoforms.fpkm.AS.FPKM <- genes_FPKM_filtering(genes.fpkm, isoforms.fpkm.AS)

stats[8,1] <- "Isoforms of genes /w more then 2 FPKM at any stage"
stats[8,2] <- nrow(unique(isoforms.fpkm.AS.FPKM[,"gene_id"]))
stats[8,3] <- nrow(unique(isoforms.fpkm.AS.FPKM[,"transcript_id"]))

###########################
#     Rmax filterezés     #
###########################

R_value <- function(isoforms.fpkm.AS.FPKM_FUN) {
  isoforms.fpkm.AS.FPKM.Rmax <- data.frame()
  ID <- unique(isoforms.fpkm.AS.FPKM_FUN$gene_id)
  for (x in 1:length(ID)) {
    ID2 <- ID[x]
    isoforms_temp <- isoforms.fpkm.AS.FPKM_FUN[isoforms.fpkm.AS.FPKM_FUN$gene_id %in% ID2,]  
    isoforms_temp$VM_Ri <- isoforms_temp$VM_FPKM / sum(isoforms_temp$VM_FPKM)
    isoforms_temp$P_Ri <- isoforms_temp$P_FPKM / sum(isoforms_temp$P_FPKM)
    isoforms_temp$YFB_Ri <- isoforms_temp$YFB_FPKM / sum(isoforms_temp$YFB_FPKM)
    isoforms_temp$FB_K_Ri <- isoforms_temp$FB_K_FPKM / sum(isoforms_temp$FB_K_FPKM)
    isoforms_temp$FB_T_Ri <- isoforms_temp$FB_T_FPKM / sum(isoforms_temp$FB_T_FPKM)
    isoforms.fpkm.AS.FPKM.Rmax <- rbind(isoforms.fpkm.AS.FPKM.Rmax, isoforms_temp)       
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL) 
  }
  isoforms.fpkm.AS.FPKM.Rmax <- isoforms.fpkm.AS.FPKM.Rmax %>%
    filter(VM_Ri > 0.1 |
             P_Ri > 0.1 |
             YFB_Ri > 0.1 |
             FB_K_Ri > 0.1 |
             FB_T_Ri > 0.1 )
  geneids <- unique(isoforms.fpkm.AS.FPKM.Rmax$gene_id)                                           # kimentjük az egyedi geneID-kat
  isoforms.fpkm.AS.FPKM.Rmax.2 <- data.frame()                                                       # megcsináljuk az üres data framet
  for (i in 1:length(geneids)) {                                                         # egy ciklussal végigmegyünk az gén ID-ken
    geneID <- geneids[i]                                                                 # kiszedjük az aktuális geneID-t
    isoforms_temp <- isoforms.fpkm.AS.FPKM.Rmax[isoforms.fpkm.AS.FPKM.Rmax$gene_id %in% geneID,]           # aktuális geneID-hoz tartozó isoformák kiszedése
    isoforms_temp$isoformcount <- nrow(isoforms_temp)                                    # az adott génhez tartozó isoformák leszámolása
    isoforms.fpkm.AS.FPKM.Rmax.2 <- rbind(isoforms.fpkm.AS.FPKM.Rmax.2, isoforms_temp)       
    pb <- txtProgressBar(min = 1, max = length(geneids), style = 3)                      # progress bar
    setTxtProgressBar(pb, i, title = NULL, label = NULL) 
  }
  isoforms.fpkm.AS.FPKM.Rmax.2 <- filter(isoforms.fpkm.AS.FPKM.Rmax.2, isoformcount != 1)
  return(isoforms.fpkm.AS.FPKM.Rmax.2)
}

isoforms.fpkm.AS.FPKM.Rmax <- R_value(isoforms.fpkm.AS.FPKM)

stats[9,1] <- "Isoforms w/ at least 10% contribution to overall gene expression"
stats[9,2] <- nrow(unique(isoforms.fpkm.AS.FPKM.Rmax[,"gene_id"]))
stats[9,3] <- nrow(unique(isoforms.fpkm.AS.FPKM.Rmax[,"transcript_id"]))

###################################
#     FOLD CHANGE KISZÁMOLÁSA     #
###################################

fold_change <- function(isoforms.fpkm.AS.FPKM.Rmax_FUN) {
  for (x in 1:length(isoforms.fpkm.AS.FPKM.Rmax_FUN$transcript_id)) {
    isoforms.fpkm.AS.FPKM.Rmax_FUN$VM_to_P1_FoldChange[x] <- max(isoforms.fpkm.AS.FPKM.Rmax_FUN[x,3:4]) / min(isoforms.fpkm.AS.FPKM.Rmax_FUN[x,3:4])
    isoforms.fpkm.AS.FPKM.Rmax_FUN$P1_to_FB_FoldChange[x] <- max(isoforms.fpkm.AS.FPKM.Rmax_FUN[x,4:7]) / min(isoforms.fpkm.AS.FPKM.Rmax_FUN[x,4:7])
    isoforms.fpkm.AS.FPKM.Rmax_FUN$VM_to_FB_FoldChange[x] <- max(isoforms.fpkm.AS.FPKM.Rmax_FUN[x,3:7]) / min(isoforms.fpkm.AS.FPKM.Rmax_FUN[x,3:7])
    pb <- txtProgressBar(min = 1, max = length(isoforms.fpkm.AS.FPKM.Rmax_FUN$transcript_id), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(isoforms.fpkm.AS.FPKM.Rmax_FUN)
}

isoforms.fpkm.AS.FPKM.Rmax.FC <- fold_change(isoforms.fpkm.AS.FPKM.Rmax)

isoforms.fpkm.AS.FPKM.Rmax.FC <- isoforms.fpkm.AS.FPKM.Rmax.FC[,c(1:7,14:16)]

results <- isoforms.fpkm.AS.FPKM.Rmax.FC

rm(genes.fpkm, isoforms.fpkm, isoforms.fpkm.AS, isoforms.fpkm.AS.FPKM, isoforms.fpkm.AS.FPKM.Rmax, isoforms.fpkm.AS.FPKM.Rmax.FC)

#################################
#     CDS RÉGIÓK VIZSGÁLATA     #
#################################

cds <- read_tsv("rmellea_transcripts.fasta.transdecoder_dir/headers.cds", col_names = F)

longest_orf <- function(cds_FUN) {
  cds.longest <- data.frame()
  ID <- unique(cds_FUN$X2)
  for (x in 1:length(ID)) {
    cds_temp <- cds_FUN[cds_FUN$X2 %in% ID[x],]
    cds.longest <- rbind(cds.longest, cds_temp[1,])
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(cds.longest)
}

cds <- longest_orf(cds) %>%
  select(transcript_id = X2, CDS_length = X5)

results.CDS <- left_join(results, cds, by = "transcript_id")

results.CDS$CDS_length <- results.CDS$CDS_length * 3

results.CDS$overall_FPKM <- results.CDS$VM_FPKM + results.CDS$P_FPKM + results.CDS$YFB_FPKM + results.CDS$FB_K_FPKM + results.CDS$FB_T_FPKM

results.CDS <- results.CDS[,c(1:2,12,3:11)]

results.CDS <- results.CDS %>%
  arrange(gene_id, desc(overall_FPKM))

CDS_length_check <- function(results.CDS_FUN){
  mydata <- data.frame()
  ID <- unique(results.CDS_FUN$gene_id)
  for (x in 1:length(ID)) {
    results_temp <- results.CDS_FUN[results.CDS_FUN$gene_id %in% ID[x],]
    results_temp$is_CDS_length_same_w_primary <- results_temp[,12] == as.integer(results_temp[1,12])
    results_temp[1,13] <- NA
    mydata <- rbind(mydata, results_temp)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(mydata)
} 

results.CDS2 <- CDS_length_check(results.CDS)

rm(cds, results, results.CDS)

########################################
#     DOMAIN ÖSSZETÉTEL VIZSGÁLATA     #
########################################

raw_interpro <- read_lines("rmellea_proteins_all.fasta.tsv")

IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)

colnames(IPR) <- "all"

IPR <- IPR %>%
  separate(all, c("transcript_id", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t")

domain_count <- function(IPR_FUN) {
  mydata <- data.frame()
  ID <- unique(IPR_FUN$transcript_id)
  for (x in 1:length(ID)) {
    IPR_temp <- IPR_FUN[IPR_FUN$transcript_id %in% ID[x],]
    mydata_temp <- data.frame()
    mydata_temp[1,1] <- ID[x]
    mydata_temp[1,2] <- nrow(IPR_temp)
    mydata <- rbind(mydata, mydata_temp)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(mydata)
}

domain.count <- domain_count(IPR)

colnames(domain.count) <- c("transcript_id", "domain_count")

results.CDS2.DC <- left_join(results.CDS2, domain.count, by = "transcript_id")

results.CDS2.DC$domain_count[is.na(results.CDS2.DC$domain_count)] <- 0

domain_same <- function(results.CDS2.DC_FUN) {
  mydata <- data.frame()
  ID <- unique(results.CDS2.DC_FUN$gene_id)
  for (x in 1:length(ID)) {
    temp <- results.CDS2.DC_FUN[results.CDS2.DC_FUN$gene_id %in% ID[x],]
    temp$no_domain_count_change <- temp[,14] == as.integer(temp[1,14])
    temp[1,15] <- NA
    mydata <- rbind(mydata, temp)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(mydata)
}

results.CDS2.DC.DS <- domain_same(results.CDS2.DC)

domain_gain <- function(results.CDS2.DC.DS_FUN) {
  mydata <- data.frame()
  ID <- unique(results.CDS2.DC.DS_FUN$gene_id)
  for (x in 1:length(ID)) {
    temp <- results.CDS2.DC.DS_FUN[results.CDS2.DC.DS_FUN$gene_id %in% ID[x],]
    temp$is_domain_gain <- temp[,14] > as.integer(temp[1,14])
    temp[1,16] <- NA
    mydata <- rbind(mydata, temp)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(mydata)
}

results.CDS2.DC.DS.DG <- domain_gain(results.CDS2.DC.DS)

domain_loss <- function(results.CDS2.DC.DS.DG_FUN) {
  mydata <- data.frame()
  ID <- unique(results.CDS2.DC.DS.DG_FUN$gene_id)
  for (x in 1:length(ID)) {
    temp <- results.CDS2.DC.DS.DG_FUN[results.CDS2.DC.DS.DG_FUN$gene_id %in% ID[x],]
    temp$is_domain_loss <- temp[,14] < as.integer(temp[1,14])
    temp[1,17] <- NA
    mydata <- rbind(mydata, temp)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(mydata)
}

results.CDS.domaincomp <- domain_loss(results.CDS2.DC.DS.DG)

rm(domain.count, results.CDS2,results.CDS2.DC, results.CDS2.DC.DS, results.CDS2.DC.DS.DG, raw_interpro)

########################
#     DOMAIN HOSSZ     #
########################

domains <- IPR %>%
  select(transcript_id, start, end) %>%
  mutate(length = as.integer(end) - as.integer(start) + 1) %>%
  select(transcript_id, length)

domain_length <- function(domains_FUN) {
  mydata <- data.frame()
  mydata_final <- data.frame()
  ID <- unique(domains_FUN$transcript_id)
  for (x in 1:length(ID)) {
    temp <- domains_FUN[domains_FUN$transcript_id %in% ID[x],]
    mydata[1,1] <- ID[x]
    mydata[1,2] <- sum(temp[,2])
    mydata_final <- rbind(mydata_final, mydata)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(mydata_final)
}

domains_collected <- domain_length(domains)

colnames(domains_collected) <- c("transcript_id", "domain_length")

results.CDS.domaincomp.length <- left_join(results.CDS.domaincomp, domains_collected, by = "transcript_id")

domain_length_compare <- function(results.CDS.domaincomp.length_FUN) {
  mydata <- data.frame()
  ID <- unique(results.CDS.domaincomp.length_FUN$gene_id)
  for (x in 1:length(ID)) {
    temp <- results.CDS.domaincomp.length_FUN[results.CDS.domaincomp.length_FUN$gene_id %in% ID[x],]
    primary_length <- temp$domain_length[1]
    temp$is_domain_damage <- (temp$domain_length - primary_length) != 0
    mydata <- rbind(mydata, temp)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(mydata)
}

df <- domain_length_compare(results.CDS.domaincomp.length)

domain_length_compare_fix <- function(df_FUN) {
  mydata <- data.frame()
  ID <- ID <- unique(df_FUN$transcript_id)
  for (x in 1:length(ID)) {
    temp <- df_FUN[df_FUN$transcript_id %in% ID[x],]
    if(is.na(temp[1,15])) {
      temp[1,19] <- NA                                         # no_domain_count_change == NA, akkor primary transzkript vagyis minden NA
    } else if (temp[1,15] == F) {
      temp[1,19] <- F                                          # no_domain_count_change == FALSE, akkor történt domainszám változás vagyis vagy gain vagy loss és nem dmg
    } else if (temp[1,15] == T & temp[1,14] == 0) {
      temp[1,19] <- F                                          # no_domain_count_change == TRUE, tehát nem történt domain szám voltozás, DE domain_count == 0, vagyis nincs domain benne akkor is FALSE
    }
    mydata <- rbind(mydata, temp)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(mydata)
}

df2 <- domain_length_compare_fix(df)

df2 <- df2[,-18]

########################
#     GO ANNOTÁCIÓ     #
########################

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
  GO_final_merged[x,2] <- paste(unique(unlist(GO_merged)), collapse = "|")
  pb <- txtProgressBar(min = 1, max = length(TID_GO), style = 3)                      # progress bar
  setTxtProgressBar(pb, x, title = NULL, label = NULL)
}

colnames(GO_final_merged) <- c("transcript_id", "GO_annotations")

df3 <- left_join(df2, GO_final_merged, by = "transcript_id")

write_tsv(df3, "rmellea_results.tsv")

write_tsv(stats, "rmellea_stats.log")
