library(topGO)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

setwd("~/Desktop/alternative_splicing/ccinerea_AmutBmut/")

isoforms.fpkm <- read_tsv("ccinerea_AmutBmut_CUFFdiff/isoforms.fpkm_tracking", col_names = T, cols(.default = col_guess(), tracking_id = col_character(), gene_id = col_character()))  # betöltjük az isoform FPKM táblát
genes.fpkm <- read_tsv("ccinerea_AmutBmut_CUFFdiff/genes.fpkm_tracking", col_names = T, cols(.default = col_guess(), tracking_id = col_character(), gene_id = col_character()))           # betöltjük az genes FPKM táblát

isoforms.fpkm <- isoforms.fpkm %>%
  select(transcript_id = tracking_id, gene_id, VM_FPKM, H_FPKM, P1_FPKM, P2_FPKM, YFB_K_FPKM, YFB_L_FPKM, YFB_T_FPKM, FB_KL_FPKM, FB_T_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat

genes.fpkm <- genes.fpkm %>%
  select(tracking_id, gene_id, VM_FPKM, H_FPKM, P1_FPKM, P2_FPKM, YFB_K_FPKM, YFB_L_FPKM, YFB_T_FPKM, FB_KL_FPKM, FB_T_FPKM)              # kiválasztjuk a minket érdeklő oszlopokat

#################################################################################
#     Azon gének kiszedése ahol az FPKM nem éri el a 2-őt egyik fázisban se     #
#################################################################################

genes_FPKM_filtering <- function(genes.fpkm_FUN, isoforms.fpkm.AS_FUN) {
  genes.fpkm.FPKM <- genes.fpkm_FUN %>% filter(VM_FPKM > 2 | 
                                                 H_FPKM > 2 | 
                                                 P1_FPKM > 2 | 
                                                 P2_FPKM > 2 | 
                                                 YFB_K_FPKM > 2 | 
                                                 YFB_L_FPKM > 2 | 
                                                 YFB_T_FPKM > 2 | 
                                                 FB_KL_FPKM > 2 | 
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
    select(gene_id, transcript_id, VM_FPKM, H_FPKM, P1_FPKM, P2_FPKM, YFB_K_FPKM, YFB_L_FPKM, YFB_T_FPKM, FB_KL_FPKM, FB_T_FPKM) 
  return(isoforms.fpkm.AS.FPKM)
}

significantTranscripts <- genes_FPKM_filtering(genes.fpkm, isoforms.fpkm)

###########################
#     AS TRANSZKRIPTEK    #
###########################

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

AS.transcripts <- alternative_splicing_associated_genes(significantTranscripts)       # lefuttatjuk az AS szűrést

###############################
#     DEVREG TRANSZKRIPTEK    #
###############################

fold_change <- function(isoforms.fpkm.AS.FPKM.Rmax_FUN) {
  for (x in 1:length(isoforms.fpkm.AS.FPKM.Rmax_FUN$transcript_id)) {
    isoforms.fpkm.AS.FPKM.Rmax_FUN$VM_to_FB_FoldChange[x] <- max(isoforms.fpkm.AS.FPKM.Rmax_FUN[x,3:11]) / min(isoforms.fpkm.AS.FPKM.Rmax_FUN[x,3:11])            # VM to FB
    pb <- txtProgressBar(min = 1, max = length(isoforms.fpkm.AS.FPKM.Rmax_FUN$transcript_id), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL)
  }
  return(isoforms.fpkm.AS.FPKM.Rmax_FUN)
}

DEVREG.transcripts <- fold_change(significantTranscripts)

DEVREG.transcripts <- DEVREG.transcripts %>%
  filter(VM_to_FB_FoldChange > 2)

##################
#     fájlok     #
##################

geneUniverse <- as_tibble(unique(genes.fpkm$gene_id))
write_tsv(geneUniverse, "ccinerea_AmutBmut_enrichment/ccinerea_AmutBmut_genelvl_enrichment_geneUniverse.tsv")

AS.transcripts <- as_tibble(unique(AS.transcripts$gene_id))
write_tsv(AS.transcripts, "ccinerea_AmutBmut_enrichment/ccinerea_AmutBmut_genelvl_enrichment_AStranscripts.tsv")

DEVREG.transcripts <- as_tibble(unique(DEVREG.transcripts$gene_id))
write_tsv(DEVREG.transcripts, "ccinerea_AmutBmut_enrichment/ccinerea_AmutBmut_genelvl_enrichment_DEVREGtranscripts.tsv")
