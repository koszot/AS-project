library(dplyr)
library(readr)
library(stringr)
library(tidyr)

setwd("~/Desktop/alternative_splicing/umaydis/")

isoforms.fpkm <- read_tsv("umaydis_CUFFdiff/isoforms.fpkm_tracking", col_names = T, cols(.default = col_guess(), tracking_id = col_character(), gene_id = col_character()))  # betöltjük az isoform FPKM táblát
genes.fpkm <- read_tsv("umaydis_CUFFdiff/genes.fpkm_tracking", col_names = T, cols(.default = col_guess(), tracking_id = col_character(), gene_id = col_character()))           # betöltjük az genes FPKM táblát

isoforms.fpkm <- isoforms.fpkm %>%
  select(transcript_id = tracking_id, gene_id, FPKM = UM_1_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat

genes.fpkm <- genes.fpkm %>%
  select(tracking_id, gene_id, FPKM = UM_1_FPKM)              # kiválasztjuk a minket érdeklő oszlopokat

#################################################################################
#     Azon gének kiszedése ahol az FPKM nem éri el a 2-őt egyik fázisban se     #
#################################################################################

genes_FPKM_filtering <- function(genes.fpkm_FUN, isoforms.fpkm.AS_FUN) {
  genes.fpkm.FPKM <- genes.fpkm_FUN %>% filter(FPKM > 2)                              # kiszedjük azokat a géneket amiknél az FPKM érték nagyobb mint 2 bármelyik fejlődési szakaszban
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
    select(gene_id, transcript_id, FPKM) 
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

##################
#     fájlok     #
##################

geneUniverse <- as_tibble(unique(genes.fpkm$gene_id))
write_tsv(geneUniverse, "umaydis_enrichment/umaydis_genelvl_enrichment_geneUniverse.tsv")

AS.transcripts <- as_tibble(unique(AS.transcripts$gene_id))
write_tsv(AS.transcripts, "umaydis_enrichment/umaydis_genelvl_enrichment_AStranscripts.tsv")
