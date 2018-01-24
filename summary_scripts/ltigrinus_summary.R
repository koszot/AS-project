library(tidyr)
library(readr)
library(dplyr)
library(stringr)

#### Fájlok betöltése ####

setwd("~/Desktop/alternative_splicing/ltigrinus/")

summary <- tibble()

genecatalog <- read_tsv("ltigrinus_genome/ltigrinus_onlygene.gtf",  col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  select(attributes) %>%
  separate(attributes, c("transcriptID_label", "transcriptID",
                         "geneID_label", "geneID"), sep = " ")

genecatalog_exons <- read_tsv("ltigrinus_genome/ltigrinus_onlyexon.gtf",  col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  select(attributes) %>%
  separate(attributes, c("transcriptID_label", "transcriptID",
                         "geneID_label", "geneID"), sep = " ") %>%
  select(transcriptID,geneID)

AS.annotation <- read_tsv("ltigrinus_genome/ltigrinus_AS_annotation.gtf",  col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  select(attributes) %>%
  separate(attributes, c("transcriptID_label", "transcriptID",
                         "geneID_label", "geneID"), sep = " ")

RRPM_transcripts <- read_tsv("ltigrinus_genome/ltigrinus_RRPM_transcripts.gtf", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand", "att2", "attributes"))

AS_count <- read_lines("ltigrinus_ASpli_binFeatures.log")

AS <- read_tsv("ltigrinus_enrichment/ltigrinus_genelvl_enrichment_AStranscripts.tsv", col_names = T, cols(value = col_character()))

DEVREG <- read_tsv("ltigrinus_enrichment/ltigrinus_genelvl_enrichment_DEVREGtranscripts.tsv", col_names = T, cols(value = col_character()))

genecatalogIPR <- read_lines("Sisbr1_GeneCatalog_proteins_20130805.aa.fasta.tsv")

#### GeneCatalog ####

summary[1,1] <- "GeneCatalog genes"
summary[1,2] <- length(unique(genecatalog$geneID))
summary[2,1] <- "GeneCatalog transcripts"
summary[2,2] <- length(unique(genecatalog$transcriptID))
summary[1:2,3] <- 0
colnames(summary) <- c("description", "N_ltigrinus", "%_ltigrinus")

#### New Annotation ####

summary[3,1] <- "New annotation genes"
summary[3,2] <- length(unique(AS.annotation$geneID))
summary[4,1] <- "New annotation transcripts"
summary[4,2] <- length(unique(AS.annotation$transcriptID))
summary[3,3] <- summary[3,2]/(summary[1,2]/100)
summary[4,3] <- summary[4,2]/(summary[2,2]/100)

#### Expressed ####

RRPM_transcripts_isoform <- RRPM_transcripts %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ") %>%
  filter(type == "transcript")

summary[5,1] <- "Expressed genes"
summary[5,2] <- length(unique(RRPM_transcripts_isoform$geneID))
summary[6,1] <- "Expressed transcripts"
summary[6,2] <- length(unique(RRPM_transcripts_isoform$transcriptID))
summary[5,3] <- summary[5,2]/(summary[3,2]/100)
summary[6,3] <- summary[6,2]/(summary[4,2]/100)

#### izoformák ####

ID <- unique(RRPM_transcripts_isoform$geneID)
isoformcount <- tibble()

for (x in 1:length(ID)) {
  annotation_temp <- RRPM_transcripts_isoform[RRPM_transcripts_isoform$geneID %in% ID[x],]
  isoformcount[x,1] <- ID[x]
  isoformcount[x,2] <- nrow(annotation_temp)
}

message("Legtöbb izoforma: ", max(isoformcount$V2))

summary[7,1] <- "Expressed genes /w 1 isoform"
summary[7,2] <- isoformcount %>% filter(V2 == 1) %>% nrow()
summary[8,1] <- "Expressed genes /w 2 isoforms"
summary[8,2] <- isoformcount %>% filter(V2 == 2) %>% nrow()
summary[9,1] <- "Expressed genes /w 3 isoforms"
summary[9,2] <- isoformcount %>% filter(V2 == 3) %>% nrow()
summary[10,1] <- "Expressed genes /w 4 isoforms"
summary[10,2] <- isoformcount %>% filter(V2 == 4) %>% nrow()
summary[11,1] <- "Expressed genes /w 5 isoforms"
summary[11,2] <- isoformcount %>% filter(V2 == 5) %>% nrow()
summary[12,1] <- "Expressed genes /w 6 isoforms"
summary[12,2] <- isoformcount %>% filter(V2 == 6) %>% nrow()
summary[7,3] <- summary[7,2]/(summary[5,2]/100)
summary[8,3] <- summary[8,2]/(summary[5,2]/100)
summary[9,3] <- summary[9,2]/(summary[5,2]/100)
summary[10,3] <- summary[10,2]/(summary[5,2]/100)
summary[11,3] <- summary[11,2]/(summary[5,2]/100)
summary[12,3] <- summary[12,2]/(summary[5,2]/100)

#### AS events ####

AS_count <- as_tibble(AS_count[-c(1,2,3,4,5,6,7,8)]) %>%
  separate(value, c("l1", "l2", "l3", "l4", "l5"), sep = "\t") %>%
  select(l2, l4) %>%
  separate(l2, c("event", "count"), sep = "=") %>%
  separate(l4, c("event2", "count2"), sep = "=") 
AS_count$count <- as.integer(AS_count$count)
AS_count$count2 <- as.integer(AS_count$count2)
AS_count[1,2] <- AS_count[1,2] + AS_count[7,4]
AS_count[2,2] <- AS_count[2,2] + AS_count[8,4]
AS_count[3,2] <- AS_count[3,2] + AS_count[9,4]
AS_count[4,2] <- AS_count[4,2] + AS_count[10,4]
AS_count <- AS_count[-(5:10),-(3:4)]

summary[13,1] <- "Alternative Splicing ES bins"
summary[13,2] <- AS_count[1,2]
summary[14,1] <- "Alternative Splicing IR bins"
summary[14,2] <- AS_count[2,2]
summary[15,1] <- "Alternative Splicing Alt5SS bins"
summary[15,2] <- AS_count[3,2]
summary[16,1] <- "Alternative Splicing Alt3SS bins"
summary[16,2] <- AS_count[4,2]
summary[13,3] <- summary[13,2]/(sum(summary[13:16,2])/100)
summary[14,3] <- summary[14,2]/(sum(summary[13:16,2])/100)
summary[15,3] <- summary[15,2]/(sum(summary[13:16,2])/100)
summary[16,3] <- summary[16,2]/(sum(summary[13:16,2])/100)

#### AS and DEVREG ####

summary[17,1] <- "Genes /w 2 FPKM and Alternative Splicing"
summary[17,2] <- AS %>% nrow()
summary[17,3] <- summary[17,2]/(summary[3,2]/100)
summary[18,1] <- "Genes /w 2 FPKM and 2 Fold Change"
summary[18,2] <- DEVREG %>% nrow()
summary[18,3] <- summary[18,2]/(summary[3,2]/100)

#### IPR and GO ####

IPR <- as_tibble(genecatalogIPR[str_detect(genecatalogIPR, "\tIPR")], header = F)
colnames(IPR) <- "all"
IPR <- IPR %>%
  separate(all, c("gene_id", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t")
GO <- IPR %>%
  filter(GO != "") 

summary[19,1] <- "GeneCatalog genes /w InterPro Domains"
summary[19,2] <- length(unique(IPR$gene_id))
summary[19,3] <- summary[19,2] / (summary[1,2]/100)
summary[20,1] <- "GeneCatalog genes /w GO annotations"
summary[20,2] <- length(unique(GO$gene_id))
summary[20,3] <- summary[20,2] / (summary[1,2]/100)

#### single és multi-exon genecatalog gének ####

ID <- unique(genecatalog_exons$geneID)

exoncount <- tibble()

for (x in 1:length(ID)) {
  temp <- genecatalog_exons[genecatalog_exons$geneID %in% ID[x],]
  exoncount[x,1] <- ID[x]
  exoncount[x,2] <- nrow(temp)
}

summary[21,1] <- "GeneCatalog single-exon genes"
summary[21,2] <- exoncount %>% filter(V2 == 1) %>% nrow()
summary[21,3] <- summary[21,2] / (summary[1,2]/100)
summary[22,1] <- "GeneCatalog multi-exon genes"
summary[22,2] <- exoncount %>% filter(V2 > 1) %>% nrow()
summary[22,3] <- summary[22,2] / (summary[1,2]/100)

#### expresszált single és multi-exon gének ####

RRPM_transcripts_geneID <- RRPM_transcripts %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ") %>%
  filter(type == "transcript") %>%
  select(geneID) %>%
  distinct()

summary[23,1] <- "Expressed single-exon genes"
summary[23,2] <- exoncount %>% filter(V2 == 1) %>% select(geneID=V1, count=V2) %>% inner_join(RRPM_transcripts_geneID) %>% nrow()
summary[23,3] <- summary[23,2] / (summary[5,2]/100)
summary[24,1] <- "Expressed multi-exon genes"
summary[24,2] <- exoncount %>% filter(V2 > 1) %>% select(geneID=V1, count=V2) %>% inner_join(RRPM_transcripts_geneID) %>% nrow()
summary[24,3] <- summary[24,2] / (summary[5,2]/100)

# fúziós géneket le kell kezelni manuálisan !!! ha a kövi FALSE
summary[23,2] + summary[24,2] == summary[5,2]

#### expresszált single és multi-exon gének AS-ben ####

isoformcount <- isoformcount %>% select(geneID=V1, count_isoforms=V2)

summary[25,1] <- "Expressed single-exon genes /w AS"
summary[25,2] <- exoncount %>% filter(V2 == 1) %>% select(geneID=V1, count=V2) %>% inner_join(RRPM_transcripts_geneID) %>% inner_join(isoformcount) %>% filter(count_isoforms > 1) %>% nrow()
summary[25,3] <- summary[25,2] / (summary[23,2]/100)
summary[26,1] <- "Expressed multi-exon genes /w AS"
summary[26,2] <- exoncount %>% filter(V2 > 1) %>% select(geneID=V1, count=V2) %>% inner_join(RRPM_transcripts_geneID) %>% inner_join(isoformcount) %>% filter(count_isoforms > 1) %>% nrow()
summary[26,3] <- summary[26,2] / (summary[24,2]/100)

summary <- summary[c(1:2,21:22,3:6,23,25,24,26,7:20),]

#### enrichment sima ####

enrichment <- read_tsv("ltigrinus_enrichment/ltigrinus_genelvl_enrichment_fisher_classic.tsv")

summary[27,1] <- "AS enriched GOs"
summary[27,2] <- enrichment %>% filter(type2 == "AS") %>% nrow()
summary[27,3] <- 0
summary[28,1] <- "DEVREG enriched GOs"
summary[28,2] <- enrichment %>% filter(type2 == "DEVREG") %>% nrow()
summary[28,3] <- 0

summary$`%_ltigrinus` <- round(summary$`%_ltigrinus`, digits = 3)

write_tsv(summary, "ltigrinus_summary.tsv")



