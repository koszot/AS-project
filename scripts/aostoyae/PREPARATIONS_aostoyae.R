#!/usr/bin/Rscript

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# beolvassuk az eredeti GFF3-ról átalakított GTF fájlt

original <- read_tsv("./genome/original.gtf", 
                  col_names = c("chr", 
                                "maker", 
                                "type", 
                                "start", 
                                "end", 
                                "att1", 
                                "strand", 
                                "att2", 
                                "attributes"))

# szétbontjuk az oszlopokat és kiszedjük az exonokat jelölő sorokat

original <- original %>%
  separate(attributes, c("tid", 
                         "tid_id", 
                         "geneid", 
                         "geneid_id", 
                         "genename", 
                         "genename_name"), 
           sep = " ") %>%
  filter(type == "exon")

# only_exon GTF fájl kiírása

original_fixed <- original[,-c(13:14)] 
original_fixed <- unite(original_fixed, attributes, 9:12, sep = " ")

write.table(original_fixed, file = "./genome/original_onlyexon.gtf", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# csak a gének határait tartalmazó GTF fájl létrehozása

exon_to_genes <- function(original_forFUN)
{
  genes <- data.frame()
  original_ID <- unique(original_forFUN$geneid_id)
  genes_list <- lapply(original_ID, function(original_ID_forFUN,original_forFUN_2) 
  {
    original_temp <- original_forFUN_2[original_forFUN_2$geneid_id %in% original_ID_forFUN,]
    genes_temp <- data.frame()
    genes_temp <- original_temp[1,1:2]
    genes_temp[1,3] <- "gene"
    genes_temp[1,4] <- original_temp[1,4]
    genes_temp[1,5] <- original_temp[length(original_temp$geneid_id),5]
    genes_temp <- cbind(genes_temp, original_temp[1,6:14])
    genes <- rbind(genes, genes_temp)
  }, 
  original_forFUN_2 = original_forFUN)
  genes <- do.call("rbind", genes_list)
}

genes <- exon_to_genes(original)
genes <- genes[,-c(13:14)]
genes <- unite(genes, attributes, 9:12, sep = " ")

write.table(genes, file = "./genome/original_onlygene.gtf", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# exon pozíciók átalakítása a RRPM-hez

exon_fix <- function(original_forFUN) 
{
  exon_fixed <- original_forFUN
  exon_fixed$chr <- exon_fixed$geneid_id
  exon_fixed$chr <- exon_fixed$chr %>% 
    str_replace("\"", "") %>%
    str_replace("\";", "")
  exon_fixed2 <- data.frame()
  chr_ID <- unique(exon_fixed$chr)
  exon_fixed_list <- lapply(chr_ID, function(ID_forFUN, exon_fixed_forFUN) 
  {
    exon_fixed_temp <- exon_fixed_forFUN[exon_fixed_forFUN$chr %in% ID_forFUN,]
    start_value <- as.integer(exon_fixed_temp[1,4])
    exon_fixed_temp[,4] <- exon_fixed_temp[,4] - start_value + 1
    exon_fixed_temp[,5] <- exon_fixed_temp[,5] - start_value + 1
    exon_fixed2 <- rbind(exon_fixed2, exon_fixed_temp) 
  }, 
  exon_fixed_forFUN = exon_fixed)
  exon_fixed3 <- do.call("rbind", exon_fixed_list)
}

exon_fixed <- exon_fix(original)
exon_fixed <- exon_fixed[,-c(13:14)]
exon_fixed <- unite(exon_fixed, attributes, 9:12, sep = " ")

write.table(exon_fixed, file = "./genome/original_fixed.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
