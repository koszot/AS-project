library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# beolvassuk az eredeti GFF annotációs fájlt és átalakítjuk GTF formtátumra ami jó a Cufflinksnek

setwd("~/Desktop/alternative_splicing/aampla/")

original <- read_tsv("aampla_genome/Auramp1_GeneCatalog_genes_20160719.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))

original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "name_id", "tid", "tid_id"), sep = " ")

original$name_id <- original$tid_id
original$name <- "transcript_id"
original$tid <- "gene_id"

original$name_id <- original$name_id %>% 
  str_replace("^", "\"") %>%
  str_replace("$", "\";")  

original$tid_id <- original$tid_id %>% 
  str_replace("^", "\"") %>%
  str_replace("$", "\";")  

original_fixed <- unite(original, attributes, 9:12, sep = " ")

write.table(original_fixed, file = "aampla_genome/aampla_onlyexon.gtf", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# csak a gének határait tartalmazó GTF fájl létrehozása

exon_to_genes <- function(original_forFUN)
{
  genes <- data.frame()
  original_ID <- unique(original_forFUN$name_id)
  genes_list <- lapply(original_ID, function(original_ID_forFUN,original_forFUN_2) 
  {
    original_temp <- original_forFUN_2[original_forFUN_2$name_id %in% original_ID_forFUN,]
    genes_temp <- data.frame()
    genes_temp <- original_temp[1,1:2]
    genes_temp[1,3] <- "gene"
    genes_temp[1,4] <- original_temp[1,4]
    genes_temp[1,5] <- original_temp[length(original_temp$name_id),5]
    genes_temp <- cbind(genes_temp, original_temp[1,6:12])
    genes <- rbind(genes, genes_temp)
  }, 
  original_forFUN_2 = original_forFUN)
  genes <- do.call("rbind", genes_list)
}

genes <- exon_to_genes(original)

genes <- unite(genes, attributes, 9:12, sep = " ")

write.table(genes, file = "./aampla_genome/aampla_onlygene.gtf", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# exon pozíciók átalakítása a RRPM-hez

exon_fix <- function(original_forFUN) 
{
  exon_fixed <- original_forFUN
  exon_fixed$chr <- exon_fixed$name_id
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

exon_fixed <- unite(exon_fixed, attributes, 9:12, sep = " ")

write.table(exon_fixed, file = "./aampla_genome/aampla_fixed.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
