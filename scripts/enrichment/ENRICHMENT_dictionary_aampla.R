library(readr)
library(dplyr)
library(tidyr)
library(stringr)

setwd("~/Desktop/MTA/AS_project/FILES_aampla/ENRICHMENT_aampla/")

# beolvassuk a headerfájlból készített szótárfájlt
fasta <- read_tsv("Auramp1_GeneCatalog_proteins_20160719.aa.fasta.headers", col_names = c("fasta_transcript_id", "gene_id"))

# beolvassuk az eredeti annotációt ami alapján az expressziós analyzis készült
original <- read_tsv("../GENOME_aampla/Auramp1_GeneCatalog_genes_20160719.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))

original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "gene_id", "tid", "annotation_transcript_id"), sep = " ") %>%
  select(gene_id, annotation_transcript_id)

original$gene_id <- original$gene_id %>% 
  str_replace("^\"", "") %>%
  str_replace("\";$", "")

original <- original[!duplicated(original), ]

final <- left_join(fasta, original, by = "gene_id")

write_tsv(final, "aampla_dictionary.tsv")
