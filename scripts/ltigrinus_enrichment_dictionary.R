library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# beolvassuk a headerfájlból készített szótárfájlt
fasta <- read_tsv("Desktop/alternative_splicing/ltigrinus/Sisbr1_GeneCatalog_proteins_20130805.aa.fasta.tsv.headers", col_names = c("fasta_transcript_id", "gene_id"))

# beolvassuk az eredeti annotációt ami alapján az expressziós analyzis készült
original <- read_tsv("Desktop/alternative_splicing/ltigrinus/ltigrinus_genome/Sisbr1_GeneCatalog_genes_20130805.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))

original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "gene_id", "tid", "annotation_transcript_id"), sep = " ") %>%
  select(gene_id, annotation_transcript_id)

original$gene_id <- original$gene_id %>% 
  str_replace("^\"", "") %>%
  str_replace("\";$", "")

original <- original[!duplicated(original), ]

final <- left_join(fasta, original, by = "gene_id")

write_tsv(final, "Desktop/alternative_splicing/ltigrinus/ltigrinus_dictionary.tsv")