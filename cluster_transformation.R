library(dplyr)
library(readr)
library(stringr)
library(tidyr)

families <- read_tsv("~/Desktop/alternative_splicing/enrichment/all_species_protein_HFX.fnodes", col_names = c("family", "gene_id"))

###################
#     AOSTOYAE    #
###################

aostoyae_families <- families[str_detect(families$gene_id, "aostoyae_"),]

####################
#     CCINEREA     #
####################

ccinerea_dictionary <- read_tsv("~/Desktop/alternative_splicing/ccinerea_AmutBmut/ccinerea_AmutBmut_dictionary.tsv", col_names = T, cols(fasta_transcript_id = col_character(),gene_id = col_character(),annotation_transcript_id = col_character()))
colnames(ccinerea_dictionary) <- c("gene_id", "gene_name", "transcript_id")

ccinerea_dictionary$gene_id <- ccinerea_dictionary$gene_id %>%
  str_replace("^", "ccinerea_")

ccinerea_dictionary$transcript_id  <- ccinerea_dictionary$transcript_id  %>%
  str_replace("^", "ccinerea_")

ccinerea_families <- families[str_detect(families$gene_id, "ccinerea_"),]

ccinerea_families <- ccinerea_families %>%
  left_join(ccinerea_dictionary, by = "gene_id") %>%
  select(family, gene_id = transcript_id)

###################
#     SCOMMUNE    #
###################

scommune_dictionary <- read_tsv("~/Desktop/alternative_splicing/scommune/scommune_dictionary.tsv", col_names = T, cols(fasta_transcript_id = col_character(),gene_id = col_character(),annotation_transcript_id = col_character()))
colnames(scommune_dictionary) <- c("gene_id", "gene_name", "transcript_id")

scommune_dictionary$gene_id <- scommune_dictionary$gene_id %>%
  str_replace("^", "scommune_")

scommune_dictionary$transcript_id  <- scommune_dictionary$transcript_id  %>%
  str_replace("^", "scommune_")

scommune_families <- families[str_detect(families$gene_id, "scommune_"),]

scommune_families <- scommune_families %>%
  left_join(scommune_dictionary, by = "gene_id") %>%
  select(family, gene_id = transcript_id)

####################
#     LTIGRINUS    #
####################

ltigrinus_dictionary <- read_tsv("~/Desktop/alternative_splicing/ltigrinus/ltigrinus_dictionary.tsv", col_names = T, cols(fasta_transcript_id = col_character(),gene_id = col_character(),annotation_transcript_id = col_character()))
colnames(ltigrinus_dictionary) <- c("gene_id", "gene_name", "transcript_id")

ltigrinus_dictionary$gene_id <- ltigrinus_dictionary$gene_id %>%
  str_replace("^", "ltigrinus_")

ltigrinus_dictionary$transcript_id  <- ltigrinus_dictionary$transcript_id  %>%
  str_replace("^", "ltigrinus_")

ltigrinus_families <- families[str_detect(families$gene_id, "ltigrinus_"),]

ltigrinus_families <- ltigrinus_families %>%
  left_join(ltigrinus_dictionary, by = "gene_id") %>%
  select(family, gene_id = transcript_id)

#################
#     AAMPLA    #
#################

aampla_dictionary <- read_tsv("~/Desktop/alternative_splicing/aampla/aampla_dictionary.tsv", col_names = T, cols(fasta_transcript_id = col_character(),gene_id = col_character(),annotation_transcript_id = col_character()))
colnames(aampla_dictionary) <- c("gene_id", "gene_name", "transcript_id")

aampla_dictionary$gene_id <- aampla_dictionary$gene_id %>%
  str_replace("^", "aampla_")

aampla_dictionary$transcript_id  <- aampla_dictionary$transcript_id  %>%
  str_replace("^", "aampla_")

aampla_families <- families[str_detect(families$gene_id, "aampla_"),]

aampla_families <- aampla_families %>%
  left_join(aampla_dictionary, by = "gene_id") %>%
  select(family, gene_id = transcript_id)

#####################
#     CNEOFORMANS   #
#####################

cneoformans_dictionary <- read_tsv("~/Desktop/alternative_splicing/cneoformans/cneoformans_dictionary.tsv", col_names = T, cols(fasta_transcript_id = col_character(),gene_id = col_character(),annotation_transcript_id = col_character()))
colnames(cneoformans_dictionary) <- c("gene_id", "gene_name", "transcript_id")

cneoformans_dictionary$gene_id <- cneoformans_dictionary$gene_id %>%
  str_replace("^", "cneoformans_")

cneoformans_dictionary$transcript_id  <- cneoformans_dictionary$transcript_id  %>%
  str_replace("^", "cneoformans_")

cneoformans_families <- families[str_detect(families$gene_id, "cneoformans_"),]

cneoformans_families <- cneoformans_families %>%
  left_join(cneoformans_dictionary, by = "gene_id") %>%
  select(family, gene_id = transcript_id)

#########################
#     PCHRYSOSPORIUM    #
#########################

pchrysosporium_dictionary <- read_tsv("~/Desktop/alternative_splicing/pchrysosporium/pchrysosporium_dictionary.tsv", col_names = T, cols(fasta_transcript_id = col_character(),gene_id = col_character(),annotation_transcript_id = col_character()))
colnames(pchrysosporium_dictionary) <- c("gene_id", "gene_name", "transcript_id")

pchrysosporium_dictionary$gene_id <- pchrysosporium_dictionary$gene_id %>%
  str_replace("^", "pchrysosporium_")

pchrysosporium_dictionary$transcript_id  <- pchrysosporium_dictionary$transcript_id  %>%
  str_replace("^", "pchrysosporium_")

pchrysosporium_families <- families[str_detect(families$gene_id, "pchrysosporium_"),]

pchrysosporium_families <- pchrysosporium_families %>%
  left_join(pchrysosporium_dictionary, by = "gene_id") %>%
  select(family, gene_id = transcript_id)

##################
#     RMELLEA    #
##################

rmellea_dictionary <- read_tsv("~/Desktop/alternative_splicing/rmellea/rmellea_dictionary.tsv", col_names = T, cols(fasta_transcript_id = col_character(),gene_id = col_character(),annotation_transcript_id = col_character()))
colnames(rmellea_dictionary) <- c("gene_id", "gene_name", "transcript_id")

rmellea_dictionary$gene_id <- rmellea_dictionary$gene_id %>%
  str_replace("^", "rmellea_")

rmellea_dictionary$transcript_id  <- rmellea_dictionary$transcript_id  %>%
  str_replace("^", "rmellea_")

rmellea_families <- families[str_detect(families$gene_id, "rmellea_"),]

rmellea_families <- rmellea_families %>%
  left_join(rmellea_dictionary, by = "gene_id") %>%
  select(family, gene_id = transcript_id)

##################
#     UMAYDIS    #
##################

umaydis_families <- families[str_detect(families$gene_id, "umaydis_"),]

################
#     MERGE    #
################

full <- aampla_families %>%
  bind_rows(aostoyae_families) %>%
  bind_rows(ccinerea_families) %>%
  bind_rows(cneoformans_families) %>%
  bind_rows(ltigrinus_families) %>%
  bind_rows(pchrysosporium_families) %>%
  bind_rows(rmellea_families) %>%
  bind_rows(scommune_families) %>%
  bind_rows(umaydis_families) %>%
  arrange(family)

write_tsv(full, "~/Desktop/alternative_splicing/enrichment/all_species_protein_HFX.fnodes.FIXED")




families %>% filter(family == "F_2_0")
full %>% filter(family == "F_2_0")
