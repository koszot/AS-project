library(readr)
library(dplyr)
library(stringr)
library(tidyr)

setwd("~/Desktop/alternative_splicing/")

###################
#     CLASSIC     #
###################

# aampla

species <- read_tsv("aampla/aampla_enrichment/aampla_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
aampla_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

# aostoyae

species <- read_tsv("aostoyae/aostoyae_enrichment/aostoyae_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
aostoyae_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

# ccinerea

species <- read_tsv("ccinerea_AmutBmut/ccinerea_AmutBmut_enrichment/ccinerea_AmutBmut_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
ccinerea_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

# cneoformans

species <- read_tsv("cneoformans/cneoformans_enrichment/cneoformans_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
cneoformans_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

# ltigrinus

species <- read_tsv("ltigrinus/ltigrinus_enrichment/ltigrinus_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
ltigrinus_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

# pchrysosporium

species <- read_tsv("pchrysosporium/pchrysosporium_enrichment/pchrysosporium_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
pchrysosporium_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

# rmellea

species <- read_tsv("rmellea/rmellea_enrichment/rmellea_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
rmellea_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

# scommune

species <- read_tsv("scommune/scommune_enrichment/scommune_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
scommune_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

# umaydis

species <- read_tsv("umaydis/umaydis_enrichment/umaydis_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
umaydis_AS_classic <- species                                                                       # átmentés
rm(species_table, species)

###############################
#     CLASSIC  SIGNIFICANT    #
###############################

# aampla

species <- read_tsv("aampla/aampla_enrichment/aampla_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
aampla_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

# aostoyae

species <- read_tsv("aostoyae/aostoyae_enrichment/aostoyae_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
aostoyae_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

# ccinerea

species <- read_tsv("ccinerea_AmutBmut/ccinerea_AmutBmut_enrichment/ccinerea_AmutBmut_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
ccinerea_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

# cneoformans

species <- read_tsv("cneoformans/cneoformans_enrichment/cneoformans_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
cneoformans_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

# ltigrinus

species <- read_tsv("ltigrinus/ltigrinus_enrichment/ltigrinus_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
ltigrinus_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

# pchrysosporium

species <- read_tsv("pchrysosporium/pchrysosporium_enrichment/pchrysosporium_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
pchrysosporium_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

# rmellea

species <- read_tsv("rmellea/rmellea_enrichment/rmellea_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
rmellea_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

# scommune

species <- read_tsv("scommune/scommune_enrichment/scommune_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
scommune_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

# umaydis

species <- read_tsv("umaydis/umaydis_enrichment/umaydis_genelvl_enrichment_fisher_classic_significant.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, families)

species <- filter(species, type2 == "AS")                                                                             # AS vagy DEVREG
species <- species$families  %>% str_split(", ")  %>% unlist() %>% str_replace("/.*$", "")
species_table <- tibble(family_id = "NA", count = 1)
for (x in 1:length(unique(species))) {
  species_table[x,1] <- species[x]
  species_table[x,2] <- length(species[species==species[x]])
}
species <- distinct(species_table, family_id, count) %>% arrange(desc(count))
umaydis_AS_classic_significant <- species                                                                       # átmentés
rm(species_table, species)

###################
#     CLASSIC     #
###################

classic <- full_join(aampla_AS_classic, aostoyae_AS_classic, by = "family_id") %>% 
  full_join(ccinerea_AS_classic, by = "family_id") %>% 
  full_join(cneoformans_AS_classic, by = "family_id") %>% 
  full_join(ltigrinus_AS_classic, by = "family_id") %>% 
  full_join(pchrysosporium_AS_classic, by = "family_id") %>% 
  full_join(rmellea_AS_classic, by = "family_id") %>% 
  full_join(scommune_AS_classic, by = "family_id") %>% 
  full_join(umaydis_AS_classic, by = "family_id") %>% 
  replace_na(list(count.x = 0)) %>%
  replace_na(list(count.y = 0)) %>%
  replace_na(list(count.x.x = 0)) %>%
  replace_na(list(count.y.y = 0)) %>%
  replace_na(list(count.x.x.x = 0)) %>%
  replace_na(list(count.y.y.y = 0)) %>%
  replace_na(list(count.x.x.x.x = 0)) %>%
  replace_na(list(count.y.y.y.y = 0)) %>%
  replace_na(list(count = 0)) %>%
  mutate(count_FULL = (count.x + count.y + count.x.x + count.y.y + count.x.x.x + count.y.y.y + count.x.x.x.x + count.y.y.y.y + count )) %>%
  arrange(desc(count_FULL))
colnames(classic) <- c("Gene_Family", "aampla_AS_classic","aostoyae_AS_classic", "ccinerea_AS_classic","cneoformans_AS_classic", 
                       "ltigrinus_AS_classic","pchrysosporium_AS_classic","rmellea_AS_classic","scommune_AS_classic","umaydis_AS_classic","ALL")

write_tsv(classic, "enrichment/family_counts_AS_classic.tsv")

###############################
#     CLASSIC  SIGNIFICANT    #
###############################

classic <- full_join(aampla_AS_classic_significant, aostoyae_AS_classic_significant, by = "family_id") %>% 
  full_join(ccinerea_AS_classic_significant, by = "family_id") %>% 
  full_join(cneoformans_AS_classic_significant, by = "family_id") %>% 
  full_join(ltigrinus_AS_classic_significant, by = "family_id") %>% 
  full_join(pchrysosporium_AS_classic_significant, by = "family_id") %>% 
  full_join(rmellea_AS_classic_significant, by = "family_id") %>% 
  full_join(scommune_AS_classic_significant, by = "family_id") %>% 
  full_join(umaydis_AS_classic_significant, by = "family_id") %>% 
  replace_na(list(count.x = 0)) %>%
  replace_na(list(count.y = 0)) %>%
  replace_na(list(count.x.x = 0)) %>%
  replace_na(list(count.y.y = 0)) %>%
  replace_na(list(count.x.x.x = 0)) %>%
  replace_na(list(count.y.y.y = 0)) %>%
  replace_na(list(count.x.x.x.x = 0)) %>%
  replace_na(list(count.y.y.y.y = 0)) %>%
  replace_na(list(count = 0)) %>%
  mutate(count_FULL = (count.x + count.y + count.x.x + count.y.y + count.x.x.x + count.y.y.y + count.x.x.x.x + count.y.y.y.y + count )) %>%
  arrange(desc(count_FULL))
colnames(classic) <- c("Gene_Family", "aampla_AS_classic_significant","aostoyae_AS_classic_significant", "ccinerea_AS_classic_significant","cneoformans_AS_classic_significant", 
                       "ltigrinus_AS_classic_significant","pchrysosporium_AS_classic_significant","rmellea_AS_classic_significant","scommune_AS_classic_significant","umaydis_AS_classic_significant","ALL")

write_tsv(classic, "enrichment/family_counts_AS_classic_significant.tsv")
