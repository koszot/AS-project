library(readr)
library(dplyr)
library(stringr)
library(tidyr)

setwd("~/Desktop/alternative_splicing/")

######################
#     GO SUMMARY     #
######################

aampla <- read_tsv("aampla/aampla_enrichment/aampla_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
aampla <- filter(aampla, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
aampla$count <- 1

aostoyae <- read_tsv("aostoyae/aostoyae_enrichment/aostoyae_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
aostoyae <- filter(aostoyae, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
aostoyae$count <- 1

ccinerea <- read_tsv("ccinerea_AmutBmut/ccinerea_AmutBmut_enrichment/ccinerea_AmutBmut_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
ccinerea <- filter(ccinerea, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
ccinerea$count <- 1

cneoformans <- read_tsv("cneoformans/cneoformans_enrichment/cneoformans_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
cneoformans <- filter(cneoformans, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
cneoformans$count <- 1

ltigrinus <- read_tsv("ltigrinus/ltigrinus_enrichment/ltigrinus_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
ltigrinus <- filter(ltigrinus, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
ltigrinus$count <- 1

pchrysosporium <- read_tsv("pchrysosporium/pchrysosporium_enrichment/pchrysosporium_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
pchrysosporium <- filter(pchrysosporium, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
pchrysosporium$count <- 1

rmellea <- read_tsv("rmellea/rmellea_enrichment/rmellea_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
rmellea <- filter(rmellea, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
rmellea$count <- 1

scommune <- read_tsv("scommune/scommune_enrichment/scommune_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
scommune <- filter(scommune, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
scommune$count <- 1

umaydis <- read_tsv("umaydis/umaydis_enrichment/umaydis_genelvl_enrichment_fisher_classic.tsv") %>%                # fájl betöltése
  select(type2, type, GO.ID, Term)
umaydis <- filter(umaydis, type2 == "AS") %>% select(GO.ID, Term)                                                # AS vagy DEVREG
umaydis$count <- 1

GO <- full_join(aampla, aostoyae, by = c("GO.ID", "Term")) %>% 
  full_join(ccinerea, by = c("GO.ID", "Term")) %>% 
  full_join(cneoformans, by = c("GO.ID", "Term")) %>% 
  full_join(ltigrinus, by = c("GO.ID", "Term")) %>% 
  full_join(pchrysosporium, by = c("GO.ID", "Term")) %>% 
  full_join(rmellea, by = c("GO.ID", "Term")) %>% 
  full_join(scommune, by = c("GO.ID", "Term")) %>% 
  full_join(umaydis, by = c("GO.ID", "Term")) %>% 
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
colnames(GO) <- c("GO_ID", "Term","aampla","aostoyae", "ccinerea","cneoformans", 
                       "ltigrinus","pchrysosporium","rmellea","scommune","umaydis","ALL")

write_tsv(GO, "enrichment/GO_counts_AS_classic.tsv")
