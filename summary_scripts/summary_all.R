library(tidyr)
library(readr)
library(dplyr)
library(stringr)

setwd("~/Desktop/alternative_splicing/")

aampla <- read_tsv("aampla/aampla_summary.tsv")
aostoyae <- read_tsv("aostoyae/aostoyae_summary.tsv")
ccinerea_AmutBmut <- read_tsv("ccinerea_AmutBmut/ccinerea_AmutBmut_summary.tsv")
cneoformans <- read_tsv("cneoformans/cneoformans_summary.tsv")
ltigrinus <- read_tsv("ltigrinus/ltigrinus_summary.tsv")
pchrysosporium <- read_tsv("pchrysosporium/pchrysosporium_summary.tsv")
rmellea <- read_tsv("rmellea/rmellea_summary.tsv")
scommune <- read_tsv("scommune/scommune_summary.tsv")
umaydis <- read_tsv("umaydis/umaydis_summary.tsv")

full <- aampla %>%
  full_join(aostoyae) %>%
  full_join(ccinerea_AmutBmut) %>%
  full_join(cneoformans) %>%
  full_join(ltigrinus) %>%
  full_join(pchrysosporium) %>%
  full_join(rmellea) %>%
  full_join(scommune) %>%
  full_join(umaydis)

full <- full[c(1:18,29,19:28),]

write_tsv(full, "~/Desktop/summary.tsv")
