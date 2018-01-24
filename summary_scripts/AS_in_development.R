library(tidyr)
library(readr)
library(dplyr)
library(stringr)
library(tibble)

setwd("~/Desktop/alternative_splicing/")

##### ALAPTÁBLA #####

AS_in_development <- as.tibble(1:9)
AS_in_development$species <- NA
AS_in_development$ALL_genes <- NA
AS_in_development$RMA <- NA
AS_in_development$VM <- NA
AS_in_development$H <- NA
AS_in_development$P <- NA
AS_in_development$P1 <- NA
AS_in_development$P2 <- NA
AS_in_development$P2_C <- NA
AS_in_development$P2_S <- NA
AS_in_development$YFB <- NA
AS_in_development$YFB_C <- NA
AS_in_development$YFB_S <- NA
AS_in_development$YFB_L <- NA
AS_in_development$FB <- NA
AS_in_development$FB_CL <- NA
AS_in_development$FB_C <- NA
AS_in_development$FB_L <- NA
AS_in_development$FB_S <- NA
AS_in_development <- AS_in_development[-1]

##### 1. AAMPLA #####

results <- read_tsv("aampla/aampla_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/aampla.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% filter(VM_FPKM > 4 | P1_FPKM > 4 | P2_FPKM > 4 | YFB_FPKM > 4 | FB_FPKM > 4)

AS_in_development$species[1] <- "AAMPLA"
AS_in_development$ALL_genes[1] <- merged %>% nrow()
AS_in_development$VM[1] <- merged %>% filter(VM_FPKM > 4) %>% nrow()
AS_in_development$P1[1] <- merged %>% filter(P1_FPKM > 4) %>% nrow()
AS_in_development$P2[1] <- merged %>% filter(P2_FPKM > 4) %>% nrow()
AS_in_development$YFB[1] <- merged %>% filter(YFB_FPKM > 4) %>% nrow()
AS_in_development$FB[1] <- merged %>% filter(FB_FPKM > 4) %>% nrow()

rm(genes.fpkm, merged, results)

##### 2. AOSTOYAE #####

results <- read_tsv("aostoyae/aostoyae_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/aostoyae.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% 
  filter(RMA_FPKM > 4 | VM_FPKM > 4 | P1_FPKM > 4 | P2_C_FPKM > 4 | 
           P2_S_FPKM > 4 | YFB_C_FPKM > 4 | YFB_S_FPKM > 4 | FB_C_FPKM > 4 | 
           FB_L_FPKM > 4 | FB_S_FPKM > 4)

AS_in_development$species[2] <- "AOSTOYAE"
AS_in_development$ALL_genes[2] <- merged %>% nrow()
AS_in_development$RMA[2] <- merged %>% filter(RMA_FPKM > 4) %>% nrow()
AS_in_development$VM[2] <- merged %>% filter(VM_FPKM > 4) %>% nrow()
AS_in_development$P1[2] <- merged %>% filter(P1_FPKM > 4) %>% nrow()
AS_in_development$P2_C[2] <- merged %>% filter(P2_C_FPKM > 4) %>% nrow()
AS_in_development$P2_S[2] <- merged %>% filter(P2_S_FPKM > 4) %>% nrow()
AS_in_development$YFB_C[2] <- merged %>% filter(YFB_C_FPKM > 4) %>% nrow()
AS_in_development$YFB_S[2] <- merged %>% filter(YFB_S_FPKM > 4) %>% nrow()
AS_in_development$FB_C[2] <- merged %>% filter(FB_C_FPKM > 4) %>% nrow()
AS_in_development$FB_L[2] <- merged %>% filter(FB_L_FPKM > 4) %>% nrow()
AS_in_development$FB_S[2] <- merged %>% filter(FB_S_FPKM > 4) %>% nrow()

rm(genes.fpkm, merged, results)

##### 3. CCINEREA #####

results <- read_tsv("ccinerea_AmutBmut/ccinerea_AmutBmut_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/ccinerea.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% 
  filter(VM_FPKM > 4 | H_FPKM > 4 | P1_FPKM > 4 | 
           P2_FPKM > 4 | YFB_K_FPKM > 4 | YFB_L_FPKM > 4 | 
           YFB_T_FPKM > 4 | FB_KL_FPKM > 4 | FB_T_FPKM > 4)

AS_in_development$species[3] <- "CCINEREA"
AS_in_development$ALL_genes[3] <- merged %>% nrow()
AS_in_development$VM[3] <- merged %>% filter(VM_FPKM > 4) %>% nrow()
AS_in_development$H[3] <- merged %>% filter(H_FPKM > 4) %>% nrow()
AS_in_development$P1[3] <- merged %>% filter(P1_FPKM > 4) %>% nrow()
AS_in_development$P2[3] <- merged %>% filter(P2_FPKM > 4) %>% nrow()
AS_in_development$YFB_C[3] <- merged %>% filter(YFB_K_FPKM > 4) %>% nrow()
AS_in_development$YFB_L[3] <- merged %>% filter(YFB_L_FPKM > 4) %>% nrow()
AS_in_development$YFB_S[3] <- merged %>% filter(YFB_T_FPKM > 4) %>% nrow()
AS_in_development$FB_CL[3] <- merged %>% filter(FB_KL_FPKM > 4) %>% nrow()
AS_in_development$FB_S[3] <- merged %>% filter(FB_T_FPKM > 4) %>% nrow()
rm(genes.fpkm, merged, results)

##### 4. CNEOFORMANS #####

results <- read_tsv("cneoformans/cneoformans_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/cneoformans.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% 
  filter(FPKM > 4)

AS_in_development$species[4] <- "CNEOFORMANS"
AS_in_development$ALL_genes[4] <- merged %>% nrow()
rm(genes.fpkm, merged, results)

##### 5. LTIGRINUS #####

results <- read_tsv("ltigrinus/ltigrinus_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/ltigrinus.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% 
  filter(VM_FPKM > 4 | P1_FPKM > 4 | P2_K_FPKM > 4 | 
           P2_T_FPKM > 4 | YFB_K_FPKM > 4 | YFB_T_FPKM > 4 | 
           FB_K_FPKM > 4 | FB_L_FPKM > 4 | FB_T_FPKM > 4)

AS_in_development$species[5] <- "LTIGRINUS"
AS_in_development$ALL_genes[5] <- merged %>% nrow()
AS_in_development$VM[5] <- merged %>% filter(VM_FPKM > 4) %>% nrow()
AS_in_development$P1[5] <- merged %>% filter(P1_FPKM > 4) %>% nrow()
AS_in_development$P2_C[5] <- merged %>% filter(P2_K_FPKM > 4) %>% nrow()
AS_in_development$P2_S[5] <- merged %>% filter(P2_T_FPKM > 4) %>% nrow()
AS_in_development$YFB_C[5] <- merged %>% filter(YFB_K_FPKM > 4) %>% nrow()
AS_in_development$YFB_S[5] <- merged %>% filter(YFB_T_FPKM > 4) %>% nrow()
AS_in_development$FB_C[5] <- merged %>% filter(FB_K_FPKM > 4) %>% nrow()
AS_in_development$FB_L[5] <- merged %>% filter(FB_L_FPKM > 4) %>% nrow()
AS_in_development$FB_S[5] <- merged %>% filter(FB_T_FPKM > 4) %>% nrow()
rm(genes.fpkm, merged, results)

##### 6. PCHRYSOSPORIUM #####

results <- read_tsv("pchrysosporium/pchrysosporium_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/pchrysosporium.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% 
  filter(VM_FPKM > 4 | YFB_FPKM > 4 | FB_FPKM > 4)

AS_in_development$species[6] <- "PCHRYSOSPORIUM"
AS_in_development$ALL_genes[6] <- merged %>% nrow()
AS_in_development$VM[6] <- merged %>% filter(VM_FPKM > 4) %>% nrow()
AS_in_development$YFB[6] <- merged %>% filter(YFB_FPKM > 4) %>% nrow()
AS_in_development$FB[6] <- merged %>% filter(FB_FPKM > 4) %>% nrow()
rm(genes.fpkm, merged, results)

##### 7. RMELLEA #####

results <- read_tsv("rmellea/rmellea_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/rmellea.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% 
  filter(VM_FPKM > 4 | P_FPKM > 4 | YFB_FPKM > 4 | 
           FB_K_FPKM > 4 | FB_T_FPKM > 4)

AS_in_development$species[7] <- "RMELLEA"
AS_in_development$ALL_genes[7] <- merged %>% nrow()
AS_in_development$VM[7] <- merged %>% filter(VM_FPKM > 4) %>% nrow()
AS_in_development$P[7] <- merged %>% filter(P_FPKM > 4) %>% nrow()
AS_in_development$YFB[7] <- merged %>% filter(YFB_FPKM > 4) %>% nrow()
AS_in_development$FB_C[7] <- merged %>% filter(FB_K_FPKM > 4) %>% nrow()
AS_in_development$FB_S[7] <- merged %>% filter(FB_T_FPKM > 4) %>% nrow()
rm(genes.fpkm, merged, results)

##### 8. SCOMMUNE #####

results <- read_tsv("scommune/scommune_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/scommune.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% 
  filter(VM_FPKM > 4 | P1_FPKM > 4 | P2_FPKM > 4 | 
           YFB_FPKM > 4 | FB_FPKM > 4)

AS_in_development$species[8] <- "SCOMMUNE"
AS_in_development$ALL_genes[8] <- merged %>% nrow()
AS_in_development$VM[8] <- merged %>% filter(VM_FPKM > 4) %>% nrow()
AS_in_development$P1[8] <- merged %>% filter(P1_FPKM > 4) %>% nrow()
AS_in_development$P2[8] <- merged %>% filter(P2_FPKM > 4) %>% nrow()
AS_in_development$YFB[8] <- merged %>% filter(YFB_FPKM > 4) %>% nrow()
AS_in_development$FB[8] <- merged %>% filter(FB_FPKM > 4) %>% nrow()
rm(genes.fpkm, merged, results)

##### 9. UMAYDIS #####

results <- read_tsv("umaydis/umaydis_results.tsv", col_names = T, cols(.default = col_guess(), transcript_id = col_character(), gene_id = col_character())) %>% select(gene_id) %>% unique()
genes.fpkm <- read_tsv("FPKM_tables/umaydis.genes.fpkm.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

merged <- left_join(results, genes.fpkm) %>% 
  filter(FPKM > 4)

AS_in_development$species[9] <- "UMAYDIS"
AS_in_development$ALL_genes[9] <- merged %>% nrow()
rm(genes.fpkm, merged, results)

##### KIÍRÁS #####

write_tsv(AS_in_development, "AS_in_development.tsv")
