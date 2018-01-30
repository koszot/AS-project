library(topGO)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

#####################
#     LOAD FILES    #
#####################

setwd("~/Desktop/alternative_splicing/cneoformans/")

geneID2GO <- readMappings(file = "cneoformans_GO_dict.tsv")  

geneUniverse <- read_tsv("cneoformans_enrichment/cneoformans_genelvl_enrichment_geneUniverse.tsv", col_names = T, cols(value = col_character()))
geneUniverse <- as.character(geneUniverse$value)

AS.transcripts <- read_tsv("cneoformans_enrichment/cneoformans_genelvl_enrichment_AStranscripts.tsv", col_names = T, cols(value = col_character()))
AS.transcripts <- as.character(AS.transcripts$value)
geneList.AS <- factor(as.integer(geneUniverse %in% AS.transcripts))
names(geneList.AS) <- geneUniverse

#####################
#     ALL VS AS     #
#####################

# Biological Process 

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList.AS,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
mysummary <- summary(attributes(resultClassic)$score <= 0.05)
numsignif <- as.integer(mysummary[[3]]) 
numsignif
allRes_all_vs_AS_BP <- GenTable(myGOdata, classicFisher = resultClassic, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes_all_vs_AS_BP
allRes_all_vs_AS_BP_significant <- allRes_all_vs_AS_BP

myterms <- allRes_all_vs_AS_BP$GO.ID
all_vs_AS_BP <- genesInTerm(myGOdata, myterms)

all_vs_AS_BP_significant <- list()
significantgenes <- as_tibble(sigGenes(myGOdata))
for (x in 1:length(all_vs_AS_BP)) {
  temp <- as_tibble(all_vs_AS_BP[[x]])
  all_vs_AS_BP_significant[[x]] <- intersect(significantgenes, temp)$value
}
names(all_vs_AS_BP_significant) <- names(all_vs_AS_BP)

# Molecular Function

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList.AS,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
mysummary <- summary(attributes(resultClassic)$score <= 0.05)
numsignif <- as.integer(mysummary[[3]])
numsignif
allRes_all_vs_AS_MF <- GenTable(myGOdata, classicFisher = resultClassic, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes_all_vs_AS_MF
allRes_all_vs_AS_MF_significant <- allRes_all_vs_AS_MF

myterms <- allRes_all_vs_AS_MF$GO.ID
all_vs_AS_MF <- genesInTerm(myGOdata, myterms)

all_vs_AS_MF_significant <- list()
significantgenes <- as_tibble(sigGenes(myGOdata))
for (x in 1:length(all_vs_AS_MF)) {
  temp <- as_tibble(all_vs_AS_MF[[x]])
  all_vs_AS_MF_significant[[x]] <- intersect(significantgenes, temp)$value
}
names(all_vs_AS_MF_significant) <- names(all_vs_AS_MF)

# Cellular Component

myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList.AS,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
mysummary <- summary(attributes(resultClassic)$score <= 0.05)
numsignif <- as.integer(mysummary[[3]])
numsignif
allRes_all_vs_AS_CC <- GenTable(myGOdata, classicFisher = resultClassic, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes_all_vs_AS_CC
allRes_all_vs_AS_CC_significant <- allRes_all_vs_AS_CC

myterms <- allRes_all_vs_AS_CC$GO.ID
all_vs_AS_CC <- genesInTerm(myGOdata, myterms)

all_vs_AS_CC_significant <- list()
significantgenes <- as_tibble(sigGenes(myGOdata))
for (x in 1:length(all_vs_AS_CC)) {
  temp <- as_tibble(all_vs_AS_CC[[x]])
  all_vs_AS_CC_significant[[x]] <- intersect(significantgenes, temp)$value
}
names(all_vs_AS_CC_significant) <- names(all_vs_AS_CC)

##################################
#     GENE FAMILITES !!!FULL!!!  #
##################################

families <- read_tsv("../enrichment/all_species_protein_HFX.fnodes.FIXED", col_names = c("family", "gene_id"))

# all_vs_AS_BP

length(all_vs_AS_BP)
for (x in 1:length(all_vs_AS_BP)) {
  all_vs_AS_BP[[x]] <- str_replace(all_vs_AS_BP[[x]], "^", "cneoformans_")
}

for (x in 1:length(allRes_all_vs_AS_BP$GO.ID)) {
  temp <- all_vs_AS_BP[allRes_all_vs_AS_BP$GO.ID[x]]
  temp <- as_tibble(unlist(temp))
  colnames(temp) <- "gene_id"
  temp <- temp %>% 
    left_join(families, by = "gene_id") %>%
    select(family)
  counted <- as_tibble(table(temp$family))
  counted <- arrange(counted, desc(n))
  counted_2 <- counted
  counted <- unite(counted, Var1, n, col = "full", sep = "/")
  counted_2 <- unite(counted_2, Var1, n, col = "full", sep = "/")
  counted <- as.list(counted)
  counted_2 <- as.list(counted_2)
  allRes_all_vs_AS_BP$families_count[[x]] <- length(counted[[1]])
  allRes_all_vs_AS_BP$families[[x]] <- paste(counted_2[[1]], collapse = ", ")
}

# all_vs_AS_MF

length(all_vs_AS_MF)
for (x in 1:length(all_vs_AS_MF)) {
  all_vs_AS_MF[[x]] <- str_replace(all_vs_AS_MF[[x]], "^", "cneoformans_")
}

for (x in 1:length(allRes_all_vs_AS_MF$GO.ID)) {
  temp <- all_vs_AS_MF[allRes_all_vs_AS_MF$GO.ID[x]]
  temp <- as_tibble(unlist(temp))
  colnames(temp) <- "gene_id"
  temp <- temp %>% 
    left_join(families, by = "gene_id") %>%
    select(family)
  counted <- as_tibble(table(temp$family))
  counted <- arrange(counted, desc(n))
  counted_2 <- counted
  counted <- unite(counted, Var1, n, col = "full", sep = "/")
  counted_2 <- unite(counted_2, Var1, n, col = "full", sep = "/")
  counted <- as.list(counted)
  counted_2 <- as.list(counted_2)
  allRes_all_vs_AS_MF$families_count[[x]] <- length(counted[[1]])
  allRes_all_vs_AS_MF$families[[x]] <- paste(counted_2[[1]], collapse = ", ")
}

# all_vs_AS_CC

length(all_vs_AS_CC)
for (x in 1:length(all_vs_AS_CC)) {
  all_vs_AS_CC[[x]] <- str_replace(all_vs_AS_CC[[x]], "^", "cneoformans_")
}

for (x in 1:length(allRes_all_vs_AS_CC$GO.ID)) {
  temp <- all_vs_AS_CC[allRes_all_vs_AS_CC$GO.ID[x]]
  temp <- as_tibble(unlist(temp))
  colnames(temp) <- "gene_id"
  temp <- temp %>% 
    left_join(families, by = "gene_id") %>%
    select(family)
  counted <- as_tibble(table(temp$family))
  counted <- arrange(counted, desc(n))
  counted_2 <- counted
  counted <- unite(counted, Var1, n, col = "full", sep = "/")
  counted_2 <- unite(counted_2, Var1, n, col = "full", sep = "/")
  counted <- as.list(counted)
  counted_2 <- as.list(counted_2)
  allRes_all_vs_AS_CC$families_count[[x]] <- length(counted[[1]])
  allRes_all_vs_AS_CC$families[[x]] <- paste(counted_2[[1]], collapse = ", ")
}

allRes_all_vs_AS_BP$type <- "Biological Process"
allRes_all_vs_AS_BP <- allRes_all_vs_AS_BP[,c(9, 1:8)]
allRes_all_vs_AS_MF$type <- "Molecular Function"
allRes_all_vs_AS_MF <- allRes_all_vs_AS_MF[,c(9, 1:8)]
allRes_all_vs_AS_CC$type <- "Cellular Component"
allRes_all_vs_AS_CC <- allRes_all_vs_AS_CC[,c(9, 1:8)]
allRes_all_vs_AS <- allRes_all_vs_AS_BP %>%
  bind_rows(allRes_all_vs_AS_MF) %>%
  bind_rows(allRes_all_vs_AS_CC)

# final table

allRes_all_vs_AS$type2 <- "AS"
allRes_all_vs_AS <- allRes_all_vs_AS[,c(10, 1:9)]

final <- allRes_all_vs_AS

write_tsv(final, "cneoformans_enrichment/cneoformans_genelvl_enrichment_fisher_classic.tsv")

rm(allRes_all_vs_AS, allRes_all_vs_AS_BP, allRes_all_vs_AS_CC, allRes_all_vs_AS_MF)
rm(all_vs_AS_BP, all_vs_AS_CC, all_vs_AS_MF, temp)

#########################################
#     GENE FAMILITES !!!SIGNIFICANT!!!  #
#########################################

# all_vs_AS_BP_significant

length(all_vs_AS_BP_significant)
for (x in 1:length(all_vs_AS_BP_significant)) {
  all_vs_AS_BP_significant[[x]] <- str_replace(all_vs_AS_BP_significant[[x]], "^", "cneoformans_")
}

for (x in 1:length(allRes_all_vs_AS_BP_significant$GO.ID)) {
  temp <- all_vs_AS_BP_significant[allRes_all_vs_AS_BP_significant$GO.ID[x]]
  temp <- as_tibble(unlist(temp))
  colnames(temp) <- "gene_id"
  temp <- temp %>% 
    left_join(families, by = "gene_id") %>%
    select(family)
  counted <- as_tibble(table(temp$family))
  counted <- arrange(counted, desc(n))
  counted_2 <- counted
  counted <- unite(counted, Var1, n, col = "full", sep = "/")
  counted_2 <- unite(counted_2, Var1, n, col = "full", sep = "/")
  counted <- as.list(counted)
  counted_2 <- as.list(counted_2)
  allRes_all_vs_AS_BP_significant$families_count[[x]] <- length(counted[[1]])
  allRes_all_vs_AS_BP_significant$families[[x]] <- paste(counted_2[[1]], collapse = ", ")
}

# all_vs_AS_MF

length(all_vs_AS_MF_significant)
for (x in 1:length(all_vs_AS_MF_significant)) {
  all_vs_AS_MF_significant[[x]] <- str_replace(all_vs_AS_MF_significant[[x]], "^", "cneoformans_")
}

for (x in 1:length(allRes_all_vs_AS_MF_significant$GO.ID)) {
  temp <- all_vs_AS_MF_significant[allRes_all_vs_AS_MF_significant$GO.ID[x]]
  temp <- as_tibble(unlist(temp))
  colnames(temp) <- "gene_id"
  temp <- temp %>% 
    left_join(families, by = "gene_id") %>%
    select(family)
  counted <- as_tibble(table(temp$family))
  counted <- arrange(counted, desc(n))
  counted_2 <- counted
  counted <- unite(counted, Var1, n, col = "full", sep = "/")
  counted_2 <- unite(counted_2, Var1, n, col = "full", sep = "/")
  counted <- as.list(counted)
  counted_2 <- as.list(counted_2)
  allRes_all_vs_AS_MF_significant$families_count[[x]] <- length(counted[[1]])
  allRes_all_vs_AS_MF_significant$families[[x]] <- paste(counted_2[[1]], collapse = ", ")
}

# all_vs_AS_CC

length(all_vs_AS_CC_significant)
for (x in 1:length(all_vs_AS_CC_significant)) {
  all_vs_AS_CC_significant[[x]] <- str_replace(all_vs_AS_CC_significant[[x]], "^", "cneoformans_")
}

for (x in 1:length(allRes_all_vs_AS_CC_significant$GO.ID)) {
  temp <- all_vs_AS_CC_significant[allRes_all_vs_AS_CC_significant$GO.ID[x]]
  temp <- as_tibble(unlist(temp))
  colnames(temp) <- "gene_id"
  temp <- temp %>% 
    left_join(families, by = "gene_id") %>%
    select(family)
  counted <- as_tibble(table(temp$family))
  counted <- arrange(counted, desc(n))
  counted_2 <- counted
  counted <- unite(counted, Var1, n, col = "full", sep = "/")
  counted_2 <- unite(counted_2, Var1, n, col = "full", sep = "/")
  counted <- as.list(counted)
  counted_2 <- as.list(counted_2)
  allRes_all_vs_AS_CC_significant$families_count[[x]] <- length(counted[[1]])
  allRes_all_vs_AS_CC_significant$families[[x]] <- paste(counted_2[[1]], collapse = ", ")
}

allRes_all_vs_AS_BP_significant$type <- "Biological Process"
allRes_all_vs_AS_BP_significant <- allRes_all_vs_AS_BP_significant[,c(9, 1:8)]
allRes_all_vs_AS_MF_significant$type <- "Molecular Function"
allRes_all_vs_AS_MF_significant <- allRes_all_vs_AS_MF_significant[,c(9, 1:8)]
allRes_all_vs_AS_CC_significant$type <- "Cellular Component"
allRes_all_vs_AS_CC_significant <- allRes_all_vs_AS_CC_significant[,c(9, 1:8)]
allRes_all_vs_AS_significant <- allRes_all_vs_AS_BP_significant %>%
  bind_rows(allRes_all_vs_AS_MF_significant) %>%
  bind_rows(allRes_all_vs_AS_CC_significant)

# final table

allRes_all_vs_AS_significant$type2 <- "AS"
allRes_all_vs_AS_significant <- allRes_all_vs_AS_significant[,c(10, 1:9)]

final_significant <- allRes_all_vs_AS_significant

write_tsv(final_significant, "cneoformans_enrichment/cneoformans_genelvl_enrichment_fisher_classic_significant.tsv")
