library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# genes input c("aostoyae", "ccinerea", "ltigrnius", "scommune")

interested_genes <- c("AROS_07204", "367962", "521662", "2572981")

##### LOAD ANNOTATION FILES #####

aostoyae <- read_tsv("~/Desktop/alternative_splicing/aostoyae/aostoyae_genome/aostoyae_AS_annotation.gtf", col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ") %>%
  select(geneID, transcriptID, strand, start, end)
aostoyae$transcriptID <- aostoyae$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")
aostoyae$geneID <- aostoyae$geneID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")

ccinerea <- read_tsv("~/Desktop/alternative_splicing/ccinerea_AmutBmut/ccinerea_AmutBmut_genome/ccinerea_AmutBmut_AS_annotation.gtf", col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ") %>%
  select(geneID, transcriptID, strand, start, end)
ccinerea$transcriptID <- ccinerea$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")
ccinerea$geneID <- ccinerea$geneID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")

ltigrinus <- read_tsv("~/Desktop/alternative_splicing/ltigrinus/ltigrinus_genome/ltigrinus_AS_annotation.gtf", col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ") %>%
  select(geneID, transcriptID, strand, start, end)
ltigrinus$transcriptID <- ltigrinus$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")
ltigrinus$geneID <- ltigrinus$geneID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")

scommune <- read_tsv("~/Desktop/alternative_splicing/scommune/scommune_genome/scommune_AS_annotation.gtf", col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ") %>%
  select(geneID, transcriptID, strand, start, end)
scommune$transcriptID <- scommune$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")
scommune$geneID <- scommune$geneID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")

##### EXONS AND INTRONS #####

# armillaria

aostoyae_gene <- aostoyae[str_detect(aostoyae$geneID, interested_genes[1]),]
aostoyae_exons <- list()
aostoyae_introns <- list()

for (x in 1:length(unique(aostoyae_gene$transcriptID))) {
  temp <- aostoyae_gene[str_detect(aostoyae_gene$transcriptID, unique(aostoyae_gene$transcriptID)[x]),]
  aostoyae_exons[[x]] <- as.integer(temp$end - temp$start + 1)
  temp_vector <- list()
  for (y in 1:length(temp$transcriptID)) {
    temp_vector[y] <- temp[y+1,"start"] - temp[y,"end"] - 1
  }
  aostoyae_introns[[x]] <- as.integer(temp_vector)
  aostoyae_introns[[x]] <- aostoyae_introns[[x]][1:(length(aostoyae_introns[[x]])-1)]
}

if (unique(aostoyae_gene$strand) == "-") {
  for (x in 1:length(aostoyae_exons)) {
    aostoyae_exons[[x]] <- aostoyae_exons[[x]][length(aostoyae_exons[[x]]):1]
    aostoyae_introns[[x]] <- aostoyae_introns[[x]][length(aostoyae_introns[[x]]):1]
  }
}

names(aostoyae_exons) <- unique(aostoyae_gene$transcriptID)
names(aostoyae_introns) <- unique(aostoyae_gene$transcriptID)

# ccinerea

ccinerea_gene <- ccinerea[str_detect(ccinerea$geneID, interested_genes[2]),] 
ccinerea_exons <- list()
ccinerea_introns <- list()

for (x in 1:length(unique(ccinerea_gene$transcriptID))) {
  temp <- ccinerea_gene[str_detect(ccinerea_gene$transcriptID, unique(ccinerea_gene$transcriptID)[x]),]
  ccinerea_exons[[x]] <- as.integer(temp$end - temp$start + 1)
  temp_vector <- list()
  for (y in 1:length(temp$transcriptID)) {
    temp_vector[y] <- temp[y+1,"start"] - temp[y,"end"] - 1
  }
  ccinerea_introns[[x]] <- as.integer(temp_vector)
  ccinerea_introns[[x]] <- ccinerea_introns[[x]][1:(length(ccinerea_introns[[x]])-1)]
}

if (unique(ccinerea_gene$strand) == "-") {
  for (x in 1:length(ccinerea_exons)) {
    ccinerea_exons[[x]] <- ccinerea_exons[[x]][length(ccinerea_exons[[x]]):1]
    ccinerea_introns[[x]] <- ccinerea_introns[[x]][length(ccinerea_introns[[x]]):1]
  }
}

names(ccinerea_exons) <- unique(ccinerea_gene$transcriptID)
names(ccinerea_introns) <- unique(ccinerea_gene$transcriptID)

# ltigrinus

ltigrinus_gene <- ltigrinus[str_detect(ltigrinus$geneID, interested_genes[3]),] 
ltigrinus_exons <- list()
ltigrinus_introns <- list()

for (x in 1:length(unique(ltigrinus_gene$transcriptID))) {
  temp <- ltigrinus_gene[str_detect(ltigrinus_gene$transcriptID, unique(ltigrinus_gene$transcriptID)[x]),]
  ltigrinus_exons[[x]] <- as.integer(temp$end - temp$start + 1)
  temp_vector <- list()
  for (y in 1:length(temp$transcriptID)) {
    temp_vector[y] <- temp[y+1,"start"] - temp[y,"end"] - 1
  }
  ltigrinus_introns[[x]] <- as.integer(temp_vector)
  ltigrinus_introns[[x]] <- ltigrinus_introns[[x]][1:(length(ltigrinus_introns[[x]])-1)]
}

if (unique(ltigrinus_gene$strand) == "-") {
  for (x in 1:length(ltigrinus_exons)) {
    ltigrinus_exons[[x]] <- ltigrinus_exons[[x]][length(ltigrinus_exons[[x]]):1]
    ltigrinus_introns[[x]] <- ltigrinus_introns[[x]][length(ltigrinus_introns[[x]]):1]
  }
}

names(ltigrinus_exons) <- unique(ltigrinus_gene$transcriptID)
names(ltigrinus_introns) <- unique(ltigrinus_gene$transcriptID)

# scommune

scommune_gene <- scommune[str_detect(scommune$geneID, interested_genes[4]),] 
scommune_exons <- list()
scommune_introns <- list()

for (x in 1:length(unique(scommune_gene$transcriptID))) {
  temp <- scommune_gene[str_detect(scommune_gene$transcriptID, unique(scommune_gene$transcriptID)[x]),]
  scommune_exons[[x]] <- as.integer(temp$end - temp$start + 1)
  temp_vector <- list()
  for (y in 1:length(temp$transcriptID)) {
    temp_vector[y] <- temp[y+1,"start"] - temp[y,"end"] - 1
  }
  scommune_introns[[x]] <- as.integer(temp_vector)
  scommune_introns[[x]] <- scommune_introns[[x]][1:(length(scommune_introns[[x]])-1)]
}

if (unique(scommune_gene$strand) == "-") {
  for (x in 1:length(scommune_exons)) {
    scommune_exons[[x]] <- scommune_exons[[x]][length(scommune_exons[[x]]):1]
    scommune_introns[[x]] <- scommune_introns[[x]][length(scommune_introns[[x]]):1]
  }
}

names(scommune_exons) <- unique(scommune_gene$transcriptID)
names(scommune_introns) <- unique(scommune_gene$transcriptID)

##### CDS REGIONS #####

# aostoyae

aostoyae_filter <- read_tsv("~/Desktop/alternative_splicing/aostoyae/aostoyae_transcripts.fasta.transdecoder_dir/filter", col_names = F)
aostoyae_CDS <- read_tsv("~/Desktop/alternative_splicing/aostoyae/aostoyae_transcripts.fasta.transdecoder_dir/headers.cds", col_names = F) %>%
  inner_join(aostoyae_filter, by = "X1") %>%
  select(X6) %>%
  separate(X6, c("transcriptID", "X6.2"), ":") %>%
  separate(X6.2, c("start", "X6.2.2"), "-") %>%
  separate(X6.2.2, c("end", "strand"), "\\(") %>%
  select(transcriptID, start, end)
aostoyae_CDS$start <- as.integer(aostoyae_CDS$start)
aostoyae_CDS$end <- as.integer(aostoyae_CDS$end)
aostoyae_CDS <- mutate(aostoyae_CDS, length = end - start + 1)

aostoyae_cdslengths <- aostoyae_exons

for (actualT in 1:length(aostoyae_exons)) {
  temp_1 <- aostoyae_CDS[str_detect(aostoyae_CDS$transcriptID, names(aostoyae_exons)[actualT]),]                              # aktuális transcriptID-hoz tartozó CDS temp tibble
  temp_2 <- aostoyae_cdslengths[[names(aostoyae_cdslengths)[actualT]]]                                                        # aktuális transcriptID-hoz tartozó exonhossz temp vector
  start_pos_exon <- "0"
  end_pos_exon <- "0"
  for (x in 0:length(temp_2)) {
    start_pos_exon[x+1] <- temp_1$start-1 <= sum(temp_2[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a kezdő érték és a logicalokat belementjük egy vectorba
  }
  start_pos_exon <- which(start_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben elkezdődik a CDS régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a CDS régió része
  
  for (x in 0:length(temp_2)) {
    end_pos_exon[x+1] <- temp_1$end <= sum(temp_2[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a záró érték és a logicalokat belementjük egy vectorba
  }
  end_pos_exon <- which(end_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben befejeződik a CDS régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a CDS régió része
  
  temp_2[start_pos_exon] <- sum(temp_2[1:start_pos_exon]) - (temp_1$start-1)          # meghatározzuk hogy a start exonba pontosan mekkora a CDS régió és bementjük a pozícióba
  
  for (x in 0:(start_pos_exon-1)) {
    temp_2[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a CDS régió kezdőpontja előtt van
  }
  
  temp_2[end_pos_exon] <- temp_2[end_pos_exon] - (sum(temp_2[1:end_pos_exon]) - temp_1$length)       # meghatározzuk hogy az utolsó exon mennyi még a CDS régió része
  
  for (x in (end_pos_exon+1):(length(temp_2)+1)) {
    temp_2[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a CDS régió kezdőpontja előtt van
  }
  
  aostoyae_cdslengths[[actualT]] <- temp_2[-(length(temp_2))]
}

# ccinerea

ccinerea_filter <- read_tsv("~/Desktop/alternative_splicing/ccinerea_AmutBmut/ccinerea_AmutBmut_transcripts.fasta.transdecoder_dir/filter", col_names = F)
ccinerea_CDS <- read_tsv("~/Desktop/alternative_splicing/ccinerea_AmutBmut/ccinerea_AmutBmut_transcripts.fasta.transdecoder_dir/headers.cds", col_names = F) %>%
  inner_join(ccinerea_filter, by = "X1") %>%
  select(X6) %>%
  separate(X6, c("transcriptID", "X6.2"), ":") %>%
  separate(X6.2, c("start", "X6.2.2"), "-") %>%
  separate(X6.2.2, c("end", "strand"), "\\(") %>%
  select(transcriptID, start, end)
ccinerea_CDS$start <- as.integer(ccinerea_CDS$start)
ccinerea_CDS$end <- as.integer(ccinerea_CDS$end)
ccinerea_CDS <- mutate(ccinerea_CDS, length = end - start + 1)

ccinerea_cdslengths <- ccinerea_exons

for (actualT in 1:length(ccinerea_exons)) {
  temp_1 <- ccinerea_CDS[str_detect(ccinerea_CDS$transcriptID, names(ccinerea_exons)[actualT]),]                              # aktuális transcriptID-hoz tartozó CDS temp tibble
  temp_2 <- ccinerea_cdslengths[[names(ccinerea_cdslengths)[actualT]]]                                                        # aktuális transcriptID-hoz tartozó exonhossz temp vector
  start_pos_exon <- "0"
  for (x in 0:length(temp_2)) {
    start_pos_exon[x+1] <- temp_1$start-1 <= sum(temp_2[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a kezdő érték és a logicalokat belementjük egy vectorba
  }
  start_pos_exon <- which(start_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben elkezdődik a CDS régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a CDS régió része
  
  end_pos_exon <- "0"
  for (x in 0:length(temp_2)) {
    end_pos_exon[x+1] <- temp_1$end <= sum(temp_2[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a záró érték és a logicalokat belementjük egy vectorba
  }
  end_pos_exon <- which(end_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben befejeződik a CDS régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a CDS régió része
  
  temp_2[start_pos_exon] <- sum(temp_2[1:start_pos_exon]) - (temp_1$start-1)          # meghatározzuk hogy a start exonba pontosan mekkora a CDS régió és bementjük a pozícióba
  
  for (x in 0:(start_pos_exon-1)) {
    temp_2[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a CDS régió kezdőpontja előtt van
  }
  
  temp_2[end_pos_exon] <- temp_2[end_pos_exon] - (sum(temp_2[1:end_pos_exon]) - temp_1$length)       # meghatározzuk hogy az utolsó exonl mennyi még a CDS régió része
  
  for (x in (end_pos_exon+1):(length(temp_2)+1)) {
    temp_2[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a CDS régió kezdőpontja előtt van
  }
  
  ccinerea_cdslengths[[actualT]] <- temp_2[-(length(temp_2))]
}

# ltigrinus

ltigrinus_filter <- read_tsv("~/Desktop/alternative_splicing/ltigrinus/ltigrinus_transcripts.fasta.transdecoder_dir/filter", col_names = F)
ltigrinus_CDS <- read_tsv("~/Desktop/alternative_splicing/ltigrinus/ltigrinus_transcripts.fasta.transdecoder_dir/headers.cds", col_names = F) %>%
  inner_join(ltigrinus_filter, by = "X1") %>%
  select(X6) %>%
  separate(X6, c("transcriptID", "X6.2"), ":") %>%
  separate(X6.2, c("start", "X6.2.2"), "-") %>%
  separate(X6.2.2, c("end", "strand"), "\\(") %>%
  select(transcriptID, start, end)
ltigrinus_CDS$start <- as.integer(ltigrinus_CDS$start)
ltigrinus_CDS$end <- as.integer(ltigrinus_CDS$end)
ltigrinus_CDS <- mutate(ltigrinus_CDS, length = end - start + 1)

ltigrinus_cdslengths <- ltigrinus_exons

for (actualT in 1:length(ltigrinus_exons)) {
  temp_1 <- ltigrinus_CDS[str_detect(ltigrinus_CDS$transcriptID, names(ltigrinus_exons)[actualT]),]                              # aktuális transcriptID-hoz tartozó CDS temp tibble
  temp_2 <- ltigrinus_cdslengths[[names(ltigrinus_cdslengths)[actualT]]]                                                        # aktuális transcriptID-hoz tartozó exonhossz temp vector
  start_pos_exon <- "0"
  for (x in 0:length(temp_2)) {
    start_pos_exon[x+1] <- temp_1$start-1 <= sum(temp_2[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a kezdő érték és a logicalokat belementjük egy vectorba
  }
  start_pos_exon <- which(start_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben elkezdődik a CDS régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a CDS régió része
  
  end_pos_exon <- "0"
  for (x in 0:length(temp_2)) {
    end_pos_exon[x+1] <- temp_1$end <= sum(temp_2[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a záró érték és a logicalokat belementjük egy vectorba
  }
  end_pos_exon <- which(end_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben befejeződik a CDS régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a CDS régió része
  
  temp_2[start_pos_exon] <- sum(temp_2[1:start_pos_exon]) - (temp_1$start-1)          # meghatározzuk hogy a start exonba pontosan mekkora a CDS régió és bementjük a pozícióba
  
  for (x in 0:(start_pos_exon-1)) {
    temp_2[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a CDS régió kezdőpontja előtt van
  }
  
  temp_2[end_pos_exon] <- temp_2[end_pos_exon] - (sum(temp_2[1:end_pos_exon]) - temp_1$length)       # meghatározzuk hogy az utolsó exonl mennyi még a CDS régió része
  
  for (x in (end_pos_exon+1):(length(temp_2)+1)) {
    temp_2[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a CDS régió kezdőpontja előtt van
  }
  
  ltigrinus_cdslengths[[actualT]] <- temp_2[-(length(temp_2))]
}


# scommune

scommune_filter <- read_tsv("~/Desktop/alternative_splicing/scommune/scommune_transcripts.fasta.transdecoder_dir/filter", col_names = F)
scommune_CDS <- read_tsv("~/Desktop/alternative_splicing/scommune/scommune_transcripts.fasta.transdecoder_dir/headers.cds", col_names = F) %>%
  inner_join(scommune_filter, by = "X1") %>%
  select(X6) %>%
  separate(X6, c("transcriptID", "X6.2"), ":") %>%
  separate(X6.2, c("start", "X6.2.2"), "-") %>%
  separate(X6.2.2, c("end", "strand"), "\\(") %>%
  select(transcriptID, start, end)
scommune_CDS$start <- as.integer(scommune_CDS$start)
scommune_CDS$end <- as.integer(scommune_CDS$end)
scommune_CDS <- mutate(scommune_CDS, length = end - start + 1)

scommune_cdslengths <- scommune_exons

for (actualT in 1:length(scommune_exons)) {
  temp_1 <- scommune_CDS[str_detect(scommune_CDS$transcriptID, names(scommune_exons)[actualT]),]                              # aktuális transcriptID-hoz tartozó CDS temp tibble
  temp_2 <- scommune_cdslengths[[names(scommune_cdslengths)[actualT]]]                                                        # aktuális transcriptID-hoz tartozó exonhossz temp vector
  start_pos_exon <- "0"
  for (x in 0:length(temp_2)) {
    start_pos_exon[x+1] <- temp_1$start-1 <= sum(temp_2[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a kezdő érték és a logicalokat belementjük egy vectorba
  }
  start_pos_exon <- which(start_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben elkezdődik a CDS régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a CDS régió része
  
  end_pos_exon <- "0"
  for (x in 0:length(temp_2)) {
    end_pos_exon[x+1] <- temp_1$end <= sum(temp_2[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a záró érték és a logicalokat belementjük egy vectorba
  }
  end_pos_exon <- which(end_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben befejeződik a CDS régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a CDS régió része
  
  temp_2[start_pos_exon] <- sum(temp_2[1:start_pos_exon]) - (temp_1$start-1)          # meghatározzuk hogy a start exonba pontosan mekkora a CDS régió és bementjük a pozícióba
  
  for (x in 0:(start_pos_exon-1)) {
    temp_2[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a CDS régió kezdőpontja előtt van
  }
  
  temp_2[end_pos_exon] <- temp_2[end_pos_exon] - (sum(temp_2[1:end_pos_exon]) - temp_1$length)       # meghatározzuk hogy az utolsó exonl mennyi még a CDS régió része
  
  for (x in (end_pos_exon+1):(length(temp_2)+1)) {
    temp_2[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a CDS régió kezdőpontja előtt van
  }
  
  scommune_cdslengths[[actualT]] <- temp_2[-(length(temp_2))]
}

##### DOMAIN COMPOSITION #####

# aostoyae

raw_interpro <- read_lines("~/Desktop/alternative_splicing/aostoyae/aostoyae_proteins_all.fasta.tsv")
aostoyae_IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)
colnames(aostoyae_IPR) <- "all"
aostoyae_IPR <- aostoyae_IPR %>%
  separate(all, c("transcriptID", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t") %>%
  select(transcriptID, start, end, IPR_domain)
aostoyae_IPR$start <- as.integer(aostoyae_IPR$start)               # átalakítjuk integerre
aostoyae_IPR$end <- as.integer(aostoyae_IPR$end)
aostoyae_IPR$start <- ((aostoyae_IPR$start - 1) * 3) + 1           # fixáljuk a protein kezdőpontot nukleotid kezdőpontra
aostoyae_IPR$end <- aostoyae_IPR$end * 3
aostoyae_IPR$length <- aostoyae_IPR$end - aostoyae_IPR$start + 1

aostoyae_CDS_forDOMAIN <- aostoyae_cdslengths

aostoyae_domains <- list()

for (y in 1:length(aostoyae_CDS_forDOMAIN )) {
  temp_1 <- aostoyae_IPR[str_detect(aostoyae_IPR$transcriptID, names(aostoyae_exons)[y]),]
  temp_2 <- aostoyae_CDS_forDOMAIN [[names(aostoyae_CDS_forDOMAIN )[y]]] 
  domain_list <- list()
  for (sub_x in 1:length(temp_1$transcriptID)) {
    temp_1_sub <- temp_1[sub_x,]
    temp_2_sub <- temp_2
    start_pos_exon <- "0"
    end_pos_exon <- "0"
    for (x in 0:length(temp_2_sub)) {
      start_pos_exon[x+1] <- temp_1_sub$start-1 <= sum(temp_2_sub[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a kezdő érték és a logicalokat belementjük egy vectorba
    }
    start_pos_exon <- which(start_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben elkezdődik a domain vagy amennyiben pontosan egyenlő akkor azt ami még nem a domain része
    for (x in 0:length(temp_2_sub)) {
      end_pos_exon[x+1] <- temp_1_sub$end <= sum(temp_2_sub[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a záró érték és a logicalokat belementjük egy vectorba
    }
    end_pos_exon <- which(end_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben befejeződik a domain régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a domain régió része
    temp_2_sub[start_pos_exon] <- sum(temp_2_sub[1:start_pos_exon]) - (temp_1_sub$start-1)          # meghatározzuk hogy a start exonba pontosan mekkora a domain és bementjük a pozícióba
    for (x in 0:(start_pos_exon-1)) {
      temp_2_sub[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a domain kezdőpontja előtt van
    }
    temp_2_sub[end_pos_exon] <- temp_2_sub[end_pos_exon] - (sum(temp_2_sub[1:end_pos_exon]) - temp_1_sub$length)       # meghatározzuk hogy az utolsó exon mennyi még a domain része
    for (x in (end_pos_exon+1):(length(temp_2_sub)+1)) {
      temp_2_sub[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a domain kezdőpontja előtt van
    }
    domain_list[[sub_x]] <- temp_2_sub[-(length(temp_2_sub))]
  }
  names(domain_list) <- temp_1$IPR_domain
  aostoyae_domains[[y]] <- domain_list
}

names(aostoyae_domains) <- names(aostoyae_CDS_forDOMAIN )

# ccinerea

raw_interpro <- read_lines("~/Desktop/alternative_splicing/ccinerea_AmutBmut/ccinerea_AmutBmut_proteins_all.fasta.tsv")
ccinerea_IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)
colnames(ccinerea_IPR) <- "all"
ccinerea_IPR <- ccinerea_IPR %>%
  separate(all, c("transcriptID", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t") %>%
  select(transcriptID, start, end, IPR_domain)
ccinerea_IPR$start <- as.integer(ccinerea_IPR$start)               # átalakítjuk integerre
ccinerea_IPR$end <- as.integer(ccinerea_IPR$end)
ccinerea_IPR$start <- ((ccinerea_IPR$start - 1) * 3) + 1           # fixáljuk a protein kezdőpontot nukleotid kezdőpontra
ccinerea_IPR$end <- ccinerea_IPR$end * 3
ccinerea_IPR$length <- ccinerea_IPR$end - ccinerea_IPR$start + 1

ccinerea_CDS_forDOMAIN <- ccinerea_cdslengths

ccinerea_domains <- list()

for (y in 1:length(ccinerea_CDS_forDOMAIN)) {
  temp_1 <- ccinerea_IPR[str_detect(ccinerea_IPR$transcriptID, names(ccinerea_exons)[y]),]
  temp_2 <- ccinerea_CDS_forDOMAIN[[names(ccinerea_CDS_forDOMAIN)[y]]] 
  domain_list <- list()
  for (sub_x in 1:length(temp_1$transcriptID)) {
    temp_1_sub <- temp_1[sub_x,]
    temp_2_sub <- temp_2
    start_pos_exon <- "0"
    end_pos_exon <- "0"
    for (x in 0:length(temp_2_sub)) {
      start_pos_exon[x+1] <- temp_1_sub$start-1 <= sum(temp_2_sub[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a kezdő érték és a logicalokat belementjük egy vectorba
    }
    start_pos_exon <- which(start_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben elkezdődik a domain vagy amennyiben pontosan egyenlő akkor azt ami még nem a domain része
    for (x in 0:length(temp_2_sub)) {
      end_pos_exon[x+1] <- temp_1_sub$end <= sum(temp_2_sub[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a záró érték és a logicalokat belementjük egy vectorba
    }
    end_pos_exon <- which(end_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben befejeződik a domain régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a domain régió része
    temp_2_sub[start_pos_exon] <- sum(temp_2_sub[1:start_pos_exon]) - (temp_1_sub$start-1)          # meghatározzuk hogy a start exonba pontosan mekkora a domain és bementjük a pozícióba
    for (x in 0:(start_pos_exon-1)) {
      temp_2_sub[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a domain kezdőpontja előtt van
    }
    temp_2_sub[end_pos_exon] <- temp_2_sub[end_pos_exon] - (sum(temp_2_sub[1:end_pos_exon]) - temp_1_sub$length)       # meghatározzuk hogy az utolsó exon mennyi még a domain része
    for (x in (end_pos_exon+1):(length(temp_2_sub)+1)) {
      temp_2_sub[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a domain kezdőpontja előtt van
    }
    domain_list[[sub_x]] <- temp_2_sub[-(length(temp_2_sub))]
  }
  names(domain_list) <- temp_1$IPR_domain
  ccinerea_domains[[y]] <- domain_list
}

names(ccinerea_domains) <- names(ccinerea_CDS_forDOMAIN)

# ltigrinus

raw_interpro <- read_lines("~/Desktop/alternative_splicing/ltigrinus/ltigrinus_proteins_all.fasta.tsv")
ltigrinus_IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)
colnames(ltigrinus_IPR) <- "all"
ltigrinus_IPR <- ltigrinus_IPR %>%
  separate(all, c("transcriptID", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t") %>%
  select(transcriptID, start, end, IPR_domain)
ltigrinus_IPR$start <- as.integer(ltigrinus_IPR$start)               # átalakítjuk integerre
ltigrinus_IPR$end <- as.integer(ltigrinus_IPR$end)
ltigrinus_IPR$start <- ((ltigrinus_IPR$start - 1) * 3) + 1           # fixáljuk a protein kezdőpontot nukleotid kezdőpontra
ltigrinus_IPR$end <- ltigrinus_IPR$end * 3
ltigrinus_IPR$length <- ltigrinus_IPR$end - ltigrinus_IPR$start + 1

ltigrinus_CDS_forDOMAIN <- ltigrinus_cdslengths

ltigrinus_domains <- list()

for (y in 1:length(ltigrinus_CDS_forDOMAIN)) {
  temp_1 <- ltigrinus_IPR[str_detect(ltigrinus_IPR$transcriptID, names(ltigrinus_exons)[y]),]
  temp_2 <- ltigrinus_CDS_forDOMAIN[[names(ltigrinus_CDS_forDOMAIN)[y]]] 
  domain_list <- list()
  for (sub_x in 1:length(temp_1$transcriptID)) {
    temp_1_sub <- temp_1[sub_x,]
    temp_2_sub <- temp_2
    start_pos_exon <- "0"
    end_pos_exon <- "0"
    for (x in 0:length(temp_2_sub)) {
      start_pos_exon[x+1] <- temp_1_sub$start-1 <= sum(temp_2_sub[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a kezdő érték és a logicalokat belementjük egy vectorba
    }
    start_pos_exon <- which(start_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben elkezdődik a domain vagy amennyiben pontosan egyenlő akkor azt ami még nem a domain része
    for (x in 0:length(temp_2_sub)) {
      end_pos_exon[x+1] <- temp_1_sub$end <= sum(temp_2_sub[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a záró érték és a logicalokat belementjük egy vectorba
    }
    end_pos_exon <- which(end_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben befejeződik a domain régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a domain régió része
    temp_2_sub[start_pos_exon] <- sum(temp_2_sub[1:start_pos_exon]) - (temp_1_sub$start-1)          # meghatározzuk hogy a start exonba pontosan mekkora a domain és bementjük a pozícióba
    for (x in 0:(start_pos_exon-1)) {
      temp_2_sub[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a domain kezdőpontja előtt van
    }
    temp_2_sub[end_pos_exon] <- temp_2_sub[end_pos_exon] - (sum(temp_2_sub[1:end_pos_exon]) - temp_1_sub$length)       # meghatározzuk hogy az utolsó exon mennyi még a domain része
    for (x in (end_pos_exon+1):(length(temp_2_sub)+1)) {
      temp_2_sub[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a domain kezdőpontja előtt van
    }
    domain_list[[sub_x]] <- temp_2_sub[-(length(temp_2_sub))]
  }
  names(domain_list) <- temp_1$IPR_domain
  ltigrinus_domains[[y]] <- domain_list
}

names(ltigrinus_domains) <- names(ltigrinus_CDS_forDOMAIN)

# scommune

raw_interpro <- read_lines("~/Desktop/alternative_splicing/scommune/scommune_proteins_all.fasta.tsv")
scommune_IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)
colnames(scommune_IPR) <- "all"
scommune_IPR <- scommune_IPR %>%
  separate(all, c("transcriptID", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t") %>%
  select(transcriptID, start, end, IPR_domain)
scommune_IPR$start <- as.integer(scommune_IPR$start)               # átalakítjuk integerre
scommune_IPR$end <- as.integer(scommune_IPR$end)
scommune_IPR$start <- ((scommune_IPR$start - 1) * 3) + 1           # fixáljuk a protein kezdőpontot nukleotid kezdőpontra
scommune_IPR$end <- scommune_IPR$end * 3
scommune_IPR$length <- scommune_IPR$end - scommune_IPR$start + 1

scommune_CDS_forDOMAIN <- scommune_cdslengths

scommune_domains <- list()

for (y in 1:length(scommune_CDS_forDOMAIN)) {
  temp_1 <- scommune_IPR[str_detect(scommune_IPR$transcriptID, names(scommune_exons)[y]),]
  temp_2 <- scommune_CDS_forDOMAIN[[names(scommune_CDS_forDOMAIN)[y]]] 
  domain_list <- list()
  for (sub_x in 1:length(temp_1$transcriptID)) {
    temp_1_sub <- temp_1[sub_x,]
    temp_2_sub <- temp_2
    start_pos_exon <- "0"
    end_pos_exon <- "0"
    for (x in 0:length(temp_2_sub)) {
      start_pos_exon[x+1] <- temp_1_sub$start-1 <= sum(temp_2_sub[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a kezdő érték és a logicalokat belementjük egy vectorba
    }
    start_pos_exon <- which(start_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben elkezdődik a domain vagy amennyiben pontosan egyenlő akkor azt ami még nem a domain része
    for (x in 0:length(temp_2_sub)) {
      end_pos_exon[x+1] <- temp_1_sub$end <= sum(temp_2_sub[1:(1+x)])                         # összeadjunk az exonokat sorba 1, 1:2, 1:3 stb. és megnézzük hogy az egyesével összeadott exonok hossza nagyobb-e mint a záró érték és a logicalokat belementjük egy vectorba
    }
    end_pos_exon <- which(end_pos_exon == TRUE)[1]                       # megkapjuk azt az exont amiben befejeződik a domain régió vagy amennyiben pontosan egyenlő akkor azt ami még nem a domain régió része
    temp_2_sub[start_pos_exon] <- sum(temp_2_sub[1:start_pos_exon]) - (temp_1_sub$start-1)          # meghatározzuk hogy a start exonba pontosan mekkora a domain és bementjük a pozícióba
    for (x in 0:(start_pos_exon-1)) {
      temp_2_sub[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a domain kezdőpontja előtt van
    }
    temp_2_sub[end_pos_exon] <- temp_2_sub[end_pos_exon] - (sum(temp_2_sub[1:end_pos_exon]) - temp_1_sub$length)       # meghatározzuk hogy az utolsó exon mennyi még a domain része
    for (x in (end_pos_exon+1):(length(temp_2_sub)+1)) {
      temp_2_sub[x] <- 0                                                                    # kinullázzuk az összes olyan exont ami a domain kezdőpontja előtt van
    }
    domain_list[[sub_x]] <- temp_2_sub[-(length(temp_2_sub))]
  }
  names(domain_list) <- temp_1$IPR_domain
  scommune_domains[[y]] <- domain_list
}

names(scommune_domains) <- names(scommune_CDS_forDOMAIN)




scommune_domains[[2]]
scommune_IPR[str_detect(scommune_IPR$transcriptID, names(scommune_exons)[1]),]

scommune_IPR[str_detect(scommune_IPR$transcriptID, names(scommune_exons)[1]),]$start - 1 + 7



