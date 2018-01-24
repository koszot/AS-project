Table of Contents
=================

   * [Alternative Splicing Analysis Pipeline](#alternative-splicing-analysis-pipeline)
      * [Preparations](#preparations)
      * [RRPM: STAR Alignment and Cufflinks Assembly](#rrpm-star-alignment-and-cufflinks-assembly)
      * [Filtering and Merge](#filtering-and-merge)
      * [Fusion correction](#fusion-correction)
      * [Alternative Splicing Statistics](#alternative-splicing-statistics)
      * [Expression: STAR Alignment, Cuffquant Assembly and Cuffdiff analysis](#expression-star-alignment-cuffquant-assembly-and-cuffdiff-analysis)
      * [ORF prediction and InterProScan analysis](#orf-prediction-and-interproscan-analysis)
      * [AS transcripts table (OUTDATED)](#as-transcripts-table-outdated)
      * [Enrichment (OUTDATED)](#enrichment-outdated)
      * [MCL Clustering a géncsaládok megkeresésére](#mcl-clustering-a-géncsaládok-megkeresésére)
      * [Summary (OUTDATED)](#summary-outdated)
      * [Silix/Hifix és enrichment (OUTDATED)](#silixhifix-és-enrichment-outdated)

# Alternative Splicing Analysis Pipeline

## Preparations

A bemeneti annotációs fájlokat át kell alakítani, hogy megfelelőek legyenek a Cufflinks számára illetve elő kell állítani azokat az annotációs és FASTA fájlokat amik az RRPM analízishez szükségesek.

### Armillaria ostoyae

```
# A GFF3 fájl átalakítása GTF formátumra
gffread ./aostoyae/genome/p3_i2_t47428_Arm_ostoy_v2.gff3 -T -o ./aostoyae/genome/original.gtf

# Az analízishez szükséges GTF fájlok elkészítése
scripts/aostoyae/PREPARATIONS_aostoyae.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./aostoyae/genome/original_onlygene.gtf > ./aostoyae/genome/original_onlygene.gtf.bed

# A géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./aostoyae/aostoyae_genome/original_genes.fasta -fi ./aostoyae/aostoyae_genome/p3_i2_t47428_Arm_ostoy_v2.scaf -bed ./aostoyae/aostoyae_genome/original_onlygene.gtf.bed

# A géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./aostoyae/aostoyae_genome/original_genes.fasta
```

### Auriculariopsis ampla

```
# Az analízishez szükséges GTF fájlok elkészítése
scripts/aampla/PREPARATIONS_aampla.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./aampla/aampla_genome/aampla_onlygene.gtf > ./aampla/aampla_genome/aampla_onlygene.gtf.bed

# A géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./aampla/aampla_genome/aampla_genes.fasta -fi ./aampla/aampla_genome/Auramp1_AssemblyScaffolds.fasta -bed ./aampla/aampla_genome/aampla_onlygene.gtf.bed

# A géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./aampla/aampla_genome/aampla_genes.fasta
```

### Coprinopsis cinerea

```
# Az analízishez szükséges GTF fájlok elkészítése
scripts/ccinerea/PREPARATIONS_ccinerea.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./ccinerea/ccinerea_AmutBmut_genome/ccinerea_AmutBmut_onlygene.gtf > ./ccinerea/ccinerea_AmutBmut_genome/ccinerea_AmutBmut_onlygene.gtf.bed

# A géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./ccinerea/ccinerea_AmutBmut_genome/ccinerea_AmutBmut_genes.fasta -fi ./ccinerea/ccinerea_AmutBmut_genome/Copci_AmutBmut1_AssemblyScaffolds.fasta -bed ./ccinerea/ccinerea_AmutBmut_genome/ccinerea_AmutBmut_onlygene.gtf.bed

# A géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./ccinerea/ccinerea_AmutBmut_genome/ccinerea_AmutBmut_genes.fasta
```

### Cryptococcus neoformans

NEM SAJÁT RNA-SEQ

Fausto Almeida, Julie M. Wolf, Thiago Aparecido da Silva, Carlos M. DeLeon-Rodriguez, Caroline Patini Rezende, André Moreira Pessoni, Fabrício Freitas Fernandes, Rafael Silva-Rocha, Roberto Martinez, Marcio L. Rodrigues, Maria Cristina Roque-Barreira & Arturo Casadevall

Galectin-3 impacts Cryptococcus neoformans infection through direct antifungal effects

RNA-seq DATA:
- R1 : RNAseq of C. noformans not treated (SRR6207454)	->	CN_R1.R1.fastq	CN_R1.R2.fastq
- R2 : RNAseq of C. noformans not treated (SRR6207453)	->	CN_R2.R1.fastq	CN_R2.R2.fastq
- R3 : RNAseq of C. noformans not treated (SRR6207452)	->	CN_R3.R1.fastq	CN_R3.R2.fastq

Annotation:
- JGI
- Cryptococcus neoformans var. grubii H99

```
# Az analízishez szükséges GTF fájlok elkészítése
scripts/cneoformans/PREPARATIONS_cneoformans.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./cneoformans_genome/cneoformans_onlygene.gtf > ./cneoformans_genome/cneoformans_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./cneoformans_genome/cneoformans_genes.fasta -fi ./cneoformans_genome/Cryptococcus_neoformans_H99.unmasked.fasta -bed ./cneoformans_genome/cneoformans_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./cneoformans_genome/cneoformans_genes.fasta
```

### Lentinus tigrinus

```
# az analízishez szükséges GTF fájlok elkészítése
scripts/ltiginus/PREPARATIONS_ltigrinus.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./ltigrinus_genome/ltigrinus_onlygene.gtf > ./ltigrinus_genome/ltigrinus_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./ltigrinus_genome/ltigrinus_genes.fasta -fi ./ltigrinus_genome/Sisbr1_AssemblyScaffolds.fasta -bed ./ltigrinus_genome/ltigrinus_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./ltigrinus_genome/ltigrinus_genes.fasta
```

### Phanerochaete chrysosporium

```
# az analízishez szükséges GTF fájlok elkészítése
scripts/pchrysosporium/PREPARATIONS_pchrysosporium.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./pchrysosporium_genome/pchrysosporium_onlygene.gtf > ./pchrysosporium_genome/pchrysosporium_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./pchrysosporium_genome/pchrysosporium_genes.fasta -fi ./pchrysosporium_genome/Phchr2_AssemblyScaffolds.fasta -bed ./pchrysosporium_genome/pchrysosporium_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./pchrysosporium_genome/pchrysosporium_genes.fasta
```

### Rickenella mellea

```
# az analízishez szükséges GTF fájlok elkészítése
scripts/rmellea/PREPARATIONS_rmellea.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./rmellea_genome/rmellea_onlygene.gtf > ./rmellea_genome/rmellea_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./rmellea_genome/rmellea_genes.fasta -fi ./rmellea_genome/Ricmel1_AssemblyScaffolds.fasta -bed ./rmellea_genome/rmellea_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./rmellea_genome/rmellea_genes.fasta
```

### Schizophyllum commune

```
# az analízishez szükséges GTF fájlok elkészítése
scripts/scommune/PREPARATIONS_scommune.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./genome/original_onlygene.gtf > ./genome/original_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./genome/original_genes.fasta -fi ./genome/Schco3_AssemblyScaffolds.fasta -bed ./genome/original_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./genome/original_genes.fasta
```

### Ustilago maydis

Marie Tollot, Daniela Assmann, Christian Becker, Janine Altmüller, Julien Y. Dutheil, Carl-Eric Wegner, Regine Kahmann 

The WOPR Protein Ros1 Is a Master Regulator of Sporogenesis and Late Effector Gene Expression in the Maize Pathogen Ustilago maydis

RNA-seq DATA:
- R1 : GSM1977379: WT_I; Ustilago maydis; RNA-Seq (SRR3038903)	->	UM_R1.R1.fastq	UM_R1.R2.fastq
- R2 : GSM1977380: WT_II; Ustilago maydis; RNA-Seq (SRR3038904)	->	UM_R2.R1.fastq	UM_R2.R2.fastq
- R3 : GSM1977380: WT_III; Ustilago maydis; RNA-Seq (SRR3038905) ->	UM_R3.R1.fastq	UM_R3.R2.fastq

Annotation:
- ftp://ftpmips.gsf.de/fungi/Ustilaginaceae/Ustilago_maydis_521/
- p3_t237631_Ust_maydi_v2GB.gff3
- p3_t237631_Ust_maydi_v2GB.scaf
- p3_t237631_Ust_maydi_v2GB.prot

```
# a GFF3 fájl átalakítása GTF formátumra
gffread ./umaydis_genome/p3_t237631_Ust_maydi_v2GB.gff3 -T -o ./umaydis_genome/original.gtf

# az analízishez szükséges GTF fájlok elkészítése
scripts/umaydis/PREPARATIONS_umaydis.R

# GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
gtf2bed < ./umaydis_genome/umaydis_onlygene.gtf > ./umaydis_genome/umaydis_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl elkészítése
bedtools getfasta -name -fo ./umaydis_genome/umaydis_genes.fasta -fi ./umaydis_genome/p3_t237631_Ust_maydi_v2GB.scaf -bed ./umaydis_genome/umaydis_onlygene.gtf.bed

# a géneket tartalmazó FASTA fájl headerjének a trimmelése
perl -pi -e 's/::.*//g' ./umaydis_genome/umaydis_genes.fasta
```

## RRPM: STAR Alignment and Cufflinks Assembly

Megnézzük milyen hosszúak a gének és az intronok a STAR és a Cufflinks beállításához.

Elkészítünk egy fájlt ami tartalmazza az összes szükséges scriptet ami az RRPM futtatásához szükséges (STAR, Cufflinks)

### Armillaria ostoyae

```
bioawk -c fastx '{ print $name, length($seq) }' < original_genes.fasta > gene_length

scripts/aostoyae/INTRON_LENGTH_aostoyae.R

# A maximális transzkriptméret: 16175 --> max-bundle-length marad 250000
# Intron min: 21 --> --alignIntronMin marad 3 
# Intron max: 6767 --> --alignIntronMax marad 30000

Cufflinks_scripts/aostoyae_STAR_CUFF_RRPM.sh
```

### Auriculariopsis ampla

```
bioawk -c fastx '{ print $name, length($seq) }' < aampla_genes.fasta > gene_length

scripts/aampla/INTRON_LENGTH_aampla.R

# A maximális transzkriptméret: 23200 --> max-bundle-length marad 250000
# Intron min: 3 --> --alignIntronMin marad 3 
# Intron max: 2000 --> --alignIntronMax marad 30000

Cufflinks_scripts/aampla_STAR_CUFF_RRPM.sh
```

### Coprinopsis cinerea

```
bioawk -c fastx '{ print $name, length($seq) }' < ccinerea_AmutBmut_genes.fasta > gene_length

scripts/ccinerea/INTRON_LENGTH_ccinerea.R

# A maximális transzkriptméret: 19640 --> max-bundle-length marad 250000
# Intron min: 11 --> --alignIntronMin marad 3 
# Intron max: 8003 --> --alignIntronMax marad 30000

Cufflinks_scripts/ccinerea_STAR_CUFF_RRPM.sh
```

### Cryptococcus neoformans

```
bioawk -c fastx '{ print $name, length($seq) }' < cneoformans_genes.fasta > gene_length

scripts/cneoformans/INTRON_LENGTH_cneoformans.R

# A maximális transzkriptméret: 15614 --> max-bundle-length marad 250000
# Intron min: 22 --> --alignIntronMin marad 3 
# Intron max: 1306 --> --alignIntronMax marad 30000

Cufflinks_scripts/cneoformans_STAR_CUFF_RRPM.sh
```

### Lentinus tigrinus

```
bioawk -c fastx '{ print $name, length($seq) }' < ltigrinus_genes.fasta > gene_length

scripts/ltigrinus/INTRON_LENGTH_ltigrinus.R

# A maximális transzkriptméret: 16440 --> max-bundle-length marad 250000
# Intron min: 11 --> --alignIntronMin marad 3 
# Intron max: 15143 --> --alignIntronMax marad 30000

Cufflinks_scripts/ltigrinus_STAR_CUFF_RRPM.sh
```

### Phanerochaete chrysosporium

```
bioawk -c fastx '{ print $name, length($seq) }' < pchrysosporium_genes.fasta > gene_length

scripts/pchrysosporium/INTRON_LENGTH_pchrysosporium.R

# A maximális transzkriptméret: 30598 --> max-bundle-length marad 250000
# Intron min: 11 --> --alignIntronMin marad 3 
# Intron max: 28051 --> --alignIntronMax marad 30000

Cufflinks_scripts/pchrysosporium_STAR_CUFF_RRPM.sh
```

### Rickenella mellea

```
bioawk -c fastx '{ print $name, length($seq) }' < rmellea_genes.fasta > gene_length

scripts/rmellea/INTRON_LENGTH_rmellea.R

# A maximális transzkriptméret: 13037 --> max-bundle-length marad 250000
# Intron min: 7 --> --alignIntronMin marad 3 
# Intron max: 1984 --> --alignIntronMax marad 30000

Cufflinks_scripts/rmellea_STAR_CUFF_RRPM.sh
```

### Schizophyllum commune

```
bioawk -c fastx '{ print $name, length($seq) }' <original_genes.fasta > gene_length

scripts/scommune/INTRON_LENGTH_scommune.R

# A maximális transzkriptméret: 51113 --> max-bundle-length marad 250000
# Intron min: 11 --> --alignIntronMin marad 3 
# Intron max: 28634 --> --alignIntronMax marad 30000

Cufflinks_scripts/scommune_STAR_CUFF_RRPM.sh
```

### Ustilago maydis

Túl sok a short read miatti unmapped MEGOLDÁS -> https://github.com/alexdobin/STAR/issues/169

```
bioawk -c fastx '{ print $name, length($seq) }' < umaydis_genes.fasta > gene_length

scripts/umaydis/INTRON_LENGTH_umaydis.R

# A maximális transzkriptméret: 16296 --> max-bundle-length marad 250000
# Intron min: 15 --> --alignIntronMin marad 3 
# Intron max: 2182 --> --alignIntronMax marad 30000

Cufflinks_scripts/umaydis_STAR_CUFF_RRPM.sh
```

## Filtering and Merge

Az RRPM outputot filterezzük, majd összeillesztjük az eredeti annotációs fájlal, hogy kitöltsük az esetleges hézagokat, ahol nem mutatott expressziós értéket az eredeti annotációban szereplő transzkript.

A GTF-ekbe átalakítja a start/stop helyeket pl.: 73e3 = 73000, ezeket kiszűrtem kézzel IGV segítségével, a beolvasásnál ha ezeket integerről stringre állítjuk akkor nem írja át, de én nem bajlodtam már vele miután lefutott.

Lépések:
- Expression Filtering: Kitörlünk minden olyan transzkriptet aminek nincsen látható expressziója
- Full Read Support Filtering: Kitörülünk minden olyan transzkriptet aminek nincs minden intron-exon junctionjére legalább 1 read
- Context Restoration: Visszaállítjuk a GTF fájl start és end pozíciójat az eredeti annotáció pozícióihoz
- Strand Filtering: Kitörlünk minden olyan transzkriptet ami nem ugyanazon a stranden van mint az alap gén (ezeket nem tudjuk biztosan megmondani strand specifikus readek nélkül)
- Annotation Merge: Összeillesztjük az eredeti annotációs fájlal

### Armillaria ostoyae

```
scripts/aostoyae/FILTERING_aostoyae.R

perl -pi -e 's/\"transcript_id \"\"/transcript_id \"/g' RRPM_transcripts.fixed.gtf
perl -pi -e 's/\"\"; gene_id \"\"/\"; gene_id \"/g' RRPM_transcripts.fixed.gtf
perl -pi -e 's/\"\";\"/\";/g' RRPM_transcripts.fixed.gtf

perl -pi -e 's/\"transcript_id \"\"/transcript_id \"/g' RRPM_transcripts.gtf
perl -pi -e 's/\"\"; gene_id \"\"/\"; gene_id \"/g' RRPM_transcripts.gtf
perl -pi -e 's/\"\";\"/\";/g' RRPM_transcripts.gtf
```

### Auriculariopsis ampla

```
scripts/aampla/FILTERING_aampla.R
```

### Coprinopsis cinerea

```
scripts/ccinerea/FILTERING_ccinerea.R
```

### Cryptococcus neoformans

```
scripts/cneoformans/FILTERING_cneoformans.R
```

### Lentinus tigrinus

```
scripts/ltigrinus/FILTERING_ltigrinus.R
```

### Phanerochaete chrysosporium

```
scripts/pchrysosporium/FILTERING_pchrysosporium.R
```

### Rickenella mellea

```
scripts/rmellea/FILTERING_rmellea.R
```

### Schizophyllum commune

```
scripts/scommune/FILTERING_scommune.R
```

### Ustilago maydis

```
scripts/umaydis/FILTERING_umaydis.R
```

## Fusion correction

Megvizsgáljuk, hogy melyik azok a gének amik fúziósak voltak az eredeti annotációba de már két külön gént alkotnak, majd ezeket a géneket manuális IGV megtekintés után átnevezzük v1, v2 stb.

### Armillaria ostoyae

```
scripts/aostoyae/FUSION_FILTER_aostoyae.R
# c("AROS_19248", "AROS_20796")
```

### Auriculariopsis ampla

```
scripts/aampla/FUSION_FILTER_aampla.R
# c("411782", "505059")
```

### Coprinopsis cinerea

```
scripts/ccinerea/FUSION_FILTER_ccinerea.R
# 371210
```

### Cryptococcus neoformans

```
scripts/cneoformans/FUSION_FILTER_cneoformans.R
# nincs
```

### Lentinus tigrinus

```
scripts/ltigrinus/FUSION_FILTER_ltigrinus.R
# c("599196", "599676")
```

### Phanerochaete chrysosporium

```
scripts/pchrysosporium/FUSION_FILTER_pchrysosporium.R
# 2916565
```

### Rickenella mellea

```
scripts/rmellea/FUSION_FILTER_rmellea.R
# c("842741", "846740")
```

### Schizophyllum commune

```
scripts/scommune/FUSION_FILTER_scommune.R
# c("2478457", "2614204", "2620485", "2644252", "2735855")
```

### Ustilago maydis

```
scripts/umaydis/FUSION_FILTER_umaydis.R
# nincs
```

## Alternative Splicing Statistics

Kinyerjük az alternative splicing statisztikákat majd az output ASpli_binFeatures.log fájlt átnevezzük xy_ASpli_binFeatures.log-ra ahol az xy a faj neve. Ezek a fájlok ezen a néven érhetőek el az adott faj könyvtárába.

### Armillaria ostoyae

```
scripts/aostoyae/AS_STATS_aostoyae.R
```

### Auriculariopsis ampla

```
scripts/aampla/AS_STATS_aampla.R
```

### Coprinopsis cinerea

```
scripts/ccinerea/AS_STATS_ccinerea.R
```

### Cryptococcus neoformans

```
scripts/cneoformans/AS_STATS_cneoformans.R
```

### Lentinus tigrinus

```
scripts/lrigrinus/AS_STATS_lrigrinus.R
```

### Phanerochaete chrysosporium

```
scripts/pchrysosporium/AS_STATS_pchrysosporium.R
```

### Rickenella mellea

```
scripts/rmellea/AS_STATS_rmellea.R
```

### Schizophyllum commune

```
scripts/scommune/AS_STATS_scommune.R
```

### Ustilago maydis

```
scripts/umaydis/AS_STATS_umaydis.R
```

## Expression: STAR Alignment, Cuffquant Assembly and Cuffdiff analysis

Futtatunk egy teljes STAR illesztést majd egy Cuffquantot, majd kinyerjük az expressziós értékeket a Cuffdiff segítségével.

### Armillaria ostoyae

```
Cufflinks_scripts/aostoyae_STAR_CUFF_expression.sh
```

### Auriculariopsis ampla

```
Cufflinks_scripts/aampla_STAR_CUFF_expression.sh
```

### Coprinopsis cinerea

```
Cufflinks_scripts/ccinerea_STAR_CUFF_expression.sh
```

### Cryptococcus neoformans

```
Cufflinks_scripts/cneoformans_STAR_CUFF_expression.sh
```

### Lentinus tigrinus

```
Cufflinks_scripts/ltigrinus_STAR_CUFF_expression.sh
```

### Phanerochaete chrysosporium

```
Cufflinks_scripts/pchrysosporium_STAR_CUFF_expression.sh
```

### Rickenella mellea

```
Cufflinks_scripts/rmellea_STAR_CUFF_expression.sh
```

### Schizophyllum commune

```
Cufflinks_scripts/scommune_STAR_CUFF_expression.sh
```

### Ustilago maydis

Túl sok a short read miatti unmapped MEGOLDÁS -> https://github.com/alexdobin/STAR/issues/169

```
Cufflinks_scripts/umaydis_STAR_CUFF_expression.sh
```

## ORF prediction and InterProScan analysis

Prediktáljuk az ORF régiókat TransDecoder segítségével, majd kiemeljük csak a headereket egy külön fájlba (headers.cds, ezt akkor még manuálisan csináltam de mostmár beleépíteném az R scriptbe).

Egy R script segítségével kiszedjük a longest ORF-hez tartalmazó headereket majd ezzel megszűrjük a fehérjeszekvenciákat tartalmazó fájlt, végül átalakítjuk úgy a FASTA-fájlt, hogy az InterProScan számára megfelelő legyen. (Elég körülményesen oldottam ezt meg de működik, nyilván mostmár máshogy írnám még ezt a részt egy R scriptbe)

### Armillaria ostoyae

```
# ORF régiók prediktálása
mkdir ~/Desktop/alternative_splicing/aostoyae/aostoyae_enrichment_analysis
cd ~/Desktop/alternative_splicing/aostoyae/aostoyae_enrichment_analysis
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ../aostoyae_genome/aostoyae_AS_annotation.gtf ../aostoyae_genome/p3_i2_t47428_Arm_ostoy_v2.scaf > aostoyae_transcripts.fasta
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ../aostoyae_genome/aostoyae_AS_annotation.gtf  > ../aostoyae_genome/aostoyae_AS_annotation.gff3
TransDecoder.LongOrfs -m 20 -S -t aostoyae_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/aostoyae/INTERPROSCAN_aostoyae.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> aostoyae_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' aostoyae_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' aostoyae_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' aostoyae_proteins_all.fasta
```

```
./interproscan.sh -i aostoyae_proteins_all.fasta -f tsv --iprlookup --goterms
```

### Auriculariopsis ampla

```
# ORF régiók prediktálása
mkdir ~/Desktop/alternative_splicing/aampla/aampla_enrichment_analysis
cd ~/Desktop/alternative_splicing/aampla/aampla_enrichment_analysis
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ../aampla_genome/aampla_AS_annotation.gtf ../aampla_genome/Auramp1_AssemblyScaffolds.fasta > aampla_transcripts.fasta
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ../aampla_genome/aampla_AS_annotation.gtf  > ../aampla_genome/aampla_AS_annotation.gff3
TransDecoder.LongOrfs -m 20 -S -t aampla_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/aampla/INTERPROSCAN_aampla.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> aampla_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' aampla_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' aampla_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' aampla_proteins_all.fasta
```

```
interproscan.sh -i aampla_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```

### Coprinopsis cinerea

```
# ORF régiók prediktálása
mkdir ~/Desktop/alternative_splicing/ccinerea_AmutBmut/ccinerea_AmutBmut_enrichment_analysis
cd ~/Desktop/alternative_splicing/ccinerea_AmutBmut/ccinerea_AmutBmut_enrichment_analysis
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ../ccinerea_AmutBmut_genome/ccinerea_AmutBmut_AS_annotation.gtf ../ccinerea_AmutBmut_genome/Copci_AmutBmut1_AssemblyScaffolds.fasta > ccinerea_AmutBmut_transcripts.fasta
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ../ccinerea_AmutBmut_genome/ccinerea_AmutBmut_AS_annotation.gtf  > ../ccinerea_AmutBmut_genome/ccinerea_AmutBmut_AS_annotation.gff3
TransDecoder.LongOrfs -m 20 -S -t ccinerea_AmutBmut_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/ccinerea/INTERPROSCAN_ccinerea.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> ccinerea_AmutBmut_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' ccinerea_AmutBmut_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' ccinerea_AmutBmut_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' ccinerea_AmutBmut_proteins_all.fasta
```

```
./interproscan.sh -i ccinerea_AmutBmut_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```

### Cryptococcus neoformans

```
# ORF régiók prediktálása
cd ~/Desktop/alternative_splicing/cneoformans
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ./cneoformans_genome/cneoformans_AS_annotation.gtf ./cneoformans_genome/Cryptococcus_neoformans_H99.unmasked.fasta > cneoformans_transcripts.fasta
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ./cneoformans_genome/cneoformans_AS_annotation.gtf  > ./cneoformans_genome/cneoformans_AS_annotation.gff3
TransDecoder.LongOrfs -m 20 -S -t cneoformans_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/cneoformans/INTERPROSCAN_cneoformans.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> cneoformans_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' cneoformans_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' cneoformans_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' cneoformans_proteins_all.fasta
```

```
interproscan.sh -i cneoformans_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```

### Lentinus tigrinus

```
# ORF régiók prediktálása
mkdir ~/Desktop/alternative_splicing/ltigrinus/ltigrinus_enrichment_analysis
cd ~/Desktop/alternative_splicing/ltigrinus/ltigrinus_enrichment_analysis
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ../ltigrinus_genome/ltigrinus_AS_annotation.gtf ../ltigrinus_genome/Sisbr1_AssemblyScaffolds.fasta > ltigrinus_transcripts.fasta
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ../ltigrinus_genome/ltigrinus_AS_annotation.gtf  > ../ltigrinus_genome/ltigrinus_AS_annotation.gff3
TransDecoder.LongOrfs -m 20 -S -t ltigrinus_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/ltigrinus/INTERPROSCAN_ltigrinus.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> ltigrinus_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' ltigrinus_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' ltigrinus_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' ltigrinus_proteins_all.fasta
```

```
interproscan.sh -i ltigrinus_proteins_all.fasta -f tsv --iprlookup --goterms
```

### Phanerochaete chrysosporium

```
# ORF régiók prediktálása
mkdir ~/Desktop/alternative_splicing/pchrysosporium/pchrysosporium_enrichment
cd ~/Desktop/alternative_splicing/pchrysosporium/pchrysosporium_enrichment
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ../pchrysosporium_genome/pchrysosporium_AS_annotation.gtf  > ../pchrysosporium_genome/pchrysosporium_AS_annotation.gff3
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ../pchrysosporium_genome/pchrysosporium_AS_annotation.gtf ../pchrysosporium_genome/Phchr2_AssemblyScaffolds.fasta > pchrysosporium_transcripts.fasta
TransDecoder.LongOrfs -m 20 -S -t pchrysosporium_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/pchrysosporium/INTERPROSCAN_pchrysosporium.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> pchrysosporium_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' pchrysosporium_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' pchrysosporium_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' pchrysosporium_proteins_all.fasta
```

```
interproscan.sh -i pchrysosporium_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```

### Rickenella mellea

```
# ORF régiók prediktálása
cd ~/Desktop/alternative_splicing/rmellea
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ./rmellea_genome/rmellea_AS_annotation.gtf ./rmellea_genome/Ricmel1_AssemblyScaffolds.fasta > rmellea_transcripts.fasta
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ./rmellea_genome/rmellea_AS_annotation.gtf  > ./rmellea_genome/rmellea_AS_annotation.gff3
TransDecoder.LongOrfs -m 20 -S -t rmellea_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/rmellea/INTERPROSCAN_rmellea.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> rmellea_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' rmellea_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' rmellea_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' rmellea_proteins_all.fasta
```

```
./interproscan.sh -i rmellea_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```

### Schizophyllum commune

```
# ORF régiók prediktálása
mkdir ~/Desktop/alternative_splicing/scommune/scommune_enrichment_analysis
cd ~/Desktop/alternative_splicing/scommune/scommune_enrichment_analysis
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ../scommune_genome/scommune_AS_annotation.gtf ../scommune_genome/Schco3_AssemblyScaffolds.fasta > scommune_transcripts.fasta
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ../scommune_genome/scommune_AS_annotation.gtf  > ../scommune_genome/scommune_AS_annotation.gff3
TransDecoder.LongOrfs -m 20 -S -t scommune_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/scommune/INTERPROSCAN_scommune.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> scommune_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' scommune_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' scommune_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' scommune_proteins_all.fasta
```

```
./interproscan.sh -i scommune_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```

### Ustilago maydis

```
# ORF régiók prediktálása
cd ~/Desktop/alternative_splicing/umaydis
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ./umaydis_genome/umaydis_AS_annotation.gtf ./umaydis_genome/p3_t237631_Ust_maydi_v2GB.scaf > umaydis_transcripts.fasta
~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ./umaydis_genome/umaydis_AS_annotation.gtf  > ./umaydis_genome/umaydis_AS_annotation.gff3
TransDecoder.LongOrfs -m 20 -S -t umaydis_transcripts.fasta

# kiszedjük a fasta headereket
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds

# R script
scripts/umaydis/INTERPROSCAN_umaydis.R

# sublime text átalakítás a longest_orfs.pep fileon, hogy egy sorba legyen a header is meg a szekvencia is
\)\n --> \)\t

# megszűrjük a filter IDjeivel a fehérje FASTA fájlt
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> umaydis_proteins_all.fasta; done <filter

# FASTA formátum helyreállítása sublime-al
\)\t --> \)\n

# *-ok törlése mivel az InterProScan nem tudja értelmezni
perl -pi -e 's/\*//g' umaydis_proteins_all.fasta

# kis átalakítás a könnyebb áttekinthetőség végett
perl -pi -e 's/::g\..*$//g' umaydis_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' umaydis_proteins_all.fasta
```

```
interproscan.sh -i umaydis_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```

## AS transcripts table (OUTDATED)

!!!!! FONTOS !!!!!

A CDS régiók, továbbá a domain összetétel összehasonlításakor primary transcriptnek az lett megadva amelyik a legmagasabb expressziós értékkel rendelkezik a transzkriptek közül összességébem ezért ez az elemzés az új primary transcript kiválasztása miatt outdated!

Ez a TSV fájl tartalmazza az AS-ben résztvevő szűrt transcripteket.

A táblázat csak olyan transzkriptt tartalmaz:
- Aminek az anyagénje részt vesz alternative splicingban
- Aminek az anyagénje legalább az egyik fejlődési fázsiban eléri a  2 FPKM-es értéket
- Ami legalább egy egyik fázisban eléri vagy meghaladja az anyagén adott fázisban lévő teljes expressziós értékének legalább 10%-át

A táblázat tartalmazza a következő adatokat:
- Expressziós értékek
- Fold Change értékek
- A CDS régiókra vonatkozó értékek és összehasonlítások
- Domain összetételek és azok összehasonlítása
- GO annotációk

### Armillaria ostoyae

```
scripts/aostoyae/AS_TRANSCRIPTS_TABLE_aostoyae.R
```

### Auriculariopsis ampla

```
scripts/aampla/AS_TRANSCRIPTS_TABLE_aampla.R
```

### Coprinopsis cinerea

```
scripts/ccinerea/AS_TRANSCRIPTS_TABLE_ccinerea.R
```

### Cryptococcus neoformans

```
scripts/cneoformans/AS_TRANSCRIPTS_TABLE_cneoformans.R
```

### Lentinus tigrinus

```
scripts/ltigrinus/AS_TRANSCRIPTS_TABLE_ltigrinus.R
```

### Phanerochaete chrysosporium

```
scripts/pchyrsosporium/AS_TRANSCRIPTS_TABLE_pchyrsosporium.R
```

### Rickenella mellea

```
scripts/rmellea/AS_TRANSCRIPTS_TABLE_rmellea.R
```

### Schizophyllum commune

```
scripts/scommune/AS_TRANSCRIPTS_TABLE_scommune.R
```

### Ustilago maydis

```
scripts/umaydis/AS_TRANSCRIPTS_TABLE_umaydis.R
```

## Enrichment (OUTDATED)

Enrichment a GeneCatalog génjeire, mivel a klaszerezés másmilyen ezért ezek a scriptek erősen megkérdőjelezhetőek, frissíteni kell a futtatásokat.

### Armillaria ostoyae

```
interproscan.sh -i p3_i2_t47428_Arm_ostoy_v2.prot -f tsv --iprlookup --goterms

aostoyae_enrichment_preparations.R

aostoyae_genelvl_enrichment_groups.R
```

### Auriculariopsis ampla

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Auramp1_GeneCatalog_proteins_20160719.aa.fasta

# >jgi\|Auramp1\| --> ""
# \| --> \t

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Auramp1_GeneCatalog_proteins_20160719.aa.fasta -f tsv --iprlookup --goterms -dp

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Auramp1\| --> ""
# \| --> \t

aampla_enrichment_dictionary.R

aampla_enrichment_preparations.R

aampla_genelvl_enrichment_groups.R
```

### Coprinopsis cinerea

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta -f tsv --iprlookup --goterms

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Copci_AmutBmut1\| --> ""
# \| --> \t

ccinerea_AmutBmut_enrichment_dictionary.R

ccinerea_AmutBmut_enrichment_preparations.R

ccinerea_AmutBmut_genelvl_enrichment_groups.R
```

### Cryptococcus neoformans

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Cryptococcus_neoformans_H99.proteins.fasta

# jgi\|Cryne_H99_1\| --> ""
# \|.*$ --> ""

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Cryptococcus_neoformans_H99.proteins.fasta -f tsv --iprlookup --goterms -dp

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Cryptococcus_neoformans_H99.proteins.fasta átalakítatlan headerjeit -> Cryptococcus_neoformans_H99.proteins.fasta.headers

# >jgi\|Cryne_H99_1\| --> ""
# \| --> \t

cneoformans_enrichment_dictionary.R

cneoformans_enrichment_preparations.R

cneoformans_genelvl_enrichment_groups.R
```

### Lentinus tigrinus

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Sisbr1_GeneCatalog_proteins_20130805.aa.fasta

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Sisbr1_GeneCatalog_proteins_20130805.aa.fasta -f tsv --iprlookup --goterms

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Sisbr1\|--> ""
# \| --> \t

ltigrinus_enrichment_dictionary.R

ltigrinus_enrichment_preparations.R

ltigrinus_genelvl_enrichment_groups.R
```

### Phanerochaete chrysosporium

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Phchr2_GeneCatalog_proteins_20131210.aa.fasta

# >jgi\|Phchr2\| --> ""
# \|.*$ --> ""

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Phchr2_GeneCatalog_proteins_20131210.aa.fasta -f tsv --iprlookup --goterms -dp

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Phchr2\| --> ""
# \| --> \t

pchrysosporium_enrichment_dictionary.R

pchrysosporium_enrichment_preparations.R

pchrysosporium_genelvl_enrichment_groups.R
```

### Rickenella mellea

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Ricmel1_GeneCatalog_proteins_20151108.aa.fasta

# jgi\|Ricmel1\| --> ""
# \|.*$ --> ""

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Ricmel1_GeneCatalog_proteins_20151108.aa.fasta -f tsv --iprlookup --goterms -dp

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Ricmel1_GeneCatalog_proteins_20151108.aa.fasta átalakítatlan headerjeit -> Ricmel1_GeneCatalog_proteins_20151108.aa.fasta.headers

# >jgi\|Ricmel1\| --> ""
# \| --> \t

rmellea_enrichment_dictionary.R

rmellea_enrichment_preparations.R

rmellea_genelvl_enrichment_groups.R
```

### Schizophyllum commune

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' Schco3_GeneCatalog_proteins_20130812.aa.fasta

# lekell futtatni az IPR scant a genecatalog génekre

interproscan.sh -i Schco3_GeneCatalog_proteins_20130812.aa.fasta -f tsv --iprlookup --goterms

# mivel nem stimmelnek az ID-k a fasta fájl proteinjei és az alap GTF fájl ID-jai között ezért meg kell egymással feleltetni őket
# ehhez először kivágjuk a Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta átalakítatlan headerjeit -> Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers

# >jgi\|Schco3\| --> ""
# \| --> \t

scommune_enrichment_dictionary.R

scommune_enrichment_preparations.R

scommune_genelvl_enrichment_groups.R
```

### Ustilago maydis

```
# *-ok törlése mivel az InterProScan nem tudja értelmezni

perl -pi -e 's/\*//g' p3_t237631_Ust_maydi_v2GB.prot

interproscan.sh -i p3_t237631_Ust_maydi_v2GB.prot -f tsv --iprlookup --goterms -dp

umaydis_enrichment_preparations.R

umaydis_genelvl_enrichment_groups.R
```

## MCL Clustering a géncsaládok megkeresésére

```
# megkeressük a géncsaládokat amikből tovább tudunk haladni a primary transcript megtalálása felé

aampla 			-> 	Auramp1_GeneCatalog_proteins_20160719.aa.fasta
aostoyae 		->	p3_i2_t47428_Arm_ostoy_v2.prot
ccinerea 		->	Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta
ltigrinus 		->	Sisbr1_GeneCatalog_proteins_20130805.aa.fasta
pchrysosporium 	->	Phchr2_GeneCatalog_proteins_20131210.aa.fasta
rmellea			->	Ricmel1_GeneCatalog_proteins_20151108.aa.fasta
scommune 		->	Schco3_GeneCatalog_proteins_20130812.aa.fasta

# orthomclAdjustFasta

cd compliantFasta
orthomclAdjustFasta aampla ../aampla.fasta 1
orthomclAdjustFasta aostoyae ../aostoyae.fasta 1
orthomclAdjustFasta ccinerea ../ccinerea.fasta 1
orthomclAdjustFasta ltigrinus ../ltigrinus.fasta 1
orthomclAdjustFasta pchrysosporium ../pchrysosporium.fasta 1
orthomclAdjustFasta rmellea ../rmellea.fasta 1
orthomclAdjustFasta scommune ../scommune.fasta 1

# orthomclFilterFasta

orthomclFilterFasta ./compliantFasta 10 20

# all VS all BLAST

mpiformatdb -N 32 -i goodProteins.fasta -o T
mv goodProteins.fasta.* ~/share/
nohup mpirun -n 32 mpiblast -d goodProteins.fasta -i goodProteins.fasta -p blastp -m 8 -o goodProteins_blasted --copy-via=none &

# orthomclBlastParser

orthomclBlastParser goodProteins_blasted ./compliantFasta >> similarSequences.txt


# mivel nem stimmelnek a transcript_id-k a protein_id-s baszakodás miatt ezért az egyes fajok dictionaryje alapján át kell alakítani a családokat tartalmazó fájlt

cluster_transformation.R
```

## Summary (OUTDATED)

```
# egyenként megcsináljuk a különböző fajok summaryjeit

summary_scripts/aampla_summary.R
summary_scripts/aostoyae_summary.R
summary_scripts/ccinerea_AmutBmut_summary.R
summary_scripts/cneoformans_summary.R
summary_scripts/ltigrinus_summary.R
summary_scripts/pchrysosporium_summary.R
summary_scripts/rmellea_summary.R
summary_scripts/scommune_summary.R
summary_scripts/umaydis_summary.R

# összesítjük a summaryket 

summary_scripts/summary_all.R
```

## Silix/Hifix és enrichment (OUTDATED)

```
# a clusterezés megkezdéséhez egy fileba kell olvasztani az összes különálló faj fasta fájlját
# hozzá kell adni az adott faj headerjéhez a faj nevét hogy meg lehessen különböztetni

aampla 			-> 	Auramp1_GeneCatalog_proteins_20160719.aa.fasta
aostoyae 		->	p3_i2_t47428_Arm_ostoy_v2.prot
ccinerea 		->	Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta
cneoformans		->	Cryptococcus_neoformans_H99.proteins.fasta
ltigrinus 		->	Sisbr1_GeneCatalog_proteins_20130805.aa.fasta
pchrysosporium 	->	Phchr2_GeneCatalog_proteins_20131210.aa.fasta
rmellea			->	Ricmel1_GeneCatalog_proteins_20151108.aa.fasta
scommune 		->	Schco3_GeneCatalog_proteins_20130812.aa.fasta
umaydis			->	p3_t237631_Ust_maydi_v2GB.prot

# mpiblast database

mpiformatdb -N 48 -i all_species_protein.fasta -o T

# áthelyezzük a fájlokat ahogy kell

mpirun -n 48 mpiblast -d all_species_protein.fasta -i all_species_protein.fasta -p blastp -m 8 -o all_species_protein_blasted --copy-via=none

# Silix a Hifix-hez

silix all_species_protein.fasta all_species_protein_blasted --net > all_species_protein_SLX.fnodes

# Hifix KITERJESZTÉSNEK .fasta-nak kell lennie!

hifix -t 48 all_species_protein.fasta all_species_protein_blasted.net all_species_protein_SLX.fnodes > all_species_protein_HFX.fnodes

# p3_i2_t47428_Arm_ostoy_v2_ --> semmire csere

# mivel nem stimmelnek a transcript_id-k a protein_id-s baszakodás miatt ezért az egyes fajok dictionaryje alapján át kell alakítani a családokat tartalmazó fájlt

cluster_transformation.R

# fajonkénti enrichment

enrichment_scripts/aampla_genelvl_enrichment_fisher_classic.R
enrichment_scripts/aostoyae_genelvl_enrichment_fisher_classic.R
enrichment_scripts/ccinerea_AmutBmut_genelvl_enrichment_fisher_classic.R
enrichment_scripts/cneoformans_genelvl_enrichment_fisher_classic.R
enrichment_scripts/ltigrinus_genelvl_enrichment_fisher_classic.R
enrichment_scripts/pchrysosporium_genelvl_enrichment_fisher_classic.R
enrichment_scripts/rmellea_genelvl_enrichment_fisher_classic.R
enrichment_scripts/scommune_genelvl_enrichment_fisher_classic.R
enrichment_scripts/umaydis_genelvl_enrichment_fisher_classic.R

# összesítés

enrichment_scripts/families_count.R
enrichment_scripts/GO_count.R
```