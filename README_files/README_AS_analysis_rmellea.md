# Rickenella mellea
## RNA-seq files
## Preparations
A bemeneti annotációs fájlokat át kell alakítani, hogy megfelelőek legyenek a Cufflinks számára illetve elő kell állítani azokat az annotációs és FASTA fájlokat amik az RRPM analízishez szükségesek. 
### Input:
- __Ricmel1_GeneCatalog_genes_20151108.gff__ : Forrás a JGI
- __Ricmel1_AssemblyScaffolds.fasta__ : Forrás a JGI
### Output:
- __rmellea_onlygene.gtf__ : Csak a géneket tartalmazó annotációs fájl
- __rmellea_onlyexon.gtf__ : Csak az exonokat tartalmazó annotációs fájl
- __rmellea_fixed.gtf__ : Pozíciófixált annotációs fájl az RRPM számára
- __rmellea_genes.fasta__ : Géneket tartalmazó fasta fájl

GTF fájlok elkészítése az RRPM számára
```
PREPARATIONS_rmellea.R
```
GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
```
gtf2bed < rmellea_onlygene.gtf > rmellea_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl elkészítése
```
bedtools getfasta -name -fo rmellea_genes.fasta -fi Ricmel1_AssemblyScaffolds.fasta -bed rmellea_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl headerjének a trimmelése
```
perl -pi -e 's/::.*//g' rmellea_genes.fasta
```
## RRPM: STAR Alignment and Cufflinks Assembly
Megnézzük milyen hosszúak a gének és az intronok a STAR és a Cufflinks beállításához.

Elkészítünk egy fájlt ami tartalmazza az összes szükséges scriptet ami az RRPM futtatásához szükséges (STAR, Cufflinks)
### Input:
- __rmellea_genes.fasta__
- __rmellea_fixed.gtf__
### Output:
- __gene_length__
- A __Cufflinks__ __outputja__ ami tartalmazza az RRPM analízis fájljait
```
bioawk -c fastx '{ print $name, length($seq) }' < rmellea_genes.fasta > gene_length
INTRON_LENGTH_rmellea.R
```
A maximális transzkriptméret: 13037 --> max-bundle-length marad 250000

Intron min: 7 --> --alignIntronMin marad 3 

Intron max: 1984 --> --alignIntronMax marad 30000

Lefuttatjuk az RRPM analízist.
```
# STAR illesztéshez szükséges indexelés

mkdir rmellea_DIR_RRPM
cp -r ./rmellea_genome ./rmellea_DIR_RRPM/rmellea_RRPM_STARindex
cd ./rmellea_DIR_RRPM/rmellea_RRPM_STARindex
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir ../rmellea_RRPM_STARindex --sjdbGTFfile ../rmellea_RRPM_STARindex/rmellea_fixed.gtf --genomeFastaFiles ../rmellea_RRPM_STARindex/rmellea_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./rmellea_DIR_RRPM/rmellea_RRPM_STARout
cd ./rmellea_DIR_RRPM/rmellea_RRPM_STARout
STAR --runThreadN 48 --genomeDir ../rmellea_RRPM_STARindex --readFilesIn ../../raw_data/R1_VM_S75_R1_001.fastq.gz,../../raw_data/R2_VM_S76_R1_001.fastq.gz,../../raw_data/R3_VM_S77_R1_001.fastq.gz,../../raw_data/KP1_S1_R1_001.fastq.gz,../../raw_data/KP3_S3_R1_001.fastq.gz,../../raw_data/P2_S5_R1_001.fastq.gz,../../raw_data/RYFB2_1_S24_R1_001.fastq.gz,../../raw_data/RYFB_4_S25_R1_001.fastq.gz,../../raw_data/RYFB_5_S26_R1_001.fastq.gz,../../raw_data/RFB_K_2_S7_R1_001.fastq.gz,../../raw_data/RFB_K_5_S8_R1_001.fastq.gz,../../raw_data/RFB_K_9_S9_R1_001.fastq.gz,../../raw_data/RFB_T_4_10_S12_R1_001.fastq.gz,../../raw_data/RFB_T_7_S13_R1_001.fastq.gz,../../raw_data/RFB_T_S10_R1_001.fastq.gz ../../raw_data/R1_VM_S75_R2_001.fastq.gz,../../raw_data/R2_VM_S76_R2_001.fastq.gz,../../raw_data/R3_VM_S77_R2_001.fastq.gz,../../raw_data/KP1_S1_R2_001.fastq.gz,../../raw_data/KP3_S3_R2_001.fastq.gz,../../raw_data/P2_S5_R2_001.fastq.gz,../../raw_data/RYFB2_1_S24_R2_001.fastq.gz,../../raw_data/RYFB_4_S25_R2_001.fastq.gz,../../raw_data/RYFB_5_S26_R2_001.fastq.gz,../../raw_data/RFB_K_2_S7_R2_001.fastq.gz,../../raw_data/RFB_K_5_S8_R2_001.fastq.gz,../../raw_data/RFB_K_9_S9_R2_001.fastq.gz,../../raw_data/RFB_T_4_10_S12_R2_001.fastq.gz,../../raw_data/RFB_T_7_S13_R2_001.fastq.gz,../../raw_data/RFB_T_S10_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 48 -g ./rmellea_DIR_RRPM/rmellea_RRPM_STARindex/rmellea_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./rmellea_DIR_RRPM/rmellea_RRPM_CUFFout ./rmellea_DIR_RRPM/rmellea_RRPM_STARout/Aligned.sortedByCoord.out.bam
```
## Filtering and Merge
Az RRPM outputot filterezzük, majd összeillesztjük az eredeti annotációs fájllal, hogy kitöltsük az esetleges hézagokat, ahol nem mutatott expressziós értéket az eredeti annotációban szereplő transzkript.

A GTF-ekbe átalakítja a start/stop helyeket pl.: 73e3 = 73000, ezeket kiszűrtem kézzel IGV segítségével, a beolvasásnál ha ezeket integerről stringre állítjuk akkor nem írja át, de én nem bajlodtam már vele miután lefutott.

Lépések:
- Expression Filtering: Kitörlünk minden olyan transzkriptet aminek nincsen látható expressziója
- Full Read Support Filtering: Kitörülünk minden olyan transzkriptet aminek nincs minden intron-exon junctionjére legalább 1 read
- Context Restoration: Visszaállítjuk a GTF fájl start és end pozíciójat az eredeti annotáció pozícióihoz
- Strand Filtering: Kitörlünk minden olyan transzkriptet ami nem ugyanazon a stranden van mint az alap gén (ezeket nem tudjuk biztosan megmondani strand specifikus readek nélkül)
- Annotation Merge: Összeillesztjük az eredeti annotációs fájlal

### Input:
- __rmellea_onlygene.gtf__
- __transcripts.gtf__ : RRPM Cufflinks assembly output file
### Output:
- __rmellea_RRPM_transcripts.gtf__ : Az újonnan felfedezett transzkripteket tartalmazó GTF
- __rmellea_AS_annotation.gtf__ : Az eredeti annotációval történő összemergeelés, hogy az esetlegesen nem detektált gének is bennelegyenek az annotációs fájlban.
- __rmellea_stats.log__ : Szűrési statisztikákat tartalmazó logfájl
```
FILTERING_rmellea.R
```
## Fusion correction
### Input:
- __isoforms.fpkm_tracking__ : Az RRPM Cufflinks output izoforms FPKM értékeket tartalmazó fájlja
- __rmellea_AS_annotation.gtf__
### Output:
- Két fúziós gént detektált a script: __842741__, __846740__

Megvizsgáljuk, hogy melyek azok a gének amik fúziósak voltak az eredeti annotációba de már két külön gént alkotnak, majd ezeket a géneket manuális IGV megtekintés után átnevezzük v1, v2 satöbbire.
```
FUSION_FILTER_rmellea.R
```
## Expression: STAR Alignment, Cuffquant Assembly and Cuffdiff analysis
### Output:
- A __CuffDiff__ __output__ fájljai az expreszziós analízisből

Futtatunk egy teljes STAR illesztést majd egy Cuffquantot, majd kinyerjük az expressziós értékeket a CuffDiff segítségével.
```
# STAR INDEXING

mkdir DIR_rmellea_expressionanalysis
cp -r ./rmellea_genome ./DIR_rmellea_expressionanalysis/rmellea_STARindex
cd ./DIR_rmellea_expressionanalysis/rmellea_STARindex
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../rmellea_STARindex/rmellea_AS_annotation.gtf --genomeFastaFiles ./Ricmel1_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# RM_VM

mkdir ./DIR_rmellea_expressionanalysis/RM_VM_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_VM_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/R1_VM_S75_R1_001.fastq.gz ../../raw_data/R1_VM_S75_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_VM_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_VM_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/R2_VM_S76_R1_001.fastq.gz ../../raw_data/R2_VM_S76_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_VM_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_VM_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/R3_VM_S77_R1_001.fastq.gz ../../raw_data/R3_VM_S77_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_VM_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_VM_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_VM_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# RM_P

mkdir ./DIR_rmellea_expressionanalysis/RM_P_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_P_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/KP1_S1_R1_001.fastq.gz ../../raw_data/KP1_S1_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_P_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_P_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/KP3_S3_R1_001.fastq.gz ../../raw_data/KP3_S3_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_P_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_P_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/P2_S5_R1_001.fastq.gz ../../raw_data/P2_S5_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_P_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_P_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_P_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_P_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_P_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_P_R3_STARout/Aligned.sortedByCoord.out.bam

# RM_YFB

mkdir ./DIR_rmellea_expressionanalysis/RM_YFB_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_YFB_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RYFB2_1_S24_R1_001.fastq.gz ../../raw_data/RYFB2_1_S24_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_YFB_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_YFB_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RYFB_4_S25_R1_001.fastq.gz ../../raw_data/RYFB_4_S25_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_YFB_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_YFB_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RYFB_5_S26_R1_001.fastq.gz ../../raw_data/RYFB_5_S26_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_YFB_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_YFB_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_YFB_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_YFB_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_YFB_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_YFB_R3_STARout/Aligned.sortedByCoord.out.bam

# RM_FB_K

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_K_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_K_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_K_2_S7_R1_001.fastq.gz ../../raw_data/RFB_K_2_S7_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_K_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_K_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_K_5_S8_R1_001.fastq.gz ../../raw_data/RFB_K_5_S8_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_K_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_K_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_K_9_S9_R1_001.fastq.gz ../../raw_data/RFB_K_9_S9_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_K_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_K_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_K_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_K_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_K_R3_STARout/Aligned.sortedByCoord.out.bam

# RM_FB_T

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_T_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_T_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_T_4_10_S12_R1_001.fastq.gz ../../raw_data/RFB_T_4_10_S12_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_T_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_T_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_T_7_S13_R1_001.fastq.gz ../../raw_data/RFB_T_7_S13_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_T_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_T_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_T_S10_R1_001.fastq.gz ../../raw_data/RFB_T_S10_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_T_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_T_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_T_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_T_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 -L VM,P,YFB,FB_K,FB_T --upper-quartile-norm ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf \
./DIR_rmellea_expressionanalysis/RM_VM_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_VM_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_VM_R3_CUFFquant/abundances.cxb \
./DIR_rmellea_expressionanalysis/RM_P_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_P_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_P_R3_CUFFquant/abundances.cxb \
./DIR_rmellea_expressionanalysis/RM_YFB_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_YFB_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_YFB_R3_CUFFquant/abundances.cxb \
./DIR_rmellea_expressionanalysis/RM_FB_K_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_FB_K_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_FB_K_R3_CUFFquant/abundances.cxb \
./DIR_rmellea_expressionanalysis/RM_FB_T_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_FB_T_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_FB_T_R3_CUFFquant/abundances.cxb
```
## ORF prediction
### Input:
- __rmellea_AS_annotation.gtf__
- __Ricmel1_AssemblyScaffolds.fasta__
### Output:
- __rmellea_transcripts.fasta__ : A TransDecoder outputja ami tartalmazza az összes transzkript szekvenciáját fasta formátumban
- __longest_orfs.cds__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF cds régióját
- __longest_orfs.pep__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF proteinjeit

Prediktáljuk az ORF régiókat TransDecoder segítségével.
```
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl rmellea_AS_annotation.gtf Ricmel1_AssemblyScaffolds.fasta > rmellea_transcripts.fasta
TransDecoder.LongOrfs -m 20 -S -t rmellea_transcripts.fasta
```
Majd kiemeljük csak a headereket a longest_orfs.cds fájlból egy külön fájlba (headers.cds, ezt akkor még manuálisan csináltam de mostmár beleépíteném az R scriptbe). Majd regexp-el átalakítjuk a headereket, hogy fel tudja őket dolgozni az R script.
```
perl -pi -e 's/>//g' headers.cds
perl -pi -e 's/ /\t/g' headers.cds
perl -pi -e 's/\ttype:.*len:/\t/g' headers.cds
perl -pi -e 's/::/\t/g' headers.cds
```
Egy R script segítségével kiszedjük a leghosszabb ORF-ekhez tartalmazó headereket. 
```
INTERPROSCAN_rmellea.R
```
Majd ezzel megszűrjük a fehérjeszekvenciákat tartalmazó fájlt, végül átalakítjuk úgy a FASTA-fájlt, hogy az InterProScan számára megfelelő legyen. (Elég körülményesen oldottam ezt meg de működik, nyilván mostmár máshogy írnám még ezt a részt egy R scriptbe)
```
# átalakítás
\)\n --> \)\t
# szűrés
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> rmellea_proteins_all.fasta; done <filter
# visszalakítás
\)\t --> \)\n
```
Töröljük a fasta fájlból a *-okat mivel az InterProScan nem tudja értelmezni
```
perl -pi -e 's/\*//g' rmellea_proteins_all.fasta
```
Egy kis regexp átalakítás az áttekinthetőségért
```
perl -pi -e 's/::g\..*$//g' rmellea_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' rmellea_proteins_all.fasta
```
## InterProScan analysis
### Input:
- __rmellea_proteins_all.fasta__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinek
### Output:
- __rmellea_proteins_all.fasta.tsv__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinhez tartozó InterProScan domainek

Lefuttatjuk at InterProScan-t.
```
interproscan.sh -i rmellea_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```