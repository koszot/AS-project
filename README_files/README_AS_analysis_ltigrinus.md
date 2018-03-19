# Lentinus tigrinus
## RNA-seq files
## Preparations
A bemeneti annotációs fájlokat át kell alakítani, hogy megfelelőek legyenek a Cufflinks számára illetve elő kell állítani azokat az annotációs és FASTA fájlokat amik az RRPM analízishez szükségesek. 
### Input:
- __Lenti6_1_GeneCatalog_genes_20130903.gff__ : Forrás a JGI
- __Lenti6_1_AssemblyScaffolds.fasta__ : Forrás a JGI
### Output:
- __ltigrinus_onlygene.gtf__ : Csak a géneket tartalmazó annotációs fájl
- __ltigrinus_onlyexon.gtf__ : Csak az exonokat tartalmazó annotációs fájl
- __ltigrinus_fixed.gtf__ : Pozíciófixált annotációs fájl az RRPM számára
- __ltigrinus_genes.fasta__ : Géneket tartalmazó fasta fájl

GTF fájlok elkészítése az RRPM számára
```
PREPARATIONS_ltigrinus.R
```
GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
```
gtf2bed < ltigrinus_onlygene.gtf > ltigrinus_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl elkészítése
```
bedtools getfasta -name -fo ltigrinus_genes.fasta -fi Lenti6_1_AssemblyScaffolds.fasta -bed ltigrinus_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl headerjének a trimmelése
```
perl -pi -e 's/::.*//g' ltigrinus_genes.fasta
```
## RRPM: STAR Alignment and Cufflinks Assembly
Megnézzük milyen hosszúak a gének és az intronok a STAR és a Cufflinks beállításához.

Elkészítünk egy fájlt ami tartalmazza az összes szükséges scriptet ami az RRPM futtatásához szükséges (STAR, Cufflinks)
### Input:
- __ltigrinus_genes.fasta__
- __ltigrinus_fixed.gtf__
### Output:
- __gene_length__
- A __Cufflinks__ __outputja__ ami tartalmazza az RRPM analízis fájljait
```
bioawk -c fastx '{ print $name, length($seq) }' < ltigrinus_genes.fasta > gene_length
INTRON_LENGTH_ltigrinus.R
```
A maximális transzkriptméret: 39394 --> max-bundle-length marad 250000

Intron min: 1 --> --alignIntronMin marad 3 

Intron max: 38009 --> --alignIntronMax marad 30000

Lefuttatjuk az RRPM analízist.
```
# STAR illesztéshez szükséges indexelés

mkdir ltigrinus_DIR_RRPM
cp -r ./ltigrinus_genome ./ltigrinus_DIR_RRPM/ltigrinus_RRPM_STARindex
cd ./ltigrinus_DIR_RRPM/ltigrinus_RRPM_STARindex
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ../ltigrinus_RRPM_STARindex --sjdbGTFfile ../ltigrinus_RRPM_STARindex/ltigrinus_fixed.gtf --genomeFastaFiles ../ltigrinus_RRPM_STARindex/ltigrinus_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./ltigrinus_DIR_RRPM/ltigrinus_RRPM_STARout
cd ./ltigrinus_DIR_RRPM/ltigrinus_RRPM_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_RRPM_STARindex --readFilesIn ../../raw_data/LT_VM_R1.R1.fastq.gz,../../raw_data/LT_VM_R2.R1.fastq.gz,../../raw_data/LT_VM_R3.all.R1.fastq.gz,../../raw_data/LT_P1_R1.all.R1.fastq.gz,../../raw_data/LT_P1_R2.all.R1.fastq.gz,../../raw_data/LT_P1_R3.all.R1.fastq.gz,../../raw_data/LT_P2_K_R1.all.R1.fastq.gz,../../raw_data/LT_P2_K_R2.all.R1.fastq.gz,../../raw_data/LT_P2_K_R3.all.R1.fastq.gz,../../raw_data/LT_P2_T_R1.all.R1.fastq.gz,../../raw_data/LT_P2_T_R2.all.R1.fastq.gz,../../raw_data/LT_P2_T_R3.all.R1.fastq.gz,../../raw_data/LT_YFB_K_R1.all.R1.fastq.gz,../../raw_data/LT_YFB_K_R2.all.R1.fastq.gz,../../raw_data/LT_YFB_K_R3.all.R1.fastq.gz,../../raw_data/LT_YFB_T_R1.all.R1.fastq.gz,../../raw_data/LT_YFB_T_R2.all.R1.fastq.gz,../../raw_data/LT_YFB_T_R3.all.R1.fastq.gz,../../raw_data/LT_FB_K_R1.all.R1.fastq.gz,../../raw_data/LT_FB_K_R3.all.R1.fastq.gz,../../raw_data/TL_FB_K_R4_R1_001.fastq.gz,../../raw_data/LT_FB_L_R1.R1.fastq.gz,../../raw_data/LT_FB_L_R2.all.R1.fastq.gz,../../raw_data/TL_FB_L_R4_R1_001.fastq.gz,../../raw_data/LT_FB_T_R2.all.R1.fastq.gz,../../raw_data/LT_FB_T_R3.all.R1.fastq.gz,../../raw_data/TL_FB_T_R4_R1_001.fastq.gz ../../raw_data/LT_VM_R1.R2.fastq.gz,../../raw_data/LT_VM_R2.R2.fastq.gz,../../raw_data/LT_VM_R3.all.R2.fastq.gz,../../raw_data/LT_P1_R1.all.R2.fastq.gz,../../raw_data/LT_P1_R2.all.R2.fastq.gz,../../raw_data/LT_P1_R3.all.R2.fastq.gz,../../raw_data/LT_P2_K_R1.all.R2.fastq.gz,../../raw_data/LT_P2_K_R2.all.R2.fastq.gz,../../raw_data/LT_P2_K_R3.all.R2.fastq.gz,../../raw_data/LT_P2_T_R1.all.R2.fastq.gz,../../raw_data/LT_P2_T_R2.all.R2.fastq.gz,../../raw_data/LT_P2_T_R3.all.R2.fastq.gz,../../raw_data/LT_YFB_K_R1.all.R2.fastq.gz,../../raw_data/LT_YFB_K_R2.all.R2.fastq.gz,../../raw_data/LT_YFB_K_R3.all.R2.fastq.gz,../../raw_data/LT_YFB_T_R1.all.R2.fastq.gz,../../raw_data/LT_YFB_T_R2.all.R2.fastq.gz,../../raw_data/LT_YFB_T_R3.all.R2.fastq.gz,../../raw_data/LT_FB_K_R1.all.R2.fastq.gz,../../raw_data/LT_FB_K_R3.all.R2.fastq.gz,../../raw_data/TL_FB_K_R4_R2_001.fastq.gz,../../raw_data/LT_FB_L_R1.R2.fastq.gz,../../raw_data/LT_FB_L_R2.all.R2.fastq.gz,../../raw_data/TL_FB_L_R4_R2_001.fastq.gz,../../raw_data/LT_FB_T_R2.all.R2.fastq.gz,../../raw_data/LT_FB_T_R3.all.R2.fastq.gz,../../raw_data/TL_FB_T_R4_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 16 -g ./ltigrinus_DIR_RRPM/ltigrinus_RRPM_STARindex/ltigrinus_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./ltigrinus_DIR_RRPM/ltigrinus_RRPM_CUFFout ./ltigrinus_DIR_RRPM/ltigrinus_RRPM_STARout/Aligned.sortedByCoord.out.bam
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
- __ltigrinus_onlygene.gtf__
- __transcripts.gtf__ : RRPM Cufflinks assembly output file
### Output:
- __ltigrinus_RRPM_transcripts.gtf__ : Az újonnan felfedezett transzkripteket tartalmazó GTF
- __ltigrinus_AS_annotation.gtf__ : Az eredeti annotációval történő összemergeelés, hogy az esetlegesen nem detektált gének is bennelegyenek az annotációs fájlban.
- __ltigrinus_stats.log__ : Szűrési statisztikákat tartalmazó logfájl
```
FILTERING_ltigrinus.R
```
## Fusion correction
### Input:
- __isoforms.fpkm_tracking__ : Az RRPM Cufflinks output izoforms FPKM értékeket tartalmazó fájlja
- __ltigrinus_AS_annotation.gtf__
### Output:
- Két fúziós gént detektált a script: __599196__, __599676__

Megvizsgáljuk, hogy melyek azok a gének amik fúziósak voltak az eredeti annotációba de már két külön gént alkotnak, majd ezeket a géneket manuális IGV megtekintés után átnevezzük v1, v2 satöbbire.
```
FUSION_FILTER_ltigrinus.R
```
## Expression: STAR Alignment, Cuffquant Assembly and Cuffdiff analysis
### Output:
- A __CuffDiff__ __output__ fájljai az expreszziós analízisből

Futtatunk egy teljes STAR illesztést majd egy Cuffquantot, majd kinyerjük az expressziós értékeket a CuffDiff segítségével.
```
# STAR INDEXING

mkdir ltigrinus_DIR_expressionanalysis
cp -r ./ltigrinus_genome ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex
cd ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --genomeFastaFiles ./Sisbr1_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# LT_VM

mkdir ./ltigrinus_DIR_expressionanalysis/LT_VM_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_VM_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_VM_R1.R1.fastq.gz ../../raw_data/LT_VM_R1.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_VM_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_VM_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_VM_R2.R1.fastq.gz ../../raw_data/LT_VM_R2.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_VM_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_VM_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_VM_R3.all.R1.fastq.gz ../../raw_data/LT_VM_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_VM_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_VM_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_VM_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_P1

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P1_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P1_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P1_R1.all.R1.fastq.gz ../../raw_data/LT_P1_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P1_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P1_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P1_R2.all.R1.fastq.gz ../../raw_data/LT_P1_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P1_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P1_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P1_R3.all.R1.fastq.gz ../../raw_data/LT_P1_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P1_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P1_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P1_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P1_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P1_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P1_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_P2_K

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_K_R1.all.R1.fastq.gz ../../raw_data/LT_P2_K_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_K_R2.all.R1.fastq.gz ../../raw_data/LT_P2_K_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_K_R3.all.R1.fastq.gz ../../raw_data/LT_P2_K_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_P2_T

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_T_R1.all.R1.fastq.gz ../../raw_data/LT_P2_T_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_T_R2.all.R1.fastq.gz ../../raw_data/LT_P2_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_T_R3.all.R1.fastq.gz ../../raw_data/LT_P2_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_YFB_K

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_K_R1.all.R1.fastq.gz ../../raw_data/LT_YFB_K_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_K_R2.all.R1.fastq.gz ../../raw_data/LT_YFB_K_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_K_R3.all.R1.fastq.gz ../../raw_data/LT_YFB_K_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_YFB_T

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_T_R1.all.R1.fastq.gz ../../raw_data/LT_YFB_T_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_T_R2.all.R1.fastq.gz ../../raw_data/LT_YFB_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_T_R3.all.R1.fastq.gz ../../raw_data/LT_YFB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_FB_K

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_K_R1.all.R1.fastq.gz ../../raw_data/LT_FB_K_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_K_R3.all.R1.fastq.gz ../../raw_data/LT_FB_K_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/TL_FB_K_R4_R1_001.fastq.gz ../../raw_data/TL_FB_K_R4_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_STARout/Aligned.sortedByCoord.out.bam

# LT_FB_L

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_L_R1.R1.fastq.gz ../../raw_data/LT_FB_L_R1.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_L_R2.all.R1.fastq.gz ../../raw_data/LT_FB_L_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/TL_FB_L_R4_R1_001.fastq.gz ../../raw_data/TL_FB_L_R4_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_STARout/Aligned.sortedByCoord.out.bam

# LT_FB_T

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_T_R2.all.R1.fastq.gz ../../raw_data/LT_FB_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_T_R3.all.R1.fastq.gz ../../raw_data/LT_FB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/TL_FB_T_R4_R1_001.fastq.gz ../../raw_data/TL_FB_T_R4_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 -L VM,P1,P2_K,P2_T,YFB_K,YFB_T,FB_K,FB_L,FB_T --upper-quartile-norm ./ltigrinus_genome/ltigrinus_AS_annotation.gtf \
./ltigrinus_DIR_expressionanalysis/LT_VM_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_VM_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_VM_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_P1_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P1_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P1_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_T_R3_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_CUFFquant/abundances.cxb
```
## ORF prediction
### Input:
- __ltigrinus_AS_annotation.gtf__
- __Sisbr1_AssemblyScaffolds.fasta__
### Output:
- __ltigrinus_transcripts.fasta__ : A TransDecoder outputja ami tartalmazza az összes transzkript szekvenciáját fasta formátumban
- __longest_orfs.cds__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF cds régióját
- __longest_orfs.pep__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF proteinjeit

Prediktáljuk az ORF régiókat TransDecoder segítségével.
```
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ltigrinus_AS_annotation.gtf Sisbr1_AssemblyScaffolds.fasta > ltigrinus_transcripts.fasta
TransDecoder.LongOrfs -m 20 -S -t ltigrinus_transcripts.fasta
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
INTERPROSCAN_ltigrinus.R
```
Majd ezzel megszűrjük a fehérjeszekvenciákat tartalmazó fájlt, végül átalakítjuk úgy a FASTA-fájlt, hogy az InterProScan számára megfelelő legyen. (Elég körülményesen oldottam ezt meg de működik, nyilván mostmár máshogy írnám még ezt a részt egy R scriptbe)
```
# átalakítás
\)\n --> \)\t
# szűrés
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> ltigrinus_proteins_all.fasta; done <filter
# visszalakítás
\)\t --> \)\n
```
Töröljük a fasta fájlból a *-okat mivel az InterProScan nem tudja értelmezni
```
perl -pi -e 's/\*//g' ltigrinus_proteins_all.fasta
```
Egy kis regexp átalakítás az áttekinthetőségért
```
perl -pi -e 's/::g\..*$//g' ltigrinus_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' ltigrinus_proteins_all.fasta
```
## InterProScan analysis
### Input:
- __ltigrinus_proteins_all.fasta__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinek
### Output:
- __ltigrinus_proteins_all.fasta.tsv__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinhez tartozó InterProScan domainek

Lefuttatjuk at InterProScan-t.
```
interproscan.sh -i ltigrinus_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```