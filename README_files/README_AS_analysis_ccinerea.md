# Coprinopsis cinerea (AmutBmut)
## RNA-seq files
## Preparations
A bemeneti annotációs fájlokat át kell alakítani, hogy megfelelőek legyenek a Cufflinks számára illetve elő kell állítani azokat az annotációs és FASTA fájlokat amik az RRPM analízishez szükségesek. 
### Input:
- __Copci_AmutBmut1_GeneCatalog_genes_20130522.gff__ : Forrás a JGI
- __Copci_AmutBmut1_AssemblyScaffolds.fasta__ : Forrás a JGI
### Output:
- __ccinerea_onlygene.gtf__ : Csak a géneket tartalmazó annotációs fájl
- __ccinerea_onlyexon.gtf__ : Csak az exonokat tartalmazó annotációs fájl
- __ccinerea_fixed.gtf__ : Pozíciófixált annotációs fájl az RRPM számára
- __ccinerea_genes.fasta__ : Géneket tartalmazó fasta fájl

GTF fájlok elkészítése az RRPM számára
```
PREPARATIONS_ccinerea.R
```
GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
```
gtf2bed < ccinerea_onlygene.gtf > ccinerea_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl elkészítése
```
bedtools getfasta -name -fo ccinerea_genes.fasta -fi Copci_AmutBmut1_AssemblyScaffolds.fasta -bed ccinerea_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl headerjének a trimmelése
```
perl -pi -e 's/::.*//g' ccinerea_genes.fasta
```
## RRPM: STAR Alignment and Cufflinks Assembly
Megnézzük milyen hosszúak a gének és az intronok a STAR és a Cufflinks beállításához.

Elkészítünk egy fájlt ami tartalmazza az összes szükséges scriptet ami az RRPM futtatásához szükséges (STAR, Cufflinks)
### Input:
- __ccinerea_genes.fasta__
- __ccinerea_fixed.gtf__
### Output:
- __gene_length__
- A __Cufflinks__ __outputja__ ami tartalmazza az RRPM analízis fájljait
```
bioawk -c fastx '{ print $name, length($seq) }' < ccinerea_genes.fasta > gene_length
INTRON_LENGTH_ccinerea.R
```
A maximális transzkriptméret: 19640 --> max-bundle-length marad 250000

Intron min: 11 --> --alignIntronMin marad 3 

Intron max: 8003 --> --alignIntronMax marad 30000

Lefuttatjuk az RRPM analízist.
```
# STAR illesztéshez szükséges indexelés

mkdir ccinerea_AmutBmut_DIR_RRPM
cp -r ./ccinerea_AmutBmut_genome ./ccinerea_AmutBmut_DIR_RRPM/ccinerea_AmutBmut_RRPM_STARindex
cd ./ccinerea_AmutBmut_DIR_RRPM/ccinerea_AmutBmut_RRPM_STARindex
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir ../ccinerea_AmutBmut_RRPM_STARindex --sjdbGTFfile ../ccinerea_AmutBmut_RRPM_STARindex/ccinerea_AmutBmut_fixed.gtf --genomeFastaFiles ../ccinerea_AmutBmut_RRPM_STARindex/ccinerea_AmutBmut_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./ccinerea_AmutBmut_DIR_RRPM/ccinerea_AmutBmut_RRPM_STARout
cd ./ccinerea_AmutBmut_DIR_RRPM/ccinerea_AmutBmut_RRPM_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_RRPM_STARindex --readFilesIn ../../raw_data/CC_H_R2.all.R1.fastq.gz,../../raw_data/CC_H_R3.all.R1.fastq.gz,../../raw_data/CC_H_R4.all.R1.fastq.gz,../../raw_data/CC_VM_R1.all.R1.fastq.gz,../../raw_data/CC_VM_R2.all.R1.fastq.gz,../../raw_data/CC_VM_R3.all.R1.fastq.gz,../../raw_data/CC_P1_R2.all.R1.fastq.gz,../../raw_data/CC_P1_R3.all.R1.fastq.gz,../../raw_data/CC_P1_R4.all.R1.fastq.gz,../../raw_data/CC_P2_R2.all.R1.fastq.gz,../../raw_data/CC_P2_R3.all.R1.fastq.gz,../../raw_data/CC_P2_R4.all.R1.fastq.gz,../../raw_data/CC_YFB_K_R1.all.R1.fastq.gz,../../raw_data/CC_YFB_K_R2.all.R1.fastq.gz,../../raw_data/CC_YFB_K_R3.all.R1.fastq.gz,../../raw_data/CC_YFB_L_R1.all.R1.fastq.gz,../../raw_data/CC_YFB_L_R2.all.R1.fastq.gz,../../raw_data/CC_YFB_L_R3.all.R1.fastq.gz,../../raw_data/CC_YFB_T_R1.all.R1.fastq.gz,../../raw_data/CC_YFB_T_R2.all.R1.fastq.gz,../../raw_data/CC_YFB_T_R3.all.R1.fastq.gz,../../raw_data/CC_FB_KL_R1.all.R1.fastq.gz,../../raw_data/CC_FB_KL_R2.all.R1.fastq.gz,../../raw_data/CC_FB_KL_R3.all.R1.fastq.gz,../../raw_data/CC_FB_T_R1.all.R1.fastq.gz,../../raw_data/CC_FB_T_R2.all.R1.fastq.gz,../../raw_data/CC_FB_T_R3.all.R1.fastq.gz ../../raw_data/CC_H_R2.all.R2.fastq.gz,../../raw_data/CC_H_R3.all.R2.fastq.gz,../../raw_data/CC_H_R4.all.R2.fastq.gz,../../raw_data/CC_VM_R1.all.R2.fastq.gz,../../raw_data/CC_VM_R2.all.R2.fastq.gz,../../raw_data/CC_VM_R3.all.R2.fastq.gz,../../raw_data/CC_P1_R2.all.R2.fastq.gz,../../raw_data/CC_P1_R3.all.R2.fastq.gz,../../raw_data/CC_P1_R4.all.R2.fastq.gz,../../raw_data/CC_P2_R2.all.R2.fastq.gz,../../raw_data/CC_P2_R3.all.R2.fastq.gz,../../raw_data/CC_P2_R4.all.R2.fastq.gz,../../raw_data/CC_YFB_K_R1.all.R2.fastq.gz,../../raw_data/CC_YFB_K_R2.all.R2.fastq.gz,../../raw_data/CC_YFB_K_R3.all.R2.fastq.gz,../../raw_data/CC_YFB_L_R1.all.R2.fastq.gz,../../raw_data/CC_YFB_L_R2.all.R2.fastq.gz,../../raw_data/CC_YFB_L_R3.all.R2.fastq.gz,../../raw_data/CC_YFB_T_R1.all.R2.fastq.gz,../../raw_data/CC_YFB_T_R2.all.R2.fastq.gz,../../raw_data/CC_YFB_T_R3.all.R2.fastq.gz,../../raw_data/CC_FB_KL_R1.all.R2.fastq.gz,../../raw_data/CC_FB_KL_R2.all.R2.fastq.gz,../../raw_data/CC_FB_KL_R3.all.R2.fastq.gz,../../raw_data/CC_FB_T_R1.all.R2.fastq.gz,../../raw_data/CC_FB_T_R2.all.R2.fastq.gz,../../raw_data/CC_FB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 48 -g ./ccinerea_AmutBmut_DIR_RRPM/ccinerea_AmutBmut_RRPM_STARindex/ccinerea_AmutBmut_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./ccinerea_AmutBmut_DIR_RRPM/ccinerea_AmutBmut_RRPM_CUFFout ./ccinerea_AmutBmut_DIR_RRPM/ccinerea_AmutBmut_RRPM_STARout/Aligned.sortedByCoord.out.bam
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
- __ccinerea_onlygene.gtf__
- __transcripts.gtf__ : RRPM Cufflinks assembly output file
### Output:
- __ccinerea_RRPM_transcripts.gtf__ : Az újonnan felfedezett transzkripteket tartalmazó GTF
- __ccinerea_AS_annotation.gtf__ : Az eredeti annotációval történő összemergeelés, hogy az esetlegesen nem detektált gének is bennelegyenek az annotációs fájlban.
- __ccinerea_stats.log__ : Szűrési statisztikákat tartalmazó logfájl
```
FILTERING_ccinerea.R
```
## Fusion correction
### Input:
- __isoforms.fpkm_tracking__ : Az RRPM Cufflinks output izoforms FPKM értékeket tartalmazó fájlja
- __ccinerea_AS_annotation.gtf__
### Output:
- Összesen 2 fúziós gén detektálva : __371210__

Megvizsgáljuk, hogy melyek azok a gének amik fúziósak voltak az eredeti annotációba de már két külön gént alkotnak, majd ezeket a géneket manuális IGV megtekintés után átnevezzük v1, v2 satöbbire.
```
FUSION_FILTER_ccinerea.R
```
## Expression: STAR Alignment, Cuffquant Assembly and Cuffdiff analysis
### Output:
- A __CuffDiff__ __output__ fájljai az expreszziós analízisből

Futtatunk egy teljes STAR illesztést majd egy Cuffquantot, majd kinyerjük az expressziós értékeket a CuffDiff segítségével.
```
# STAR INDEXING

mkdir DIR_ccinerea_AmutBmut_expressionanalysis
cp -r ./ccinerea_AmutBmut_genome ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --genomeFastaFiles ./Copci_AmutBmut1_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# CC_VM

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_VM_R1.all.R1.fastq.gz ../../raw_data/CC_VM_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_VM_R2.all.R1.fastq.gz ../../raw_data/CC_VM_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_VM_R3.all.R1.fastq.gz ../../raw_data/CC_VM_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_H

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_H_R2.all.R1.fastq.gz ../../raw_data/CC_H_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_H_R3.all.R1.fastq.gz ../../raw_data/CC_H_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_H_R4.all.R1.fastq.gz ../../raw_data/CC_H_R4.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_STARout/Aligned.sortedByCoord.out.bam

# CC_P1

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P1_R2.all.R1.fastq.gz ../../raw_data/CC_P1_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P1_R3.all.R1.fastq.gz ../../raw_data/CC_P1_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P1_R4.all.R1.fastq.gz ../../raw_data/CC_P1_R4.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_STARout/Aligned.sortedByCoord.out.bam

# CC_P2

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P2_R2.all.R1.fastq.gz ../../raw_data/CC_P2_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P2_R3.all.R1.fastq.gz ../../raw_data/CC_P2_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P2_R4.all.R1.fastq.gz ../../raw_data/CC_P2_R4.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_STARout/Aligned.sortedByCoord.out.bam

# CC_YFB_K

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_K_R1.all.R1.fastq.gz ../../raw_data/CC_YFB_K_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_K_R2.all.R1.fastq.gz ../../raw_data/CC_YFB_K_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_K_R3.all.R1.fastq.gz ../../raw_data/CC_YFB_K_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_YFB_L

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_L_R1.all.R1.fastq.gz ../../raw_data/CC_YFB_L_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_L_R2.all.R1.fastq.gz ../../raw_data/CC_YFB_L_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_L_R3.all.R1.fastq.gz ../../raw_data/CC_YFB_L_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_YFB_T

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_T_R1.all.R1.fastq.gz ../../raw_data/CC_YFB_T_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_T_R2.all.R1.fastq.gz ../../raw_data/CC_YFB_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_T_R3.all.R1.fastq.gz ../../raw_data/CC_YFB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_FB_KL

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_KL_R1.all.R1.fastq.gz ../../raw_data/CC_FB_KL_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_KL_R2.all.R1.fastq.gz ../../raw_data/CC_FB_KL_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_KL_R3.all.R1.fastq.gz ../../raw_data/CC_FB_KL_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_FB_T

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_T_R1.all.R1.fastq.gz ../../raw_data/CC_FB_T_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_T_R2.all.R1.fastq.gz ../../raw_data/CC_FB_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_T_R3.all.R1.fastq.gz ../../raw_data/CC_FB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 -L VM,H,P1,P2,YFB_K,YFB_L,YFB_T,FB_KL,FB_T --upper-quartile-norm ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_CUFFquant/abundances.cxb
```
## ORF prediction
### Input:
- __ccinerea_AS_annotation.gtf__
- __Copci_AmutBmut1_AssemblyScaffolds.fasta__
### Output:
- __ccinerea_transcripts.fasta__ : A TransDecoder outputja ami tartalmazza az összes transzkript szekvenciáját fasta formátumban
- __longest_orfs.cds__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF cds régióját
- __longest_orfs.pep__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF proteinjeit

Prediktáljuk az ORF régiókat TransDecoder segítségével.
```
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ccinerea_AS_annotation.gtf Copci_AmutBmut1_AssemblyScaffolds.fasta > ccinerea_transcripts.fasta
TransDecoder.LongOrfs -m 20 -S -t ccinerea_transcripts.fasta
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
INTERPROSCAN_ccinerea.R
```
Majd ezzel megszűrjük a fehérjeszekvenciákat tartalmazó fájlt, végül átalakítjuk úgy a FASTA-fájlt, hogy az InterProScan számára megfelelő legyen. (Elég körülményesen oldottam ezt meg de működik, nyilván mostmár máshogy írnám még ezt a részt egy R scriptbe)
```
# átalakítás
\)\n --> \)\t
# szűrés
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> ccinerea_proteins_all.fasta; done <filter
# visszalakítás
\)\t --> \)\n
```
Töröljük a fasta fájlból a *-okat mivel az InterProScan nem tudja értelmezni
```
perl -pi -e 's/\*//g' ccinerea_proteins_all.fasta
```
Egy kis regexp átalakítás az áttekinthetőségért
```
perl -pi -e 's/::g\..*$//g' ccinerea_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' ccinerea_proteins_all.fasta
```
## InterProScan analysis
### Input:
- __ccinerea_proteins_all.fasta__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinek
### Output:
- __ccinerea_proteins_all.fasta.tsv__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinhez tartozó InterProScan domainek

Lefuttatjuk at InterProScan-t.
```
interproscan.sh -i ccinerea_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```