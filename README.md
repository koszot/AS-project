# Armillaria ostoyae
## RNA-seq files
## Preparations
FONTOS: A scriptek nem lettek frissítve csak az output fájlok elnevezései mivel ez volt az egyik első leelemzett genome, és az elnevezési metodika a későbbi fajoknál jött be, ezért a scriptek újrafuttatás esetén átnézésre/átalakításra szorulnak, hogy az újabb fajok scriptjeihez hasonlóak legyenek.

A bemeneti annotációs fájlokat át kell alakítani, hogy megfelelőek legyenek a Cufflinks számára illetve elő kell állítani azokat az annotációs és FASTA fájlokat amik az RRPM analízishez szükségesek. 
### Input:
- __p3_i2_t47428_Arm_ostoy_v2.gff3__ : saját annotációs fájl
- __p3_i2_t47428_Arm_ostoy_v2.scaf__ : saját scaffoldokat tartalmazó FASTA
### Output:
- __aostoyae.gtf__ : A GTF-é alakított GFF3 fájl
- __aostoyae_onlygene.gtf__ : Csak a géneket tartalmazó annotációs fájl
- __aostoyae_onlyexon.gtf__ : Csak az exonokat tartalmazó annotációs fájl
- __aostoyae_fixed.gtf__ : Pozíciófixált annotációs fájl az RRPM számára
- __aostoyae_genes.fasta__ : Géneket tartalmazó fasta fájl

A GFF3 fájlt átalakítjuk GTF formátumra a további analízishez
```
gffread p3_i2_t47428_Arm_ostoy_v2.gff3 -T -o aostoyae.gtf
```
GTF fájlok elkészítése az RRPM számára
```
PREPARATIONS_aostoyae.R
```
GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
```
gtf2bed < aostoyae_onlygene.gtf > aostoyae_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl elkészítése
```
bedtools getfasta -name -fo aostoyae_genes.fasta -fi p3_i2_t47428_Arm_ostoy_v2.scaf -bed aostoyae_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl headerjének a trimmelése
```
perl -pi -e 's/::.*//g' aostoyae_genes.fasta
```
## RRPM: STAR Alignment and Cufflinks Assembly
Megnézzük milyen hosszúak a gének és az intronok a STAR és a Cufflinks beállításához.

Elkészítünk egy fájlt ami tartalmazza az összes szükséges scriptet ami az RRPM futtatásához szükséges (STAR, Cufflinks)
### Input:
- __aostoyae_genes.fasta__
- __aostoyae_fixed.gtf__
### Output:
- __gene_length__
- A __Cufflinks__ __outputja__ ami tartalmazza az RRPM analízis fájljait
```
bioawk -c fastx '{ print $name, length($seq) }' < aostoyae_genes.fasta > gene_length
INTRON_LENGTH_aostoyae.R
```
A maximális transzkriptméret: 16175 --> max-bundle-length marad 250000

Intron min: 21 --> --alignIntronMin marad 3 

Intron max: 6767 --> --alignIntronMax marad 30000

Lefuttatjuk az RRPM analízist.
```
mkdir DIR_RRPM
cp -r ./genome ./DIR_RRPM/RRPM_STARindex
cd ./DIR_RRPM/RRPM_STARindex
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ../RRPM_STARindex --sjdbGTFfile ../RRPM_STARindex/original_fixed.gtf --genomeFastaFiles ../RRPM_STARindex/original_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

mkdir ./DIR_RRPM/RRPM_STARout
cd ./DIR_RRPM/RRPM_STARout
STAR --runThreadN 16 --genomeDir ../RRPM_STARindex --readFilesIn ../../raw_data/ARO_FB_C_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_C_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_C_R3.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R3.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R3.all.R1.fastq.gz,../../raw_data/ARO_P1_R1.all.R1.fastq.gz,../../raw_data/ARO_P1_R2.all.R1.fastq.gz,../../raw_data/ARO_P1_R3.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R1.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R2.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R3.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R1.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R2.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R3.all.R1.fastq.gz,../../raw_data/ARO_RMA_R1.all.R1.fastq.gz,../../raw_data/ARO_RMA_R2.all.R1.fastq.gz,../../raw_data/ARO_RMA_R3.all.R1.fastq.gz,../../raw_data/ARO_VM_R1.all.R1.fastq.gz,../../raw_data/ARO_VM_R2.all.R1.fastq.gz,../../raw_data/ARO_VM_R3.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R1.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R2.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R3.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R1.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R2.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_C_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_C_R3.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R3.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R3.all.R2.fastq.gz,../../raw_data/ARO_P1_R1.all.R2.fastq.gz,../../raw_data/ARO_P1_R2.all.R2.fastq.gz,../../raw_data/ARO_P1_R3.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R1.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R2.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R3.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R1.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R2.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R3.all.R2.fastq.gz,../../raw_data/ARO_RMA_R1.all.R2.fastq.gz,../../raw_data/ARO_RMA_R2.all.R2.fastq.gz,../../raw_data/ARO_RMA_R3.all.R2.fastq.gz,../../raw_data/ARO_VM_R1.all.R2.fastq.gz,../../raw_data/ARO_VM_R2.all.R2.fastq.gz,../../raw_data/ARO_VM_R3.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R1.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R2.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R3.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R1.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R2.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cufflinks -p 16 -g ./DIR_RRPM/RRPM_STARindex/original_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./DIR_RRPM/RRPM_CUFFout ./DIR_RRPM/RRPM_STARout/Aligned.sortedByCoord.out.bam
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
- __aostoyae_onlygene.gtf__
- __transcripts.gtf__ : RRPM Cufflinks assembly output file
### Output:
- __aostoyae_RRPM_transcripts.gtf__ : Az újonnan felfedezett transzkripteket tartalmazó GTF
- __aostoyae_AS_annotation.gtf__ : Az eredeti annotációval történő összemergeelés, hogy az esetlegesen nem detektált gének is bennelegyenek az annotációs fájlban.
- __aostoyae_stats.log__ : Szűrési statisztikákat tartalmazó logfájl
```
FILTERING_aostoyae.R
```
## Fusion correction
### Input:
- __isoforms.fpkm_tracking__ : Az RRPM Cufflinks output izoforms FPKM értékeket tartalmazó fájlja
- __aostoyae_AS_annotation.gtf__
### Output:
- Összesen 2 fúziós gén detektálva : __AROS_19248__, __AROS_20796__

Megvizsgáljuk, hogy melyek azok a gének amik fúziósak voltak az eredeti annotációba de már két külön gént alkotnak, majd ezeket a géneket manuális IGV megtekintés után átnevezzük v1, v2 satöbbire.
```
FUSION_FILTER_aostoyae.R
```
## Expression: STAR Alignment, Cuffquant Assembly and Cuffdiff analysis
### Output:
- A __CuffDiff__ __output__ fájljai az expreszziós analízisből

Futtatunk egy teljes STAR illesztést majd egy Cuffquantot, majd kinyerjük az expressziós értékeket a CuffDiff segítségével.
```
# STAR INDEXING

cp -r ./genome ./DIR_expressionanalysis/STARindex
cd ./DIR_expressionanalysis/STARindex
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../STARindex/RRPM_transcripts.fixed.gtf --genomeFastaFiles ./p3_i2_t47428_Arm_ostoy_v2.scaf --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# ARO_RMA

mkdir ./DIR_expressionanalysis/ARO_RMA_R1_STARout
cd ./DIR_expressionanalysis/ARO_RMA_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_RMA_R1.all.R1.fastq.gz ../../raw_data/ARO_RMA_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_RMA_R2_STARout
cd ./DIR_expressionanalysis/ARO_RMA_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_RMA_R2.all.R1.fastq.gz ../../raw_data/ARO_RMA_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_RMA_R3_STARout
cd ./DIR_expressionanalysis/ARO_RMA_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_RMA_R3.all.R1.fastq.gz ../../raw_data/ARO_RMA_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_VM

mkdir ./DIR_expressionanalysis/ARO_VM_R1_STARout
cd ./DIR_expressionanalysis/ARO_VM_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_VM_R1.all.R1.fastq.gz ../../raw_data/ARO_VM_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_VM_R2_STARout
cd ./DIR_expressionanalysis/ARO_VM_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_VM_R2.all.R1.fastq.gz ../../raw_data/ARO_VM_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_VM_R3_STARout
cd ./DIR_expressionanalysis/ARO_VM_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_VM_R3.all.R1.fastq.gz ../../raw_data/ARO_VM_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_P1

mkdir ./DIR_expressionanalysis/ARO_P1_R1_STARout
cd ./DIR_expressionanalysis/ARO_P1_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P1_R1.all.R1.fastq.gz ../../raw_data/ARO_P1_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P1_R2_STARout
cd ./DIR_expressionanalysis/ARO_P1_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P1_R2.all.R1.fastq.gz ../../raw_data/ARO_P1_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P1_R3_STARout
cd ./DIR_expressionanalysis/ARO_P1_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P1_R3.all.R1.fastq.gz ../../raw_data/ARO_P1_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_P2_C

mkdir ./DIR_expressionanalysis/ARO_P2_C_R1_STARout
cd ./DIR_expressionanalysis/ARO_P2_C_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_C_R1.all.R1.fastq.gz ../../raw_data/ARO_P2_C_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P2_C_R2_STARout
cd ./DIR_expressionanalysis/ARO_P2_C_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_C_R2.all.R1.fastq.gz ../../raw_data/ARO_P2_C_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P2_C_R3_STARout
cd ./DIR_expressionanalysis/ARO_P2_C_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_C_R3.all.R1.fastq.gz ../../raw_data/ARO_P2_C_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_P2_S

mkdir ./DIR_expressionanalysis/ARO_P2_S_R1_STARout
cd ./DIR_expressionanalysis/ARO_P2_S_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_S_R1.all.R1.fastq.gz ../../raw_data/ARO_P2_S_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P2_S_R2_STARout
cd ./DIR_expressionanalysis/ARO_P2_S_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_S_R2.all.R1.fastq.gz ../../raw_data/ARO_P2_S_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P2_S_R3_STARout
cd ./DIR_expressionanalysis/ARO_P2_S_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_S_R3.all.R1.fastq.gz ../../raw_data/ARO_P2_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_YFB_C

mkdir ./DIR_expressionanalysis/ARO_YFB_C_R1_STARout
cd ./DIR_expressionanalysis/ARO_YFB_C_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_C_R1.all.R1.fastq.gz ../../raw_data/ARO_YFB_C_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_YFB_C_R2_STARout
cd ./DIR_expressionanalysis/ARO_YFB_C_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_C_R2.all.R1.fastq.gz ../../raw_data/ARO_YFB_C_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_YFB_C_R3_STARout
cd ./DIR_expressionanalysis/ARO_YFB_C_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_C_R3.all.R1.fastq.gz ../../raw_data/ARO_YFB_C_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_YFB_S

mkdir ./DIR_expressionanalysis/ARO_YFB_S_R1_STARout
cd ./DIR_expressionanalysis/ARO_YFB_S_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_S_R1.all.R1.fastq.gz ../../raw_data/ARO_YFB_S_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_YFB_S_R2_STARout
cd ./DIR_expressionanalysis/ARO_YFB_S_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_S_R2.all.R1.fastq.gz ../../raw_data/ARO_YFB_S_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_YFB_S_R3_STARout
cd ./DIR_expressionanalysis/ARO_YFB_S_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_S_R3.all.R1.fastq.gz ../../raw_data/ARO_YFB_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_FB_C

mkdir ./DIR_expressionanalysis/ARO_FB_C_R1_STARout
cd ./DIR_expressionanalysis/ARO_FB_C_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_C_R1.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_C_R2_STARout
cd ./DIR_expressionanalysis/ARO_FB_C_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_C_R2.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_C_R3_STARout
cd ./DIR_expressionanalysis/ARO_FB_C_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_C_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_FB_L

mkdir ./DIR_expressionanalysis/ARO_FB_L_R1_STARout
cd ./DIR_expressionanalysis/ARO_FB_L_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_L_R1.all.R1.fastq.gz ../../raw_data/ARO_FB_L_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_L_R2_STARout
cd ./DIR_expressionanalysis/ARO_FB_L_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_L_R2.all.R1.fastq.gz ../../raw_data/ARO_FB_L_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_L_R3_STARout
cd ./DIR_expressionanalysis/ARO_FB_L_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_L_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_L_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_FB_S

mkdir ./DIR_expressionanalysis/ARO_FB_S_R1_STARout
cd ./DIR_expressionanalysis/ARO_FB_S_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_S_R1.all.R1.fastq.gz ../../raw_data/ARO_FB_S_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_S_R2_STARout
cd ./DIR_expressionanalysis/ARO_FB_S_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_S_R2.all.R1.fastq.gz ../../raw_data/ARO_FB_S_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_S_R3_STARout
cd ./DIR_expressionanalysis/ARO_FB_S_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_S_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# CUFFQUANT

# ARO_RMA

cuffquant -o ./DIR_expressionanalysis/ARO_RMA_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_RMA_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_RMA_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_RMA_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_RMA_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_RMA_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_VM

cuffquant -o ./DIR_expressionanalysis/ARO_VM_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_VM_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_VM_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_P1

cuffquant -o ./DIR_expressionanalysis/ARO_P1_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P1_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P1_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P1_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P1_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P1_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_P2_C

cuffquant -o ./DIR_expressionanalysis/ARO_P2_C_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_C_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P2_C_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_C_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P2_C_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_C_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_P2_S

cuffquant -o ./DIR_expressionanalysis/ARO_P2_S_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_S_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P2_S_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_S_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P2_S_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_S_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_YFB_C

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_C_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_C_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_C_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_C_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_C_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_C_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_YFB_S

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_S_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_S_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_S_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_S_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_S_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_S_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_FB_C

cuffquant -o ./DIR_expressionanalysis/ARO_FB_C_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_C_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_C_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_C_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_C_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_C_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_FB_L

cuffquant -o ./DIR_expressionanalysis/ARO_FB_L_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_L_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_L_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_L_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_L_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_L_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_FB_S

cuffquant -o ./DIR_expressionanalysis/ARO_FB_S_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_S_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_S_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_S_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_S_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_S_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -o CUFFdiff -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 96 -L RMA,VM,P1,P2_C,P2_S,YFB_C,YFB_S,FB_C,FB_L,FB_S --upper-quartile-norm ./DIR_RRPM/RRPM/RRPM_transcripts.fixed.gtf \
./DIR_expressionanalysis/ARO_RMA_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_RMA_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_RMA_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_VM_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_VM_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_VM_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_P1_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P1_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P1_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_P2_C_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P2_C_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P2_C_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_P2_S_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P2_S_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P2_S_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_YFB_C_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_YFB_C_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_YFB_C_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_YFB_S_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_YFB_S_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_YFB_S_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_FB_C_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_C_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_C_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_FB_L_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_L_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_L_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_FB_S_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_S_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_S_R3_CUFFquant/abundances.cxb
```
## ORF prediction
### Input:
- __aostoyae_AS_annotation.gtf__
- __p3_i2_t47428_Arm_ostoy_v2.scaf__
### Output:
- __aostoyae_transcripts.fasta__ : A TransDecoder outputja ami tartalmazza az összes transzkript szekvenciáját fasta formátumban
- __longest_orfs.cds__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF cds régióját
- __longest_orfs.pep__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF proteinjeit

Prediktáljuk az ORF régiókat TransDecoder segítségével.
```
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl aostoyae_AS_annotation.gtf p3_i2_t47428_Arm_ostoy_v2.scaf > aostoyae_transcripts.fasta
TransDecoder.LongOrfs -m 20 -S -t aostoyae_transcripts.fasta
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
INTERPROSCAN_aostoyae.R
```
Majd ezzel megszűrjük a fehérjeszekvenciákat tartalmazó fájlt, végül átalakítjuk úgy a FASTA-fájlt, hogy az InterProScan számára megfelelő legyen. (Elég körülményesen oldottam ezt meg de működik, nyilván mostmár máshogy írnám még ezt a részt egy R scriptbe)
```
# átalakítás
\)\n --> \)\t
# szűrés
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> aostoyae_proteins_all.fasta; done <filter
# visszalakítás
\)\t --> \)\n
```
Töröljük a fasta fájlból a *-okat mivel az InterProScan nem tudja értelmezni
```
perl -pi -e 's/\*//g' aostoyae_proteins_all.fasta
```
Egy kis regexp átalakítás az áttekinthetőségért
```
perl -pi -e 's/::g\..*$//g' aostoyae_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' aostoyae_proteins_all.fasta
```

## InterProScan analysis
### Input:
- __aostoyae_proteins_all.fasta__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinek
### Output:
- __aostoyae_proteins_all.fasta.tsv__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinhez tartozó InterProScan domainek

Lefuttatjuk at InterProScan-t.
```
interproscan.sh -i aostoyae_proteins_all.fasta -f tsv --iprlookup --goterms
```
## Enrichment
### Input:
- __p3_i2_t47428_Arm_ostoy_v2.prot__
- __isoforms.fpkm_tracking__ : A CuffDiff analízis output fájlja ami az izoformák expressziós értékeit tartalmazza
- __genes.fpkm_tracking__ : A CuffDiff analízis output fájlja ami a gének expressziós értékeit tartalmazza
### Output:
- __p3_i2_t47428_Arm_ostoy_v2.prot.tsv__ : Az InterProScan domaineket tartalmazó output fájlja
- __aostoyae_GO_dict.tsv__ : A géneket és a GO-kat tartalmazó táblázat

Lefuttatjuk az InterProScan-t a GeneCatalog génjeire.
```
interproscan.sh -i p3_i2_t47428_Arm_ostoy_v2.prot -f tsv --iprlookup --goterms
```
Készítünk egy táblázatot ami tartalmazza a génekhez tartozó GO azonosítókat az enrichment számára.
```
ENRICHMENT_preparations_aostoyae.R
```
FONTOS: A következő két scriptet újra kell futtatni mert megváltozott a definíciója a DEVREG-nek és az AS-nek is.
```
ENRICHMENT_groups_aostoyae.R
ENRICHMENT_fisher_classic_aostoyae.R
```
