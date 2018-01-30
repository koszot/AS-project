# Rickenella mellea
## RNA-seq files
Marie Tollot, Daniela Assmann, Christian Becker, Janine Altmüller, Julien Y. Dutheil, Carl-Eric Wegner, Regine Kahmann 

The WOPR Protein Ros1 Is a Master Regulator of Sporogenesis and Late Effector Gene Expression in the Maize Pathogen Ustilago maydis

RNA-seq DATA:
- R1 : GSM1977379: WT_I; Ustilago maydis; RNA-Seq (SRR3038903)	->	UM_R1.R1.fastq	UM_R1.R2.fastq
- R2 : GSM1977380: WT_II; Ustilago maydis; RNA-Seq (SRR3038904)	->	UM_R2.R1.fastq	UM_R2.R2.fastq
- R3 : GSM1977380: WT_III; Ustilago maydis; RNA-Seq (SRR3038905) ->	UM_R3.R1.fastq	UM_R3.R2.fastq
## Preparations
A bemeneti annotációs fájlokat át kell alakítani, hogy megfelelőek legyenek a Cufflinks számára illetve elő kell állítani azokat az annotációs és FASTA fájlokat amik az RRPM analízishez szükségesek. 
### Input:
- __p3_t237631_Ust_maydi_v2GB.gff3__ : ftp://ftpmips.gsf.de/fungi/Ustilaginaceae/Ustilago_maydis_521/
- __p3_t237631_Ust_maydi_v2GB.scaf__ : ftp://ftpmips.gsf.de/fungi/Ustilaginaceae/Ustilago_maydis_521/
### Output:
- __umaydis_onlygene.gtf__ : Csak a géneket tartalmazó annotációs fájl
- __umaydis_onlyexon.gtf__ : Csak az exonokat tartalmazó annotációs fájl
- __umaydis_fixed.gtf__ : Pozíciófixált annotációs fájl az RRPM számára
- __umaydis_genes.fasta__ : Géneket tartalmazó fasta fájl

A GFF3 fájl átalakítása GTF formátumra
```
gffread p3_t237631_Ust_maydi_v2GB.gff3 -T -o original.gtf
```
GTF fájlok elkészítése az RRPM számára
```
PREPARATIONS_umaydis.R
```
GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
```
gtf2bed < umaydis_onlygene.gtf > umaydis_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl elkészítése
```
bedtools getfasta -name -fo umaydis_genes.fasta -fi p3_t237631_Ust_maydi_v2GB.scaf -bed umaydis_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl headerjének a trimmelése
```
perl -pi -e 's/::.*//g' umaydis_genes.fasta
```
## RRPM: STAR Alignment and Cufflinks Assembly
Megnézzük milyen hosszúak a gének és az intronok a STAR és a Cufflinks beállításához.

Elkészítünk egy fájlt ami tartalmazza az összes szükséges scriptet ami az RRPM futtatásához szükséges (STAR, Cufflinks)
### Input:
- __umaydis_genes.fasta__
- __umaydis_fixed.gtf__
### Output:
- __gene_length__
- A __Cufflinks__ __outputja__ ami tartalmazza az RRPM analízis fájljait
```
bioawk -c fastx '{ print $name, length($seq) }' < scommune_genes.fasta > gene_length
INTRON_LENGTH_scommune.R
```
A maximális transzkriptméret: 16296 --> max-bundle-length marad 250000

Intron min: 15 --> --alignIntronMin marad 3 

Intron max: 2182 --> --alignIntronMax marad 30000

Lefuttatjuk az RRPM analízist.

Túl sok a short read miatti unmapped MEGOLDÁS -> https://github.com/alexdobin/STAR/issues/169
```
# STAR illesztéshez szükséges indexelés

mkdir umaydis_DIR_RRPM
cp -r ./umaydis_genome ./umaydis_DIR_RRPM/umaydis_RRPM_STARindex
cd ./umaydis_DIR_RRPM/umaydis_RRPM_STARindex
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ../umaydis_RRPM_STARindex --sjdbGTFfile ../umaydis_RRPM_STARindex/umaydis_fixed.gtf --genomeFastaFiles ../umaydis_RRPM_STARindex/umaydis_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./umaydis_DIR_RRPM/umaydis_RRPM_STARout
cd ./umaydis_DIR_RRPM/umaydis_RRPM_STARout
STAR --runThreadN 32 --genomeDir ../umaydis_RRPM_STARindex --readFilesIn ../../raw_data/UM_R1.R1.fastq.gz,../../raw_data/UM_R2.R1.fastq.gz,../../raw_data/UM_R3.R1.fastq.gz ../../raw_data/UM_R1.R2.fastq.gz,../../raw_data/UM_R2.R2.fastq.gz,../../raw_data/UM_R3.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2
cd ../..

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 32 -g ./umaydis_DIR_RRPM/umaydis_RRPM_STARindex/umaydis_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./umaydis_DIR_RRPM/umaydis_RRPM_CUFFout_v2 ./umaydis_DIR_RRPM/umaydis_RRPM_STARout/Aligned.sortedByCoord.out.bam
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
- __umaydis_onlygene.gtf__
- __transcripts.gtf__ : RRPM Cufflinks assembly output file
### Output:
- __umaydis_RRPM_transcripts.gtf__ : Az újonnan felfedezett transzkripteket tartalmazó GTF
- __umaydis_AS_annotation.gtf__ : Az eredeti annotációval történő összemergeelés, hogy az esetlegesen nem detektált gének is bennelegyenek az annotációs fájlban.
- __umaydis_stats.log__ : Szűrési statisztikákat tartalmazó logfájl
```
FILTERING_umaydis.R
```
## Fusion correction
### Input:
- __isoforms.fpkm_tracking__ : Az RRPM Cufflinks output izoforms FPKM értékeket tartalmazó fájlja
- __umaydis_AS_annotation.gtf__
### Output:
- Nulla fúziós gént detektált a script

Megvizsgáljuk, hogy melyek azok a gének amik fúziósak voltak az eredeti annotációba de már két külön gént alkotnak, majd ezeket a géneket manuális IGV megtekintés után átnevezzük v1, v2 satöbbire.
```
FUSION_FILTER_umaydis.R
```
## Expression: STAR Alignment, Cuffquant Assembly and Cuffdiff analysis
### Output:
- A __CuffDiff__ __output__ fájljai az expreszziós analízisből

Túl sok a short read miatti unmapped MEGOLDÁS -> https://github.com/alexdobin/STAR/issues/169

Futtatunk egy teljes STAR illesztést majd egy Cuffquantot, majd kinyerjük az expressziós értékeket a CuffDiff segítségével.
```
# STAR INDEXING

mkdir DIR_umaydis_expressionanalysis
cp -r ./umaydis_genome ./DIR_umaydis_expressionanalysis/umaydis_STARindex
cd ./DIR_umaydis_expressionanalysis/umaydis_STARindex
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../umaydis_STARindex/umaydis_AS_annotation.gtf --genomeFastaFiles ./p3_t237631_Ust_maydi_v2GB.scaf --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# UM

mkdir ./DIR_umaydis_expressionanalysis/UM_R1_STARout
cd ./DIR_umaydis_expressionanalysis/UM_R1_STARout
STAR --runThreadN 32 --genomeDir ../umaydis_STARindex --readFilesIn ../../raw_data/UM_R1.R1.fastq.gz ../../raw_data/UM_R1.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2
cd ../..

mkdir ./DIR_umaydis_expressionanalysis/UM_R2_STARout
cd ./DIR_umaydis_expressionanalysis/UM_R2_STARout
STAR --runThreadN 32 --genomeDir ../umaydis_STARindex --readFilesIn ../../raw_data/UM_R2.R1.fastq.gz ../../raw_data/UM_R2.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2
cd ../..

mkdir ./DIR_umaydis_expressionanalysis/UM_R3_STARout
cd ./DIR_umaydis_expressionanalysis/UM_R3_STARout
STAR --runThreadN 32 --genomeDir ../umaydis_STARindex --readFilesIn ../../raw_data/UM_R3.R1.fastq.gz ../../raw_data/UM_R3.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2
cd ../..

cuffquant -q -o ./DIR_umaydis_expressionanalysis/UM_R1_CUFFquant -b ./DIR_umaydis_expressionanalysis/umaydis_STARindex/p3_t237631_Ust_maydi_v2GB.scaf -p 32 ./DIR_umaydis_expressionanalysis/umaydis_STARindex/umaydis_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_umaydis_expressionanalysis/UM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_umaydis_expressionanalysis/UM_R2_CUFFquant -b ./DIR_umaydis_expressionanalysis/umaydis_STARindex/p3_t237631_Ust_maydi_v2GB.scaf -p 32 ./DIR_umaydis_expressionanalysis/umaydis_STARindex/umaydis_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_umaydis_expressionanalysis/UM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_umaydis_expressionanalysis/UM_R3_CUFFquant -b ./DIR_umaydis_expressionanalysis/umaydis_STARindex/p3_t237631_Ust_maydi_v2GB.scaf -p 32 ./DIR_umaydis_expressionanalysis/umaydis_STARindex/umaydis_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_umaydis_expressionanalysis/UM_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./DIR_umaydis_expressionanalysis/umaydis_STARindex/p3_t237631_Ust_maydi_v2GB.scaf -p 32 -L UM_1,UM_2 --upper-quartile-norm ./DIR_umaydis_expressionanalysis/umaydis_STARindex/umaydis_AS_annotation.gtf \
./DIR_umaydis_expressionanalysis/UM_R1_CUFFquant/abundances.cxb,./DIR_umaydis_expressionanalysis/UM_R2_CUFFquant/abundances.cxb,./DIR_umaydis_expressionanalysis/UM_R3_CUFFquant/abundances.cxb \
./DIR_umaydis_expressionanalysis/UM_R1_CUFFquant/abundances.cxb,./DIR_umaydis_expressionanalysis/UM_R2_CUFFquant/abundances.cxb,./DIR_umaydis_expressionanalysis/UM_R3_CUFFquant/abundances.cxb
```
## ORF prediction
### Input:
- __umaydis_AS_annotation.gtf__
- __p3_t237631_Ust_maydi_v2GB.scaf__
### Output:
- __umaydis_transcripts.fasta__ : A TransDecoder outputja ami tartalmazza az összes transzkript szekvenciáját fasta formátumban
- __longest_orfs.cds__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF cds régióját
- __longest_orfs.pep__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF proteinjeit

Prediktáljuk az ORF régiókat TransDecoder segítségével.
```
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl umaydis_AS_annotation.gtf p3_t237631_Ust_maydi_v2GB.scaf > umaydis_transcripts.fasta
TransDecoder.LongOrfs -m 20 -S -t umaydis_transcripts.fasta
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
INTERPROSCAN_umaydis.R
```
Majd ezzel megszűrjük a fehérjeszekvenciákat tartalmazó fájlt, végül átalakítjuk úgy a FASTA-fájlt, hogy az InterProScan számára megfelelő legyen. (Elég körülményesen oldottam ezt meg de működik, nyilván mostmár máshogy írnám még ezt a részt egy R scriptbe)
```
# átalakítás
\)\n --> \)\t
# szűrés
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> umaydis_proteins_all.fasta; done <filter
# visszalakítás
\)\t --> \)\n
```
Töröljük a fasta fájlból a *-okat mivel az InterProScan nem tudja értelmezni
```
perl -pi -e 's/\*//g' umaydis_proteins_all.fasta
```
Egy kis regexp átalakítás az áttekinthetőségért
```
perl -pi -e 's/::g\..*$//g' umaydis_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' umaydis_proteins_all.fasta
```
## InterProScan analysis
### Input:
- __umaydis_proteins_all.fasta__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinek
### Output:
- __umaydis_proteins_all.fasta.tsv__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinhez tartozó InterProScan domainek

Lefuttatjuk at InterProScan-t.
```
interproscan.sh -i umaydis_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```