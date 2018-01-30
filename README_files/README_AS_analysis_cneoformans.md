# Cryptococcus neoformans
## RNA-seq files
Fausto Almeida, Julie M. Wolf, Thiago Aparecido da Silva, Carlos M. DeLeon-Rodriguez, Caroline Patini Rezende, André Moreira Pessoni, Fabrício Freitas Fernandes, Rafael Silva-Rocha, Roberto Martinez, Marcio L. Rodrigues, Maria Cristina Roque-Barreira & Arturo Casadevall

Galectin-3 impacts Cryptococcus neoformans infection through direct antifungal effects

RNA-seq DATA:
- R1 : RNAseq of C. noformans not treated (SRR6207454)	->	CN_R1.R1.fastq	CN_R1.R2.fastq
- R2 : RNAseq of C. noformans not treated (SRR6207453)	->	CN_R2.R1.fastq	CN_R2.R2.fastq
- R3 : RNAseq of C. noformans not treated (SRR6207452)	->	CN_R3.R1.fastq	CN_R3.R2.fastq
## Preparations
A bemeneti annotációs fájlokat át kell alakítani, hogy megfelelőek legyenek a Cufflinks számára illetve elő kell állítani azokat az annotációs és FASTA fájlokat amik az RRPM analízishez szükségesek. 
### Input:
- __Cryptococcus_neoformans_H99.genes.gff__ : Forrás a JGI
- __Cryptococcus_neoformans_H99.unmasked.fasta__ : Forrás a JGI
### Output:
- __cneoformans_onlygene.gtf__ : Csak a géneket tartalmazó annotációs fájl
- __cneoformans_onlyexon.gtf__ : Csak az exonokat tartalmazó annotációs fájl
- __cneoformans_fixed.gtf__ : Pozíciófixált annotációs fájl az RRPM számára
- __cneoformans_genes.fasta__ : Géneket tartalmazó fasta fájl

GTF fájlok elkészítése az RRPM számára
```
PREPARATIONS_cneoformans.R
```
GTF fájl bed formátumra alakítása a géneket tartalmazó FASTA fájl elkészítéséhez
```
gtf2bed < cneoformans_onlygene.gtf > cneoformans_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl elkészítése
```
bedtools getfasta -name -fo cneoformans_genes.fasta -fi Cryptococcus_neoformans_H99.unmasked.fasta -bed cneoformans_onlygene.gtf.bed
```
A géneket tartalmazó FASTA fájl headerjének a trimmelése
```
perl -pi -e 's/::.*//g' cneoformans_genes.fasta
```
## RRPM: STAR Alignment and Cufflinks Assembly
Megnézzük milyen hosszúak a gének és az intronok a STAR és a Cufflinks beállításához.

Elkészítünk egy fájlt ami tartalmazza az összes szükséges scriptet ami az RRPM futtatásához szükséges (STAR, Cufflinks)
### Input:
- __cneoformans_genes.fasta__
- __cneoformans_fixed.gtf__
### Output:
- __gene_length__
- A __Cufflinks__ __outputja__ ami tartalmazza az RRPM analízis fájljait
```
bioawk -c fastx '{ print $name, length($seq) }' < cneoformans_genes.fasta > gene_length
INTRON_LENGTH_cneoformans.R
```
A maximális transzkriptméret: 15614 --> max-bundle-length marad 250000

Intron min: 22 --> --alignIntronMin marad 3 

Intron max: 1306 --> --alignIntronMax marad 30000

Lefuttatjuk az RRPM analízist.
```
# STAR illesztéshez szükséges indexelés

mkdir cneoformans_DIR_RRPM
cp -r ./cneoformans_genome ./cneoformans_DIR_RRPM/cneoformans_RRPM_STARindex
cd ./cneoformans_DIR_RRPM/cneoformans_RRPM_STARindex
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir ../cneoformans_RRPM_STARindex --sjdbGTFfile ../cneoformans_RRPM_STARindex/cneoformans_fixed.gtf --genomeFastaFiles ../cneoformans_RRPM_STARindex/cneoformans_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./cneoformans_DIR_RRPM/cneoformans_RRPM_STARout
cd ./cneoformans_DIR_RRPM/cneoformans_RRPM_STARout
STAR --runThreadN 48 --genomeDir ../cneoformans_RRPM_STARindex --readFilesIn ../../raw_data/CN_R1.R1.fastq.gz,../../raw_data/CN_R2.R1.fastq.gz,../../raw_data/CN_R3.R1.fastq.gz ../../raw_data/CN_R1.R2.fastq.gz,../../raw_data/CN_R2.R2.fastq.gz,../../raw_data/CN_R3.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 48 -g ./cneoformans_DIR_RRPM/cneoformans_RRPM_STARindex/cneoformans_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./cneoformans_DIR_RRPM/cneoformans_RRPM_CUFFout ./cneoformans_DIR_RRPM/cneoformans_RRPM_STARout/Aligned.sortedByCoord.out.bam
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
- __cneoformans_onlygene.gtf__
- __transcripts.gtf__ : RRPM Cufflinks assembly output file
### Output:
- __cneoformans_RRPM_transcripts.gtf__ : Az újonnan felfedezett transzkripteket tartalmazó GTF
- __cneoformans_AS_annotation.gtf__ : Az eredeti annotációval történő összemergeelés, hogy az esetlegesen nem detektált gének is bennelegyenek az annotációs fájlban.
- __cneoformans_stats.log__ : Szűrési statisztikákat tartalmazó logfájl
```
FILTERING_cneoformans.R
```
## Fusion correction
### Input:
- __isoforms.fpkm_tracking__ : Az RRPM Cufflinks output izoforms FPKM értékeket tartalmazó fájlja
- __cneoformans_AS_annotation.gtf__
### Output:
- Nincs fúziós gén detektálva

Megvizsgáljuk, hogy melyek azok a gének amik fúziósak voltak az eredeti annotációba de már két külön gént alkotnak, majd ezeket a géneket manuális IGV megtekintés után átnevezzük v1, v2 satöbbire.
```
FUSION_FILTER_cneoformans.R
```
## Expression: STAR Alignment, Cuffquant Assembly and Cuffdiff analysis
### Output:
- A __CuffDiff__ __output__ fájljai az expreszziós analízisből

Futtatunk egy teljes STAR illesztést majd egy Cuffquantot, majd kinyerjük az expressziós értékeket a CuffDiff segítségével.
```
# STAR INDEXING

mkdir DIR_cneoformans_expressionanalysis
cp -r ./cneoformans_genome ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex
cd ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../cneoformans_STARindex/cneoformans_AS_annotation.gtf --genomeFastaFiles ./Cryptococcus_neoformans_H99.unmasked.fasta --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# CN

mkdir ./DIR_cneoformans_expressionanalysis/CN_R1_STARout
cd ./DIR_cneoformans_expressionanalysis/CN_R1_STARout
STAR --runThreadN 32 --genomeDir ../cneoformans_STARindex --readFilesIn ../../raw_data/CN_R1.R1.fastq.gz ../../raw_data/CN_R1.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_cneoformans_expressionanalysis/CN_R2_STARout
cd ./DIR_cneoformans_expressionanalysis/CN_R2_STARout
STAR --runThreadN 32 --genomeDir ../cneoformans_STARindex --readFilesIn ../../raw_data/CN_R2.R1.fastq.gz ../../raw_data/CN_R2.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_cneoformans_expressionanalysis/CN_R3_STARout
cd ./DIR_cneoformans_expressionanalysis/CN_R3_STARout
STAR --runThreadN 32 --genomeDir ../cneoformans_STARindex --readFilesIn ../../raw_data/CN_R3.R1.fastq.gz ../../raw_data/CN_R3.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_cneoformans_expressionanalysis/CN_R1_CUFFquant -b ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex/Cryptococcus_neoformans_H99.unmasked.fasta -p 32 ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex/cneoformans_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_cneoformans_expressionanalysis/CN_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_cneoformans_expressionanalysis/CN_R2_CUFFquant -b ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex/Cryptococcus_neoformans_H99.unmasked.fasta -p 32 ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex/cneoformans_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_cneoformans_expressionanalysis/CN_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_cneoformans_expressionanalysis/CN_R3_CUFFquant -b ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex/Cryptococcus_neoformans_H99.unmasked.fasta -p 32 ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex/cneoformans_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_cneoformans_expressionanalysis/CN_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex/Cryptococcus_neoformans_H99.unmasked.fasta -p 32 -L CN_1,CN_2 --upper-quartile-norm ./DIR_cneoformans_expressionanalysis/cneoformans_STARindex/cneoformans_AS_annotation.gtf \
./DIR_cneoformans_expressionanalysis/CN_R1_CUFFquant/abundances.cxb,./DIR_cneoformans_expressionanalysis/CN_R2_CUFFquant/abundances.cxb,./DIR_cneoformans_expressionanalysis/CN_R3_CUFFquant/abundances.cxb \
./DIR_cneoformans_expressionanalysis/CN_R1_CUFFquant/abundances.cxb,./DIR_cneoformans_expressionanalysis/CN_R2_CUFFquant/abundances.cxb,./DIR_cneoformans_expressionanalysis/CN_R3_CUFFquant/abundances.cxb
```
## ORF prediction
### Input:
- __cneoformans_AS_annotation.gtf__
- __Cryptococcus_neoformans_H99.unmasked.fasta__
### Output:
- __cneoformans_transcripts.fasta__ : A TransDecoder outputja ami tartalmazza az összes transzkript szekvenciáját fasta formátumban
- __longest_orfs.cds__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF cds régióját
- __longest_orfs.pep__ : A TransDecoder outputja ami tartalmazza az összes prediktált ORF proteinjeit

Prediktáljuk az ORF régiókat TransDecoder segítségével.
```
~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl cneoformans_AS_annotation.gtf Cryptococcus_neoformans_H99.unmasked.fasta > cneoformans_transcripts.fasta
TransDecoder.LongOrfs -m 20 -S -t cneoformans_transcripts.fasta
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
INTERPROSCAN_cneoformans.R
```
Majd ezzel megszűrjük a fehérjeszekvenciákat tartalmazó fájlt, végül átalakítjuk úgy a FASTA-fájlt, hogy az InterProScan számára megfelelő legyen. (Elég körülményesen oldottam ezt meg de működik, nyilván mostmár máshogy írnám még ezt a részt egy R scriptbe)
```
# átalakítás
\)\n --> \)\t
# szűrés
while read l; do perl -n -e "print if /^>$l::.*/" longest_orfs.pep >> cneoformans_proteins_all.fasta; done <filter
# visszalakítás
\)\t --> \)\n
```
Töröljük a fasta fájlból a *-okat mivel az InterProScan nem tudja értelmezni
```
perl -pi -e 's/\*//g' cneoformans_proteins_all.fasta
```
Egy kis regexp átalakítás az áttekinthetőségért
```
perl -pi -e 's/::g\..*$//g' cneoformans_proteins_all.fasta
perl -pi -e 's/^>.*::/>/g' cneoformans_proteins_all.fasta
```
## InterProScan analysis
### Input:
- __cneoformans_proteins_all.fasta__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinek
### Output:
- __cneoformans_proteins_all.fasta.tsv__ : Minden transzkriptre a leghosszabb ORF-ek alapján prediktált proteinhez tartozó InterProScan domainek

Lefuttatjuk at InterProScan-t.
```
interproscan.sh -i cneoformans_proteins_all.fasta -f tsv --iprlookup --goterms -dp
```