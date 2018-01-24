# STAR illesztéshez szükséges indexelés

mkdir DIR_RRPM
cp -r ./genome ./DIR_RRPM/RRPM_STARindex
cd ./DIR_RRPM/RRPM_STARindex
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ../RRPM_STARindex --sjdbGTFfile ../RRPM_STARindex/original_fixed.gtf --genomeFastaFiles ../RRPM_STARindex/original_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# kiszámoljuk az intron hosszakat a GFF3 fájlból az aostoyae_intron_length.pl segítségével
# min: 21
# max: 6767
# nem változtatom meg a beállításokat MERT a legkisebb maradjon egy értelmes reading frame, a maximumnál meg a gének határa úgyis limitál az RRPM alatt

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./DIR_RRPM/RRPM_STARout
cd ./DIR_RRPM/RRPM_STARout
STAR --runThreadN 16 --genomeDir ../RRPM_STARindex --readFilesIn ../../raw_data/ARO_FB_C_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_C_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_C_R3.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R3.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R3.all.R1.fastq.gz,../../raw_data/ARO_P1_R1.all.R1.fastq.gz,../../raw_data/ARO_P1_R2.all.R1.fastq.gz,../../raw_data/ARO_P1_R3.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R1.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R2.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R3.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R1.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R2.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R3.all.R1.fastq.gz,../../raw_data/ARO_RMA_R1.all.R1.fastq.gz,../../raw_data/ARO_RMA_R2.all.R1.fastq.gz,../../raw_data/ARO_RMA_R3.all.R1.fastq.gz,../../raw_data/ARO_VM_R1.all.R1.fastq.gz,../../raw_data/ARO_VM_R2.all.R1.fastq.gz,../../raw_data/ARO_VM_R3.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R1.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R2.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R3.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R1.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R2.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_C_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_C_R3.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R3.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R3.all.R2.fastq.gz,../../raw_data/ARO_P1_R1.all.R2.fastq.gz,../../raw_data/ARO_P1_R2.all.R2.fastq.gz,../../raw_data/ARO_P1_R3.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R1.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R2.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R3.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R1.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R2.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R3.all.R2.fastq.gz,../../raw_data/ARO_RMA_R1.all.R2.fastq.gz,../../raw_data/ARO_RMA_R2.all.R2.fastq.gz,../../raw_data/ARO_RMA_R3.all.R2.fastq.gz,../../raw_data/ARO_VM_R1.all.R2.fastq.gz,../../raw_data/ARO_VM_R2.all.R2.fastq.gz,../../raw_data/ARO_VM_R3.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R1.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R2.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R3.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R1.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R2.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# STAR illesztés az RRPM-hez !!!HA AZ SORTOLÁST NEM TUDJA MEGCSINÁLNI A STAR!!!

mkdir ./DIR_RRPM/RRPM_STARout
cd ./DIR_RRPM/RRPM_STARout
STAR --runThreadN 12 --genomeDir ../RRPM_STARindex --readFilesIn ../../raw_data/ARO_FB_C_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_C_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_C_R3.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_L_R3.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R1.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R2.all.R1.fastq.gz,../../raw_data/ARO_FB_S_R3.all.R1.fastq.gz,../../raw_data/ARO_P1_R1.all.R1.fastq.gz,../../raw_data/ARO_P1_R2.all.R1.fastq.gz,../../raw_data/ARO_P1_R3.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R1.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R2.all.R1.fastq.gz,../../raw_data/ARO_P2_C_R3.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R1.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R2.all.R1.fastq.gz,../../raw_data/ARO_P2_S_R3.all.R1.fastq.gz,../../raw_data/ARO_RMA_R1.all.R1.fastq.gz,../../raw_data/ARO_RMA_R2.all.R1.fastq.gz,../../raw_data/ARO_RMA_R3.all.R1.fastq.gz,../../raw_data/ARO_VM_R1.all.R1.fastq.gz,../../raw_data/ARO_VM_R2.all.R1.fastq.gz,../../raw_data/ARO_VM_R3.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R1.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R2.all.R1.fastq.gz,../../raw_data/ARO_YFB_C_R3.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R1.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R2.all.R1.fastq.gz,../../raw_data/ARO_YFB_S_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_C_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_C_R3.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_L_R3.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R1.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R2.all.R2.fastq.gz,../../raw_data/ARO_FB_S_R3.all.R2.fastq.gz,../../raw_data/ARO_P1_R1.all.R2.fastq.gz,../../raw_data/ARO_P1_R2.all.R2.fastq.gz,../../raw_data/ARO_P1_R3.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R1.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R2.all.R2.fastq.gz,../../raw_data/ARO_P2_C_R3.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R1.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R2.all.R2.fastq.gz,../../raw_data/ARO_P2_S_R3.all.R2.fastq.gz,../../raw_data/ARO_RMA_R1.all.R2.fastq.gz,../../raw_data/ARO_RMA_R2.all.R2.fastq.gz,../../raw_data/ARO_RMA_R3.all.R2.fastq.gz,../../raw_data/ARO_VM_R1.all.R2.fastq.gz,../../raw_data/ARO_VM_R2.all.R2.fastq.gz,../../raw_data/ARO_VM_R3.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R1.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R2.all.R2.fastq.gz,../../raw_data/ARO_YFB_C_R3.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R1.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R2.all.R2.fastq.gz,../../raw_data/ARO_YFB_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated  --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# a BAM file sortolása !!!HA AZ SORTOLÁST NEM TUDJA MEGCSINÁLNI A STAR!!!

sambamba sort -m 10GB -t 16 -o ./DIR_RRPM/RRPM_STARout/Aligned.sortedByCoord.out.bam -p ./DIR_RRPM/RRPM_STARout/Aligned.out.bam

# kiszámoljuk a gének maximális mértetét hogy azt tudjuk megadni max-bundle-length -nek

#/Users/genomelab3/Downloads/bioawk-master/bioawk -c fastx '{ print $name, length($seq) }' <original_genes.fasta > gene_length

# max méret : 16175
# ezt se változatatom meg egyenlőre

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 16 -g ./DIR_RRPM/RRPM_STARindex/original_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./DIR_RRPM/RRPM_CUFFout ./DIR_RRPM/RRPM_STARout/Aligned.sortedByCoord.out.bam
