# STAR illesztéshez szükséges indexelés

mkdir DIR_RRPM
cp -r ./genome ./DIR_RRPM/RRPM_STARindex
cd ./DIR_RRPM/RRPM_STARindex
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir ../RRPM_STARindex --sjdbGTFfile ../RRPM_STARindex/original_fixed.gtf --genomeFastaFiles ../RRPM_STARindex/original_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./DIR_RRPM/RRPM_STARout
cd ./DIR_RRPM/RRPM_STARout
STAR --runThreadN 16 --genomeDir ../RRPM_STARindex --readFilesIn ../../raw_data/SC_VM_R1.all.R2.fastq.gz,../../raw_data/SC_VM_R2.all.R2.fastq.gz,../../raw_data/SC_VM_R3.all.R2.fastq.gz,../../raw_data/SC_P1_R1.all.R2.fastq.gz,../../raw_data/SC_P1_R2.all.R2.fastq.gz,../../raw_data/SC_P1_R3.all.R2.fastq.gz,../../raw_data/SC_P2_R1.all.R2.fastq.gz,../../raw_data/SC_P2_R2.all.R2.fastq.gz,../../raw_data/SC_P2_R3.all.R2.fastq.gz,../../raw_data/SC_YFB_R1.all.R2.fastq.gz,../../raw_data/SC_YFB_R2.all.R2.fastq.gz,../../raw_data/SC_YFB_R3.all.R2.fastq.gz,../../raw_data/SC_FB_R1.all.R2.fastq.gz,../../raw_data/SC_FB_R2.all.R2.fastq.gz,../../raw_data/SC_FB_R3.all.R2.fastq.gz ../../raw_data/SC_VM_R1.all.R1.fastq.gz,../../raw_data/SC_VM_R2.all.R1.fastq.gz,../../raw_data/SC_VM_R3.all.R1.fastq.gz,../../raw_data/SC_P1_R1.all.R1.fastq.gz,../../raw_data/SC_P1_R2.all.R1.fastq.gz,../../raw_data/SC_P1_R3.all.R1.fastq.gz,../../raw_data/SC_P2_R1.all.R1.fastq.gz,../../raw_data/SC_P2_R2.all.R1.fastq.gz,../../raw_data/SC_P2_R3.all.R1.fastq.gz,../../raw_data/SC_YFB_R1.all.R1.fastq.gz,../../raw_data/SC_YFB_R2.all.R1.fastq.gz,../../raw_data/SC_YFB_R3.all.R1.fastq.gz,../../raw_data/SC_FB_R1.all.R1.fastq.gz,../../raw_data/SC_FB_R2.all.R1.fastq.gz,../../raw_data/SC_FB_R3.all.R1.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# STAR illesztés az RRPM-hez !!!HA AZ SORTOLÁST NEM TUDJA MEGCSINÁLNI A STAR!!!

mkdir ./DIR_RRPM/RRPM_STARout
cd ./DIR_RRPM/RRPM_STARout
STAR --runThreadN 12 --genomeDir ../RRPM_STARindex --readFilesIn ../../raw_data/SC_VM_R1.all.R2.fastq.gz,../../raw_data/SC_VM_R2.all.R2.fastq.gz,../../raw_data/SC_VM_R3.all.R2.fastq.gz,../../raw_data/SC_P1_R1.all.R2.fastq.gz,../../raw_data/SC_P1_R2.all.R2.fastq.gz,../../raw_data/SC_P1_R3.all.R2.fastq.gz,../../raw_data/SC_P2_R1.all.R2.fastq.gz,../../raw_data/SC_P2_R2.all.R2.fastq.gz,../../raw_data/SC_P2_R3.all.R2.fastq.gz,../../raw_data/SC_YFB_R1.all.R2.fastq.gz,../../raw_data/SC_YFB_R2.all.R2.fastq.gz,../../raw_data/SC_YFB_R3.all.R2.fastq.gz,../../raw_data/SC_FB_R1.all.R2.fastq.gz,../../raw_data/SC_FB_R2.all.R2.fastq.gz,../../raw_data/SC_FB_R3.all.R2.fastq.gz ../../raw_data/SC_VM_R1.all.R1.fastq.gz,../../raw_data/SC_VM_R2.all.R1.fastq.gz,../../raw_data/SC_VM_R3.all.R1.fastq.gz,../../raw_data/SC_P1_R1.all.R1.fastq.gz,../../raw_data/SC_P1_R2.all.R1.fastq.gz,../../raw_data/SC_P1_R3.all.R1.fastq.gz,../../raw_data/SC_P2_R1.all.R1.fastq.gz,../../raw_data/SC_P2_R2.all.R1.fastq.gz,../../raw_data/SC_P2_R3.all.R1.fastq.gz,../../raw_data/SC_YFB_R1.all.R1.fastq.gz,../../raw_data/SC_YFB_R2.all.R1.fastq.gz,../../raw_data/SC_YFB_R3.all.R1.fastq.gz,../../raw_data/SC_FB_R1.all.R1.fastq.gz,../../raw_data/SC_FB_R2.all.R1.fastq.gz,../../raw_data/SC_FB_R3.all.R1.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated  --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# a BAM file sortolása !!!HA AZ SORTOLÁST NEM TUDJA MEGCSINÁLNI A STAR!!!

sambamba sort -m 10GB -t 16 -o ./DIR_RRPM/RRPM_STARout/Aligned.sortedByCoord.out.bam -p ./DIR_RRPM/RRPM_STARout/Aligned.out.bam

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 16 -g ./DIR_RRPM/RRPM_STARindex/original_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./DIR_RRPM/RRPM_CUFFout ./DIR_RRPM/RRPM_STARout/Aligned.sortedByCoord.out.bam

