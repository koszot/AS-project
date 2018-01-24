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

