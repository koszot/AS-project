# STAR illesztéshez szükséges indexelés

mkdir aampla_DIR_RRPM
cp -r ./aampla_genome ./aampla_DIR_RRPM/aampla_RRPM_STARindex
cd ./aampla_DIR_RRPM/aampla_RRPM_STARindex
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir ../aampla_RRPM_STARindex --sjdbGTFfile ../aampla_RRPM_STARindex/aampla_fixed.gtf --genomeFastaFiles ../aampla_RRPM_STARindex/aampla_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./aampla_DIR_RRPM/aampla_RRPM_STARout
cd ./aampla_DIR_RRPM/aampla_RRPM_STARout
STAR --runThreadN 48 --genomeDir ../aampla_RRPM_STARindex --readFilesIn ../../raw_data/AUA_VM_R1.all.R1.fastq.gz,../../raw_data/AUA_VM_R2.all.R1.fastq.gz,../../raw_data/AUA_VM_R3.all.R1.fastq.gz,../../raw_data/AUA_P1_R1.all.R1.fastq.gz,../../raw_data/AUA_P1_R2.all.R1.fastq.gz,../../raw_data/AUA_P1_R3.all.R1.fastq.gz,../../raw_data/AUA_P2_R1.all.R1.fastq.gz,../../raw_data/AUA_P2_R2.all.R1.fastq.gz,../../raw_data/AUA_P2_R3.all.R1.fastq.gz,../../raw_data/AUA_YFB_R1.all.R1.fastq.gz,../../raw_data/AUA_YFB_R2.all.R1.fastq.gz,../../raw_data/AUA_YFB_R3.all.R1.fastq.gz,../../raw_data/AUA_FB_R1.all.R1.fastq.gz,../../raw_data/AUA_FB_R2.all.R1.fastq.gz,../../raw_data/AUA_FB_R3.all.R1.fastq.gz ../../raw_data/AUA_VM_R1.all.R2.fastq.gz,../../raw_data/AUA_VM_R2.all.R2.fastq.gz,../../raw_data/AUA_VM_R3.all.R2.fastq.gz,../../raw_data/AUA_P1_R1.all.R2.fastq.gz,../../raw_data/AUA_P1_R2.all.R2.fastq.gz,../../raw_data/AUA_P1_R3.all.R2.fastq.gz,../../raw_data/AUA_P2_R1.all.R2.fastq.gz,../../raw_data/AUA_P2_R2.all.R2.fastq.gz,../../raw_data/AUA_P2_R3.all.R2.fastq.gz,../../raw_data/AUA_YFB_R1.all.R2.fastq.gz,../../raw_data/AUA_YFB_R2.all.R2.fastq.gz,../../raw_data/AUA_YFB_R3.all.R2.fastq.gz,../../raw_data/AUA_FB_R1.all.R2.fastq.gz,../../raw_data/AUA_FB_R2.all.R2.fastq.gz,../../raw_data/AUA_FB_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 48 -g ./aampla_DIR_RRPM/aampla_RRPM_STARindex/aampla_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./aampla_DIR_RRPM/aampla_RRPM_CUFFout ./aampla_DIR_RRPM/aampla_RRPM_STARout/Aligned.sortedByCoord.out.bam




