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

