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

