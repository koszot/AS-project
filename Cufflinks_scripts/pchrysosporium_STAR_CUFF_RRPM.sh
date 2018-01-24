# STAR illesztéshez szükséges indexelés

mkdir pchrysosporium_DIR_RRPM
cp -r ./pchrysosporium_genome ./pchrysosporium_DIR_RRPM/pchrysosporium_RRPM_STARindex
cd ./pchrysosporium_DIR_RRPM/pchrysosporium_RRPM_STARindex
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ../pchrysosporium_RRPM_STARindex --sjdbGTFfile ../pchrysosporium_RRPM_STARindex/pchrysosporium_fixed.gtf --genomeFastaFiles ../pchrysosporium_RRPM_STARindex/pchrysosporium_genes.fasta --sjdbOverhang 100 --genomeChrBinNbits 12
cd ../..

# STAR illesztés az RRPM-hez  !!!HA AZ SORTOLÁST IS MEG TUDJA CSINÁLNI A STAR!!!

mkdir ./pchrysosporium_DIR_RRPM/pchrysosporium_RRPM_STARout
cd ./pchrysosporium_DIR_RRPM/pchrysosporium_RRPM_STARout
STAR --runThreadN 32 --genomeDir ../pchrysosporium_RRPM_STARindex --readFilesIn ../../raw_data/PH_VM_2_S18_R1_001.fastq.gz,../../raw_data/PH_VM_5_S19_R1_001.fastq.gz,../../raw_data/PH_VM_6_S20_R1_001.fastq.gz,../../raw_data/PH_YFB_R21_S22_R1_001.fastq.gz,../../raw_data/PH_YFB_R24_S23_R1_001.fastq.gz,../../raw_data/PH_FB_10_S74_R1_001.fastq.gz,../../raw_data/PH_FB_12_S17_R1_001.fastq.gz ../../raw_data/PH_VM_2_S18_R2_001.fastq.gz,../../raw_data/PH_VM_5_S19_R2_001.fastq.gz,../../raw_data/PH_VM_6_S20_R2_001.fastq.gz,../../raw_data/PH_YFB_R21_S22_R2_001.fastq.gz,../../raw_data/PH_YFB_R24_S23_R2_001.fastq.gz,../../raw_data/PH_FB_10_S74_R2_001.fastq.gz,../../raw_data/PH_FB_12_S17_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# Cufflinks assembly az alternatív splice variánsok meghatározására

cufflinks -p 32 -g ./pchrysosporium_DIR_RRPM/pchrysosporium_RRPM_STARindex/pchrysosporium_fixed.gtf --max-intron-length 30000 --min-intron-length 3 --overlap-radius 25 --max-bundle-length 250000 -o ./pchrysosporium_DIR_RRPM/pchrysosporium_RRPM_CUFFout ./pchrysosporium_DIR_RRPM/pchrysosporium_RRPM_STARout/Aligned.sortedByCoord.out.bam




