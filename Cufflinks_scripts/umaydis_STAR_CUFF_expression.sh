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










