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










