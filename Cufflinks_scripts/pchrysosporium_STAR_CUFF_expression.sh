# STAR INDEXING

mkdir DIR_pchrysosporium_expressionanalysis
cp -r ./pchrysosporium_genome ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex
cd ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf --genomeFastaFiles ./Phchr2_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# PH_VM

mkdir ./DIR_pchrysosporium_expressionanalysis/PH_VM_R1_STARout
cd ./DIR_pchrysosporium_expressionanalysis/PH_VM_R1_STARout
STAR --runThreadN 32 --genomeDir ../pchrysosporium_STARindex --readFilesIn ../../raw_data/PH_VM_2_S18_R1_001.fastq.gz ../../raw_data/PH_VM_2_S18_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_pchrysosporium_expressionanalysis/PH_VM_R2_STARout
cd ./DIR_pchrysosporium_expressionanalysis/PH_VM_R2_STARout
STAR --runThreadN 32 --genomeDir ../pchrysosporium_STARindex --readFilesIn ../../raw_data/PH_VM_5_S19_R1_001.fastq.gz ../../raw_data/PH_VM_5_S19_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_pchrysosporium_expressionanalysis/PH_VM_R3_STARout
cd ./DIR_pchrysosporium_expressionanalysis/PH_VM_R3_STARout
STAR --runThreadN 32 --genomeDir ../pchrysosporium_STARindex --readFilesIn ../../raw_data/PH_VM_6_S20_R1_001.fastq.gz ../../raw_data/PH_VM_6_S20_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_pchrysosporium_expressionanalysis/PH_VM_R1_CUFFquant -b ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/Phchr2_AssemblyScaffolds.fasta -p 32 ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_pchrysosporium_expressionanalysis/PH_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_pchrysosporium_expressionanalysis/PH_VM_R2_CUFFquant -b ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/Phchr2_AssemblyScaffolds.fasta -p 32 ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_pchrysosporium_expressionanalysis/PH_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_pchrysosporium_expressionanalysis/PH_VM_R3_CUFFquant -b ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/Phchr2_AssemblyScaffolds.fasta -p 32 ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_pchrysosporium_expressionanalysis/PH_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# PH_YFB

mkdir ./DIR_pchrysosporium_expressionanalysis/PH_YFB_R1_STARout
cd ./DIR_pchrysosporium_expressionanalysis/PH_YFB_R1_STARout
STAR --runThreadN 32 --genomeDir ../pchrysosporium_STARindex --readFilesIn ../../raw_data/PH_YFB_R21_S22_R1_001.fastq.gz ../../raw_data/PH_YFB_R21_S22_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_pchrysosporium_expressionanalysis/PH_YFB_R2_STARout
cd ./DIR_pchrysosporium_expressionanalysis/PH_YFB_R2_STARout
STAR --runThreadN 32 --genomeDir ../pchrysosporium_STARindex --readFilesIn ../../raw_data/PH_YFB_R24_S23_R1_001.fastq.gz ../../raw_data/PH_YFB_R24_S23_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_pchrysosporium_expressionanalysis/PH_YFB_R1_CUFFquant -b ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/Phchr2_AssemblyScaffolds.fasta -p 32 ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_pchrysosporium_expressionanalysis/PH_YFB_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_pchrysosporium_expressionanalysis/PH_YFB_R2_CUFFquant -b ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/Phchr2_AssemblyScaffolds.fasta -p 32 ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_pchrysosporium_expressionanalysis/PH_YFB_R2_STARout/Aligned.sortedByCoord.out.bam

# PH_FB

mkdir ./DIR_pchrysosporium_expressionanalysis/PH_FB_R1_STARout
cd ./DIR_pchrysosporium_expressionanalysis/PH_FB_R1_STARout
STAR --runThreadN 32 --genomeDir ../pchrysosporium_STARindex --readFilesIn ../../raw_data/PH_FB_10_S74_R1_001.fastq.gz ../../raw_data/PH_FB_10_S74_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_pchrysosporium_expressionanalysis/PH_FB_R2_STARout
cd ./DIR_pchrysosporium_expressionanalysis/PH_FB_R2_STARout
STAR --runThreadN 32 --genomeDir ../pchrysosporium_STARindex --readFilesIn ../../raw_data/PH_FB_12_S17_R1_001.fastq.gz ../../raw_data/PH_FB_12_S17_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_pchrysosporium_expressionanalysis/PH_FB_R1_CUFFquant -b ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/Phchr2_AssemblyScaffolds.fasta -p 32 ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_pchrysosporium_expressionanalysis/PH_FB_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_pchrysosporium_expressionanalysis/PH_FB_R2_CUFFquant -b ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/Phchr2_AssemblyScaffolds.fasta -p 32 ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_pchrysosporium_expressionanalysis/PH_FB_R2_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/Phchr2_AssemblyScaffolds.fasta -p 32 -L VM,YFB,FB --upper-quartile-norm ./DIR_pchrysosporium_expressionanalysis/pchrysosporium_STARindex/pchrysosporium_AS_annotation.gtf \
./DIR_pchrysosporium_expressionanalysis/PH_VM_R1_CUFFquant/abundances.cxb,./DIR_pchrysosporium_expressionanalysis/PH_VM_R2_CUFFquant/abundances.cxb,./DIR_pchrysosporium_expressionanalysis/PH_VM_R3_CUFFquant/abundances.cxb \
./DIR_pchrysosporium_expressionanalysis/PH_YFB_R1_CUFFquant/abundances.cxb,./DIR_pchrysosporium_expressionanalysis/PH_YFB_R2_CUFFquant/abundances.cxb \
./DIR_pchrysosporium_expressionanalysis/PH_FB_R1_CUFFquant/abundances.cxb,./DIR_pchrysosporium_expressionanalysis/PH_FB_R2_CUFFquant/abundances.cxb
