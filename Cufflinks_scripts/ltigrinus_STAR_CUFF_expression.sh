# STAR INDEXING

mkdir ltigrinus_DIR_expressionanalysis
cp -r ./ltigrinus_genome ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex
cd ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --genomeFastaFiles ./Sisbr1_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# LT_VM

mkdir ./ltigrinus_DIR_expressionanalysis/LT_VM_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_VM_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_VM_R1.R1.fastq.gz ../../raw_data/LT_VM_R1.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_VM_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_VM_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_VM_R2.R1.fastq.gz ../../raw_data/LT_VM_R2.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_VM_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_VM_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_VM_R3.all.R1.fastq.gz ../../raw_data/LT_VM_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_VM_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_VM_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_VM_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_P1

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P1_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P1_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P1_R1.all.R1.fastq.gz ../../raw_data/LT_P1_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P1_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P1_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P1_R2.all.R1.fastq.gz ../../raw_data/LT_P1_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P1_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P1_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P1_R3.all.R1.fastq.gz ../../raw_data/LT_P1_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P1_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P1_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P1_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P1_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P1_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P1_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_P2_K

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_K_R1.all.R1.fastq.gz ../../raw_data/LT_P2_K_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_K_R2.all.R1.fastq.gz ../../raw_data/LT_P2_K_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_K_R3.all.R1.fastq.gz ../../raw_data/LT_P2_K_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_P2_T

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_T_R1.all.R1.fastq.gz ../../raw_data/LT_P2_T_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_T_R2.all.R1.fastq.gz ../../raw_data/LT_P2_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_P2_T_R3.all.R1.fastq.gz ../../raw_data/LT_P2_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_YFB_K

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_K_R1.all.R1.fastq.gz ../../raw_data/LT_YFB_K_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_K_R2.all.R1.fastq.gz ../../raw_data/LT_YFB_K_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_K_R3.all.R1.fastq.gz ../../raw_data/LT_YFB_K_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_YFB_T

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_T_R1.all.R1.fastq.gz ../../raw_data/LT_YFB_T_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_T_R2.all.R1.fastq.gz ../../raw_data/LT_YFB_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_YFB_T_R3.all.R1.fastq.gz ../../raw_data/LT_YFB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_STARout/Aligned.sortedByCoord.out.bam

# LT_FB_K

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_K_R1.all.R1.fastq.gz ../../raw_data/LT_FB_K_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_K_R3.all.R1.fastq.gz ../../raw_data/LT_FB_K_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/TL_FB_K_R4_R1_001.fastq.gz ../../raw_data/TL_FB_K_R4_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_STARout/Aligned.sortedByCoord.out.bam

# LT_FB_L

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_L_R1.R1.fastq.gz ../../raw_data/LT_FB_L_R1.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_L_R2.all.R1.fastq.gz ../../raw_data/LT_FB_L_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/TL_FB_L_R4_R1_001.fastq.gz ../../raw_data/TL_FB_L_R4_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_STARout/Aligned.sortedByCoord.out.bam

# LT_FB_T

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_T_R2.all.R1.fastq.gz ../../raw_data/LT_FB_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R3_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R3_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/LT_FB_T_R3.all.R1.fastq.gz ../../raw_data/LT_FB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_STARout
cd ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_STARout
STAR --runThreadN 16 --genomeDir ../ltigrinus_STARindex --readFilesIn ../../raw_data/TL_FB_T_R4_R1_001.fastq.gz ../../raw_data/TL_FB_T_R4_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R1_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_CUFFquant -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/ltigrinus_AS_annotation.gtf --max-bundle-frags 1000000 ./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./ltigrinus_DIR_expressionanalysis/ltigrinus_STARindex/Sisbr1_AssemblyScaffolds.fasta -p 16 -L VM,P1,P2_K,P2_T,YFB_K,YFB_T,FB_K,FB_L,FB_T --upper-quartile-norm ./ltigrinus_genome/ltigrinus_AS_annotation.gtf \
./ltigrinus_DIR_expressionanalysis/LT_VM_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_VM_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_VM_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_P1_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P1_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P1_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_P2_K_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P2_K_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P2_K_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_P2_T_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P2_T_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_P2_T_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_YFB_K_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_YFB_T_R3_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_FB_K_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_K_R3_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_K_R4_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_FB_L_R1_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_L_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_L_R4_CUFFquant/abundances.cxb \
./ltigrinus_DIR_expressionanalysis/LT_FB_T_R2_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_T_R3_CUFFquant/abundances.cxb,./ltigrinus_DIR_expressionanalysis/LT_FB_T_R4_CUFFquant/abundances.cxb
