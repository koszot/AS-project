# STAR INDEXING

mkdir DIR_aampla_expressionanalysis
cp -r ./aampla_genome ./DIR_aampla_expressionanalysis/aampla_STARindex
cd ./DIR_aampla_expressionanalysis/aampla_STARindex
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../aampla_STARindex/aampla_AS_annotation.gtf --genomeFastaFiles ./Auramp1_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# AUA_VM

mkdir ./DIR_aampla_expressionanalysis/AUA_VM_R1_STARout
cd ./DIR_aampla_expressionanalysis/AUA_VM_R1_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_VM_R1.all.R1.fastq.gz ../../raw_data/AUA_VM_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_VM_R2_STARout
cd ./DIR_aampla_expressionanalysis/AUA_VM_R2_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_VM_R2.all.R1.fastq.gz ../../raw_data/AUA_VM_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_VM_R3_STARout
cd ./DIR_aampla_expressionanalysis/AUA_VM_R3_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_VM_R3.all.R1.fastq.gz ../../raw_data/AUA_VM_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_VM_R1_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_VM_R2_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_VM_R3_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# AUA_P1

mkdir ./DIR_aampla_expressionanalysis/AUA_P1_R1_STARout
cd ./DIR_aampla_expressionanalysis/AUA_P1_R1_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_P1_R1.all.R1.fastq.gz ../../raw_data/AUA_P1_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_P1_R2_STARout
cd ./DIR_aampla_expressionanalysis/AUA_P1_R2_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_P1_R2.all.R1.fastq.gz ../../raw_data/AUA_P1_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_P1_R3_STARout
cd ./DIR_aampla_expressionanalysis/AUA_P1_R3_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_P1_R3.all.R1.fastq.gz ../../raw_data/AUA_P1_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_P1_R1_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_P1_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_P1_R2_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_P1_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_P1_R3_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_P1_R3_STARout/Aligned.sortedByCoord.out.bam

# AUA_P2

mkdir ./DIR_aampla_expressionanalysis/AUA_P2_R1_STARout
cd ./DIR_aampla_expressionanalysis/AUA_P2_R1_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_P2_R1.all.R1.fastq.gz ../../raw_data/AUA_P2_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_P2_R2_STARout
cd ./DIR_aampla_expressionanalysis/AUA_P2_R2_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_P2_R2.all.R1.fastq.gz ../../raw_data/AUA_P2_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_P2_R3_STARout
cd ./DIR_aampla_expressionanalysis/AUA_P2_R3_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_P2_R3.all.R1.fastq.gz ../../raw_data/AUA_P2_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_P2_R1_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_P2_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_P2_R2_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_P2_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_P2_R3_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_P2_R3_STARout/Aligned.sortedByCoord.out.bam

# AUA_YFB

mkdir ./DIR_aampla_expressionanalysis/AUA_YFB_R1_STARout
cd ./DIR_aampla_expressionanalysis/AUA_YFB_R1_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_YFB_R1.all.R1.fastq.gz ../../raw_data/AUA_YFB_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_YFB_R2_STARout
cd ./DIR_aampla_expressionanalysis/AUA_YFB_R2_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_YFB_R2.all.R1.fastq.gz ../../raw_data/AUA_YFB_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_YFB_R3_STARout
cd ./DIR_aampla_expressionanalysis/AUA_YFB_R3_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_YFB_R3.all.R1.fastq.gz ../../raw_data/AUA_YFB_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_YFB_R1_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_YFB_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_YFB_R2_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_YFB_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_YFB_R3_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_YFB_R3_STARout/Aligned.sortedByCoord.out.bam

# AUA_FB

mkdir ./DIR_aampla_expressionanalysis/AUA_FB_R1_STARout
cd ./DIR_aampla_expressionanalysis/AUA_FB_R1_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_FB_R1.all.R1.fastq.gz ../../raw_data/AUA_FB_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_FB_R2_STARout
cd ./DIR_aampla_expressionanalysis/AUA_FB_R2_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_FB_R2.all.R1.fastq.gz ../../raw_data/AUA_FB_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_aampla_expressionanalysis/AUA_FB_R3_STARout
cd ./DIR_aampla_expressionanalysis/AUA_FB_R3_STARout
STAR --runThreadN 32 --genomeDir ../aampla_STARindex --readFilesIn ../../raw_data/AUA_FB_R3.all.R1.fastq.gz ../../raw_data/AUA_FB_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_FB_R1_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_FB_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_FB_R2_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_FB_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_aampla_expressionanalysis/AUA_FB_R3_CUFFquant -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_aampla_expressionanalysis/AUA_FB_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./DIR_aampla_expressionanalysis/aampla_STARindex/Auramp1_AssemblyScaffolds.fasta -p 32 -L VM,P1,P2,YFB,FB --upper-quartile-norm ./DIR_aampla_expressionanalysis/aampla_STARindex/aampla_AS_annotation.gtf \
./DIR_aampla_expressionanalysis/AUA_VM_R1_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_VM_R2_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_VM_R3_CUFFquant/abundances.cxb \
./DIR_aampla_expressionanalysis/AUA_P1_R1_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_P1_R2_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_P1_R3_CUFFquant/abundances.cxb \
./DIR_aampla_expressionanalysis/AUA_P2_R1_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_P2_R2_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_P2_R3_CUFFquant/abundances.cxb \
./DIR_aampla_expressionanalysis/AUA_YFB_R1_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_YFB_R2_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_YFB_R3_CUFFquant/abundances.cxb \
./DIR_aampla_expressionanalysis/AUA_FB_R1_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_FB_R2_CUFFquant/abundances.cxb,./DIR_aampla_expressionanalysis/AUA_FB_R3_CUFFquant/abundances.cxb
