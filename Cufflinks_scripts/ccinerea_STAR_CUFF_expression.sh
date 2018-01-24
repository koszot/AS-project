# STAR INDEXING

mkdir DIR_ccinerea_AmutBmut_expressionanalysis
cp -r ./ccinerea_AmutBmut_genome ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --genomeFastaFiles ./Copci_AmutBmut1_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# CC_VM

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_VM_R1.all.R1.fastq.gz ../../raw_data/CC_VM_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_VM_R2.all.R1.fastq.gz ../../raw_data/CC_VM_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_VM_R3.all.R1.fastq.gz ../../raw_data/CC_VM_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_H

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_H_R2.all.R1.fastq.gz ../../raw_data/CC_H_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_H_R3.all.R1.fastq.gz ../../raw_data/CC_H_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_H_R4.all.R1.fastq.gz ../../raw_data/CC_H_R4.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_STARout/Aligned.sortedByCoord.out.bam

# CC_P1

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P1_R2.all.R1.fastq.gz ../../raw_data/CC_P1_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P1_R3.all.R1.fastq.gz ../../raw_data/CC_P1_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P1_R4.all.R1.fastq.gz ../../raw_data/CC_P1_R4.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_STARout/Aligned.sortedByCoord.out.bam

# CC_P2

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P2_R2.all.R1.fastq.gz ../../raw_data/CC_P2_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P2_R3.all.R1.fastq.gz ../../raw_data/CC_P2_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_P2_R4.all.R1.fastq.gz ../../raw_data/CC_P2_R4.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_STARout/Aligned.sortedByCoord.out.bam

# CC_YFB_K

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_K_R1.all.R1.fastq.gz ../../raw_data/CC_YFB_K_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_K_R2.all.R1.fastq.gz ../../raw_data/CC_YFB_K_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_K_R3.all.R1.fastq.gz ../../raw_data/CC_YFB_K_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_YFB_L

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_L_R1.all.R1.fastq.gz ../../raw_data/CC_YFB_L_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_L_R2.all.R1.fastq.gz ../../raw_data/CC_YFB_L_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_L_R3.all.R1.fastq.gz ../../raw_data/CC_YFB_L_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_YFB_T

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_T_R1.all.R1.fastq.gz ../../raw_data/CC_YFB_T_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_T_R2.all.R1.fastq.gz ../../raw_data/CC_YFB_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_YFB_T_R3.all.R1.fastq.gz ../../raw_data/CC_YFB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_FB_KL

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_KL_R1.all.R1.fastq.gz ../../raw_data/CC_FB_KL_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_KL_R2.all.R1.fastq.gz ../../raw_data/CC_FB_KL_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_KL_R3.all.R1.fastq.gz ../../raw_data/CC_FB_KL_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_STARout/Aligned.sortedByCoord.out.bam

# CC_FB_T

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_T_R1.all.R1.fastq.gz ../../raw_data/CC_FB_T_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_T_R2.all.R1.fastq.gz ../../raw_data/CC_FB_T_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_STARout
cd ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_STARout
STAR --runThreadN 48 --genomeDir ../ccinerea_AmutBmut_STARindex --readFilesIn ../../raw_data/CC_FB_T_R3.all.R1.fastq.gz ../../raw_data/CC_FB_T_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_CUFFquant -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/Copci_AmutBmut1_AssemblyScaffolds.fasta -p 48 -L VM,H,P1,P2,YFB_K,YFB_L,YFB_T,FB_KL,FB_T --upper-quartile-norm ./DIR_ccinerea_AmutBmut_expressionanalysis/ccinerea_AmutBmut_STARindex/ccinerea_AmutBmut_AS_annotation.gtf \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_VM_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R3_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_H_R4_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R3_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P1_R4_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R3_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_P2_R4_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_K_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_L_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_YFB_T_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_KL_R3_CUFFquant/abundances.cxb \
./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R1_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R2_CUFFquant/abundances.cxb,./DIR_ccinerea_AmutBmut_expressionanalysis/CC_FB_T_R3_CUFFquant/abundances.cxb
