# STAR INDEXING

cp -r ./genome ./DIR_expressionanalysis/STARindex
cd ./DIR_expressionanalysis/STARindex
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../STARindex/RRPM_transcripts.fixed.gtf --genomeFastaFiles ./p3_i2_t47428_Arm_ostoy_v2.scaf --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# ARO_RMA

mkdir ./DIR_expressionanalysis/ARO_RMA_R1_STARout
cd ./DIR_expressionanalysis/ARO_RMA_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_RMA_R1.all.R1.fastq.gz ../../raw_data/ARO_RMA_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_RMA_R2_STARout
cd ./DIR_expressionanalysis/ARO_RMA_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_RMA_R2.all.R1.fastq.gz ../../raw_data/ARO_RMA_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_RMA_R3_STARout
cd ./DIR_expressionanalysis/ARO_RMA_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_RMA_R3.all.R1.fastq.gz ../../raw_data/ARO_RMA_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_VM

mkdir ./DIR_expressionanalysis/ARO_VM_R1_STARout
cd ./DIR_expressionanalysis/ARO_VM_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_VM_R1.all.R1.fastq.gz ../../raw_data/ARO_VM_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_VM_R2_STARout
cd ./DIR_expressionanalysis/ARO_VM_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_VM_R2.all.R1.fastq.gz ../../raw_data/ARO_VM_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_VM_R3_STARout
cd ./DIR_expressionanalysis/ARO_VM_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_VM_R3.all.R1.fastq.gz ../../raw_data/ARO_VM_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_P1

mkdir ./DIR_expressionanalysis/ARO_P1_R1_STARout
cd ./DIR_expressionanalysis/ARO_P1_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P1_R1.all.R1.fastq.gz ../../raw_data/ARO_P1_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P1_R2_STARout
cd ./DIR_expressionanalysis/ARO_P1_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P1_R2.all.R1.fastq.gz ../../raw_data/ARO_P1_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P1_R3_STARout
cd ./DIR_expressionanalysis/ARO_P1_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P1_R3.all.R1.fastq.gz ../../raw_data/ARO_P1_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_P2_C

mkdir ./DIR_expressionanalysis/ARO_P2_C_R1_STARout
cd ./DIR_expressionanalysis/ARO_P2_C_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_C_R1.all.R1.fastq.gz ../../raw_data/ARO_P2_C_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P2_C_R2_STARout
cd ./DIR_expressionanalysis/ARO_P2_C_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_C_R2.all.R1.fastq.gz ../../raw_data/ARO_P2_C_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P2_C_R3_STARout
cd ./DIR_expressionanalysis/ARO_P2_C_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_C_R3.all.R1.fastq.gz ../../raw_data/ARO_P2_C_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_P2_S

mkdir ./DIR_expressionanalysis/ARO_P2_S_R1_STARout
cd ./DIR_expressionanalysis/ARO_P2_S_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_S_R1.all.R1.fastq.gz ../../raw_data/ARO_P2_S_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P2_S_R2_STARout
cd ./DIR_expressionanalysis/ARO_P2_S_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_S_R2.all.R1.fastq.gz ../../raw_data/ARO_P2_S_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_P2_S_R3_STARout
cd ./DIR_expressionanalysis/ARO_P2_S_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_P2_S_R3.all.R1.fastq.gz ../../raw_data/ARO_P2_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_YFB_C

mkdir ./DIR_expressionanalysis/ARO_YFB_C_R1_STARout
cd ./DIR_expressionanalysis/ARO_YFB_C_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_C_R1.all.R1.fastq.gz ../../raw_data/ARO_YFB_C_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_YFB_C_R2_STARout
cd ./DIR_expressionanalysis/ARO_YFB_C_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_C_R2.all.R1.fastq.gz ../../raw_data/ARO_YFB_C_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_YFB_C_R3_STARout
cd ./DIR_expressionanalysis/ARO_YFB_C_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_C_R3.all.R1.fastq.gz ../../raw_data/ARO_YFB_C_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_YFB_S

mkdir ./DIR_expressionanalysis/ARO_YFB_S_R1_STARout
cd ./DIR_expressionanalysis/ARO_YFB_S_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_S_R1.all.R1.fastq.gz ../../raw_data/ARO_YFB_S_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_YFB_S_R2_STARout
cd ./DIR_expressionanalysis/ARO_YFB_S_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_S_R2.all.R1.fastq.gz ../../raw_data/ARO_YFB_S_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_YFB_S_R3_STARout
cd ./DIR_expressionanalysis/ARO_YFB_S_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_YFB_S_R3.all.R1.fastq.gz ../../raw_data/ARO_YFB_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_FB_C

mkdir ./DIR_expressionanalysis/ARO_FB_C_R1_STARout
cd ./DIR_expressionanalysis/ARO_FB_C_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_C_R1.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_C_R2_STARout
cd ./DIR_expressionanalysis/ARO_FB_C_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_C_R2.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_C_R3_STARout
cd ./DIR_expressionanalysis/ARO_FB_C_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_C_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_C_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_FB_L

mkdir ./DIR_expressionanalysis/ARO_FB_L_R1_STARout
cd ./DIR_expressionanalysis/ARO_FB_L_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_L_R1.all.R1.fastq.gz ../../raw_data/ARO_FB_L_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_L_R2_STARout
cd ./DIR_expressionanalysis/ARO_FB_L_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_L_R2.all.R1.fastq.gz ../../raw_data/ARO_FB_L_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_L_R3_STARout
cd ./DIR_expressionanalysis/ARO_FB_L_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_L_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_L_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# ARO_FB_S

mkdir ./DIR_expressionanalysis/ARO_FB_S_R1_STARout
cd ./DIR_expressionanalysis/ARO_FB_S_R1_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_S_R1.all.R1.fastq.gz ../../raw_data/ARO_FB_S_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_S_R2_STARout
cd ./DIR_expressionanalysis/ARO_FB_S_R2_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_S_R2.all.R1.fastq.gz ../../raw_data/ARO_FB_S_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/ARO_FB_S_R3_STARout
cd ./DIR_expressionanalysis/ARO_FB_S_R3_STARout
STAR --runThreadN 16 --genomeDir ../STARindex --readFilesIn ../../raw_data/ARO_FB_S_R3.all.R1.fastq.gz ../../raw_data/ARO_FB_S_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# CUFFQUANT

# ARO_RMA

cuffquant -o ./DIR_expressionanalysis/ARO_RMA_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_RMA_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_RMA_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_RMA_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_RMA_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_RMA_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_VM

cuffquant -o ./DIR_expressionanalysis/ARO_VM_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_VM_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_VM_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_P1

cuffquant -o ./DIR_expressionanalysis/ARO_P1_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P1_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P1_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P1_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P1_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P1_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_P2_C

cuffquant -o ./DIR_expressionanalysis/ARO_P2_C_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_C_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P2_C_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_C_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P2_C_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_C_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_P2_S

cuffquant -o ./DIR_expressionanalysis/ARO_P2_S_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_S_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P2_S_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_S_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_P2_S_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_P2_S_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_YFB_C

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_C_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_C_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_C_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_C_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_C_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_C_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_YFB_S

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_S_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_S_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_S_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_S_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_YFB_S_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_YFB_S_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_FB_C

cuffquant -o ./DIR_expressionanalysis/ARO_FB_C_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_C_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_C_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_C_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_C_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_C_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_FB_L

cuffquant -o ./DIR_expressionanalysis/ARO_FB_L_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_L_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_L_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_L_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_L_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_L_R3_STARout/Aligned.sortedByCoord.out.bam

# ARO_FB_S

cuffquant -o ./DIR_expressionanalysis/ARO_FB_S_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_S_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_S_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_S_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/ARO_FB_S_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/p3_i2_t47428_Arm_ostoy_v2.scaf -p 16 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/ARO_FB_S_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -o CUFFdiff -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 96 -L RMA,VM,P1,P2_C,P2_S,YFB_C,YFB_S,FB_C,FB_L,FB_S --upper-quartile-norm ./DIR_RRPM/RRPM/RRPM_transcripts.fixed.gtf \
./DIR_expressionanalysis/ARO_RMA_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_RMA_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_RMA_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_VM_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_VM_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_VM_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_P1_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P1_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P1_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_P2_C_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P2_C_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P2_C_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_P2_S_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P2_S_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_P2_S_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_YFB_C_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_YFB_C_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_YFB_C_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_YFB_S_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_YFB_S_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_YFB_S_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_FB_C_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_C_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_C_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_FB_L_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_L_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_L_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/ARO_FB_S_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_S_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/ARO_FB_S_R3_CUFFquant/abundances.cxb
