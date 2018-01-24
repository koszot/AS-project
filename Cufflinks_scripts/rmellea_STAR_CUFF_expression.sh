# STAR INDEXING

mkdir DIR_rmellea_expressionanalysis
cp -r ./rmellea_genome ./DIR_rmellea_expressionanalysis/rmellea_STARindex
cd ./DIR_rmellea_expressionanalysis/rmellea_STARindex
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../rmellea_STARindex/rmellea_AS_annotation.gtf --genomeFastaFiles ./Ricmel1_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# RM_VM

mkdir ./DIR_rmellea_expressionanalysis/RM_VM_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_VM_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/R1_VM_S75_R1_001.fastq.gz ../../raw_data/R1_VM_S75_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_VM_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_VM_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/R2_VM_S76_R1_001.fastq.gz ../../raw_data/R2_VM_S76_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_VM_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_VM_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/R3_VM_S77_R1_001.fastq.gz ../../raw_data/R3_VM_S77_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_VM_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_VM_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_VM_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# RM_P

mkdir ./DIR_rmellea_expressionanalysis/RM_P_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_P_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/KP1_S1_R1_001.fastq.gz ../../raw_data/KP1_S1_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_P_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_P_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/KP3_S3_R1_001.fastq.gz ../../raw_data/KP3_S3_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_P_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_P_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/P2_S5_R1_001.fastq.gz ../../raw_data/P2_S5_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_P_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_P_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_P_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_P_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_P_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_P_R3_STARout/Aligned.sortedByCoord.out.bam

# RM_YFB

mkdir ./DIR_rmellea_expressionanalysis/RM_YFB_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_YFB_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RYFB2_1_S24_R1_001.fastq.gz ../../raw_data/RYFB2_1_S24_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_YFB_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_YFB_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RYFB_4_S25_R1_001.fastq.gz ../../raw_data/RYFB_4_S25_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_YFB_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_YFB_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RYFB_5_S26_R1_001.fastq.gz ../../raw_data/RYFB_5_S26_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_YFB_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_YFB_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_YFB_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_YFB_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_YFB_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_YFB_R3_STARout/Aligned.sortedByCoord.out.bam

# RM_FB_K

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_K_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_K_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_K_2_S7_R1_001.fastq.gz ../../raw_data/RFB_K_2_S7_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_K_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_K_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_K_5_S8_R1_001.fastq.gz ../../raw_data/RFB_K_5_S8_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_K_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_K_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_K_9_S9_R1_001.fastq.gz ../../raw_data/RFB_K_9_S9_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_K_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_K_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_K_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_K_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_K_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_K_R3_STARout/Aligned.sortedByCoord.out.bam

# RM_FB_T

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_T_R1_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_T_R1_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_T_4_10_S12_R1_001.fastq.gz ../../raw_data/RFB_T_4_10_S12_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_T_R2_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_T_R2_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_T_7_S13_R1_001.fastq.gz ../../raw_data/RFB_T_7_S13_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_rmellea_expressionanalysis/RM_FB_T_R3_STARout
cd ./DIR_rmellea_expressionanalysis/RM_FB_T_R3_STARout
STAR --runThreadN 32 --genomeDir ../rmellea_STARindex --readFilesIn ../../raw_data/RFB_T_S10_R1_001.fastq.gz ../../raw_data/RFB_T_S10_R2_001.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_T_R1_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_T_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_T_R2_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_T_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -q -o ./DIR_rmellea_expressionanalysis/RM_FB_T_R3_CUFFquant -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf --max-bundle-frags 1000000 ./DIR_rmellea_expressionanalysis/RM_FB_T_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -q -o CUFFdiff -b ./DIR_rmellea_expressionanalysis/rmellea_STARindex/Ricmel1_AssemblyScaffolds.fasta -p 32 -L VM,P,YFB,FB_K,FB_T --upper-quartile-norm ./DIR_rmellea_expressionanalysis/rmellea_STARindex/rmellea_AS_annotation.gtf \
./DIR_rmellea_expressionanalysis/RM_VM_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_VM_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_VM_R3_CUFFquant/abundances.cxb \
./DIR_rmellea_expressionanalysis/RM_P_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_P_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_P_R3_CUFFquant/abundances.cxb \
./DIR_rmellea_expressionanalysis/RM_YFB_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_YFB_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_YFB_R3_CUFFquant/abundances.cxb \
./DIR_rmellea_expressionanalysis/RM_FB_K_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_FB_K_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_FB_K_R3_CUFFquant/abundances.cxb \
./DIR_rmellea_expressionanalysis/RM_FB_T_R1_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_FB_T_R2_CUFFquant/abundances.cxb,./DIR_rmellea_expressionanalysis/RM_FB_T_R3_CUFFquant/abundances.cxb












