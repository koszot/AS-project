# STAR INDEXING

mkdir DIR_expressionanalysis
cp -r ./genome ./DIR_expressionanalysis/STARindex
cd ./DIR_expressionanalysis/STARindex
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ../STARindex/RRPM_transcripts.fixed.gtf --genomeFastaFiles ./Schco3_AssemblyScaffolds.fasta --sjdbOverhang 100
cd ../..

# STAR ALIGNMENT

# SC_VM

mkdir ./DIR_expressionanalysis/SC_VM_R1_STARout
cd ./DIR_expressionanalysis/SC_VM_R1_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_VM_R1.all.R1.fastq.gz ../../raw_data/SC_VM_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_VM_R2_STARout
cd ./DIR_expressionanalysis/SC_VM_R2_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_VM_R2.all.R1.fastq.gz ../../raw_data/SC_VM_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_VM_R3_STARout
cd ./DIR_expressionanalysis/SC_VM_R3_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_VM_R3.all.R1.fastq.gz ../../raw_data/SC_VM_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# SC_P1

mkdir ./DIR_expressionanalysis/SC_P1_R1_STARout
cd ./DIR_expressionanalysis/SC_P1_R1_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_P1_R1.all.R1.fastq.gz ../../raw_data/SC_P1_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_P1_R2_STARout
cd ./DIR_expressionanalysis/SC_P1_R2_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_P1_R2.all.R1.fastq.gz ../../raw_data/SC_P1_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_P1_R3_STARout
cd ./DIR_expressionanalysis/SC_P1_R3_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_P1_R3.all.R1.fastq.gz ../../raw_data/SC_P1_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# SC_P2

mkdir ./DIR_expressionanalysis/SC_P2_R1_STARout
cd ./DIR_expressionanalysis/SC_P2_R1_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_P2_R1.all.R1.fastq.gz ../../raw_data/SC_P2_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_P2_R2_STARout
cd ./DIR_expressionanalysis/SC_P2_R2_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_P2_R2.all.R1.fastq.gz ../../raw_data/SC_P2_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_P2_R3_STARout
cd ./DIR_expressionanalysis/SC_P2_R3_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_P2_R3.all.R1.fastq.gz ../../raw_data/SC_P2_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# SC_YFB

mkdir ./DIR_expressionanalysis/SC_YFB_R1_STARout
cd ./DIR_expressionanalysis/SC_YFB_R1_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_YFB_R1.all.R1.fastq.gz ../../raw_data/SC_YFB_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_YFB_R2_STARout
cd ./DIR_expressionanalysis/SC_YFB_R2_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_YFB_R2.all.R1.fastq.gz ../../raw_data/SC_YFB_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_YFB_R3_STARout
cd ./DIR_expressionanalysis/SC_YFB_R3_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_YFB_R3.all.R1.fastq.gz ../../raw_data/SC_YFB_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# SC_FB

mkdir ./DIR_expressionanalysis/SC_FB_R1_STARout
cd ./DIR_expressionanalysis/SC_FB_R1_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_FB_R1.all.R1.fastq.gz ../../raw_data/SC_FB_R1.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_FB_R2_STARout
cd ./DIR_expressionanalysis/SC_FB_R2_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_FB_R2.all.R1.fastq.gz ../../raw_data/SC_FB_R2.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

mkdir ./DIR_expressionanalysis/SC_FB_R3_STARout
cd ./DIR_expressionanalysis/SC_FB_R3_STARout
STAR --runThreadN 12 --genomeDir ../STARindex --readFilesIn ../../raw_data/SC_FB_R3.all.R1.fastq.gz ../../raw_data/SC_FB_R3.all.R2.fastq.gz --readFilesCommand zcat --genomeChrBinNbits 12 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --twopassMode Basic --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitBAMsortRAM 80000000000 --alignIntronMax 30000 --alignIntronMin 3
cd ../..

# CUFFQUANT

# SC_VM

cuffquant -o ./DIR_expressionanalysis/SC_VM_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_VM_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_VM_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_VM_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_VM_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_VM_R3_STARout/Aligned.sortedByCoord.out.bam

# SC_P1

cuffquant -o ./DIR_expressionanalysis/SC_P1_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_P1_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_P1_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_P1_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_P1_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_P1_R3_STARout/Aligned.sortedByCoord.out.bam

# SC_P2

cuffquant -o ./DIR_expressionanalysis/SC_P2_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_P2_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_P2_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_P2_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_P2_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_P2_R3_STARout/Aligned.sortedByCoord.out.bam

# SC_YFB

cuffquant -o ./DIR_expressionanalysis/SC_YFB_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_YFB_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_YFB_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_YFB_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_YFB_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_YFB_R3_STARout/Aligned.sortedByCoord.out.bam

# SC_FB

cuffquant -o ./DIR_expressionanalysis/SC_FB_R1_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_FB_R1_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_FB_R2_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_FB_R2_STARout/Aligned.sortedByCoord.out.bam

cuffquant -o ./DIR_expressionanalysis/SC_FB_R3_CUFFquant -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 ./DIR_expressionanalysis/STARindex/RRPM_transcripts.fixed.gtf --max-bundle-frags 1000000 ./DIR_expressionanalysis/SC_FB_R3_STARout/Aligned.sortedByCoord.out.bam

# CUFFDIFF

cuffdiff -o CUFFdiff -b ./DIR_expressionanalysis/STARindex/Schco3_AssemblyScaffolds.fasta -p 12 -L VM,P1,P2,YFB,FB --upper-quartile-norm ./genome/RRPM_transcripts.fixed.gtf \
./DIR_expressionanalysis/SC_VM_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_VM_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_VM_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/SC_P1_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_P1_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_P1_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/SC_P2_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_P2_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_P2_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/SC_YFB_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_YFB_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_YFB_R3_CUFFquant/abundances.cxb \
./DIR_expressionanalysis/SC_FB_R1_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_FB_R2_CUFFquant/abundances.cxb,./DIR_expressionanalysis/SC_FB_R3_CUFFquant/abundances.cxb




