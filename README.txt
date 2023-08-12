This is a workflow to process fastq reads into enrichment scores
Steps
1. Merge reads with BBmerge
# bbmerge.sh in=reads.fq out=merged.fq outu=unmerged.fq ihist=ihist.txt
2. Filter reads 
3. Trim reads
4. Align reads
5. Convert into counts
6. Use MAGECK to conduct statistics

Nextflow will run each module in its own docker container with all dependencies installed
It will have two profiles for running locally or globally


# test case bb merge. Must specify read 2 otherwise it goes interleaved
bbmerge.sh in=test_data/CJ019-C1-8hshort_R1_001.fastq in2=test_data/CJ019-C1-8hshort_R2_001.fastq out=testmerge.fq outu=unmergedtest.fq ihist=ihisttest.txt

making more test files
head -n 20000 /Volumes/Jacob_Drive/78k_ngs_andrew/Reads/R1-1_R1_001.fastq > test_data/R1-1_short_R1_001.fastq
head -n 20000 /Volumes/Jacob_Drive/78k_ngs_andrew/Reads/R1-1_R2_001.fastq > test_data/R1-1_short_R2_001.fastq
head -n 20000 /Volumes/Jacob_Drive/78k_ngs_andrew/Reads/R1-2_R1_001.fastq > test_data/R1-2_short_R1_001.fastq
head -n 20000 /Volumes/Jacob_Drive/78k_ngs_andrew/Reads/R1-2_R2_001.fastq > test_data/R1-2_short_R2_001.fastq