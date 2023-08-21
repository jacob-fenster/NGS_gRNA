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

testing out bbmerge parameters
docker_containers/bbmerge/software/bbmap/bbmerge.sh trimnonoverlapping=t in=output/R1-1_reallyshort_R1.fastq in2=output/R1-1_reallyshort_R2.fastq out=output/R1-1_merged.fq outu1=output/unmerg1.fq outu2=output/unmerg2.f1q
testing allowed overlap errors
docker_containers/bbmerge/software/bbmap/bbmerge.sh trimnonoverlapping=t in=test_data/R1-1_short_R1_001.fastq in2=test_data/R1-1_short_R2_001.fastq out=output/R1-1merged.fq efilter=6 ihist=hist_efil_6.txt 
THis is a test merge that can be more easily made with each
docker_containers/bbmerge/software/bbmap/bbmerge.sh trimnonoverlapping=t in=/Volumes/Jacob_Drive/78k_ngs_andrew/Reads/R1-1_R1_001.fastq in2=/Volumes/Jacob_Drive/78k_ngs_andrew/Reads/R1-1_R2_001.fastq ihist=R1-1_hist.txt 


process test_merge {
    publishDir 'output', pattern: '*.txt', mode: 'copy'
    input:
    //the Channel.fromFilePairs sends in tuples of [sampleID, [file[0], file[1]]
        tuple val(sampleID), file(read)
    //this will output two channels in order with a list emitted
    output:
        path '*_merged.fq'
        path '*.txt'
    script:
    // here I could have also used 'each' in inputs to iterate over parameters
    """
    bbmerge.sh in=${read[0]} in2=${read[1]} out=${sampleID}_strict_merged.fq ihist=${sampleID}_strict.txt trimnonoverlapping=t strict=t
    bbmerge.sh in=${read[0]} in2=${read[1]} out=${sampleID}_std_merged.fq ihist=${sampleID}_std.txt trimnonoverlapping=t 
    bbmerge.sh in=${read[0]} in2=${read[1]} out=${sampleID}_vloose_merged.fq ihist=${sampleID}_vloose.txt trimnonoverlapping=t veryloose=t
    bbmerge.sh in=${read[0]} in2=${read[1]} out=${sampleID}_mloose_merged.fq ihist=${sampleID}_mloose.txt trimnonoverlapping=t maxloose=t
    """
}


