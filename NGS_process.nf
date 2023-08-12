
process mergefastq {
    // this is to make copies of this process's files for validation
    publishDir 'output', mode: 'copy'
    input:
    //the Channel.fromFilePairs sends in tuples of [sampleID, [file[0], file[1]]
        tuple val(sampleID), file(read)
    output:
        path '*_merged.fq'
    script:
    """
    bbmerge.sh in=${read[0]} in2=${read[1]} out=${sampleID}_merged.fq
    """
}

workflow {
    //this makes a channel of fastq file pairs with R1 and R2 identifiers
    //and sends them in as tuples of [sampleID, [file[0], file[1]]
    Channel.fromFilePairs('test_data/*_R{1,2}_*.fastq') | mergefastq
}