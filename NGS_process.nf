params.gRNA_database = "$baseDir/databases/Lib3_dCas9/Lib3_dCas9"

process Mergefastq {
    // Final merge process, trained with test_merge
    publishDir 'output', pattern: '*_hist.txt', mode: 'copy'
    input:
    //the Channel.fromFilePairs sends in tuples of [sampleID, [file[0], file[1]]
        tuple val(sampleID), file(read)
        each mode
    output:
        path '*_merged.fq'
        path '*_hist.txt'
    script:
    """
    bbmerge.sh in=${read[0]} in2=${read[1]} out=${sampleID}_${mode.replaceAll(/=/, "_")}_merged.fq ihist=${sampleID}_${mode.replaceAll(/=/, "_")}_hist.txt trimnonoverlapping=t ${mode}
    """
}


process Fastq_to_Fasta {
    // converts Fastq files to Fasta
    input:
        path fastq
    output:
        path "*.fasta"
    
    script:
    """
    awk 'NR%4==1 {printf ">%s\\n", substr(\$0, 2)} NR%4==2 {print}' "${fastq}" > "${fastq.getBaseName()}.fasta"
    """
}

process Blastn {
    input:
        path fasta
    output:
        path "*.txt"
    script:
    """
    blastn -query ${fasta} -db ${params.gRNA_database} -strand both -task blastn -max_target_seqs 1 -outfmt 6 -out ${fasta.getBaseName()}.txt
    """
}

workflow {
    //this makes a channel of fastq file pairs with R1 and R2 identifiers
    //and sends them in as tuples of [sampleID, [file[0], file[1]]
    bbmerge_strictness = ['strict=t', 'loose=t', '']
    
    fastq_input = Channel.fromFilePairs('test_data/*_R{1,2}_*.fastq')
    (merged_strictness, hist_strict) = Mergefastq(fastq_input, bbmerge_strictness)
   
    align_strict = Fastq_to_Fasta(merged_strictness) | Blastn
    
    //if we don't flatten this will give a list for each merged read
    //this can now be sent into blast alignments

}
workflow {
    bbmerge_errors = ['efilter=20', 'efilter=2']
    fastq_input = Channel.fromFilePairs('test_data/*_R{1,2}_*.fastq')
    (merged_errors, hist_errors) = Mergefastq(fastq_input, bbmerge_errors)
    align_errors = Fastq_to_Fasta(merged_errors) | Blastn
}