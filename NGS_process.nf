process bbmerge {
    publishDir "${params.publish_dir}", pattern: '*.txt', mode: 'copy'
    //publishDir "${params.publish_dir}", pattern: '*_merged.fastq', mode: 'symlink'
    input:
        tuple val(sampleID), file(read)
    output:
        path '*_merged.fastq'
        path '*_hist.txt'
    script:
    """
    bbmerge.sh in=${read[0]} in2=${read[1]} out=${sampleID}_merged.fastq ihist=${sampleID}_hist.txt trimnonoverlapping=t 
    """
}

process exact_seed {
    //This inputs a list of merged fastq file paths and calculates the number of exact matching seed reads
    //outputs a counts per seed .csv and a filtered reads stats .csv file 
    publishDir "${params.publish_dir}", pattern: '*-counts.tsv', mode: 'copy'
    publishDir "${params.publish_dir}", pattern: '*-filread_stats.tsv', mode: 'copy'
    input: 
        path merged_fastq
    output:
        path "*-counts.tsv"
        path "*-read_stats.tsv"
    script:
    """
    python3 $projectDir/scripts/exact_seed_reads_df.py ${params.seed_db} ${params.exp_tag}-counts.csv ${params.exp_tag}-filread_stats.csv ${merged_fastq}
    """
}

params.input_reads_dir = "/Volumes/Extreme_SSD/78K_CRISPRi/Expanded"
// params.input_design_matricies = "/Users/jacobfenster/Documents/NGS_gRNA/design_matricies"
params.publish_dir = "/Users/jacobfenster/Documents/NGS_gRNA/output/20231126_exact_counts"
params.seed_db = "/Users/jacobfenster/Documents/NGS_gRNA/databases/Lib3_dCas9_78K_seeds.csv"
params.exp_tag = "20231126_carbon_sources"

workflow {
    raw_fastq = Channel.fromFilePairs("${params.input_reads_dir}/*_R{1,2}_*.fastq")
    (merged_fastq, hist) = bbmerge(raw_fastq)
    (reads, fil_readstats) = exact_seed(merged_fastq.collect())
}