//work in progress to go from DB search -> MSA -> filtering -> pMSA for various methods
process query_first {
//this process just puts the desired query based on params.query to first of fasta list
    input:
        path fasta
    output:
        path "*_rearranged.fasta"
    when:
    def fastaFile = new File("${params.inputdir}/${fasta}")
        def entryCount = 0
        fastaFile.eachLine { line ->
            if (line.startsWith(">")) {
                entryCount++
            }
        }
        return entryCount > 1
    script:
    """
    python3 $projectDir/scripts/rearrange_fasta_query.py "${fasta}" "${params.query}" "${fasta.getBaseName()}_rearranged.fasta" "${fasta.getBaseName()}_query.fasta"
    """
}
process clustalo_MSA {
//starting with default parameters. This accepts an unaligned fasta list and aligns with Clustal-W
//truncates descriptions at 170? characters 
    input: 
        path fasta
    output:
        path "*_clustalo.fasta"
    script:
    """
    clustalo -i ${fasta} -o ${fasta.getBaseName()}_clustalo.fasta --output-order=input-order --outfmt=fasta
    """
}

process clustalo_stats_filter {
//This accepts a clustalo fasta MSA, calculates summary stats in a df, filters based on params
//then outputs the filtered fasta MSAs along with histograms of the stats of the filtered MSAs
    publishDir "${params.outputdir}/stats", pattern: '*_filtered.csv', mode: 'copy'
    input: 
        path fasta_MSA
    output:
        path "*_filtered.fasta"
        path "*_filtered.csv"
    script:
    """
    python3 $projectDir/scripts/clustalo_fasta_stats_filter.py "${fasta_MSA}" ${params.percent_ID} ${params.coverage} ${params.gap_ratio}
    """
}

process hist_csv_stats {
//this inputs a .tsv from a blastp, parses desired columns, and outputs histograms of each column
    publishDir "${params.outputdir}/hist", pattern: '*.png', mode: 'copy'
    input:
        path MSA_csv
        each columns
    output:
        path "*.png"
    script:
    """
    python3 $projectDir/scripts/hist_normpercent.py "${MSA_csv}" "${columns}"
    """
}

process fasta_to_a3m {
    //publishDir 'testoutput', mode: 'copy'
    input:
        path MSA
    output:
        path "*.a3m"
    """
    $projectDir/scripts/reformat.pl -v 0 fas a3m "${MSA}" "${MSA.getBaseName()}.a3m"
    """
}


process pMSA_taxID {
//this pairs alignments based on TaxID (initially trained on OrthoDB databases)
    publishDir "${params.outputdir}/pMSAs", pattern:'*.a3m',  mode: 'copy'
    input:
        path MSA
    output:
        path '*.a3m'
    script:
    """
    python3 $projectDir/scripts/pMSA_taxID.py ${MSA}
    """
}


// %%%%%%%%%%%%%%%%%% WORKFLOW %%%%%%%%%%%%%%%%%% 
params.query = "ASFV Georgia 2007/1(FR682468)" //this is the search term that will be used to set the first query
params.coverage = 0
params.percent_ID = 100
params.gap_ratio = 0
columns = ["%ID", "coverage"]
params.outputdir = 'output/FastaBank_clustalo'
params.inputdir = 'data/FastaBank'

workflow {
    //this sends a list of all paths of .fasta in the given Dir in a list at once
    input = Channel.fromPath("${params.inputdir}/*.fasta")
    rearranged = query_first(input)
    MSAs = clustalo_MSA(rearranged)
    (filtered_MSAs, stats) = clustalo_stats_filter(MSAs)
    images = hist_csv_stats(stats, columns)
    fasta_to_a3m(filtered_MSAs).collect() | pMSA_taxID
}

// %%%%%%%%%%%%%%%%%% WORKFLOW %%%%%%%%%%%%%%%%%% 



//not as useful processes
process BLASTp_filter {
//this process pulls the query (params.query) from a fasta file and makes it the first line. 
//then it builds a BLASTp database from it and searches the query against
//last, it filters the original fasta file and outputs the fasta and the .tsv of the seqeunces
    input: 
        path raw_fasta
    output:
        path '*_filtered.fasta'
        path '*_filtered.tsv'

    script:
    """
    python3 $projectDir/scripts/rearrange_fasta_query.py ${raw_fasta} ${params.query} ${raw_fasta.getBaseName()}_rearranged.fasta ${raw_fasta.getBaseName()}_query.fasta
    makeblastdb -in ${raw_fasta.getBaseName()}_rearranged.fasta -dbtype prot -out ${raw_fasta.getBaseName()}
    blastp -query ${raw_fasta.getBaseName()}_query.fasta -db ${raw_fasta.getBaseName()} -out ${raw_fasta.getBaseName()}.tsv -evalue 1 -max_target_seqs 1000000 -outfmt "7 delim=    std qcovs qcovhsp ppos"
    python3 $projectDir/scripts/blastp-tsv_filter.py ${raw_fasta.getBaseName()}_rearranged.fasta ${raw_fasta.getBaseName()}.tsv "${params.percent_ID}" "${params.coverage}"
    """
}
process pHMMER_MSA {
//this process inputs a single query fasta and a database fasta and does a MSA with pHMMER
//converts the stockholm 1.0 to a fasta output. Haven't looked at how it does naming
    input:
        tuple path(db) path(query)
    output:
        path "*.fasta"
    script:
    """
    phmmer -A ${db.getBaseName}.sto ${query} ${db}
    python3 $projectDir/scripts/stockholm_to_fasta.py ${db.getBaseName}.sto
    """
}
process pMSA {
    publishDir 'output/pMSAs', pattern:'*.fasta',  mode: 'copy'
    input: 
    // list of paths for all .fasta files. Useful for internal Fastabank/MSAbank
        path MSA
    output: 
        path "*.fasta"
    script:
    """
    python3 $projectDir/scripts/pMSA.py ${MSA}
    """
}