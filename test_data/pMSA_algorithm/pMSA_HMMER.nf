
process HMMER_MSA {
    input: 
        path fasta
    output:
        path '*-hmmalign.a3m'
    when:  //only allow this process to run with fastas with more than one entry
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
    phmmer --notextw -A ${fasta.getBaseName()}-phmmer.sto ${fasta.getBaseName()}_query.fasta ${fasta.getBaseName()}_rearranged.fasta
    python3 $projectDir/scripts/phmmer_seed.py ${fasta.getBaseName()}-phmmer.sto ${fasta.getBaseName()}_rearranged.fasta ${fasta.getBaseName()}-phmmerseed.sto
    hmmbuild ${fasta.getBaseName()}.hmm ${fasta.getBaseName()}-phmmerseed.sto > ${fasta.getBaseName()}-hmmbuild.out
    hmmalign -o ${fasta.getBaseName()}-hmmalign.sto ${fasta.getBaseName()}.hmm ${fasta.getBaseName()}_rearranged.fasta
    python3 $projectDir/scripts/hmmalign_cleanup.py ${fasta.getBaseName()}-hmmalign.sto 0.5 ${fasta.getBaseName()}-hmmalign.fasta
    $projectDir/scripts/reformat.pl -v 0 fas a3m "${fasta.getBaseName()}-hmmalign.fasta" "${fasta.getBaseName()}-hmmalign.a3m"
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


params.query = "ASFV Georgia 2007/1(FR682468)" //this is the search term that will be used to set the first query
params.outputdir = 'output/hmmer-Euk-complex'
params.inputdir = 'data/superkingdom'

workflow {
OG_fastas = Channel.fromPath("${params.inputdir}/*.fasta")
MSAs = HMMER_MSA(OG_fastas)
MSAs.collect() | pMSA_taxID
}