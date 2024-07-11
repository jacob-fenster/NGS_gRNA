//This is to search the NCBI collection of all viral genomes with jackhmmer, filter, and assemble into MSAs
// input proteome must be in format of >
process split_proteome {
    input:
        path proteome
    output:
        path "*.fasta"
    """
    python3 $projectDir/scripts/split_proteome.py ${proteome}
    """
}

process jackhmmer_viral_proteome {
    publishDir "${params.outputdir}/jhmmer_sto", pattern:'*.sto',  mode: 'copy'
    input:
        path protein
    output:
        path "*.a3m"
        path "*.sto"
    script:
    """
    jackhmmer -N 3 --cpu ${params.cpus} --notextw -A ${protein.getBaseName()}_jhmmer.sto ${protein} ${params.seqdb}
    python3 $projectDir/scripts/jackhmmer_MSA.py ${protein.getBaseName()}_jhmmer.sto ${params.gap_ratio} ${protein.getBaseName()}_jhmmer.fasta
    $projectDir/scripts/reformat.pl -v 0 fas a3m ${protein.getBaseName()}_jhmmer.fasta ${protein.getBaseName()}.a3m
    """
}

process reuse_jackhmmer_sto {
    //this is a workaround to avoid re-running jackhmmer
    input:
        path jhmmer_sto
    output:
        path "*.a3m"
    script:
    """
    python3 $projectDir/scripts/jackhmmer_MSA.py ${jhmmer_sto} ${params.gap_ratio} ${jhmmer_sto.getBaseName()}.fasta
    $projectDir/scripts/reformat.pl -v 0 fas a3m ${jhmmer_sto.getBaseName()}.fasta ${jhmmer_sto.getBaseName()[0..-8]}.a3m
    """
}

process pMSA {
//this pairs alignments based on TaxID (initially trained on OrthoDB databases)
// pMSA_taxID has been updated to pair on genome NCBI fasta_cds_aa from nuccore
    //publishDir "${params.outputdir}/pMSAs", pattern:'*.a3m',  mode: 'copy'
    input:
        path MSA
    output:
        path '*.a3m'
    script:
    """
    python3 $projectDir/scripts/pMSA.py ${MSA}
    """
}

process hhfilter_MSA {
    //filters using hhfilter based on 100% id. outputs .a3m and .fasta format (for downstream stats)
    // have not containerized hh-suite yet
    publishDir "${params.outputdir}/pMSAs", pattern:'*fil.a3m',  mode: 'copy'
    //publishDir "${params.outputdir}/pMSAs_temp", pattern:'*fil.fasta',  mode: 'copy'
    //publishDir "${params.outputdir}/pMSAs_temp", pattern:'*_fil80.a3m',  mode: 'copy'
    input:
        path pMSA
    output:
        path '*fil.a3m'
        path '*fil.fasta'
        path '*_fil80.a3m'
    script:
    """
    export PATH=$PATH:/Users/jacobfenster/Documents/ASFV_Postdoc/hhsuite/bin
    hhfilter -v 0 -id ${params.fil_id} -i ${pMSA} -o ${pMSA.getBaseName()}fil.a3m
    hhfilter -v 0 -id 80 -i ${pMSA.getBaseName()}fil.a3m -o ${pMSA.getBaseName()}_fil80.a3m
    $projectDir/scripts/reformat.pl a3m fas ${pMSA.getBaseName()}fil.a3m ${pMSA.getBaseName()}fil.fasta
    """
}

process filter_MSA_exact {
    //filters using hhfilter based on 100% id. outputs .a3m and .fasta format (for downstream stats)
    // have not containerized hh-suite yet
    publishDir "${params.outputdir}/pMSAs", pattern:'*filexact.a3m',  mode: 'copy'
    publishDir "${params.outputdir}/pMSAs_temp", pattern:'*filexact.fasta',  mode: 'copy'
    publishDir "${params.outputdir}/pMSAs_temp", pattern:'*_fil80.a3m',  mode: 'copy'
    input:
        path pMSA
    output:
        path '*filexact.a3m'
        path '*filexact.fasta'
        path '*_fil80.a3m'
    script:
    """
    export PATH=$PATH:/Users/jacobfenster/Documents/ASFV_Postdoc/hhsuite/bin
    python3 $projectDir/scripts/MSA_cluster_exact_match.py ${pMSA} ${pMSA.getBaseName()}filexact.a3m
    hhfilter -v 0 -id 80 -i ${pMSA} -o ${pMSA.getBaseName()}_fil80.a3m
    $projectDir/scripts/reformat.pl a3m fas ${pMSA.getBaseName()}filexact.a3m ${pMSA.getBaseName()}filexact.fasta
    """
}

process pMSA_stats {
    publishDir "${params.outputdir}/stats/figures", pattern: '*.png',  mode: 'move'
    publishDir "${params.outputdir}/stats", pattern: '*summary.csv',  mode: 'copy'
    input:
        path pMSA
        path pMSA_80_a3m
    output:
        path '*summary.csv'
        path '*.png'
    script:
    """
    python3 $projectDir/scripts/hhfilter_cluster_workaround.py temp_Nseq.csv ${pMSA_80_a3m}
    python3 $projectDir/scripts/fasta_MSA_stats_hhfilter_workaround.py ${params.exp_tag} temp_Nseq.csv ${pMSA}
    """

}

//parameters
params.seqdb = "/Users/jacobfenster/Local_Documents/Databases/proteins_of_all_viral_genomes/CompleteGenome_virus_CDS_NCBI_v2_cleanup.fasta"
params.outputdir = "/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/clustered_at_99"
params.proteome = "/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/data/jackhmmer_test_proteome/Georgia-2007-LR743116_trans-rename.fasta"
params.cpus = 10 //for jackhmmer searching
params.gap_ratio = 0.5 //alignments below this will be removed from MSAs
params.exp_tag = "ASFV_Georgia2007"
params.sto_dir = "/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/jhmmer_sto"
params.fil_id = 99
params.workflow='reuse_jh_hhfilter'

//This is a combined workflow script that can be customized for runnning jackhmmer or inputting previously run jackhmmer .sto files
//You can also cluster pMSAs at 100% with custom python script or any %ID by setting params.fil_id
workflow reuse_jh_hhfilter {
    //This is to reuse existing .sto files from jackhmmer and filter pMSAs at params.fil_id %ID
    jhmmer_sto = Channel.fromPath("${params.sto_dir}/*.sto")
    MSAs = reuse_jackhmmer_sto(jhmmer_sto)
    pMSAs = pMSA(MSAs.collect()).flatten()
    (pMSAs_hhfil_a3m, pMSAs_hhfil_fasta, pMSAs_80_a3m) = hhfilter_MSA(pMSAs)
    //pMSA_stats(pMSAs_hhfil_fasta.collect(), pMSAs_80_a3m.collect())

}

workflow reuse_jh_filexact {
    //This is to reuse existing .sto files from jackhmmer and filter pMSAs at 100% id 
    jhmmer_sto = Channel.fromPath("${params.sto_dir}/*.sto")
    MSAs = reuse_jackhmmer_sto(jhmmer_sto)
    pMSAs = pMSA(MSAs.collect()).flatten()
    (pMSAs_filexact_a3m, pMSAs_filexact_fasta, pMSAs_80_a3m) = filter_MSA_exact(pMSAs)
    pMSA_stats(pMSAs_filexact_fasta.collect(), pMSAs_80_a3m.collect())
} 


workflow jh_hhfilter {
    //this is to generate new jackhmmer .sto files by searching each protein in params.proteome against params.seqdb
    proteome = Channel.fromPath(params.proteome)
    proteins = split_proteome(proteome).flatten()
    (jhmmer_sto, MSAs) = jackhmmer_viral_proteome(proteins)
    pMSAs = pMSA(MSAs.collect()).flatten()
    (pMSAs_hhfil_a3m, pMSAs_hhfil_fasta, pMSAs_80_a3m) = hhfilter_MSA(pMSAs)
    pMSA_stats(pMSAs_hhfil_fasta.collect(), pMSAs_80_a3m.collect())
}

workflow jh_exact {
    //this is to generate new jackhmmer .sto files by searching each protein in params.proteome against params.seqdb
    proteome = Channel.fromPath(params.proteome)
    proteins = split_proteome(proteome).flatten()
    (jhmmer_sto, MSAs) = jackhmmer_viral_proteome(proteins)
    pMSAs = pMSA(MSAs.collect()).flatten()
    (pMSAs_filexact_a3m, pMSAs_filexact_fasta, pMSAs_80_a3m) = filter_MSA_exact(pMSAs)
    pMSA_stats(pMSAs_filexact_fasta.collect(), pMSAs_80_a3m.collect())
}
workflow {
    if ( params.workflow == 'reuse_jh_filexact' )
        reuse_jh_filexact()
    if ( params.workflow == 'reuse_jh_hhfilter' )
        reuse_jh_hhfilter()
    if ( params.workflow == 'jh_hhfilter' )
        jh_hhfilter()
    if (params.workflow == 'jh_exact' )
        jh_exact()
}