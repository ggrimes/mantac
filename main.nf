
// set nextflow to use dsl2
nextflow.enable.dsl = 2

// add parameters
params.bam="*.{bam,bai}"
params.config="configManta.py.ini"
params.region="region.bed"
params.reference="hg38.fa"

// add help message
def helpMessage (  ) {
    log.info """
    Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --bam '*.{bam,bai}' --config 'configManta.py.ini' --region 'region.bed' --reference 'hg38.fa'
    Mandatory arguments:
        --bam                          Path to input bam file.  This must be surrounded by quotes.
        --config                       Path to configManta.py.ini file.  This must be surrounded by quotes.
        --region                       Path to region.bed file.  This must be surrounded by quotes.
        --reference                    Path to reference genome.  This must be surrounded by quotes.
    Other options:
        --help                         Shows this help message
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// add log
log.info """\
        Running MANTA SV caller
        bam: ${params.bam}
        config: ${params.config}
        region: ${params.region}
        reference: ${params.reference}
        """.stripIndent()

// write process to run MANTA SV caller
process MANTA {
    publishDir "results/${sample}", mode: 'copy'
    tag "${sample}"
    cpus 16
    memory '32 GB'
    container 'depot.galaxyproject.org/singularity/manta:1.6.0--py27h2d50403_0'
    
    input:
        tuple val(sample), path(bam)
        each path(config)
        each path(region)
    output:
        path "manta/results" into manta_results
    
    script:
    """
    configManta.py \
    --tumorBam ${bam[0]} \
    --reference ${params.reference} \
    --callRegions ${region} \
    --exome \
    --runDir manta \
    --generateEvidenceBam \
    --config ${config}
    
    manta/runWorkflow.py -j ${ task.cpus}
    """
}
// setup channels
channel.fromFilePairs(params.bam,checkIfExists: true)
    .set { bam_files }
channel.fromPath(params.config,checkIfExists: true)
    .set { config_files }    

/*
* Run MANTA SV caller
*/    
workflow {
    MANTA(bam_files, config_files)

}   

workflow.onComplete {
        log.info "Pipeline completed at: $workflow.complete"
        log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}    
