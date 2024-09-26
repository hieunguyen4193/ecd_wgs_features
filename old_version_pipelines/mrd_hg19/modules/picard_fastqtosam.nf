process PICARD_FASTQTOSAM {
    tag "$sampleid"
    label 'process_medium'
    label 'no_publish'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.5--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.5--hdfd78af_0' }"

    input:
    tuple val(sampleid), path(reads)

    output:
    tuple val(sampleid), path("*.unmapped.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'TMP_DIR="tmp"'
    def prefix = task.ext.prefix ?: "${sampleid}"
    if (!task.memory) {
        log.warn '[Picard FastqToSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072
    def sample_name = args.contains("--SAMPLE_NAME") || args.contains("-SM") ? "" : "SAMPLE_NAME=${sampleid}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        FastqToSam \\
        ${args} \\
        FASTQ=${reads[0]} FASTQ2=${reads[1]} \\
        ${sample_name} \\
        OUTPUT=${prefix}.unmapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard FastqToSam --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}