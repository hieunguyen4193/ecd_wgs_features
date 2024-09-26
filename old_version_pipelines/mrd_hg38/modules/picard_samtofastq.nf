process PICARD_SAMTOFASTQ {
    tag "$sampleid"
    label "process_single"
    label 'no_publish'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.5--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.5--hdfd78af_0' }"

    input:
    tuple val(sampleid), path(bam)  
    
    output:
    tuple val(sampleid), path("*.fastq.gz") , emit: fastq
    path("versions.yml")                       , emit: versions

    when:
    task.ext.when == null || task.ext.when    

    script:
    def args = task.ext.args ?: 'TMP_DIR="tmp" INTERLEAVE=true'
    def prefix = task.ext.prefix ?: "${sampleid}.noUMI"
    if (!task.memory) {
        log.warn '[Picard SamToFastq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072
    """
    picard \\
        -Xmx${avail_mem}M \\
        SamToFastq \\
        ${args} \\
        I=$bam \\
        F=${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SamToFastq --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}