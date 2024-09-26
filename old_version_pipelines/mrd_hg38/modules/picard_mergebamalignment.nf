process PICARD_MERGEBAMALIGNMENT {
    tag "$sampleid"
    label "process_low"
    label 'no_publish'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.5--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.5--hdfd78af_0' }"

    input:
    tuple val(sampleid), file(unmapped), file(mapped)
    path(fasta)
    path(dict)

    output:
    tuple val(sampleid), file("*.merge.bam") , emit: bam
    path('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'TMP_DIR="tmp" SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true'
    def prefix = task.ext.prefix ?: "${sampleid}"
    if (!task.memory) {
        log.warn '[Picard MergeBamAlignment] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072
    """
    picard \\
        -Xmx${avail_mem}M \\
        MergeBamAlignment \\
        ${args} \\
        UNMAPPED=$unmapped \\
        ALIGNED=$mapped \\
        O=${prefix}.merge.bam \\
        R=$fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MergeBamAlignment --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}