process FGBIO_FILTERCONSENSUSREADS {
    tag "$sampleid"
    label 'process_medium'
    label 'no_publish'  

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:1.4.0--hdfd78af_1' :
        'quay.io/biocontainers/fgbio:1.4.0--hdfd78af_1' }"

    input:
    tuple val(sampleid), path(bam)
    path(fasta)
    path(dict)

    output:
    tuple val(sampleid), path("*.consensus.mapped.filtered.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--min-reads=1 --max-read-error-rate=0.05 --max-base-error-rate=0.1 --min-base-quality=2 --max-no-call-fraction=0.2'
    def prefix = task.ext.prefix ?: "${sampleid}"
    def avail_mem = (task.memory.giga).intValue()
    """
    fgbio \\
        -Djava.io.tmpdir="tmp" -Xmx${avail_mem}g \\
        FilterConsensusReads \\
        --input $bam \\
        $args \\
        --ref $fasta \\
        --output ${prefix}.consensus.mapped.filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}