process FGBIO_CALLMOLECULARCONSENSUSREADS {
    tag "$sampleid"
    label 'process_medium'
    label 'no_publish'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:1.4.0--hdfd78af_1' :
        'quay.io/biocontainers/fgbio:1.4.0--hdfd78af_1' }"

    input:
    tuple val(sampleid), path(bam)

    output:
    tuple val(sampleid), path("*.consensus.unmapped.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--error-rate-pre-umi=45 --error-rate-post-umi=30 --min-consensus-base-quality=2 --min-input-base-quality=10 --min-reads=1'
    def prefix = task.ext.prefix ?: "${sampleid}"
    def avail_mem = (task.memory.giga).intValue()
    """
    fgbio \\
        -Djava.io.tmpdir="tmp" -Xmx${avail_mem}g \\
        CallMolecularConsensusReads \\
        --input $bam \\
        $args \\
        --output ${prefix}.consensus.unmapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}