process FGBIO_GROUPREADSBYUMI {
    tag "$sampleid"
    label 'process_medium'
    label 'no_publish'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:1.4.0--hdfd78af_1' :
        'quay.io/biocontainers/fgbio:1.4.0--hdfd78af_1' }"

    input:
    tuple val(sampleid), path(bam)
    val(strategy)

    output:
    tuple val(sampleid), path("*_umi-grouped.bam")  , emit: bam
    path("*_umi_histogram.txt"), emit: histogram
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--edits=0 --min-map-q=20'
    def prefix = task.ext.prefix ?: "${sampleid}"
    def avail_mem = (task.memory.giga).intValue()
    """

    fgbio \\
        -Djava.io.tmpdir="tmp" -Xmx${avail_mem}g \\
        GroupReadsByUmi \\
        -s $strategy \\
        $args \\
        -i $bam \\
        -o ${prefix}_umi-grouped.bam \\
        -f ${prefix}_umi_histogram.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}