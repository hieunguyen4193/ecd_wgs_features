process FGBIO_FASTQTOBAM {
    tag "$sampleid"
    label 'process_low'
    label 'no_publish'  

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:1.4.0--hdfd78af_1' :
        'quay.io/biocontainers/fgbio:1.4.0--hdfd78af_1' }"

    input:
    tuple val(sampleid), path(reads)

    output:
    tuple val(sampleid), path("*.bam") , emit: bam , optional: true
    tuple val(sampleid), path("*.cram"), emit: cram, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--read-structure=8M92T 8M92T'
    def prefix = task.ext.prefix ?: "${sampleid}.unmapped.withUMI"
    def suffix = task.ext.suffix ?: "bam"
    def sample_name = args.contains("--sample") ? "" : "--sample ${sampleid}"
    def library_name = args.contains("--library") ? "" : "--library ${sampleid}"
    def avail_mem = (task.memory.giga).intValue()
    """

    fgbio \\
        -Djava.io.tmpdir="tmp" -Xmx${avail_mem}g \\
        FastqToBam \\
        ${args} \\
        --input ${reads} \\
        --output ${prefix}.${suffix} \\
        ${sample_name} \\
        ${library_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}