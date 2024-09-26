process BBMAP_REPAIR {
    tag "$sampleid"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(sampleid), path(reads)

    output:
    tuple val(sampleid), path('*repaired*.fastq.gz'), emit: reads
    path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sampleid}"
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    repair.sh \\
        -Xmx\$maxmem \\
        in1=${reads[0]} in2=${reads[1]} \\
        out1=${prefix}_repaired_R1.fastq.gz out2=${prefix}_repaired_R2.fastq.gz \\
        outsingle=${prefix}_single.fastq.gz \\
        threads=$task.cpus \\
        $args \\
        &> ${prefix}.bbduk.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}