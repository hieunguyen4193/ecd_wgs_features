process Dist_bed {
    tag "$sampleid"
    label 'no_publish'

    container "staphb/samtools:latest"
    
    input:
        tuple val(sampleid), path(data)
    output:
        tuple val(sampleid), path("${sampleid}.dist.bed"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    awk -v OFS='\\t' '{\$4=\$3-\$2; print \$0}' $data > ${sampleid}.dist.bed
    """
}