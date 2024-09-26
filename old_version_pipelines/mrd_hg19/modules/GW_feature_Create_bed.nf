process Create_bed {
    tag "$sampleid"
    label 'no_publish'

    container "staphb/samtools:latest"
    
    input:
        tuple val(sampleid), path(bam), path(bai)
    output:
        tuple val(sampleid), path("${sampleid}.bed"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    samtools view -@ 8 ${bam} | awk '{if (\$9 > 0){print \$0}}'  | cut -f3,4 > ${sampleid}.bed
    """
}