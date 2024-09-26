process Create_full_bed {
    tag "$sampleid"
    label 'no_publish'

    container "staphb/samtools:latest"
    
    input:
        tuple val(sampleid), path(data)
    output:
        tuple val(sampleid), path("${sampleid}.full.bed"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    awk -v OFS='\\t' '{\$3=\$2+1; print \$0}' $data > ${sampleid}.full.bed
    """
}