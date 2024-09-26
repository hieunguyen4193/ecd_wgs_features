process QC_Prepare_Mpileup {
    tag "$sampleid"
    label 'no_publish'

    container "staphb/samtools:latest"
    
    input:
        tuple val(sampleid), path(bam)
    output:
        path("${sampleid}.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    a=\$(samtools mpileup ${bam} | awk -v X="1" '\$4>=X' | wc -l)
    echo "$sampleid,\$a"  > ${sampleid}.csv
    """
}