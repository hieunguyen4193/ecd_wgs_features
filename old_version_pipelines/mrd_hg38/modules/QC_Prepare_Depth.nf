process QC_Prepare_Depth {
    tag "$sampleid"
    label 'no_publish'

    container "gene110/samtools_python:v4"
    
    input:
        tuple val(sampleid), path(bam)
    output:
        path("${sampleid}.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    cal_depth.py $sampleid ${bam} . 4
    """
}