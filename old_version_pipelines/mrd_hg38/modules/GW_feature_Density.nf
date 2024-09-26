process Density {
    tag "$sampleid"
    label 'no_publish'

    container "gene110/samtools_python:v4"
    
    input:
        tuple val(sampleid), path(data)
    output:
        path("${sampleid}.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    count_density.py $data .
    """
}