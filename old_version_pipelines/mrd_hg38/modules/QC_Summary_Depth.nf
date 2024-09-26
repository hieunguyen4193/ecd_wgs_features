process Summary_Depth {
    tag "$sampleid"
    label 'process_low'

    container "gene110/samtools_python:v4"
    
    input:
        path(datas)
    output:
        path("Depth.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    summary_depth.py . .
    """
}