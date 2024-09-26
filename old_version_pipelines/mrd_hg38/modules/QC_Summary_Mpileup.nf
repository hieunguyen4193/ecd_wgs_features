process Summary_Mpileup {
    tag "$sampleid"
    label 'process_low'

    container "gene110/samtools_python:v4"
    
    input:
        path(datas)
    output:
        path("Mpileup.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    summary_mpileup.py . .
    """
}