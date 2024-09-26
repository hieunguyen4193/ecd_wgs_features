process SUMMARY_ICHORCNA {
    tag "$sampleid"
    label 'process_low'

    container "gene110/samtools_python:v4"
    
    input:
        path(datas)
    output:
        path("IchorCNA.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    get_tumor_fraction.py . .
    """
}