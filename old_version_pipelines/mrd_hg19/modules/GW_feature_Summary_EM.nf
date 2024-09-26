process Summary_EM {
    tag "$sampleid"
    label 'process_low'

    container "gene110/samtools_python:v4"
    
    input:
        path(datas)
        path(example_samples)
    output:
        path("EM.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    summary_feature.py . . EM $example_samples
    """
}