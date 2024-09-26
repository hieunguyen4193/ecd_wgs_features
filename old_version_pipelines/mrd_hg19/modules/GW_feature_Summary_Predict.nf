process Summary_Predict {
    tag "$sampleid"
    label 'process_medium'

    container "gene110/samtools_python:v5"
    
    input:
        path(nucleosome)
        path(em)
        path(flen)
        path(ichorcna)
        path(model)
        path(example_samples)
    output:
        path("output.xlsx"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    predict.py . $model . $example_samples
    """
}