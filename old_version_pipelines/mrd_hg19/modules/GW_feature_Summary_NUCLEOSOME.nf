process Summary_NUCLEOSOME {
    tag "$sampleid"
    label 'process_low'

    container "gene110/samtools_python:v4"
    
    input:
        path(datas)
        path(example_samples)
    output:
        path("NUCLEOSOME.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    merge_feature_distribution.py . . $example_samples
    """
}