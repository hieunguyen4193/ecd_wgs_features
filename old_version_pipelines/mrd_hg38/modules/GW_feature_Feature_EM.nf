process Feature_EM {
    tag "$sampleid"
    label 'process_low'
    label 'no_publish'

    container "tronghieunguyen/ecd_features:latest"
    
    input:
        tuple val(sampleid), path(data)
    output:
        path("${sampleid}_GWfeature_EM.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    features_from_GW_data.R \
        --feature_type EM \
        --sample $sampleid \
        --input . \
        --output .
    """
}