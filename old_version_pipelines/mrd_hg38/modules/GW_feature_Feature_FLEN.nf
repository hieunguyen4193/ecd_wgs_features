process Feature_FLEN {
    tag "$sampleid"
    label 'no_publish'

    container "tronghieunguyen/ecd_features:latest"
    
    input:
        tuple val(sampleid), path(data)
    output:
        path("${sampleid}_GWfeature_FLEN.csv"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    features_from_GW_data.R \
        --feature_type FLEN \
        --sample $sampleid \
        --input . \
        --output .
    """
}