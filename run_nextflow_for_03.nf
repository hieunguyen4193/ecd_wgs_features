// ----- ----- ----- CHANNEL ----- ----- -----
params.input_file= "$params.input/*.tsv"
params.feature_version= ""
params.old_nuc= ""
params.generate_feature=""
src=file(params.src)

motif_order=file(params.motif_order)
Channel
    .fromPath( params.input_file )
    .ifEmpty { error "Cannot find any reads matching: ${params.input_file}"  }
    .view()
    .set {input_ch}

process processing_bam_file_to_image_matrix {
    cache "deep";
    publishDir "$params.output/output", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 3

    input:
        file(input_file) from input_ch
        file src
        file motif_order
    output:
        file("*") into output_ch
    script:
    """
    python $src --input $input_file --output . --motif_order $motif_order --feature_version $params.feature_version --old_nuc $params.old_nuc --generate_feature $params.generate_feature
    """
}

// example command:
//  nextflow run run_nextflow_for_03.nf \
//  --input /mnt/archiving/DATA_HIEUNHO/ecd_wgs_features_trong_hieu/ready-to-use \
//  --output /datassd/hieunguyen/2024/outdir/ecd_wgs_features/images \
//  --src /datassd/hieunguyen/2024/src/ecd_wgs_features/03_generate_WGS_features.py \
//  --motif_order /datassd/hieunguyen/2024/src/ecd_wgs_features/motif_order.csv \
//  --feature_version "old" --generate_feature "image_only" --old_nuc "none" -resume