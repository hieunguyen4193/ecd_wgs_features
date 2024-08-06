// ----- ----- ----- CHANNEL ----- ----- -----
params.input_file= "$params.input/*.tsv"

Channel
    .fromPath( params.input_file )
    .ifEmpty { error "Cannot find any reads matching: ${params.INPUT_PAIRS}"  }
    // .view()
    .into {input_ch}

process processing_bam_file_to_image_matrix {
    cache "deep"; tag "$sample_id"
    publishDir "$params.output/output", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 10

    input:
        file(input_file) from input_ch
    output:
        file("*") into output_ch
    script:
    """
    python generate_image_matrix --input $input_file --output .
    """
}