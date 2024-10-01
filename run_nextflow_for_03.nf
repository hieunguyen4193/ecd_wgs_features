// ----- ----- ----- CHANNEL ----- ----- -----
params.input_file= "$params.input/*.tsv"
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
    maxForks 10

    input:
        file(input_file) from input_ch
        file src
        file motif_order
    output:
        file("*") into output_ch
    script:
    """
    python $src --input $input_file --output . --motif_order $motif_order
    """
}