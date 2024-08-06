// ----- ----- ----- CHANNEL ----- ----- -----
params.input_file= "$params.input/*.tsv"
src=file(params.src)
Channel
    .fromPath( params.input_file )
    .ifEmpty { error "Cannot find any reads matching: ${params.input_file}"  }
    // .view()
    .into {input_ch}

process processing_bam_file_to_image_matrix {
    cache "deep";
    publishDir "$params.output/output", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 20

    input:
        file(input_file) from input_ch
        file src
    output:
        file("*") into output_ch
    script:
    """
    python $src --input $input_file --output .
    """
}