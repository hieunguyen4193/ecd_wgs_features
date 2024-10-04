// ----- ----- ----- CHANNEL ----- ----- -----
params.input_file= "$params.input/*.tsv"
hg19=file(params.hg19)
nucleosome_ref=file(params.nucleosome_ref)
src=file(params.src)

Channel
    .fromPath( params.input_file )
    .ifEmpty { error "Cannot find any reads matching: ${params.input_file}"  }
    .view()
    .set {input_ch}

process process_02_generate_EM_FLEN_NUC_features {  
    cache "deep";
    publishDir "$params.output/02_EM_FLEN_NUC_features", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 25

    input:
        file(input_file) from input_ch
        file src
    output:
        file("*") into output_ch
    script:
    """
    bash $src -i ${input_file} -o . -f ${hg19} -r ${nucleosome_ref} -c ${params.clean_up}
    """
}

// nextflow run run_pipeline_FinaleDB.nf \
// --input_file /mnt/NAS_PROJECT/vol_ECDteam/hieunho/data/finaledb_extract \
// --output /datassd/hieunguyen/2024/outdir/FinaleDB \
// --hg19 /datassd/hieunguyen/2024/resources/hg19.fa \
// --nucleosome_ref /datassd/hieunguyen/2024/resources/rpr_map_EXP0779.sorted.bed -resume

