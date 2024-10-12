// ----- ----- ----- CHANNEL ----- ----- -----
params.input_file= "$params.input/*.tsv"
hg19=file(params.hg19)
nucleosome_ref=file(params.nucleosome_ref)
prep_src=file(params.prep_src)
src=file(params.src)
convert_bed=file(params.src_convert_bed)

Channel
    .fromPath( params.input_file )
    .ifEmpty { error "Cannot find any reads matching: ${params.input_file}"  }
    // .view()
    .set {input_ch}

process preprocess_input {  
    cache "deep";
    publishDir "$params.output/01_preprocess_FinaleDB_files", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 25

    input:
        file(input_file) from input_ch
        file prep_src
    output:
        file("*.tsv") into prep_ch
    script:
    """
    bash 
    bash $prep_src -i ${input_file} -o . 
    """
}

process process_02_generate_EM_FLEN_NUC_features {  
    cache "deep";
    publishDir "$params.output/02_EM_FLEN_NUC_features", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 25

    input:
        file(input_file) from prep_ch
        file src
        file convert_bed
    output:
        file("*_prep.frag.tsv") into output_ch
    script:
    """
    bash $src -i ${input_file} -o . -f ${hg19} -r ${nucleosome_ref} -c ${params.clean_up}
    """
}

// nextflow run run_pipeline_FinaleDB.nf \
// --input /mnt/NAS_PROJECT/vol_ECDteam/hieunho/data/finaledb_extract \
// --output /datassd/DATA_HIEUNGUYEN/2024/outdir/FinaleDB \
// --hg19 /datassd/DATA_HIEUNGUYEN/2024/resources/hg19.fa \
// --nucleosome_ref /datassd/DATA_HIEUNGUYEN/2024/resources/rpr_map_EXP0779.sorted.bed \
// --src /datassd/DATA_HIEUNGUYEN/2024/src/ecd_wgs_features/02_calculate_EM_FLEN_NUC_from_FRAG.sh \
// --clean_up true \
// --src_convert_bed /datassd/DATA_HIEUNGUYEN/2024/src/ecd_wgs_features/convert_full_bed_nucleosome.py \
// --prep_src /datassd/DATA_HIEUNGUYEN/2024/src/ecd_wgs_features/preprocess_FinaleDB_input.sh
// -resume

