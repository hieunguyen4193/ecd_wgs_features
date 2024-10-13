// ----- ----- ----- CHANNEL ----- ----- -----
params.SampleSheet=""
hg19=file(params.hg19)
nucleosome_ref=file(params.nucleosome_ref)
prep_src=file(params.prep_src)
src=file(params.src)
convert_bed=file(params.src_convert_bed)
params.maxForks=""
Channel
    .fromPath( params.SampleSheet )
    .splitCsv( header:true )
    .map { row -> tuple(row.SampleID,  file(row.path))}  
    .view()
    .set { Input_ch }


process process_02_generate_EM_FLEN_NUC_features {  
    cache "deep";
    publishDir "$params.output/FinaleDB_features", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks "${params.maxForks}"

    input:
        tuple sample_id, file(input_frag_file) from Input_ch
        file src
        file prep_src
        file convert_bed
    output:
        file("*${sample_id}*") into output_ch
    script:
    """
    bash ${prep_src} -i ${input_frag_file} -o .
    bash ${src} -i ${sample_id}_prep.frag.tsv -o . -f ${hg19} -r ${nucleosome_ref} -c ${params.clean_up}
    """
}

// nextflow run run_pipeline_FinaleDB.nf \
// --SampleSheet /datassd/DATA_HIEUNGUYEN/2024/src/ecd_wgs_features/FinaleDB_samplesheet.csv \
// --output /datassd/DATA_HIEUNGUYEN/2024/outdir/FinaleDB \
// --hg19 /datassd/DATA_HIEUNGUYEN/2024/resources/hg19.fa \
// --nucleosome_ref /datassd/DATA_HIEUNGUYEN/2024/resources/rpr_map_EXP0779.sorted.bed \
// --src /datassd/DATA_HIEUNGUYEN/2024/src/ecd_wgs_features/02_calculate_EM_FLEN_NUC_from_FRAG.sh \
// --clean_up true \
// --src_convert_bed /datassd/DATA_HIEUNGUYEN/2024/src/ecd_wgs_features/convert_full_bed_nucleosome.py \
// --prep_src /datassd/DATA_HIEUNGUYEN/2024/src/ecd_wgs_features/preprocess_FinaleDB_input.sh \
// --maxForks 40 \
// -resume

