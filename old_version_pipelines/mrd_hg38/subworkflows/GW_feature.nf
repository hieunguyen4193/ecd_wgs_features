//
// Collect mapping metrics
//
include { Sorted_bam_and_index_for_feature } from "../modules/GW_feature_Sorted_bam_and_index_for_feature"
include { Prepare_EM } from "../modules/GW_feature_Prepare_EM"
include { Feature_EM } from "../modules/GW_feature_Feature_EM"
include { Summary_EM } from "../modules/GW_feature_Summary_EM"
include { Prepare_FLEN } from "../modules/GW_feature_Prepare_FLEN"
include { Feature_FLEN } from "../modules/GW_feature_Feature_FLEN"
include { Summary_FLEN } from "../modules/GW_feature_Summary_FLEN"
include { Create_bed } from "../modules/GW_feature_Create_bed"
include { Create_full_bed } from "../modules/GW_feature_Create_full_bed"
include { Bedtool_closest } from "../modules/GW_feature_Bedtool_closest"
include { Dist_bed } from "../modules/GW_feature_Dist_bed"
include { Density } from "../modules/GW_feature_Density"
include { Summary_NUCLEOSOME } from "../modules/GW_feature_Summary_NUCLEOSOME"

workflow GW_FEATURE {
    take:
        ch_bam
    main:

        ch_preprocessing_script = file("$projectDir/bin/preprocessing_script.awk")
        ch_fa = params.fasta ? file(params.fasta) : []
        ch_ref_nucleosome = file("$projectDir/bin/hglft_genome_rpr_map_Budhraja_STM2023_formated12.bed")
        ch_example_samples = file("$projectDir/bin/rerun_samples")

        /////////////////////// PREPROCESSING FOR AFTER UMI BAM ///////////////////////
        Sorted_bam_and_index_for_feature ( ch_bam, ch_preprocessing_script )


        /////////////////////// EM feature ///////////////////////
        Prepare_EM ( Sorted_bam_and_index_for_feature.out.output, ch_fa )
        Feature_EM ( Prepare_EM.out.output )
        Summary_EM ( Feature_EM.out.output.collect(), ch_example_samples )


        /////////////////////// FLEN feature ///////////////////////
        Prepare_FLEN ( Sorted_bam_and_index_for_feature.out.output )
        Feature_FLEN ( Prepare_FLEN.out.output )
        Summary_FLEN (Feature_FLEN.out.output.collect(), ch_example_samples )


        /////////////////////// NUCLEOSOME feature /////////////////
        Create_bed ( Sorted_bam_and_index_for_feature.out.output )
        Create_full_bed ( Create_bed.out.output )
        Bedtool_closest ( Create_full_bed.out.output, ch_ref_nucleosome )
        Dist_bed ( Bedtool_closest.out.output)
        Density (Dist_bed.out.output )
        Summary_NUCLEOSOME ( Density.out.output.collect(), ch_example_samples )

    emit:
        nucleosome = Summary_NUCLEOSOME.out.output
        em = Summary_EM.out.output
        flen = Summary_FLEN.out.output
}