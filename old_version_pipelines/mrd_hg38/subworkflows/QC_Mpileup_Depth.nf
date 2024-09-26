//
// Collect mapping metrics
//
include { QC_Prepare_Depth } from "../modules/QC_Prepare_Depth"
include { QC_Prepare_Mpileup } from "../modules/QC_Prepare_Mpileup"
include { Summary_Depth } from "../modules/QC_Summary_Depth"
include { Summary_Mpileup } from "../modules/QC_Summary_Mpileup"

workflow QC_MPILEUP_DEPTH {
    take:
        ch_bam
    main:
        /////////////////////// EM feature ///////////////////////
        QC_Prepare_Depth ( ch_bam )
        QC_Prepare_Mpileup ( ch_bam )


        Summary_Depth ( QC_Prepare_Depth.out.output.collect() )
        Summary_Mpileup ( QC_Prepare_Mpileup.out.output.collect() )
}