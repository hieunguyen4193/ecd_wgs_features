//
// Run Facets
//

include { HMMCOPY_READCOUNTER } from "../modules/hmmcopy_readcounter"
include { ICHORCNA_RUN } from "../modules/ichorcna_run"
include { SUMMARY_ICHORCNA } from "../modules/ichorcna_summary"

workflow ICHORCNA_SWF {
    take:
        ch_input // sample, bam, bai
    main:
        
        // get input channel
        ch_gc_wig = params.gc_wig ? file(params.gc_wig) : []
        ch_map_wig = params.map_wig ? file(params.map_wig) : []
        ch_ichorcna_pon = params.ichorcna_pon ? file(params.ichorcna_pon) : []
        ch_centromere = params.centromere ? file(params.centromere) : []

        // HMMcopy readcounter
        HMMCOPY_READCOUNTER ( ch_input )

        // ichorCNA
        ICHORCNA_RUN ( 
            HMMCOPY_READCOUNTER.out.wig,
            ch_gc_wig,
            ch_map_wig,
            ch_ichorcna_pon,
            ch_centromere
            )

        SUMMARY_ICHORCNA( ICHORCNA_RUN.out.ichorcna_tumor_fraction.collect() )
    emit:
        ichorcna = SUMMARY_ICHORCNA.out.output








        

}