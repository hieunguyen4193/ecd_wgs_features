//
// Collect mapping metrics
//
include { Summary_Predict } from "../modules/GW_feature_Summary_Predict"

workflow GW_FEATURE_PREDICT {
    take:
        nucleosome
        em
        flen
        ichorcna
    main:

        ch_model = file("$projectDir/bin/model_files")
        ch_example_samples = file("$projectDir/bin/rerun_samples")

        /////////////////////// PREPROCESSING FOR AFTER UMI BAM ///////////////////////
        Summary_Predict ( nucleosome, em, flen, ichorcna, ch_model, ch_example_samples)
}