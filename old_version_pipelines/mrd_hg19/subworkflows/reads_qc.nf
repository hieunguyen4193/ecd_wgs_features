//
// Perform reads quality control.
//

include { FASTQC } from '../modules/fastqc'
include { MULTIQC } from "../modules/multiqc"

workflow READS_QC_SWF {
    take:
        ch_fastq
    main:
        
        ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        multiqc_dirs = Channel.empty()
        ch_versions =  Channel.empty()
        
        // Run fastqc
        FASTQC ( ch_fastq )
        multiqc_dirs = multiqc_dirs.mix( 
            FASTQC.out.fastqc_zip)
        ch_versions = ch_versions.mix(FASTQC.out.versions)

        // Run multiqc
        MULTIQC ( 
            multiqc_dirs.collect(),
            ch_multiqc_config,
            '-m fastqc',
            'Reads_QC_report.html'
        )

    emit:
        versions = ch_versions
        multiqc = multiqc_dirs

}