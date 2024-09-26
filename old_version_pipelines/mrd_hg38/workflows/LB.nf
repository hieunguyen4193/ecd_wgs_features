/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW
//
include { READS_QC_SWF } from "../subworkflows/reads_qc"
include { MAPPING_SWF } from "../subworkflows/mapping"
include { GW_FEATURE } from "../subworkflows/GW_feature"
include { ICHORCNA_SWF } from "../subworkflows/ichorCNA"
include { GW_FEATURE_PREDICT } from "../subworkflows/GW_feature_predict"
include { QC_MPILEUP_DEPTH as QC_MPILEUP_DEPTH_PRE  } from "../subworkflows/QC_Mpileup_Depth"
include { QC_MPILEUP_DEPTH as QC_MPILEUP_DEPTH_POST  } from "../subworkflows/QC_Mpileup_Depth"
//
// MODULES
//
include { MULTIQC as MULTIQC_BRIEF } from "../modules/multiqc"
include { DUMPSOFTWAREVERSIONS } from "../modules/versions/versions"
include { SAMTOOLS_SORT } from "../modules/samtools_sort"





workflow LB_WF {

    // Initialize ChannelDump
    ChannelDump.initialise(workflow, params.outdir)

    // create fastq channel
    ch_fastq = Channel
        .fromFilePairs("$params.fqdir/*_R{1,2}*")
        .ifEmpty { error "Cannot find any reads matching: ${params.fqdir}/*_R{1,2}*" }

    // multiqc channels
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    multiqc_dirs = Channel.empty()
    ch_versions = Channel.empty()

    // Reads QC
    READS_QC_SWF ( ch_fastq )
    multiqc_dirs = multiqc_dirs.mix(READS_QC_SWF.out.multiqc)
    ch_versions = ch_versions.mix(READS_QC_SWF.out.versions)

    // Mapping
    MAPPING_SWF ( ch_fastq )
    multiqc_dirs = multiqc_dirs.mix(MAPPING_SWF.out.multiqc)
    ch_versions = ch_versions.mix(MAPPING_SWF.out.versions)

    MAPPING_SWF.out.bam
        .collect(flat: false)
        .subscribe { ChannelDump.writeCsv(it, ["sample", "bam", "bai"], 'after_UMI.csv') }

    QC_MPILEUP_DEPTH_POST(MAPPING_SWF.out.bam)

    GW_FEATURE(MAPPING_SWF.out.bam)

    SAMTOOLS_SORT(MAPPING_SWF.out.bam_before)

    ICHORCNA_SWF(SAMTOOLS_SORT.out.bam)

    QC_MPILEUP_DEPTH_PRE(SAMTOOLS_SORT.out.bam)

    SAMTOOLS_SORT.out.bam
        .collect(flat: false)
        .subscribe { ChannelDump.writeCsv(it, ["sample", "bam", "bai"], 'before_UMI.csv') }

    GW_FEATURE_PREDICT(GW_FEATURE.out.nucleosome, GW_FEATURE.out.em, GW_FEATURE.out.flen, ICHORCNA_SWF.out.ichorcna)

    

    //
    // MODULE: Pipeline reporting
    //
    DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    multiqc_dirs = multiqc_dirs.mix(DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    
    //
    // MODULE: MultiQC
    //

    MULTIQC_BRIEF ( 
        multiqc_dirs.collect(),
        ch_multiqc_config,
        '-m fastqc -m picard -m custom_content',
        'Metrics.html'
    )


}