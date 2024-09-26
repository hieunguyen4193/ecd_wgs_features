//
// mapping reads by bwa, sort bam and umi process by fgbio
//
include { PICARD_FASTQTOSAM } from "../modules/picard_fastqtosam"
include { FGBIO_EXTRACTUMISFROMBAM } from "../modules/fgbio_extractumisfrombam"
include { PICARD_SAMTOFASTQ as PICARD_SAMTOFASTQ_PRE } from "../modules/picard_samtofastq"
include { BWA_MEM as BWA_MEM_PRE  } from "../modules/bwa_mem"
include { PICARD_MERGEBAMALIGNMENT as PICARD_MERGEBAMALIGNMENT_PRE  } from "../modules/picard_mergebamalignment"
include { FGBIO_CORRECTUMIS } from "../modules/fgbio_correctumi"
include { FGBIO_GROUPREADSBYUMI } from "../modules/fgbio_groupreadsbyumi"
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from "../modules/fgbio_callmolecularconsensusreads"
include { PICARD_SAMTOFASTQ as PICARD_SAMTOFASTQ_POST } from "../modules/picard_samtofastq"
include { BWA_MEM as BWA_MEM_POST  } from "../modules/bwa_mem"
include { PICARD_MERGEBAMALIGNMENT as PICARD_MERGEBAMALIGNMENT_POST  } from "../modules/picard_mergebamalignment"
include { FGBIO_FILTERCONSENSUSREADS } from "../modules/fgbio_filterconsensusreads"
include { SAMTOOLS_SORT } from "../modules/samtools_sort"
include { BBMAP_REFORMAT } from "../modules/bbmap_reformat"
include { MULTIQC } from "../modules/multiqc"

workflow MAPPING_SWF {
    take:
        ch_fastq
    main:
        
        ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        multiqc_dirs = Channel.empty()
        ch_versions =  Channel.empty()
        
        // get input channel
        ch_index = params.index ? file(params.index) : []
        ch_fasta = params.fasta ? file(params.fasta) : []
        ch_fai = params.fai ? file(params.fai) : []
        ch_dict = params.dict ? file(params.dict) : []

        // repair fastq files
        BBMAP_REFORMAT (ch_fastq)
        ch_versions = ch_versions.mix(BBMAP_REFORMAT.out.versions)

        // Convert fastq to ummapped bam and extract UMI from reads
        PICARD_FASTQTOSAM( BBMAP_REFORMAT.out.reads )
        ch_versions = ch_versions.mix(PICARD_FASTQTOSAM.out.versions)
        
        FGBIO_EXTRACTUMISFROMBAM ( PICARD_FASTQTOSAM.out.bam )
        ch_versions = ch_versions.mix(FGBIO_EXTRACTUMISFROMBAM.out.versions)

        // Convert ummaped no UMI bam to fastq
        PICARD_SAMTOFASTQ_PRE ( FGBIO_EXTRACTUMISFROMBAM.out.bam )
        ch_versions = ch_versions.mix(PICARD_SAMTOFASTQ_PRE.out.versions)

        // Align reads first time
        BWA_MEM_PRE ( PICARD_SAMTOFASTQ_PRE.out.fastq, ch_index, false )
        ch_versions = ch_versions.mix(BWA_MEM_PRE.out.versions)

        // Merge umapped with mapped bam
        PICARD_MERGEBAMALIGNMENT_PRE (
            FGBIO_EXTRACTUMISFROMBAM.out.bam.join(BWA_MEM_PRE.out.bam),
            ch_fasta,
            ch_dict
        )
        ch_versions = ch_versions.mix(PICARD_MERGEBAMALIGNMENT_PRE.out.versions)

        // Correct UMI
        FGBIO_CORRECTUMIS ( PICARD_MERGEBAMALIGNMENT_PRE.out.bam )
        multiqc_dirs = multiqc_dirs.mix( FGBIO_CORRECTUMIS.out.metrics )
        ch_versions = ch_versions.mix(FGBIO_CORRECTUMIS.out.versions)
        
        // Group reads by UMI
        FGBIO_GROUPREADSBYUMI ( FGBIO_CORRECTUMIS.out.bam, "paired" )
        multiqc_dirs = multiqc_dirs.mix( FGBIO_GROUPREADSBYUMI.out.histogram )
        ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions)

        // Call molecular consensus
        FGBIO_CALLMOLECULARCONSENSUSREADS ( FGBIO_GROUPREADSBYUMI.out.bam )
        ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions)

        // Convert ummaped no UMI bam to fastq
        PICARD_SAMTOFASTQ_POST ( FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam )
        ch_versions = ch_versions.mix(PICARD_SAMTOFASTQ_POST.out.versions)

        // Align reads first time
        BWA_MEM_POST ( PICARD_SAMTOFASTQ_POST.out.fastq, ch_index, false )
        ch_versions = ch_versions.mix(BWA_MEM_POST.out.versions)

        // Merge umapped with mapped bam
        PICARD_MERGEBAMALIGNMENT_POST (
            FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam.join(BWA_MEM_POST.out.bam),
            ch_fasta,
            ch_dict
        )
        ch_versions = ch_versions.mix(PICARD_MERGEBAMALIGNMENT_POST.out.versions)

        // filter consensus reads
        FGBIO_FILTERCONSENSUSREADS ( PICARD_MERGEBAMALIGNMENT_POST.out.bam, ch_fasta, ch_dict )
        ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions)

        // sort bam
        SAMTOOLS_SORT ( FGBIO_FILTERCONSENSUSREADS.out.bam )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

        // Run multiqc
        MULTIQC ( 
            multiqc_dirs.collect(),
            ch_multiqc_config,
            '-m fgbio',
            'Mapping_report.html'
        )

    emit:
        bam = SAMTOOLS_SORT.out.bam
        bam_before = BWA_MEM_PRE.out.bam
        versions = ch_versions
        multiqc = multiqc_dirs

}