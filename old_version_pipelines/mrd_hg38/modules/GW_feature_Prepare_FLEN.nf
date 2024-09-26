process Prepare_FLEN {
    tag "$sampleid"
    label 'no_publish'

    container "staphb/samtools:latest"
    
    input:
        tuple val(sampleid), path(bam), path(bai)
    output:
        tuple val(sampleid), path("${sampleid}.FLEN.txt"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    samtools view ${bam} | cut -f9 > ${sampleid}.FLEN.txt;
    """
}