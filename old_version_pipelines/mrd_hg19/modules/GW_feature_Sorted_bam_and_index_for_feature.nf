process Sorted_bam_and_index_for_feature {
    tag "$sampleid"
    label 'process_low'
    label 'no_publish'

    container "staphb/samtools:latest"
    
    input:
        tuple val(sampleid), path(bam)
        path(preprocessing_script)
    output:
        tuple val(sampleid), path("${sampleid}.sorted.bam.sorted.bam"), path("${sampleid}.sorted.bam.sorted.bam.bai"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    samtools view -h -f 3 ${bam} | samtools sort -n -@ 4 -o tmp.bam;
    samtools view -h tmp.bam | awk -f $preprocessing_script - > tmp.sam;
    samtools sort -@ 4 -O BAM -o ${sampleid}.sorted.bam.sorted.bam tmp.sam;
    rm -rf tmp.sam;
    rm -rf tmp.bam
    samtools index -@ 4 ${sampleid}.sorted.bam.sorted.bam
    """
}