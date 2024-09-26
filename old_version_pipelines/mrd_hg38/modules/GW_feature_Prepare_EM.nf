process Prepare_EM {
    tag "$sampleid"
    label 'process_low'
    label 'no_publish'

    container "gene110/samtools_bedtools:latest"
    
    input:
        tuple val(sampleid), path(bam), path(bai)
        path(fa)
    output:
        tuple val(sampleid), path("${sampleid}.EM.txt"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    samtools view -@ 4 -q 30 -f 2 ${bam} | awk '{if(\$9 > 0) {start=\$4 - 1; end= \$4 - 1 + 4; name= \$1"_"\$9; strand = "+"} else {start=\$8 - \$9 - 4 - 1; end= \$8 - \$9 -1; name= \$1"_"\$9; strand = "-"}; print \$3 "\\t" start "\\t" end "\\t" name "\\t" "1" "\\t" strand}' > ${sampleid}.endcoord4bp.bed;
    bedtools getfasta -s -name -tab -fi ${fa} -bed ${sampleid}.endcoord4bp.bed  > ${sampleid}.EM.txt
    rm -rf ${sampleid}.endcoord4bp.bed;
    """
}