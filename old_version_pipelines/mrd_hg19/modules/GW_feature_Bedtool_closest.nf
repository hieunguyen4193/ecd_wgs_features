process Bedtool_closest {
    tag "$sampleid"
    label 'no_publish'

    container "gene110/samtools_bedtools:latest"
    
    input:
        tuple val(sampleid), path(data)
        path(ch_ref_nucleosome)
    output:
        tuple val(sampleid), path("${sampleid}.rpr_map.bed"), emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
    """
    bedtools closest -a $data -b $ch_ref_nucleosome | cut -f1,2,10 > ${sampleid}.rpr_map.bed
    """
}