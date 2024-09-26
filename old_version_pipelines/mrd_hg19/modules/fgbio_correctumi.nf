process FGBIO_CORRECTUMIS {
    tag "$sampleid"
    label "process_medium"
    label 'no_publish'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:1.4.0--hdfd78af_1' :
        'quay.io/biocontainers/fgbio:1.4.0--hdfd78af_1' }"

    input:
    tuple val(sampleid), path(bam)  

    output:
    tuple val(sampleid), path("*.fixedumi.bam")  , emit: bam
    path("*.metrics.txt"), emit: metrics
    path('versions.yml')    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--max-mismatches=3 --min-distance=1 -t RX -u GAGACGAT TTCCAAGG CGCATGAT ACGGAACA CGGCTAAT GCTATCCT TGGACTCT ATCCAGAG CTTAGGAC GTGCCATA TCGCTGTT TTCGTTGG AAGCACTG GTCGAAGA ACCACGAT GATTACCG GCACAACT GCGTCATT GAAGGAAG ACTGAGGT TGAAGACG GTTACGCA AGCGTGTT GATCGAGT TTGCGAAG CTGTTGAC GATGTGTG ACGTTCAG TTGCAGAC CAATGTGG ACGACTTG ACTAGGAG'
    def prefix = task.ext.prefix ?: "${sampleid}"
    def avail_mem = (task.memory.giga).intValue()

    """
    fgbio \\
        -Djava.io.tmpdir="tmp" -Xmx${avail_mem}g \\
        CorrectUmis \\
        ${args} \\
        -i ${bam} \\
        -o ${prefix}.fixedumi.bam \\
        -M ${prefix}.metrics.txt \\
        -r ${prefix}.rejected.bam
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}