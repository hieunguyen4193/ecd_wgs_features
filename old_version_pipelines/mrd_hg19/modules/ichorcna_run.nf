process ICHORCNA_RUN {
    tag "$sampleid"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ichorcna:0.3.2--pl5321r42hdfd78af_2' :
        'quay.io/biocontainers/r-ichorcna:0.3.2--pl5321r42hdfd78af_2' }"

    input:
    tuple val(sampleid), path(wig)
    path gc_wig
    path map_wig
    path panel_of_normals
    path centromere

    output:
    tuple val(sampleid), path("*.cna.seg")    , emit: cna_seg
    tuple val(sampleid), path("*.params.txt") , emit: ichorcna_params
    path("*.params.txt") , emit: ichorcna_tumor_fraction
    path "*genomeWide.pdf"                , emit: genome_plot
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sampleid}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    runIchorCNA.R \\
        $args \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        --outDir .

    cp */*genomeWide.pdf .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
}