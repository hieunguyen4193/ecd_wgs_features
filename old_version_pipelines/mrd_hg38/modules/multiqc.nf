process MULTIQC {
    label 'process_single'
    
    // conda "bioconda::multiqc=1.14"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
        path  multiqc_dirs, stageAs: "?/*"
        path  config
        val module
        val filename
    
    output:
        path "*.html", emit: report
        path "*_data/*.txt", emit: data
        path "versions.yml"        , emit: versions
    
    when:
    task.ext.when == null || task.ext.when


    script:
    def filename_command = ""
    if (filename){
        filename_command = "--filename $filename"
    } else if (task.ext.prefix){
        filename_command = "--filename ${task.ext.prefix}"
    }
    def args = task.ext.args ?: ''

    """
    multiqc $module $args --config $config $filename_command $multiqc_dirs 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}