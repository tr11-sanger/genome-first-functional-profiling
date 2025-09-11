process CONCATENATE {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), val(out_fn), path(files, stageAs: "input_files/?/*")

    output:
    tuple val(meta), path("concatenated/${out_fn}"), emit: concatenated_file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def out_dir = "concatenated"

    """
    mkdir -p $out_dir
    find input_files/*/* -exec cat {} + > ${out_dir}/${out_fn}
    """

    stub:
    def args = task.ext.args ?: ''
    def out_dir = "concatenated"

    """
    mkdir -p $out_dir
    touch ${out_dir}/${out_fn}
    """
}
