process SOURMASH2FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(sourmash_results)

    output:
    path "*.fasta", emit: fasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash2fasta: \$(sourmash2fasta --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash2fasta: \$(sourmash2fasta --version)
    END_VERSIONS
    """
}
