process SEQSTATS {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("stats.json"), emit: stats
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    def in_cmd = "-i \"${seqs}\""

    def script = file("${moduleDir}/bin/seq_stats.py")

    """
    python ${script} ${args} ${in_cmd} -o "stats.json" $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqstats: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    """
    touch "stats.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqstats: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
