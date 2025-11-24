process SEQKIT_SHUFFLE_FASTX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fastx, stageAs: 'input/*')
    val subsample_n
    val shuffle

    output:
    tuple val(meta), path("subsampled/*"), emit: fastx
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ? "${task.ext.prefix}.${meta.id}" : "${meta.id}"
    def call_gzip = fastx.toString().endsWith('.gz') ? "| gzip -c ${args2}" : ''
    def head_cmd = shuffle ? "shuf -n ${subsample_n}" : "head -n ${subsample_n}"

    if (meta.single_end) {
    """
        seqkit seq --name --only-id ${fastx} \\
        | ${head_cmd} \\
        > seq_ids

        mkdir -p subsampled/input

        seqkit faidx --region-file seq_ids --full-head ${fastx} \\
        ${call_gzip} \\
        > subsampled/${fastx.name}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | cut -d' ' -f2)
        END_VERSIONS
        """
    } else {
        """
        seqkit seq --name --only-id ${fastx[0]} \\
        | ${head_cmd} \\
        > seq_ids
    
        mkdir -p subsampled/input

        seqkit faidx --region-file seq_ids --full-head ${fastx[0]} \\
        ${call_gzip} \\
        > subsampled/${fastx[0].name}

        seqkit faidx --region-file seq_ids --full-head ${fastx[1]} \\
        ${call_gzip} \\
        > subsampled/${fastx[1].name}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | cut -d' ' -f2)
        END_VERSIONS
        """
    }
    stub:
    prefix = task.ext.prefix ? "${task.ext.prefix}.${meta.id}" : "${meta.id}"
    if (meta.single_end) {
        """
        mkdir -p subsampled/input
        touch subsampled/${fastx.name}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | cut -d' ' -f2)
        END_VERSIONS
        """
    } else {
        """
        mkdir -p subsampled/input
        touch subsampled/${fastx[0].name}
        touch subsampled/${fastx[1].name}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | cut -d' ' -f2)
        END_VERSIONS
        """

    }
}