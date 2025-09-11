process CHUNKFASTX {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("chunked/*"), emit: chunked_reads
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    def prefix = reads.getName().tokenize('.')[0]
    def extension = reads.getName().tokenize('.')[1..-1].join('.')
    if (extension.endsWith('.gz')) {
        extension = extension.tokenize('.')[0..-2].join('.')
    }
    def out_fn = "${prefix}.${extension}"

    def reads_cmd = "-i \"${reads}\""

    def script = file("${moduleDir}/bin/chunk_fastx.py")

    """
    mkdir -p chunked
    python ${script} ${args} ${reads_cmd} -o "chunked/${out_fn}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chunkfastx: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def prefix = reads.getName().tokenize('.')[0]
    def extension = reads.getName().tokenize('.')[1..-1].join('.')
    def out_fn = "${prefix}.${extension}"
    """
    mkdir -p chunked
    touch "chunked/${prefix}.${extension}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        renamepairedfastxheaders: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
