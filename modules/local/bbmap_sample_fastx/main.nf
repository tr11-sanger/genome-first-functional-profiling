process BBMAP_SAMPLE_FASTX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(fastx, stageAs: 'input/*')
    val subsample_n
    val shuffle

    output:
    tuple val(meta), path("subsampled/input/*"), emit: fastx
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    single_file = (fastx instanceof Collection) ? (fastx.size() == 1) : true
    in_reads  = single_file ? "in=${fastx}" : "in=${fastx[0]} in2=${fastx[1]}"
    out_reads = meta.single_end ? "out=subsampled/${fastx.name}" : "out=subsampled/${fastx[0].name} out2=subsampled/${fastx[1].name} outs=${prefix}_singleton"

    """
    mkdir -p subsampled/input
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    reformat.sh \\
        -Xmx\$maxmem \\
        $in_reads \\
        $out_reads \\
        threads=${task.cpus} \\
        allowidenticalnames=t \\
        trimreaddescription=t \\
        samplereadstarget=${subsample_n} \\
        sampleseed=42 \\
        ${args} \\
        &> ${prefix}.reformat.sh.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > subsampled/${fastx[0].name}
    echo "" | gzip > subsampled/${fastx[1].name}
    touch ${prefix}.reformat.sh.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}