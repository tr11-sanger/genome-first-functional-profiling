process BAM2CSV {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.22.1--py39hdd5828d_3'
        : 'biocontainers/pysam:0.22.1--py39hdd5828d_3'}"

    publishDir "${params.databases.cache_path}", mode: 'copy'
    errorStrategy 'retry'

    input:
    tuple val(meta), path(bam)
    path(script)
    val delete_bam

    output:
    tuple val(meta), path("*.csv.gz"), emit: csv

    script:
    def rm_cmd = delete_bam ? "rm -r ${bam}/" : "" 
    """
    python ${script} -i ${bam} -o "${meta.id}.csv" 
    gzip "${meta.id}.csv"
    ${rm_cmd}
    """

    stub:
    """
    touch "${meta.id}.csv.gz"
    """
}
