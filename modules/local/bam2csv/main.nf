process BAM2CSV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.databases.cache_path}", mode: 'copy'
    errorStrategy 'retry'

    input:
    tuple val(meta), path(bam)
    val delete_bam

    output:
    tuple val(meta), path("*.csv.gz"), emit: csv

    script:
    def rm_cmd = delete_bam ? "rm -r ${bam}/" : "" 
    def script = "${moduleDir}/bin/bam2csv.py"
    """
    bin/python ${script} -i ${bam} -o "${meta.id}.csv" 
    gzip "${meta.id}.csv"
    ${rm_cmd}
    """

    stub:
    """
    touch "${meta.id}.csv.gz"
    """
}
