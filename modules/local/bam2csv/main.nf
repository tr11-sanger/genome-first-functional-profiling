process BAM2CSV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.databases.cache_path}", mode: 'copy'
    errorStrategy 'retry'

    input:
    tuple val(meta), path(bam), path(genome2cds), path(genome_contig_mapping), path(genome_species), path(refs)
    val delete_bam

    output:
    tuple val(meta), path("output/*species_coverage.tsv"), emit: species_profile
    tuple val(meta), path("output/*genome_coverage.tsv"), emit: genome_profile
    tuple val(meta), path("output/*species_cds_coverage.tsv"), path("output/*species_index.txt"), emit: species_cds_profile
    tuple val(meta), path("output/*mapping_statistics.json"), emit: mapping_statistics

    script:
    def prefix = task.ext.prefix ? task.ext.prefix : meta.id 
    def rm_cmd = delete_bam ? "rm -r ${bam}/" : "" 
    def script = "${moduleDir}/bin/bam2csv.py"
    """
    python ${script} \
        --genome_cds_filepaths "${genome2cds}" \\
        --genome_contig_mapping "${genome_contig_mapping}" \\
        --genome_species "${genome_species}" \\
        --bam "${bam}" \\
        --refs "${refs}" \\
        --output_prefix "${prefix}" \\
        --output_dir output
    ${rm_cmd}
    """

    stub:
    def prefix = task.ext.prefix ? task.ext.prefix : meta.id 
    """
    mkdir output
    touch "output/${prefix}_species_coverage.tsv"
    touch "output/${prefix}_genome_coverage.tsv"
    touch "output/${prefix}_species_cds_coverage.tsv"
    touch "output/${prefix}_species_index.txt"
    touch "output/${prefix}_mapping_statistics.json"
    """
}
