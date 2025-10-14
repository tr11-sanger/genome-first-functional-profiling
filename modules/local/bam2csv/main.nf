process BAM2CSV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(genome2cds), path(genome_contig_mapping), path(genome_species), path(refs)
    val delete_bam

    output:
    tuple val(meta), path("output/*species_coverage.tsv.gz"), emit: species_profile
    tuple val(meta), path("output/*genome_coverage.tsv.gz"), emit: genome_profile
    tuple val(meta), path("output/*species_cds_coverage.tsv.gz"), path("output/*species_index.txt.gz"), emit: species_cds_profile

    script:
    def prefix = task.ext.prefix ? task.ext.prefix : meta.id 
    def rm_cmd = delete_bam ? "rm -r ${bam}/" : "" 
    def script = "${moduleDir}/bin/bam2csv_sqlite.py"
    """
    python ${script} \\
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
    touch "output/${prefix}_species_coverage.tsv.gz"
    touch "output/${prefix}_genome_coverage.tsv.gz"
    touch "output/${prefix}_species_cds_coverage.tsv.gz"
    touch "output/${prefix}_species_index.txt.gz"
    """
}
