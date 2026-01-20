process SQLITE2PROFILE_TOP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(sqlite), path(genome2cds)

    output:
    tuple val(meta), path("output/*species_coverage.tsv.gz"), emit: species_profile
    tuple val(meta), path("output/*genome_coverage.tsv.gz"), emit: genome_profile
    tuple val(meta), path("output/*species_cds_coverage.tsv.gz"), path("output/*species_index.txt.gz"), emit: species_cds_profile

    script:
    def prefix = task.ext.prefix ? task.ext.prefix : meta.id 
    """
    sqlite2profile_top.py \\
        --sqlite "${sqlite}" \\
        --genome_cds_filepaths "${genome2cds}" \\
        --output_prefix "${prefix}" \\
        --min_ani "${params.min_ani}" \\
        --min_coverage_ratio "${params.min_coverage_ratio}" \\
        --output_dir output
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
