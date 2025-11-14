process BAM2CSV_GREEDY {
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
    """
    samtools view -@ 1 ${bam} | bam2csv_sqlite_greedy.py \\
        --genome_cds_filepaths "${genome2cds}" \\
        --genome_contig_mapping "${genome_contig_mapping}" \\
        --genome_species "${genome_species}" \\
        --refs "${refs}" \\
        --output_prefix "${prefix}" \\
        --min_align_length "${params.min_align_length}" \\
        --max_align_length "${params.max_align_length}" \\
        --min_ani "${params.min_ani}" \\
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
