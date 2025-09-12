process SOURMASH2FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph:0.7.0--h919a2d8_0' :
        'biocontainers/sylph:0.7.0--h919a2d8_0' }"

    input:
    tuple val(meta), path(sourmash_results)
    path(lookup)

    output:
    path "genomes_table.csv" , emit: genomes
    path "contig_mapping.csv", emit: contigs
    path "genomes.fna"       , emit: fasta
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    lookup="/data/pam/team162/tr11/scratch/sourmash/genome_fp_lookup_sorted.txt"
    gunzip -c $sourmash_results | cut -d"," -f10 | tail -n +2 | sort > data.tmp
    join -t',' -1 1 -2 1 -o 2.1,2.2 data.tmp $lookup > genomes_table.csv
    rm data.tmp
    cut -d',' -f2 genomes_table.csv | xargs cat > genomes.fna
    
    : > contig_mapping.csv
    while read line; do
      IFS=, read g fp <<< "\$line"
      if [[ \$fp =~ \.gz$ ]]; then
        gunzip -c \$fp | grep '>' | {
        while read c; do
          c_=(\$c)
          echo "\$g,\${c_[0]:1}" >> contig_mapping.csv
        done
        }
      else
        cat \$fp | grep '>' | {
        while read c; do
          c_=(\$c)
          echo "\$g,\${c_[0]:1}" >> contig_mapping.csv
        done
        }
      fi
    done < genomes_table.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch genomes_table.csv
    touch contig_mapping.csv
    touch genomes.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """
}
