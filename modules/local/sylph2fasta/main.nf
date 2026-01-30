process SYLPH2FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(sylph_results)
    path(lookup)

    output:
    tuple val(meta), path("genomes_table.csv") , emit: genomes
    tuple val(meta), path("contig_mapping.csv"), emit: contigs
    tuple val(meta), path("genomes.fna")       , emit: fasta
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # gunzip -c $sylph_results | cut -d"," -f10 | tail -n +2 | sort > data.tmp
    sourmash2fasta.py -s $sylph_results -o data.tmp
    cat data.tmp | sort > data.tmp2
    join -t',' -1 1 -2 1 -o 2.1,2.2 data.tmp2 $lookup > genomes_table.csv
    rm data.tmp
    rm data.tmp2
    
    echo -n '' > contig_mapping.csv
    echo -n '' > genomes.fna
    while read line; do
      IFS=, read g fp <<< "\$line"
      if [[ \$fp =~ \\.gz\$ ]]; then
        gunzip -c \$fp >> genomes.fna
        gunzip -c \$fp | grep '>' | {
        while read c; do
          c_=(\$c)
          echo "\$g,\${c_[0]:1}" >> contig_mapping.csv
        done
        }
      else
        cat \$fp >> genomes.fna
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
