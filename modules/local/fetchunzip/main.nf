process FETCHUNZIP {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    publishDir "${params.databases.cache_path}", mode: 'copy'
    errorStrategy 'retry'

    input:
    tuple val(meta), val(dir_name), path(fp)

    output:
    tuple val(meta), path(dir_name)

    script:
    if (fp.name[-7..-1] == '.tar.gz') {
        """
        #!/bin/bash
        mkdir -p "${dir_name}"
        tar -xvzf "\$(readlink ${fp.name})" -C "${dir_name}"
        exit 0
        """
    }
    else {
        if (fp.name[-3..-1] == '.gz') {
            """
        #!/bin/bash
        mkdir -p "${dir_name}"
        gunzip -c "\$(readlink ${fp.name})" > "${dir_name}/${fp.name[0..-4]}"
        exit 0
        """
        }
        else {
            """
        #!/bin/bash
        mkdir -p "${dir_name}"
        cp "\$(readlink ${fp.name})" "${dir_name}/${fp.name}"
        exit 0
        """
        }
    }

    stub:
    db_files = meta.files.collect { _k, v -> v }
    db_files_cmd = db_files
        .collect { fn ->
            def new_fp = "${dir_name}/${meta.base_dir}/${fn}"
            return "mkdir -p \"\$(dirname \"${new_fp}\")\" && touch \"${new_fp}\""
        }
        .join('\n')
    if (fp.name[-7..-1] == '.tar.gz') {
        """
        #!/bin/bash
        mkdir -p "${dir_name}/${meta.base_dir}"
        ${db_files_cmd}
        """
    }
    else {
        """
        #!/bin/bash
        mkdir -p "${dir_name}"
        touch "${dir_name}/${fp.name}"
        """
    }
}
