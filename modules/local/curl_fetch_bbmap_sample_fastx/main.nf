process CURL_FETCH_BBMAP_SAMPLE_FASTX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/bbmap:39.15--h92535d8_0' :
    //     'biocontainers/bbmap:39.15--h92535d8_0' }"

    input:
    tuple val(meta), val(urls)
    val   subsample_n
    val   max_retries
    val   wait_retry
    val   timeout
    val   resume
    val   soft_fail
    val   delete_original

    output:
    tuple val(meta), path('*.sampled.fastq.gz'), env(SUCCESS), emit: fastx
    tuple val(meta), path('*.download.log')                  , emit: download_log
    tuple val(meta), path('*.status')                        , emit: status
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh 2>&1 | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions
    tuple val("${task.process}"), val('curl'),  eval('curl --version 2>&1 | head -1 | sed -e "s/curl //;s/ .*//"'), emit: versions_curl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_single = meta.single_end ? 'true' : 'false'

    // Build curl flags
    def curl_resume = resume ? '-C -' : ''
    def curl_flags = "--retry ${max_retries} --retry-delay ${wait_retry} --connect-timeout ${timeout} --retry-connrefused ${curl_resume} -f -L -sS"

    // Determine URL list and download targets
    def url_list = urls instanceof List ? urls : [urls]
    def url1 = url_list[0]
    def url2 = url_list.size() > 1 ? url_list[1] : ''

    // Download filenames
    def dl_file1 = is_single == 'true' ? "${prefix}.fastq.gz" : "${prefix}_R1.fastq.gz"
    def dl_file2 = "${prefix}_R2.fastq.gz"

    // Soft fail empty file creation
    def empty_fastx = ''
    if (meta.single_end) {
        empty_fastx = "echo '' | gzip > ${prefix}.sampled.fastq.gz"
    } else {
        empty_fastx = "echo '' | gzip > ${prefix}_R1.sampled.fastq.gz ; echo '' | gzip > ${prefix}_R2.sampled.fastq.gz"
    }

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Helper: fetch a URL or copy a local file
    fetch_url() {
        local url="\$1"
        local output="\$2"
        local log_file="\$3"

        if [[ "\$url" == /* ]] || [[ "\$url" == file://* ]]; then
            local local_path="\${url#file://}"
            echo "Copying local file \${local_path} -> \${output}" >> "\$log_file"
            cp "\${local_path}" "\${output}" >> "\$log_file" 2>&1
        else
            echo "Downloading \${url} -> \${output}" >> "\$log_file"
            curl \${CURL_FLAGS} -o "\${output}" "\${url}" >> "\$log_file" 2>&1
        fi
    }

    PREFIX="${prefix}"
    MAX_RETRIES="${max_retries}"
    WAIT_RETRY="${wait_retry}"
    IS_SINGLE="${is_single}"
    SOFT_FAIL="${soft_fail}"
    DELETE_ORIGINAL="${delete_original}"
    URL1="${url1}"
    URL2="${url2}"
    DL_FILE1="${dl_file1}"
    DL_FILE2="${dl_file2}"
    CURL_FLAGS="${curl_flags}"
    DOWNLOAD_LOG="\${PREFIX}.download.log"

    SUBSAMPLE_N="${subsample_n}"
    BBMAP_ARGS="${args}"

    SUCCESS=false

    : > "\${DOWNLOAD_LOG}"

    for attempt in \$(seq 1 \${MAX_RETRIES}); do
        echo "=== Attempt \${attempt}/\${MAX_RETRIES} ===" >> "\${DOWNLOAD_LOG}"
        echo "Started: \$(date -Iseconds)" >> "\${DOWNLOAD_LOG}"

        # --- Download ---
        DOWNLOAD_OK=true

        if ! fetch_url "\${URL1}" "\${DL_FILE1}" "\${DOWNLOAD_LOG}"; then
            echo "Download failed for \${URL1}" >> "\${DOWNLOAD_LOG}"
            DOWNLOAD_OK=false
        fi

        if [ "\${DOWNLOAD_OK}" = "true" ] && [ -n "\${URL2}" ]; then
            if ! fetch_url "\${URL2}" "\${DL_FILE2}" "\${DOWNLOAD_LOG}"; then
                echo "Download failed for \${URL2}" >> "\${DOWNLOAD_LOG}"
                DOWNLOAD_OK=false
            fi
        fi

        if [ "\${DOWNLOAD_OK}" = "false" ]; then
            echo "Download failed on attempt \${attempt}" >> "\${DOWNLOAD_LOG}"
            rm -f "\${DL_FILE1}" "\${DL_FILE2}"
            if [ "\${attempt}" -lt "\${MAX_RETRIES}" ]; then
                echo "Sleeping \${WAIT_RETRY}s before next attempt..." >> "\${DOWNLOAD_LOG}"
                sleep "\${WAIT_RETRY}"
            fi
            continue
        fi

        echo "Download succeeded, running reformat.sh..." >> "\${DOWNLOAD_LOG}"

        # --- Run reformat.sh for subsampling ---
        TOOL_OK=true

        if [ "\${IS_SINGLE}" = "true" ]; then
            if ! reformat.sh \\
                in="\${DL_FILE1}" \\
                out="\${PREFIX}.sampled.fastq.gz" \\
                samplereadstarget=\${SUBSAMPLE_N} \\
                sampleseed=42 \\
                \${BBMAP_ARGS} \\
                >> "\${DOWNLOAD_LOG}" 2>&1; then
                TOOL_OK=false
            fi
        else
            if ! reformat.sh \\
                in="\${DL_FILE1}" \\
                in2="\${DL_FILE2}" \\
                out="\${PREFIX}_R1.sampled.fastq.gz" \\
                out2="\${PREFIX}_R2.sampled.fastq.gz" \\
                samplereadstarget=\${SUBSAMPLE_N} \\
                sampleseed=42 \\
                \${BBMAP_ARGS} \\
                >> "\${DOWNLOAD_LOG}" 2>&1; then
                TOOL_OK=false
            fi
        fi

        if [ "\${TOOL_OK}" = "true" ]; then
            echo "reformat.sh succeeded on attempt \${attempt}" >> "\${DOWNLOAD_LOG}"
            if [ "\${DELETE_ORIGINAL}" = "true" ]; then
                echo "Deleting original downloaded files" >> "\${DOWNLOAD_LOG}"
                rm -f "\${DL_FILE1}" "\${DL_FILE2}"
            fi
            SUCCESS=true
            break
        else
            echo "reformat.sh failed on attempt \${attempt}" >> "\${DOWNLOAD_LOG}"
            rm -f "\${PREFIX}.sampled.fastq.gz" "\${PREFIX}_R1.sampled.fastq.gz" "\${PREFIX}_R2.sampled.fastq.gz"
            rm -f "\${DL_FILE1}" "\${DL_FILE2}"
            if [ "\${attempt}" -lt "\${MAX_RETRIES}" ]; then
                echo "Sleeping \${WAIT_RETRY}s before next attempt..." >> "\${DOWNLOAD_LOG}"
                sleep "\${WAIT_RETRY}"
            fi
        fi
    done

    # Write status
    echo "\${SUCCESS}" > "\${PREFIX}.status"

    if [ "\${SUCCESS}" = "false" ]; then
        echo "All \${MAX_RETRIES} attempts failed" >> "\${DOWNLOAD_LOG}"
        if [ "\${SOFT_FAIL}" = "true" ]; then
            echo "Soft fail mode: creating empty output files" >> "\${DOWNLOAD_LOG}"
            ${empty_fastx}
        else
            exit 1
        fi
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_fastx = meta.single_end ? "echo '' | gzip > ${prefix}.sampled.fastq.gz" : "echo '' | gzip > ${prefix}_R1.sampled.fastq.gz ; echo '' | gzip > ${prefix}_R2.sampled.fastq.gz"
    """
    SUCCESS=true
    $touch_fastx
    touch "${prefix}.download.log"
    echo "true" > ${prefix}.status
    """
}
