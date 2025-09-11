include { FETCHUNZIP } from '../../../modules/local/fetchunzip/main'

workflow FETCHDB {
    take:
    fetch_ch // channel: val(meta)
    cache_path // String

    main:
    local_ch = fetch_ch
        .filter { meta -> meta.local_path }
        .map { meta -> [meta, file(meta.local_path, checkIfExists: true)] }

    cache_path_ch = fetch_ch
        .filter { meta -> ((!meta.local_path) && meta.remote_path) }
        .map { meta -> [meta, file("${cache_path}/${meta.id}")] }

    download_ch = cache_path_ch
        .filter { _meta, cache_fp -> (cache_path.isEmpty() || (!cache_fp.exists()) || params.force_download_dbs==true) }
        .map { meta, _cache_fp -> [meta, meta.id, file(meta.remote_path, checkIfExists: true)] }

    cache_ch = cache_path_ch
       .filter { _meta, cache_fp -> ((!cache_path.isEmpty()) && cache_fp.exists() && (params.force_download_dbs==false)) }
       .map { meta, cache_fp -> [meta, cache_fp] }

    FETCHUNZIP(download_ch)
    downloaded_ch = FETCHUNZIP.out

    emit:
    dbs = local_ch.mix(cache_ch).mix(downloaded_ch)
}
