/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { FETCHDB } from '../subworkflows/local/fetchdb/main'
include { BBMAP_SAMPLE_FASTX } from '../modules/local/bbmap_sample_fastx/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { SYLPH_PROFILE } from '../modules/nf-core/sylph/profile/main'
include { SYLPH_QUERY } from '../modules/local/sylph/query/main'
include { SOURMASH_GATHER } from '../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH } from '../modules/nf-core/sourmash/sketch/main'
include { SOURMASH2FASTA } from '../modules/local/sourmash2fasta/main'
include { SYLPH2FASTA } from '../modules/local/sylph2fasta/main'
include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN_BAM2SQLITE } from '../modules/local/bowtie2/align_bam2sqlite/main'
include { SQLITE2PROFILE_TOP } from '../modules/local/sqlite2profile_top/main'
include { SQLITE2PROFILE_GREEDY } from '../modules/local/sqlite2profile_greedy/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GFFP {
    main:
    ch_versions = channel.empty()

    // Parse samplesheet and fetch reads
    reads_ch = channel.fromList(
            samplesheetToList(
                params.samplesheet, 
                "${workflow.projectDir}/assets/schema_input.json"
            )
            .withIndex().collect{ elem, idx -> [idx] + elem }
        )
        .map {
            idx, sample, reads1, reads2, single_end ->
            return [
                ['id': sample, 'idx': idx, 'single_end': single_end=='true'],
                (reads2 == []) ? [file(reads1)] : [file(reads1), file(reads2)],
            ]
        }
    
    if (params.reads_subsampling != -1) {
        BBMAP_SAMPLE_FASTX(reads_ch, params.reads_subsampling, true)
        reads_ch = BBMAP_SAMPLE_FASTX.out.fastx
    }

    // Fetch databases
    db_ch = channel
        .from(
            params.databases.collect { k, v ->
                if (v instanceof Map) {
                    if (v.containsKey('chunked') && v['chunked']) {
                        v.collect { k_, v_ ->
                            if (v_ instanceof Map) {
                                if (v_.containsKey('files')) {
                                    return [id: k, chunk_id: k_] + v_
                                }
                            }
                        }
                    } else if (v.containsKey('files')) {
                        return [id: k] + v
                    }
                }
            }.flatten()
        )
        .filter { it -> it }

    FETCHDB(db_ch, "${launchDir}/${params.databases.cache_path}")
    dbs_path_ch = FETCHDB.out.dbs

    dbs_path_ch
        .branch { meta, _fp ->
            genome_catalogue: meta.id == 'genome_catalogue'
        }
        .set { dbs }
    
    if (!params.skip_qc) {
        FASTP(
            reads_ch.map{ meta, reads -> [meta, reads, []]},
            false,
            false,
            false,
        )
        ch_versions = ch_versions.mix(FASTP.out.versions_fastp)
        reads_ch = FASTP.out.reads
        qc_stats = FASTP.out.fastp_json
    } else {
        qc_stats = channel.empty()
    }

    // Run fast genome profiling
    sylph_db = dbs.genome_catalogue
        .map { meta, fp ->
            file("${fp}/${meta.files.sylph}")
        }
        .first()
    SYLPH_PROFILE(reads_ch, sylph_db)
    // SYLPH_QUERY(reads_ch, sylph_db)

    genome_fp_lookup_table = dbs.genome_catalogue
        .map { meta, fp ->
            file("${fp}/${meta.files.filepath_lookup}")
        }
        .first()
    genomes_ch = channel.empty()
    contigs_ch = channel.empty()
    if (params.sourmash_genome_selector) {
        sourmash_db = dbs.genome_catalogue
            .map { meta, fp ->
                file("${fp}/${meta.files.sourmash}")
            }
            .first()
        SOURMASH_SKETCH(reads_ch)
        SOURMASH_GATHER(SOURMASH_SKETCH.out.signatures, sourmash_db, false, false, false, false)
    
    
        // Create bowtie2 index and align reads
        sourmash2fasta_ch = SOURMASH_GATHER.out.result
            .map{ meta,fp -> [meta, fp, file(params.genome_species)] }
    
        SOURMASH2FASTA(
            sourmash2fasta_ch,
            genome_fp_lookup_table
        )
    
        // filter out if no genomes/contigs
        genomes_ch = genomes_ch.mix(
            SOURMASH2FASTA.out.fasta
                .filter { _meta, fp -> fp.exists() & (fp.readLines().size() > 0) }
        )
        contigs_ch = contigs_ch.mix(
            SOURMASH2FASTA.out.contigs
                .filter { _meta, fp -> fp.exists() & (fp.readLines().size() > 0) }
        )
    }
    if (params.sylph_genome_selector) {
        SYLPH2FASTA(
            SYLPH_PROFILE.out.profile_out,
            genome_fp_lookup_table
        )

        // filter out if no genomes/contigs
        genomes_ch = genomes_ch.mix(
            SYLPH2FASTA.out.fasta
                .filter { _meta, fp -> fp.exists() & (fp.readLines().size() > 0) }
        )
        contigs_ch = contigs_ch.mix(
            SYLPH2FASTA.out.contigs
                .filter { _meta, fp -> fp.exists() & (fp.readLines().size() > 0) }
        )
    }

    BOWTIE2_BUILD(genomes_ch)

    align_in_ch = reads_ch
        .join(BOWTIE2_BUILD.out.index, remainder: true)
        .filter { _meta, _reads, index -> index }
        .join(contigs_ch)
        .join(genomes_ch)
        .multiMap{ meta, reads, index, genome_contigs, fasta -> 
            reads: [meta, reads]
            index: [meta, index]
            fasta: [meta, fasta]
            genome_contigs: [meta, genome_contigs]
            genome_species: [meta, file(params.genome_species)]
        }
    BOWTIE2_ALIGN_BAM2SQLITE(
        align_in_ch.reads, 
        align_in_ch.index, 
        align_in_ch.fasta, 
        align_in_ch.genome_contigs, 
        align_in_ch.genome_species, 
    )

    profile_ch = BOWTIE2_ALIGN_BAM2SQLITE.out.sqlite
        .map { meta, sqlite -> 
            [meta, sqlite, file(params.genome2cds)] 
        }
    if (params.greedy_profile) {
        SQLITE2PROFILE_GREEDY(profile_ch)
        // taxonomic_profile = SQLITE2PROFILE_GREEDY.out.species_profile
        // functional_profile = SQLITE2PROFILE_GREEDY.out.species_cds_profile
    }
    if (params.top_profile) {
        SQLITE2PROFILE_TOP(profile_ch)
        // taxonomic_profile = SQLITE2PROFILE_TOP.out.species_profile
        // functional_profile = SQLITE2PROFILE_TOP.out.species_cds_profile
    }

    MULTIQC(
        qc_stats,
        [],
        [],
        [],
        [],
        [],
    )

    emit:
    // sylph_profile = SYLPH_PROFILE.out.profile_out
    // sylph_query = SYLPH_QUERY.out.profile_out
    // sourmash_profile = SOURMASH_GATHER.out.result
    // taxonomic_profile = taxonomic_profile
    // functional_profile = functional_profile
    versions = ch_versions                         // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
