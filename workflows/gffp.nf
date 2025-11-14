/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { FETCHDB } from '../subworkflows/local/fetchdb/main'
include { SYLPH_PROFILE } from '../modules/nf-core/sylph/profile/main'
include { SYLPH_QUERY } from '../modules/local/sylph/query/main'
include { SOURMASH_GATHER } from '../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH } from '../modules/nf-core/sourmash/sketch/main'
include { SOURMASH2FASTA } from '../modules/local/sourmash2fasta/main'
include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN } from '../modules/nf-core/bowtie2/align/main'
include { BAM2CSV_TOP } from '../modules/local/bam2csv_top/main'
include { BAM2CSV_GREEDY } from '../modules/local/bam2csv_greedy/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GFFP {
    main:
    ch_versions = Channel.empty()

    // Parse samplesheet and fetch reads
    reads_ch = Channel.fromList(
            samplesheetToList(
                params.samplesheet, 
                "${workflow.projectDir}/assets/schema_input.json"
            )
            .withIndex().collect{ elem, idx -> [idx] + elem }
        )
        .map {
            idx, sample, reads1, reads2, single_end ->
            return [
                ['id': sample, 'idx': idx, 'single_end': single_end],
                (reads2 == []) ? [file(reads1)] : [file(reads1), file(reads2)],
            ]
        }

    // Fetch databases
    db_ch = Channel
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
        .filter { it }

    FETCHDB(db_ch, "${launchDir}/${params.databases.cache_path}")
    dbs_path_ch = FETCHDB.out.dbs

    dbs_path_ch
        .branch { meta, _fp ->
            genome_catalogue: meta.id == 'genome_catalogue'
        }
        .set { dbs }
    

    // Run fast genome profiling
    sylph_db = dbs.genome_catalogue
        .map { meta, fp ->
            file("${fp}/${meta.files.sylph}")
        }
        .first()
    SYLPH_PROFILE(reads_ch, sylph_db)
    SYLPH_QUERY(reads_ch, sylph_db)

    sourmash_db = dbs.genome_catalogue
        .map { meta, fp ->
            file("${fp}/${meta.files.sourmash}")
        }
        .first()
    SOURMASH_SKETCH(reads_ch)
    SOURMASH_GATHER(SOURMASH_SKETCH.out.signatures, sourmash_db, false, false, false, false)


    // Create bowtie2 index and align reads
    genome_fp_lookup_table = dbs.genome_catalogue
        .map { meta, fp ->
            file("${fp}/${meta.files.filepath_lookup}")
        }
        .first()
    sourmash2fasta_ch = SOURMASH_GATHER.out.result
        .map{ meta,fp -> [meta, fp, file(params.genome_species)] }

    SOURMASH2FASTA(
        sourmash2fasta_ch,
        genome_fp_lookup_table
    )

    // filter out if no genomes
    genomes_ch = SOURMASH2FASTA.out.fasta
        .filter { _meta, fp -> fp.exists() & (fp.readLines().size() > 0) }

    BOWTIE2_BUILD(genomes_ch)

    align_in_ch = reads_ch
        .join(BOWTIE2_BUILD.out.index, remainder: true)
        .filter { _meta, _reads, index -> index }
        .join(genomes_ch)
        .multiMap{ meta, reads, index, fasta -> 
            reads: [meta, reads]
            index: [meta, index]
            fasta: [meta, fasta]
        }
    BOWTIE2_ALIGN(align_in_ch.reads, align_in_ch.index, align_in_ch.fasta, false, false)

    contigs_ch = SOURMASH2FASTA.out.contigs
        .filter { _meta, fp -> fp.exists() & (fp.readLines().size() > 0) }
    profile_ch = BOWTIE2_ALIGN.out.bam
        .join(contigs_ch)
        .join(genomes_ch)    
        .map { meta, bam, genome_contigs, refs -> 
            [meta, bam, file(params.genome2cds), genome_contigs, file(params.genome_species), refs] 
        }
    if (params.greedy_read_reassignment) {
        BAM2CSV_GREEDY(profile_ch, false)
        taxonomic_profile = BAM2CSV_GREEDY.out.species_profile
        functional_profile = BAM2CSV_GREEDY.out.species_cds_profile
    } else {
        BAM2CSV_TOP(profile_ch, false)
        taxonomic_profile = BAM2CSV_TOP.out.species_profile
        functional_profile = BAM2CSV_TOP.out.species_cds_profile
    }

    emit:
    sylph_profile = SYLPH_PROFILE.out.profile_out
    sylph_query = SYLPH_QUERY.out.profile_out
    sourmash_profile = SOURMASH_GATHER.out.result
    taxonomic_profile = taxonomic_profile
    functional_profile = functional_profile
    versions = ch_versions                         // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
