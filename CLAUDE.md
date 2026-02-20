# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Genome-first functional profiling (GFFP) pipeline built with **Nextflow DSL2** on the **nf-core** framework. It maps whole genome sequencing reads to genomes and then to genes for fast, accurate functional profiling of metagenomes.

**Version:** 0.1 (early-stage development)

## Running the Pipeline

```bash
# Local with Docker
nextflow run main.nf -profile local,singularity --samplesheet samplesheet.csv --outdir results

# LSF farm cluster
nextflow run main.nf -profile farm,singularity --samplesheet samplesheet.csv --outdir results

# Test profile (uses minimal test data)
nextflow run main.nf -profile test,singularity --outdir test_results
```

Required: Nextflow >= 24.10.5 and a container runtime (Singularity/Docker).

## Running Tests

Uses **nf-test** framework:
```bash
nf-test test                        # all tests
nf-test test tests/default.nf.test  # specific test
```

nf-core module tests are excluded via `nf-test.config`. Tests use the `test` profile.

## Architecture

### Pipeline Flow (workflows/gffp.nf)

1. **Input parsing** — CSV samplesheet validated against `assets/schema_input.json`
2. **Read fetching & preprocessing** — Two paths depending on `params.reads_subsampling`:
   - **Subsampling enabled** (`!= -1`): CURL_FETCH_BBMAP_SAMPLE_FASTX fetches remote/local reads and subsamples in one step, then optionally BBMAP_REFORMAT_STANDARDISE + BBMAP_REPAIR
   - **No subsampling** (default): CURL_FETCH_BBMAP_REFORMAT_STANDARDISE fetches and standardises reads in one step, then optionally BBMAP_REPAIR
   - Both CURL_FETCH modules support retry logic and soft-fail mode for unreliable remote sources
   - Standardisation/repair controlled by `params.skip_standardise`
3. **Database fetching** — FETCHDB subworkflow handles local/remote/cached databases
4. **QC** — FASTP (skippable via `params.skip_qc`)
5. **Fast genome profiling** — SYLPH_PROFILE against genome catalogue; optionally SOURMASH_GATHER
6. **Genome selection** — SYLPH2FASTA or SOURMASH2FASTA selects relevant genomes per sample
7. **Alignment** — BOWTIE2_BUILD + BOWTIE2_ALIGN_BAM2SQLITE (custom module combining alignment with BAM-to-SQLite conversion)
8. **Functional profiling** — SQLITE2PROFILE_TOP (top read assignment strategy)
9. **Reporting** — MULTIQC aggregates QC results

### Module Organization

- `modules/nf-core/` — Standard bioinformatics tools from nf-core/modules (bowtie2, fastp, sourmash, sylph, multiqc, etc.)
- `modules/local/` — Custom modules specific to this pipeline:
  - `curl_fetch_bbmap_reformat_standardise/` — Combined fetch + standardise with retry/soft-fail
  - `curl_fetch_bbmap_sample_fastx/` — Combined fetch + subsample with retry/soft-fail
  - `bowtie2/align_bam2sqlite/` — Combined alignment + SQLite DB creation
  - `sqlite2profile_top/`, `sqlite2profile_greedy/` — Coverage profiling from SQLite
  - `sylph2fasta/`, `sourmash2fasta/` — Genome selection from profiling results
  - `bbmap/reformat_standardise/`, `bbmap_sample_fastx/` — BBMap wrappers (standalone)
  - `bam2csv_top/`, `bam2csv_greedy/` — Legacy BAM-to-CSV profiling (not used in main workflow)
  - `fetchunzip/`, `chunkfastx/`, `concatenate/`, `seqstats/` — Utilities
- `subworkflows/local/fetchdb/` — Database download/caching subworkflow
- `bin/` — Python scripts used by local modules (bam2sqlite.py, sqlite2profile_*.py, sylph2fasta.py, sourmash2fasta.py)

### Key Design Patterns

- **Meta maps**: All channels carry `[meta, data]` tuples where meta contains `id`, `idx`, `single_end`
- **Output directory structure**: `${outdir}/${(int)(meta.idx/1000)}/${meta.id}/...` — samples are bucketed into subdirectories of 1000
- **Two assignment strategies**: "top" (winner-takes-all) and "greedy" for read-to-genome assignment
- **Genome selectors**: Controlled by `params.sylph_genome_selector` (default: true) and `params.sourmash_genome_selector` (default: false)
- **Version reporting**: Newer modules (CURL_FETCH_*, FASTP) use `topic: versions` for automatic collection; older modules emit `path "versions.yml"` collected via `ch_versions`
- **Soft-fail pattern**: CURL_FETCH modules support `params.ftp_fetch.soft_fail` — on failure, they create empty output files and a `.status` file instead of terminating

### Configuration Cascade

`nextflow.config` → `conf/modules.config` (per-process args, publishDir, resources) → profile-specific config (`conf/farm.config` or `conf/local.config`)

Strict mode is enabled (`nextflow.enable.strict = true`). Container registry: quay.io.

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bowtie2_num_secondary_mappings` | 8 | Secondary alignments to report |
| `min_ani` | 0.95 | Minimum average nucleotide identity |
| `skip_qc` | false | Skip FASTP QC |
| `skip_standardise` | true | Skip read standardization |
| `sylph_genome_selector` | true | Use Sylph for genome selection |
| `sourmash_genome_selector` | false | Use Sourmash for genome selection |
| `run_profiling` | true | Enable functional profiling |
| `reads_subsampling` | -1 | Subsample reads (-1 = disabled) |
| `remove_paired_suffix` | 1 | Remove paired-end suffix from read names (1=yes, 0=no) |
| `genome_species` | '' | Path to genome-to-species mapping file |
| `genome2cds` | '' | Path to genome-to-CDS mapping file list |
| `force_download_dbs` | false | Force re-download of databases |
| `download_dbs` | false | Enable database downloading |
| `ftp_fetch.max_retries` | 5 | Max download retry attempts |
| `ftp_fetch.wait_retry` | 60 | Seconds between retries |
| `ftp_fetch.timeout` | 60 | Connection timeout in seconds |
| `ftp_fetch.resume` | true | Enable curl resume for partial downloads |
| `ftp_fetch.soft_fail` | true | Create empty outputs on download failure instead of failing |
| `ftp_fetch.delete_original` | false | Delete original downloaded files after successful processing |

### Plugins

- `nf-schema@2.2.0` — Parameter/samplesheet validation
- `nf-amazon` — AWS execution support

### Remaining Tech Debt

- **MULTIQC runs per-sample** rather than aggregating across all samples (issue #10)
- **Incomplete version collection** — versions from BOWTIE2_BUILD, BOWTIE2_ALIGN_BAM2SQLITE, SYLPH_PROFILE, SYLPH2FASTA, SOURMASH2FASTA, SQLITE2PROFILE_TOP are not collected (issue #11)
- **Unused nf-core modules** installed but not referenced: pyrodigal, fastqc, minimap2, hmmer, sourmash/index, bbmap/bbmerge (issue #13)
- **Duplicated preprocessing logic** — the subsampling/non-subsampling branches share standardise+repair code (issue #15)
- **Legacy scripts in `bin/`** — multiple generations of profiling scripts (bam2csv.py, sqlite2profile_top.py, sqlite2profile_top1.py) are no longer used (issue #18)
- **`params.min_ani` hardcoded in BOWTIE2_ALIGN_BAM2SQLITE** script block instead of `task.ext.args` (issue #9)
