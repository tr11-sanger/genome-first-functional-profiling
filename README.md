# Genome-First Functional Profiling (GFFP)

A Nextflow DSL2 pipeline for fast, accurate functional profiling of metagenomes. GFFP maps whole-genome sequencing reads to reference genomes and then to genes, enabling simultaneous taxonomic and functional characterisation of microbial communities.

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Quick Start](#quick-start)
- [Input](#input)
- [Databases](#databases)
- [Parameters](#parameters)
- [Output](#output)
- [Pipeline Steps](#pipeline-steps)
- [Project Structure](#project-structure)
- [Tools and Citations](#tools-and-citations)
- [Testing](#testing)
- [Credits and Authorship](#credits-and-authorship)
- [License](#license)

## Pipeline Overview

```
Reads ──> Fetch & Preprocess ──> QC ──> Fast Genome Profiling ──> Genome Selection
                                              (Sylph / Sourmash)        │
                                                                        v
           Functional Profiles <── SQLite Profiling <── Bowtie2 Alignment
```

1. **Read fetching and preprocessing** — Fetch reads from local or remote sources with automatic retry, optional subsampling, header standardisation, and read repair.
2. **Quality control** — Adapter trimming and quality filtering with fastp.
3. **Fast genome profiling** — Rapid taxonomic screening with Sylph (and optionally Sourmash) to identify which reference genomes are present in each sample.
4. **Genome selection** — Extract relevant reference genome sequences per sample from a pre-built genome catalogue.
5. **Competitive alignment** — Align reads against only the selected genomes using Bowtie2, with results piped directly into a SQLite database.
6. **Functional profiling** — Assign reads to species and coding sequences using winner-takes-all (top) or greedy multi-mapping strategies.
7. **Reporting** — Aggregate QC metrics with MultiQC.

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) >= 24.10.5
- A container runtime: [Singularity](https://sylabs.io/docs/) or [Docker](https://www.docker.com/)

### Run

```bash
# Local execution with Docker
nextflow run main.nf \
    -profile local,singularity \
    --samplesheet samplesheet.csv \
    --outdir results

# LSF farm cluster
nextflow run main.nf \
    -profile farm,singularity \
    --samplesheet samplesheet.csv \
    --outdir results

# Test profile (minimal test data)
nextflow run main.nf \
    -profile test,singularity \
    --outdir test_results
```

## Input

### Samplesheet

A CSV file passed via `--samplesheet` with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Sample identifier (no spaces, min 3 characters) |
| `fastq1` | Yes | Path or URL to forward reads (FASTQ/FASTA, optionally gzipped) |
| `fastq2` | No | Path or URL to reverse reads (paired-end only; `.fq.gz` or `.fastq.gz`) |
| `single_end` | Yes | `true` for single-end, `false` for paired-end |

Example:

```csv
sample,fastq1,fastq2,single_end
SAMPLE_A,/data/reads/A_R1.fastq.gz,/data/reads/A_R2.fastq.gz,false
SAMPLE_B,ftp://server/B.fastq.gz,,true
```

Read paths can be local filesystem paths or remote URLs (HTTP/FTP). Remote reads are fetched with curl, with configurable retry logic and optional soft-fail mode.

## Databases

The pipeline requires a **genome catalogue** containing pre-built indices for genome profiling and selection. The catalogue is configured under `params.databases.genome_catalogue`:

```groovy
params {
    genome_species = '/path/to/genome_species.tsv'   // genome-to-species mapping
    genome2cds     = '/path/to/genome2cds.txt'       // genome-to-CDS mapping file list

    databases {
        cache_path = 'download_cache/databases'

        genome_catalogue {
            remote_path = 'https://example.com/catalogue.tar.gz'  // or leave empty for local
            local_path  = '/path/to/local/catalogue'              // use if pre-downloaded
            files {
                sylph                  = 'catalogue.syldb'
                sourmash               = 'catalogue.sbt.zip'
                sylph_filepath_lookup  = 'sylph_lookup.csv'
                sourmash_filepath_lookup = 'sourmash_lookup.csv'
            }
        }
    }
}
```

### Required database files

| File | Description |
|------|-------------|
| `sylph` | Pre-sketched Sylph database for fast genome profiling |
| `sylph_filepath_lookup` | CSV mapping genome accessions to FASTA file paths |
| `genome_species` | TSV mapping genome accessions to species names |
| `genome2cds` | Text file listing paths to genome-to-CDS mapping TSV files |

### Optional database files (Sourmash pathway)

| File | Description |
|------|-------------|
| `sourmash` | Sourmash SBT index for genome profiling |
| `sourmash_filepath_lookup` | CSV mapping genome accessions to FASTA file paths |

Databases can be provided locally (`local_path`) or downloaded automatically from a remote source (`remote_path`). Downloaded databases are cached in `params.databases.cache_path` and reused on subsequent runs unless `--force_download_dbs` is set.

## Parameters

### Core parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | (required) | Path to input CSV samplesheet |
| `--outdir` | `results` | Output directory |
| `--run_profiling` | `true` | Enable functional profiling |
| `--skip_qc` | `false` | Skip fastp quality control |
| `--skip_standardise` | `true` | Skip read header standardisation and repair |
| `--reads_subsampling` | `-1` | Subsample to N reads per sample (`-1` = disabled) |

### Alignment and profiling

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bowtie2_num_secondary_mappings` | `8` | Number of secondary alignments to report (`-k`) |
| `--min_ani` | `0.95` | Minimum average nucleotide identity for read assignment |
| `--remove_paired_suffix` | `1` | Remove paired-end suffix from read names (`1` = yes, `0` = no) |

### Genome selector

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sylph_genome_selector` | `true` | Use Sylph for genome selection |
| `--sourmash_genome_selector` | `false` | Use Sourmash for genome selection |

Both selectors can be enabled simultaneously; their selected genomes are merged.

### Remote read fetching

Parameters under `--ftp_fetch.*` control the download behaviour for remote read URLs:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ftp_fetch.max_retries` | `5` | Maximum download retry attempts |
| `--ftp_fetch.wait_retry` | `60` | Seconds to wait between retries |
| `--ftp_fetch.timeout` | `60` | Connection timeout in seconds |
| `--ftp_fetch.resume` | `true` | Resume partial downloads with curl |
| `--ftp_fetch.soft_fail` | `true` | Create empty outputs on failure instead of terminating |

### Database management

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--force_download_dbs` | `false` | Force re-download of all databases |
| `--download_dbs` | `false` | Enable database downloading |
| `--genome_species` | `''` | Path to genome-to-species mapping TSV |
| `--genome2cds` | `''` | Path to genome-to-CDS mapping file list |

## Output

Results are organised per sample in bucketed subdirectories:

```
results/
├── 0/                              # bucket for samples 0-999
│   └── SAMPLE_A/
│       ├── qc/
│       │   ├── SAMPLE_A_qc.fastp.json
│       │   └── SAMPLE_A_multiqc_report.html
│       ├── sylph/
│       │   └── SAMPLE_A_profile.tsv
│       ├── bowtie2/
│       │   └── SAMPLE_A.sqlite.gz
│       └── profiles/
│           └── top/
│               ├── species_coverage.tsv.gz
│               ├── genome_coverage.tsv.gz
│               ├── species_cds_coverage.tsv.gz
│               └── species_index.tsv.gz
├── 1/                              # bucket for samples 1000-1999
│   └── ...
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.html
```

### Key output files

| File | Description |
|------|-------------|
| `*_qc.fastp.json` | Per-sample QC statistics from fastp |
| `*_profile.tsv` | Sylph genome profiling results |
| `*.sqlite.gz` | Gzipped SQLite database of read-to-genome alignments |
| `species_coverage.tsv.gz` | Per-species coverage depth and breadth |
| `genome_coverage.tsv.gz` | Per-genome coverage depth and breadth |
| `species_cds_coverage.tsv.gz` | CDS-level functional coverage per species |
| `species_index.tsv.gz` | Species index for CDS profiles |

## Pipeline Steps

### 1. Read fetching and preprocessing

Reads are fetched from local paths or remote URLs using composite modules that combine download with processing in a single job:

- **CURL_FETCH_BBMAP_SAMPLE_FASTX** — Fetch + subsample (when `--reads_subsampling` is set)
- **CURL_FETCH_BBMAP_REFORMAT_STANDARDISE** — Fetch + standardise headers (default path)

Both modules include automatic retry logic for unreliable remote sources and a soft-fail mode that produces empty output files instead of terminating the pipeline.

When `--skip_standardise` is `false`, additional steps run:
- **BBMAP_REFORMAT_STANDARDISE** — Standardise read headers, de-interleave paired reads
- **BBMAP_REPAIR** — Remove unpaired reads from paired-end data

### 2. Quality control

**fastp** trims adapters and filters low-quality reads (skip with `--skip_qc`). Per-sample QC reports are aggregated by **MultiQC**.

### 3. Fast genome profiling

**Sylph** performs rapid k-mer-based profiling against a pre-sketched genome catalogue to estimate which reference genomes are present in each sample. Optionally, **Sourmash** can be used as an alternative or complementary selector.

### 4. Genome selection and alignment

Selected genome accessions are looked up in a filepath table, and their FASTA sequences are concatenated into a per-sample reference. **Bowtie2** builds an index and aligns reads against only the relevant genomes, reporting multiple secondary alignments (`-k`). The SAM output is piped directly into **bam2sqlite.py**, which parses alignments, computes ANI, and writes a gzipped SQLite database.

### 5. Functional profiling

**sqlite2profile_top2.py** reads the SQLite alignment database and generates coverage profiles using a winner-takes-all (top) read assignment strategy with iterative coverage-ratio filtering:

- **Species coverage** — depth and breadth per species
- **Genome coverage** — depth and breadth per genome
- **CDS coverage** — functional coverage at the coding sequence level

A greedy (multi-mapping) assignment strategy is also available via `SQLITE2PROFILE_GREEDY`.

## Project Structure

```
.
├── main.nf                         # Entry point
├── nextflow.config                 # Global configuration
├── workflows/
│   └── gffp.nf                     # Main workflow logic
├── modules/
│   ├── nf-core/                    # Standard nf-core modules
│   │   ├── bowtie2/                #   (build, align)
│   │   ├── fastp/
│   │   ├── multiqc/
│   │   ├── sylph/profile/
│   │   ├── sourmash/               #   (sketch, gather)
│   │   └── bbmap/repair/
│   └── local/                      # Custom pipeline modules
│       ├── bowtie2/align_bam2sqlite/
│       ├── curl_fetch_bbmap_reformat_standardise/
│       ├── curl_fetch_bbmap_sample_fastx/
│       ├── sqlite2profile_top/
│       ├── sqlite2profile_greedy/
│       ├── sylph2fasta/
│       ├── sourmash2fasta/
│       └── ...
├── subworkflows/
│   └── local/fetchdb/              # Database download/caching
├── bin/                            # Python scripts
│   ├── bam2sqlite.py               # SAM-to-SQLite converter
│   ├── sqlite2profile_top2.py      # Top assignment profiler
│   ├── sqlite2profile_greedy.py    # Greedy assignment profiler
│   ├── sylph2fasta.py              # Genome selector (Sylph)
│   └── sourmash2fasta.py           # Genome selector (Sourmash)
├── conf/
│   ├── modules.config              # Per-process args, publishDir, resources
│   ├── farm.config                 # LSF cluster profile
│   └── local.config                # Local Docker profile
├── assets/
│   └── schema_input.json           # Samplesheet validation schema
└── tests/
    └── default.nf.test             # Pipeline-level nf-test
```

## Tools and Citations

### Core tools

- **Sylph** — Ultra-fast genome profiling of metagenomes via k-mer sketching.
  > Sun J, Liao S, Shi M, et al. Sylph: taxonomic profiling with ANI through k-mer sketching. *Nat Biotechnol.* 2024. doi: [10.1038/s41587-024-02412-y](https://doi.org/10.1038/s41587-024-02412-y)

- **Bowtie2** — Fast and sensitive read alignment.
  > Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nat Methods.* 2012;9(4):357-359. doi: [10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)

- **Sourmash** — k-mer-based sequence analysis and comparison.
  > Brown CT, Irber L. sourmash: a library for MinHash sketching of DNA. *J Open Source Softw.* 2016;1(5):27. doi: [10.21105/joss.00027](https://doi.org/10.21105/joss.00027)

- **fastp** — Fast all-in-one preprocessing for FASTQ files.
  > Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics.* 2018;34(17):i884-i890. doi: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

- **BBMap/BBTools** — Read preprocessing, standardisation, repair, and subsampling.
  > Bushnell B. BBTools software package. 2014. Available at: [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)

- **MultiQC** — Aggregate QC reports across samples.
  > Ewels P, Magnusson M, Lundin S, Kaller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics.* 2016;32(19):3047-3048. doi: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

### Framework

- **Nextflow** — Workflow management.
  > Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol.* 2017;35(4):316-319. doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

- **nf-core** — Community-curated pipeline framework.
  > Ewels PA, Peltzer A, Fillinger S, et al. The nf-core framework for community-curated bioinformatics pipelines. *Nat Biotechnol.* 2020;38(3):276-278. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)

## Testing

The pipeline uses the [nf-test](https://www.nf-test.com/) framework:

```bash
# Run all tests
nf-test test

# Run a specific test
nf-test test tests/default.nf.test
```

Tests use the `test` profile with minimal test data. nf-core module tests are excluded via `nf-test.config`.

## Credits and Authorship

**Author:** Timothy J. Rozday ([@timrozday](https://github.com/timrozday))

Built with the [nf-core](https://nf-co.re) framework.

## License

This project is licensed under the [MIT License](LICENSE).
